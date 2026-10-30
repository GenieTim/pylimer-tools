#include "MetalForceBalanceSolver.h"
#include "../utils/utilityMacros.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

namespace pylimer_tools::sim::mehp {

static const char* metalForceBalanceShaders = R"(
#include <metal_stdlib>
using namespace metal;

kernel void csr_spmv(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* x [[buffer(3)]],
    device float* y [[buffer(4)]],
    constant uint& numRows [[buffer(5)]],
    uint row [[thread_position_in_grid]])
{
    if (row >= numRows) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float sum = 0.0f;
    for (uint k = start; k < end; ++k) {
        sum += values[k] * x[col_indices[k]];
    }
    y[row] = sum;
}

// Compute nextX = x - alpha * gradient, and deltaX = -alpha * gradient
kernel void bb_update_x(
    device const float* x [[buffer(0)]],
    device const float* gradient [[buffer(1)]],
    device float* nextX [[buffer(2)]],
    device float* deltaX [[buffer(3)]],
    constant float& alpha [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (id >= n) return;
    float dx = -alpha * gradient[id];
    deltaX[id] = dx;
    nextX[id] = x[id] + dx;
}

// Compute deltaG = nextG - gradient
kernel void bb_calc_delta_g(
    device const float* nextG [[buffer(0)]],
    device const float* gradient [[buffer(1)]],
    device float* deltaG [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    uint id [[thread_position_in_grid]])
{
    if (id >= n) return;
    deltaG[id] = nextG[id] - gradient[id];
}

// Compute gradient = A*x - b
kernel void bb_calc_gradient(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* x [[buffer(3)]],
    device const float* b [[buffer(4)]],
    device float* gradient [[buffer(5)]],
    constant uint& n [[buffer(6)]],
    uint row [[thread_position_in_grid]])
{
    if (row >= n) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float sum = 0.0f;
    for (uint k = start; k < end; ++k) {
        sum += values[k] * x[col_indices[k]];
    }
    gradient[row] = sum - b[row];
}

// CG update: x += alpha * p, r -= alpha * Ap
kernel void cg_update_x_r(
    device float* x [[buffer(0)]],
    device float* r [[buffer(1)]],
    device const float* p [[buffer(2)]],
    device const float* Ap [[buffer(3)]],
    constant float& alpha [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (id >= n) return;
    x[id] += alpha * p[id];
    r[id] -= alpha * Ap[id];
}

// CG update: p = r + beta * p
kernel void cg_update_p(
    device float* p [[buffer(0)]],
    device const float* r [[buffer(1)]],
    constant float& beta [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    uint id [[thread_position_in_grid]])
{
    if (id >= n) return;
    p[id] = r[id] + beta * p[id];
}

// Dot product with parallel reduction
kernel void dot_product(
    device const float* a [[buffer(0)]],
    device const float* b [[buffer(1)]],
    device float* partialSums [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    float val = (gid < n) ? a[gid] * b[gid] : 0.0f;
    sharedData[tid] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    for (uint s = tg_size / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sharedData[tid] += sharedData[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    
    if (tid == 0) {
        partialSums[tg_id] = sharedData[0];
    }
}
)";

class MetalFBSBackend
{
public:
  static MetalFBSBackend& instance()
  {
    static MetalFBSBackend backend;
    return backend;
  }

  bool isAvailable() const { return available; }
  id<MTLDevice> getDevice() const { return device; }
  id<MTLCommandQueue> getQueue() const { return queue; }

  id<MTLComputePipelineState> getSpmvPipe() const { return pipeSpmv; }
  id<MTLComputePipelineState> getDotPipe() const { return pipeDot; }
  id<MTLComputePipelineState> getBbUpdateXPipe() const { return pipeBbUpdateX; }
  id<MTLComputePipelineState> getBbCalcGradPipe() const
  {
    return pipeBbCalcGrad;
  }
  id<MTLComputePipelineState> getBbCalcDeltaGPipe() const
  {
    return pipeBbCalcDeltaG;
  }
  id<MTLComputePipelineState> getCgUpdateXRPipe() const
  {
    return pipeCgUpdateXR;
  }
  id<MTLComputePipelineState> getCgUpdatePPipe() const { return pipeCgUpdateP; }

private:
  bool available = false;
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> pipeSpmv = nil;
  id<MTLComputePipelineState> pipeDot = nil;
  id<MTLComputePipelineState> pipeBbUpdateX = nil;
  id<MTLComputePipelineState> pipeBbCalcGrad = nil;
  id<MTLComputePipelineState> pipeBbCalcDeltaG = nil;
  id<MTLComputePipelineState> pipeCgUpdateXR = nil;
  id<MTLComputePipelineState> pipeCgUpdateP = nil;

  MetalFBSBackend()
  {
    @autoreleasepool {
      device = MTLCreateSystemDefaultDevice();
      if (!device)
        return;

      queue = [device newCommandQueue];
      if (!queue)
        return;

      NSError* error = nil;
      id<MTLLibrary> library = [device
        newLibraryWithSource:[NSString
                               stringWithUTF8String:metalForceBalanceShaders]
                     options:nil
                       error:&error];
      if (!library) {
        std::cerr << "Failed to compile Metal force balance shaders: " <<
          [[error localizedDescription] UTF8String] << std::endl;
        return;
      }

      pipeSpmv = [device
        newComputePipelineStateWithFunction:[library
                                              newFunctionWithName:@"csr_spmv"]
                                      error:&error];
      pipeDot = [device newComputePipelineStateWithFunction:
                          [library newFunctionWithName:@"dot_product"]
                                                      error:&error];
      pipeBbUpdateX = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"bb_update_x"]
                                                            error:&error];
      pipeBbCalcGrad =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_gradient"]
                                              error:&error];
      pipeBbCalcDeltaG =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_delta_g"]
                                              error:&error];
      pipeCgUpdateXR = [device newComputePipelineStateWithFunction:
                                 [library newFunctionWithName:@"cg_update_x_r"]
                                                             error:&error];
      pipeCgUpdateP = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"cg_update_p"]
                                                            error:&error];

      if (pipeSpmv && pipeDot && pipeBbUpdateX && pipeBbCalcGrad &&
          pipeBbCalcDeltaG && pipeCgUpdateXR && pipeCgUpdateP) {
        available = true;
      }
    }
  }
};

bool
isMetalAvailable()
{
  return MetalFBSBackend::instance().isAvailable();
}

Eigen::VectorXd
solveSLEMetalBarzilaiBorwein(const Eigen::SparseMatrix<double>& A,
                             const Eigen::VectorXd& b,
                             double learningRate,
                             double tolerance,
                             int maxIterations,
                             int& iterationsDone,
                             bool shortAlpha,
                             const Eigen::VectorXd& initialX)
{
  Eigen::Index n = b.size();
  if (n == 0) {
    iterationsDone = 0;
    return Eigen::VectorXd::Zero(0);
  }
  if (b.isZero(1e-9)) {
    iterationsDone = 1;
    return Eigen::VectorXd::Zero(n);
  }

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLDevice> device = backend.getDevice();
    id<MTLCommandQueue> queue = backend.getQueue();
    id<MTLComputePipelineState> pUpdateX = backend.getBbUpdateXPipe();
    id<MTLComputePipelineState> pCalcGrad = backend.getBbCalcGradPipe();
    id<MTLComputePipelineState> pCalcDeltaG = backend.getBbCalcDeltaGPipe();
    id<MTLComputePipelineState> pDot = backend.getDotPipe();

    uint32_t N = static_cast<uint32_t>(n);
    uint32_t nVal = N;

    // Convert sparse matrix A to CSR in float
    std::vector<uint32_t> row_ptrs(N + 1);
    std::vector<uint32_t> col_indices(static_cast<size_t>(A.nonZeros()));
    std::vector<float> values(static_cast<size_t>(A.nonZeros()));
    for (uint32_t i = 0; i <= N; ++i)
      row_ptrs[i] = static_cast<uint32_t>(A.outerIndexPtr()[i]);
    for (Eigen::Index i = 0; i < A.nonZeros(); ++i) {
      col_indices[static_cast<size_t>(i)] =
        static_cast<uint32_t>(A.innerIndexPtr()[i]);
      values[static_cast<size_t>(i)] = static_cast<float>(A.valuePtr()[i]);
    }

    std::vector<float> b_float(N);
    for (uint32_t i = 0; i < N; ++i)
      b_float[i] = static_cast<float>(b[i]);

    id<MTLBuffer> bufRowPtr =
      [device newBufferWithBytes:row_ptrs.data()
                          length:sizeof(uint32_t) * row_ptrs.size()
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufColIdx =
      [device newBufferWithBytes:col_indices.data()
                          length:sizeof(uint32_t) * col_indices.size()
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufValues =
      [device newBufferWithBytes:values.data()
                          length:sizeof(float) * values.size()
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufB =
      [device newBufferWithBytes:b_float.data()
                          length:sizeof(float) * N
                         options:MTLResourceStorageModeShared];

    id<MTLBuffer> bufX =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufNextX =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufG =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufNextG =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufDeltaX =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufDeltaG =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];

    float* xPtr = static_cast<float*>([bufX contents]);
    if (initialX.size() == n) {
      for (uint32_t i = 0; i < N; ++i)
        xPtr[i] = static_cast<float>(initialX[i]);
    } else {
      std::memset(xPtr, 0, sizeof(float) * N);
    }

    uint32_t tg_size = 256;
    uint32_t num_tgs = (N + tg_size - 1) / tg_size;
    id<MTLBuffer> bufPartial =
      [device newBufferWithLength:sizeof(float) * num_tgs
                          options:MTLResourceStorageModeShared];

    auto gpu_dot = [&](id<MTLBuffer> buf1, id<MTLBuffer> buf2) -> float {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:pDot];
      [enc setBuffer:buf1 offset:0 atIndex:0];
      [enc setBuffer:buf2 offset:0 atIndex:1];
      [enc setBuffer:bufPartial offset:0 atIndex:2];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float* parts = static_cast<float*>([bufPartial contents]);
      float sum = 0.0f;
      for (uint32_t i = 0; i < num_tgs; ++i)
        sum += parts[i];
      return sum;
    };

    // Compute initial gradient: g = A*x - b
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:pCalcGrad];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufG offset:0 atIndex:5];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:6];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    float initGNormSq = gpu_dot(bufG, bufG);
    float initGNorm = std::sqrt(initGNormSq);
    if (initGNorm < 1e-12f) {
      iterationsDone = 0;
      return initialX.size() == n ? initialX : Eigen::VectorXd::Zero(n);
    }

    float alpha = static_cast<float>(learningRate);
    int iter = 0;
    float curGNorm = initGNorm;

    while (iter < maxIterations &&
           (curGNorm / initGNorm) > static_cast<float>(tolerance)) {
      // Encode:
      // 1. nextX = x - alpha * g, deltaX = -alpha * g
      // 2. nextG = A * nextX - b
      // 3. deltaG = nextG - g
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

      [enc setComputePipelineState:pUpdateX];
      [enc setBuffer:bufX offset:0 atIndex:0];
      [enc setBuffer:bufG offset:0 atIndex:1];
      [enc setBuffer:bufNextX offset:0 atIndex:2];
      [enc setBuffer:bufDeltaX offset:0 atIndex:3];
      [enc setBytes:&alpha length:sizeof(float) atIndex:4];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:5];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pCalcGrad];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufNextX offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufNextG offset:0 atIndex:5];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:6];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pCalcDeltaG];
      [enc setBuffer:bufNextG offset:0 atIndex:0];
      [enc setBuffer:bufG offset:0 atIndex:1];
      [enc setBuffer:bufDeltaG offset:0 atIndex:2];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float sTy = gpu_dot(bufDeltaX, bufDeltaG);
      float yTy = gpu_dot(bufDeltaG, bufDeltaG);
      float sTs = gpu_dot(bufDeltaX, bufDeltaX);
      float nextGNormSq = gpu_dot(bufNextG, bufNextG);
      curGNorm = std::sqrt(nextGNormSq);

      float alphaShort = (yTy > 1e-20f && std::isfinite(sTy / yTy))
                           ? (sTy / yTy)
                           : static_cast<float>(learningRate);
      float alphaLong = (std::abs(sTy) > 1e-20f && std::isfinite(sTs / sTy))
                          ? (sTs / sTy)
                          : static_cast<float>(learningRate);
      alpha = shortAlpha ? alphaShort : alphaLong;
      if (!std::isfinite(alpha)) {
        alpha = static_cast<float>(learningRate);
      }

      std::swap(bufX, bufNextX);
      std::swap(bufG, bufNextG);
      iter++;

      // Periodically recompute gradient to prevent float drift
      if (iter % 20 == 0) {
        id<MTLCommandBuffer> reCmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> reEnc = [reCmd computeCommandEncoder];
        [reEnc setComputePipelineState:pCalcGrad];
        [reEnc setBuffer:bufRowPtr offset:0 atIndex:0];
        [reEnc setBuffer:bufColIdx offset:0 atIndex:1];
        [reEnc setBuffer:bufValues offset:0 atIndex:2];
        [reEnc setBuffer:bufX offset:0 atIndex:3];
        [reEnc setBuffer:bufB offset:0 atIndex:4];
        [reEnc setBuffer:bufG offset:0 atIndex:5];
        [reEnc setBytes:&nVal length:sizeof(uint32_t) atIndex:6];
        [reEnc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [reEnc endEncoding];
        [reCmd commit];
        [reCmd waitUntilCompleted];
        curGNorm = std::sqrt(gpu_dot(bufG, bufG));
      }
    }

    iterationsDone = iter;
    float* finalXPtr = static_cast<float*>([bufX contents]);
    Eigen::VectorXd res(n);
    for (uint32_t i = 0; i < N; ++i) {
      res[i] = static_cast<double>(finalXPtr[i]);
    }
    return res;
  }
}

Eigen::VectorXd
solveSLEMetalCG(const Eigen::SparseMatrix<double>& A,
                const Eigen::VectorXd& b,
                double tolerance,
                int maxIterations,
                int& iterationsDone,
                const Eigen::VectorXd& initialX)
{
  Eigen::Index n = b.size();
  if (n == 0) {
    iterationsDone = 0;
    return Eigen::VectorXd::Zero(0);
  }
  if (b.isZero(1e-9)) {
    iterationsDone = 1;
    return Eigen::VectorXd::Zero(n);
  }

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLDevice> device = backend.getDevice();
    id<MTLCommandQueue> queue = backend.getQueue();
    id<MTLComputePipelineState> pSpmv = backend.getSpmvPipe();
    id<MTLComputePipelineState> pDot = backend.getDotPipe();
    id<MTLComputePipelineState> pUpdateXR = backend.getCgUpdateXRPipe();
    id<MTLComputePipelineState> pUpdateP = backend.getCgUpdatePPipe();

    uint32_t N = static_cast<uint32_t>(n);
    uint32_t nVal = N;

    std::vector<uint32_t> row_ptrs(N + 1);
    std::vector<uint32_t> col_indices(static_cast<size_t>(A.nonZeros()));
    std::vector<float> values(static_cast<size_t>(A.nonZeros()));
    for (uint32_t i = 0; i <= N; ++i)
      row_ptrs[i] = static_cast<uint32_t>(A.outerIndexPtr()[i]);
    for (Eigen::Index i = 0; i < A.nonZeros(); ++i) {
      col_indices[static_cast<size_t>(i)] =
        static_cast<uint32_t>(A.innerIndexPtr()[i]);
      values[static_cast<size_t>(i)] = static_cast<float>(A.valuePtr()[i]);
    }

    std::vector<float> b_float(N);
    for (uint32_t i = 0; i < N; ++i)
      b_float[i] = static_cast<float>(b[i]);

    id<MTLBuffer> bufRowPtr =
      [device newBufferWithBytes:row_ptrs.data()
                          length:sizeof(uint32_t) * row_ptrs.size()
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufColIdx =
      [device newBufferWithBytes:col_indices.data()
                          length:sizeof(uint32_t) * col_indices.size()
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufValues =
      [device newBufferWithBytes:values.data()
                          length:sizeof(float) * values.size()
                         options:MTLResourceStorageModeShared];

    id<MTLBuffer> bufX =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufR =
      [device newBufferWithBytes:b_float.data()
                          length:sizeof(float) * N
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufP =
      [device newBufferWithBytes:b_float.data()
                          length:sizeof(float) * N
                         options:MTLResourceStorageModeShared];
    id<MTLBuffer> bufAp =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];

    float* xPtr = static_cast<float*>([bufX contents]);
    if (initialX.size() == n) {
      for (uint32_t i = 0; i < N; ++i)
        xPtr[i] = static_cast<float>(initialX[i]);
      // If non-zero initialX, r = b - A*x
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:pSpmv];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufAp offset:0 atIndex:4];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:5];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float* rPtr = static_cast<float*>([bufR contents]);
      float* pPtr = static_cast<float*>([bufP contents]);
      float* apPtr = static_cast<float*>([bufAp contents]);
      for (uint32_t i = 0; i < N; ++i) {
        float val = b_float[i] - apPtr[i];
        rPtr[i] = val;
        pPtr[i] = val;
      }
    } else {
      std::memset(xPtr, 0, sizeof(float) * N);
    }

    uint32_t tg_size = 256;
    uint32_t num_tgs = (N + tg_size - 1) / tg_size;
    id<MTLBuffer> bufPartial =
      [device newBufferWithLength:sizeof(float) * num_tgs
                          options:MTLResourceStorageModeShared];

    auto gpu_dot = [&](id<MTLBuffer> buf1, id<MTLBuffer> buf2) -> float {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:pDot];
      [enc setBuffer:buf1 offset:0 atIndex:0];
      [enc setBuffer:buf2 offset:0 atIndex:1];
      [enc setBuffer:bufPartial offset:0 atIndex:2];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float* parts = static_cast<float*>([bufPartial contents]);
      float sum = 0.0f;
      for (uint32_t i = 0; i < num_tgs; ++i)
        sum += parts[i];
      return sum;
    };

    float rho = gpu_dot(bufR, bufR);
    float initial_rho = rho;
    if (initial_rho < 1e-24f) {
      iterationsDone = 0;
      return initialX.size() == n ? initialX : Eigen::VectorXd::Zero(n);
    }
    int iter = 0;
    float tolSq = static_cast<float>(tolerance * tolerance);

    for (; iter < maxIterations; ++iter) {
      if (rho / initial_rho < tolSq)
        break;

      // 1. Ap = A * p
      {
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pSpmv];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:bufP offset:0 atIndex:3];
        [enc setBuffer:bufAp offset:0 atIndex:4];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
      }

      float pAp = gpu_dot(bufP, bufAp);
      if (std::abs(pAp) <= 1e-20f || !std::isfinite(pAp)) {
        break;
      }
      float alpha = rho / pAp;

      // 2. x += alpha * p, r -= alpha * Ap, and dot(r, r)
      {
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pUpdateXR];
        [enc setBuffer:bufX offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufP offset:0 atIndex:2];
        [enc setBuffer:bufAp offset:0 atIndex:3];
        [enc setBytes:&alpha length:sizeof(float) atIndex:4];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        [enc setComputePipelineState:pDot];
        [enc setBuffer:bufR offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
      }

      float* parts = static_cast<float*>([bufPartial contents]);
      float rho_new = 0.0f;
      for (uint32_t i = 0; i < num_tgs; ++i)
        rho_new += parts[i];

      float beta = rho_new / rho;
      rho = rho_new;

      // 3. p = r + beta * p
      {
        id<MTLCommandBuffer> cmd = [queue commandBuffer];
        id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
        [enc setComputePipelineState:pUpdateP];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBytes:&beta length:sizeof(float) atIndex:2];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc endEncoding];
        [cmd commit];
        [cmd waitUntilCompleted];
      }
    }

    iterationsDone = iter;
    float* finalXPtr = static_cast<float*>([bufX contents]);
    Eigen::VectorXd res(n);
    for (uint32_t i = 0; i < N; ++i) {
      res[i] = static_cast<double>(finalXPtr[i]);
    }
    return res;
  }
}

} // namespace pylimer_tools::sim::mehp
