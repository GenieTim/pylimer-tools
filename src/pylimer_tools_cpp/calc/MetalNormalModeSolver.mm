#include "MetalNormalModeSolver.h"
#include "../utils/utilityMacros.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433
#endif

namespace pylimer_tools::calc {

static const char* metalShaders = R"(
#include <metal_stdlib>
using namespace metal;

// MARK: - Jacobi Eigenvalue Kernels

kernel void jacobi_rotate(
    device float* A [[buffer(0)]],
    device float* V [[buffer(1)]],
    constant uint& n [[buffer(2)]],
    constant uint& p [[buffer(3)]],
    constant uint& q [[buffer(4)]],
    constant float& c [[buffer(5)]],
    constant float& s [[buffer(6)]],
    constant uint& computeV [[buffer(7)]],
    uint k [[thread_position_in_grid]])
{
    if (k >= n) return;
    
    if (k != p && k != q) {
        float akp = A[k * n + p];
        float akq = A[k * n + q];
        float new_p = c * akp - s * akq;
        float new_q = s * akp + c * akq;
        
        A[k * n + p] = new_p;
        A[p * n + k] = new_p;
        A[k * n + q] = new_q;
        A[q * n + k] = new_q;
    }
    
    if (computeV) {
        float vkp = V[k * n + p];
        float vkq = V[k * n + q];
        V[k * n + p] = c * vkp - s * vkq;
        V[k * n + q] = s * vkp + c * vkq;
    }
}

// MARK: - KPM (Kernel Polynomial Method) Kernels

kernel void kpm_step_first(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* u0 [[buffer(3)]],
    device float* u1 [[buffer(4)]],
    constant float& c [[buffer(5)]],
    constant float& w [[buffer(6)]],
    constant uint& n [[buffer(7)]],
    uint row [[thread_position_in_grid]])
{
    if (row >= n) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float y = 0.0f;
    for (uint k = start; k < end; ++k) {
        y += values[k] * u0[col_indices[k]];
    }
    u1[row] = (y - c * u0[row]) / w;
}

kernel void kpm_step_recurrence(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* u_curr [[buffer(3)]],
    device const float* u_prev [[buffer(4)]],
    device float* u_next [[buffer(5)]],
    constant float& c [[buffer(6)]],
    constant float& w [[buffer(7)]],
    constant uint& n [[buffer(8)]],
    uint row [[thread_position_in_grid]])
{
    if (row >= n) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float y = 0.0f;
    for (uint k = start; k < end; ++k) {
        y += values[k] * u_curr[col_indices[k]];
    }
    u_next[row] = (2.0f / w) * (y - c * u_curr[row]) - u_prev[row];
}

kernel void kpm_dot(
    device const float* v [[buffer(0)]],
    device const float* u [[buffer(1)]],
    device float* partialSums [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    float val = (gid < n) ? v[gid] * u[gid] : 0.0f;
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

class MetalBackend
{
public:
  static MetalBackend& instance()
  {
    static MetalBackend backend;
    return backend;
  }

  bool isAvailable() const { return available; }
  id<MTLDevice> getDevice() const { return device; }
  id<MTLCommandQueue> getQueue() const { return queue; }

  id<MTLComputePipelineState> getJacobiPipe() const { return pipeJacobi; }
  id<MTLComputePipelineState> getKpmFirstPipe() const { return pipeKpmFirst; }
  id<MTLComputePipelineState> getKpmRecPipe() const { return pipeKpmRec; }
  id<MTLComputePipelineState> getKpmDotPipe() const { return pipeKpmDot; }

private:
  bool available = false;
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;
  id<MTLComputePipelineState> pipeJacobi = nil;
  id<MTLComputePipelineState> pipeKpmFirst = nil;
  id<MTLComputePipelineState> pipeKpmRec = nil;
  id<MTLComputePipelineState> pipeKpmDot = nil;

  MetalBackend()
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
        newLibraryWithSource:[NSString stringWithUTF8String:metalShaders]
                     options:nil
                       error:&error];
      if (!library) {
        std::cerr << "Failed to compile Metal shaders: " <<
          [[error localizedDescription] UTF8String] << std::endl;
        return;
      }

      pipeJacobi = [device newComputePipelineStateWithFunction:
                             [library newFunctionWithName:@"jacobi_rotate"]
                                                         error:&error];
      pipeKpmFirst = [device newComputePipelineStateWithFunction:
                               [library newFunctionWithName:@"kpm_step_first"]
                                                           error:&error];
      pipeKpmRec =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"kpm_step_recurrence"]
                                              error:&error];
      pipeKpmDot = [device
        newComputePipelineStateWithFunction:[library
                                              newFunctionWithName:@"kpm_dot"]
                                      error:&error];

      if (pipeJacobi && pipeKpmFirst && pipeKpmRec && pipeKpmDot) {
        available = true;
      }
    }
  }
};

bool
isMetalAvailable()
{
  return MetalBackend::instance().isAvailable();
}

void
computeAllEigenvaluesMetal(const Eigen::SparseMatrix<double>& mat,
                           const bool includeEigenvectors,
                           Eigen::VectorXd& outEvals,
                           Eigen::MatrixXd& outEvecs)
{
  Eigen::Index n = mat.rows();
  if (n == 0) {
    outEvals = Eigen::VectorXd::Zero(0);
    outEvecs = Eigen::MatrixXd::Zero(0, 0);
    return;
  }
  if (n == 1) {
    outEvals = Eigen::VectorXd::Zero(1);
    outEvals[0] = mat.coeff(0, 0);
    if (includeEigenvectors) {
      outEvecs = Eigen::MatrixXd::Identity(1, 1);
    }
    return;
  }

  MetalBackend& backend = MetalBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLDevice> device = backend.getDevice();
    id<MTLCommandQueue> queue = backend.getQueue();
    id<MTLComputePipelineState> pipe = backend.getJacobiPipe();

    uint32_t N = static_cast<uint32_t>(n);

    // Convert sparse matrix to dense float format in shared GPU buffer
    Eigen::MatrixXd denseMat = Eigen::MatrixXd(mat);
    id<MTLBuffer> bufA =
      [device newBufferWithLength:sizeof(float) * N * N
                          options:MTLResourceStorageModeShared];
    float* A = (float*)[bufA contents];
    for (uint32_t i = 0; i < N; ++i) {
      for (uint32_t j = 0; j < N; ++j) {
        A[i * N + j] = static_cast<float>(denseMat(i, j));
      }
    }

    id<MTLBuffer> bufV =
      [device newBufferWithLength:sizeof(float) * N * N
                          options:MTLResourceStorageModeShared];
    float* V = (float*)[bufV contents];
    std::memset(V, 0, sizeof(float) * N * N);
    for (uint32_t i = 0; i < N; ++i) {
      V[i * N + i] = 1.0f;
    }

    uint32_t computeV = includeEigenvectors ? 1 : 0;
    const int maxSweeps = 50;
    const float tol = 1e-6f;

    // Threshold cyclic Jacobi:
    // Compute initial off-diagonal norm
    float sumOffDiag = 0.0f;
    for (uint32_t p = 0; p < N; ++p) {
      for (uint32_t q = p + 1; q < N; ++q) {
        sumOffDiag += std::abs(A[p * N + q]);
      }
    }

    for (int sweep = 0; sweep < maxSweeps; ++sweep) {
      float threshold =
        (sweep < 4) ? (0.2f * sumOffDiag / static_cast<float>(N * (N - 1) / 2))
                    : 0.0f;
      float maxOff = 0.0f;

      for (uint32_t p = 0; p < N; ++p) {
        for (uint32_t q = p + 1; q < N; ++q) {
          float apq = A[p * N + q];
          float absApq = std::abs(apq);
          maxOff = std::max(maxOff, absApq);

          if (absApq <= threshold && absApq <= tol)
            continue;
          if (absApq <= tol)
            continue;

          float app = A[p * N + p];
          float aqq = A[q * N + q];
          float tau = (aqq - app) / (2.0f * apq);
          float t = (tau >= 0.0f ? 1.0f : -1.0f) /
                    (std::abs(tau) + std::sqrt(1.0f + tau * tau));
          float c = 1.0f / std::sqrt(1.0f + t * t);
          float s = t * c;

          // Update diagonals
          A[p * N + p] = app - t * apq;
          A[q * N + q] = aqq + t * apq;
          A[p * N + q] = 0.0f;
          A[q * N + p] = 0.0f;

          // Dispatch GPU kernel to rotate columns/rows in parallel
          id<MTLCommandBuffer> cmd = [queue commandBuffer];
          id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
          [enc setComputePipelineState:pipe];
          [enc setBuffer:bufA offset:0 atIndex:0];
          [enc setBuffer:bufV offset:0 atIndex:1];
          [enc setBytes:&N length:sizeof(uint32_t) atIndex:2];
          [enc setBytes:&p length:sizeof(uint32_t) atIndex:3];
          [enc setBytes:&q length:sizeof(uint32_t) atIndex:4];
          [enc setBytes:&c length:sizeof(float) atIndex:5];
          [enc setBytes:&s length:sizeof(float) atIndex:6];
          [enc setBytes:&computeV length:sizeof(uint32_t) atIndex:7];
          [enc dispatchThreads:MTLSizeMake(N, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
          [enc endEncoding];
          [cmd commit];
          [cmd waitUntilCompleted];
        }
      }

      if (maxOff <= tol)
        break;
    }

    // Extract eigenvalues and index permutation
    std::vector<std::pair<double, size_t>> sortedEvals(N);
    for (uint32_t i = 0; i < N; ++i) {
      sortedEvals[i] = { static_cast<double>(A[i * N + i]),
                         static_cast<size_t>(i) };
    }
    std::sort(sortedEvals.begin(), sortedEvals.end());

    outEvals.resize(n);
    for (uint32_t i = 0; i < N; ++i) {
      outEvals[i] = sortedEvals[i].first;
    }

    if (includeEigenvectors) {
      outEvecs.resize(n, n);
      for (uint32_t col = 0; col < N; ++col) {
        size_t origCol = sortedEvals[col].second;
        for (uint32_t row = 0; row < N; ++row) {
          outEvecs(row, col) = static_cast<double>(V[row * N + origCol]);
        }
      }
    }
  }
}

Eigen::VectorXd
computeChebyshevMomentsKPMMetal(const Eigen::SparseMatrix<double>& mat,
                                size_t numMoments,
                                size_t numRandomVectors,
                                double c,
                                double w)
{
  Eigen::Index n = mat.rows();
  if (n == 0 || numMoments == 0) {
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));
  }

  MetalBackend& backend = MetalBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLDevice> device = backend.getDevice();
    id<MTLCommandQueue> queue = backend.getQueue();
    id<MTLComputePipelineState> pStepFirst = backend.getKpmFirstPipe();
    id<MTLComputePipelineState> pStepRec = backend.getKpmRecPipe();
    id<MTLComputePipelineState> pDot = backend.getKpmDotPipe();

    uint32_t N = static_cast<uint32_t>(n);
    uint32_t nVal = N;
    float c_val = static_cast<float>(c);
    float w_val = static_cast<float>(w);

    // Convert sparse matrix to CSR format in float
    std::vector<uint32_t> row_ptrs(N + 1);
    std::vector<uint32_t> col_indices(mat.nonZeros());
    std::vector<float> values(mat.nonZeros());
    for (uint32_t i = 0; i <= N; ++i)
      row_ptrs[i] = static_cast<uint32_t>(mat.outerIndexPtr()[i]);
    for (int i = 0; i < mat.nonZeros(); ++i) {
      col_indices[i] = static_cast<uint32_t>(mat.innerIndexPtr()[i]);
      values[i] = static_cast<float>(mat.valuePtr()[i]);
    }

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

    id<MTLBuffer> bufU[3];
    for (int k = 0; k < 3; ++k) {
      bufU[k] = [device newBufferWithLength:sizeof(float) * N
                                    options:MTLResourceStorageModeShared];
    }
    id<MTLBuffer> bufV =
      [device newBufferWithLength:sizeof(float) * N
                          options:MTLResourceStorageModeShared];

    uint32_t tg_size = 256;
    uint32_t num_tgs = (N + tg_size - 1) / tg_size;
    id<MTLBuffer> bufPartials =
      [device newBufferWithLength:sizeof(float) * numMoments * num_tgs
                          options:MTLResourceStorageModeShared];

    Eigen::VectorXd avgMoments =
      Eigen::VectorXd::Zero(static_cast<Eigen::Index>(numMoments));

    std::mt19937_64 rng(42);
    std::uniform_int_distribution<int> coin(0, 1);
    double scale = 1.0 / std::sqrt(static_cast<double>(n));

    for (size_t r = 0; r < numRandomVectors; ++r) {
      std::vector<float> v_float(N);
      for (uint32_t i = 0; i < N; ++i) {
        v_float[i] = static_cast<float>(coin(rng) ? scale : -scale);
      }
      std::memcpy([bufV contents], v_float.data(), sizeof(float) * N);
      std::memcpy([bufU[0] contents], v_float.data(), sizeof(float) * N);

      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

      // Moment 0: dot(v, u0)
      [enc setComputePipelineState:pDot];
      [enc setBuffer:bufV offset:0 atIndex:0];
      [enc setBuffer:bufU[0] offset:0 atIndex:1];
      [enc setBuffer:bufPartials offset:0 atIndex:2];
      [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      if (numMoments > 1) {
        // Moment 1: u1 = (A*u0 - c*u0)/w
        [enc setComputePipelineState:pStepFirst];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:bufU[0] offset:0 atIndex:3];
        [enc setBuffer:bufU[1] offset:0 atIndex:4];
        [enc setBytes:&c_val length:sizeof(float) atIndex:5];
        [enc setBytes:&w_val length:sizeof(float) atIndex:6];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:7];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

        [enc setComputePipelineState:pDot];
        [enc setBuffer:bufV offset:0 atIndex:0];
        [enc setBuffer:bufU[1] offset:0 atIndex:1];
        [enc setBuffer:bufPartials
                offset:1 * num_tgs * sizeof(float)
               atIndex:2];
        [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

        int prev_idx = 0;
        int curr_idx = 1;
        int next_idx = 2;

        for (size_t m = 1; m < numMoments - 1; ++m) {
          [enc setComputePipelineState:pStepRec];
          [enc setBuffer:bufRowPtr offset:0 atIndex:0];
          [enc setBuffer:bufColIdx offset:0 atIndex:1];
          [enc setBuffer:bufValues offset:0 atIndex:2];
          [enc setBuffer:bufU[curr_idx] offset:0 atIndex:3];
          [enc setBuffer:bufU[prev_idx] offset:0 atIndex:4];
          [enc setBuffer:bufU[next_idx] offset:0 atIndex:5];
          [enc setBytes:&c_val length:sizeof(float) atIndex:6];
          [enc setBytes:&w_val length:sizeof(float) atIndex:7];
          [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:8];
          [enc dispatchThreads:MTLSizeMake(N, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];

          [enc setComputePipelineState:pDot];
          [enc setBuffer:bufV offset:0 atIndex:0];
          [enc setBuffer:bufU[next_idx] offset:0 atIndex:1];
          [enc setBuffer:bufPartials
                  offset:(m + 1) * num_tgs * sizeof(float)
                 atIndex:2];
          [enc setBytes:&nVal length:sizeof(uint32_t) atIndex:3];
          [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
          [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

          int tmp = prev_idx;
          prev_idx = curr_idx;
          curr_idx = next_idx;
          next_idx = tmp;
        }
      }

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float* parts = (float*)[bufPartials contents];
      for (size_t m = 0; m < numMoments; ++m) {
        float sum = 0.0f;
        for (uint32_t tg = 0; tg < num_tgs; ++tg) {
          sum += parts[m * num_tgs + tg];
        }
        avgMoments[static_cast<Eigen::Index>(m)] += static_cast<double>(sum);
      }
    }

    avgMoments /= static_cast<double>(numRandomVectors);

    // Apply Jackson damping kernel
    for (size_t m = 0; m < numMoments; ++m) {
      double m_d = static_cast<double>(m);
      double M_d = static_cast<double>(numMoments);
      double gm =
        ((M_d - m_d + 1.0) * std::cos(M_PI * m_d / (M_d + 1.0)) +
         std::sin(M_PI * m_d / (M_d + 1.0)) / std::tan(M_PI / (M_d + 1.0))) /
        (M_d + 1.0);
      avgMoments[static_cast<Eigen::Index>(m)] *= gm;
    }

    return avgMoments;
  }
}

} // namespace pylimer_tools::calc
