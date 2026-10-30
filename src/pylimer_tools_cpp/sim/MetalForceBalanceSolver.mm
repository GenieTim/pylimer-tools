#include "MetalForceBalanceSolver.h"
#include "../utils/utilityMacros.h"

#import <Foundation/Foundation.h>
#import <Metal/Metal.h>
#import <simd/simd.h>

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iostream>
#include <vector>

namespace pylimer_tools::sim::mehp {

#pragma pack(push, 4)
struct CGState
{
  float rho;
  float rho_new;
  float pAp;
  float alpha;
  float beta;
  float initial_rho;
  float tolSq;
  uint32_t converged;
  uint32_t iter;
  uint32_t maxIters;
};

struct BBState
{
  float initGNormSq;
  float curGNormSq;
  float alpha;
  float learningRate;
  float tolSq;
  uint32_t shortAlpha;
  uint32_t converged;
  uint32_t iter;
  uint32_t maxIters;
};
#pragma pack(pop)

static const char* metalForceBalanceShaders = R"(
#include <metal_stdlib>
using namespace metal;

struct CGState {
    float rho;
    float rho_new;
    float pAp;
    float alpha;
    float beta;
    float initial_rho;
    float tolSq;
    uint converged;
    uint iter;
    uint maxIters;
};

struct BBState {
    float initGNormSq;
    float curGNormSq;
    float alpha;
    float learningRate;
    float tolSq;
    uint shortAlpha;
    uint converged;
    uint iter;
    uint maxIters;
};

// 1. Scalar SpMV (baseline/fallback)
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

// 2. SIMDgroup-Cooperative SpMV (8 threads per row, coalesced memory access, simd_shuffle_down)
kernel void csr_spmv_coop8(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* x [[buffer(3)]],
    device float* y [[buffer(4)]],
    device const CGState* state [[buffer(5)]],
    constant uint& numRows [[buffer(6)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= numRows) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float sum = 0.0f;
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * x[col_indices[k]];
    }
    sum += simd_shuffle_down(sum, 4);
    sum += simd_shuffle_down(sum, 2);
    sum += simd_shuffle_down(sum, 1);
    if (lane == 0) {
        y[row] = sum;
    }
}

// SIMDgroup-Cooperative SpMV (stateless)
kernel void csr_spmv_coop8_stateless(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* x [[buffer(3)]],
    device float* y [[buffer(4)]],
    constant uint& numRows [[buffer(5)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= numRows) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float sum = 0.0f;
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * x[col_indices[k]];
    }
    sum += simd_shuffle_down(sum, 4);
    sum += simd_shuffle_down(sum, 2);
    sum += simd_shuffle_down(sum, 1);
    if (lane == 0) {
        y[row] = sum;
    }
}

// Dot product with parallel reduction
kernel void dot_product(
    device const float* a [[buffer(0)]],
    device const float* b [[buffer(1)]],
    device float* partialSums [[buffer(2)]],
    device const CGState* state [[buffer(3)]],
    constant uint& n [[buffer(4)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state != nullptr && state->converged != 0) return;
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

// CG Step 1: Compute alpha = rho / pAp on GPU
kernel void cg_compute_alpha(
    device const float* partialSums [[buffer(0)]],
    device CGState* state [[buffer(1)]],
    constant uint& num_tgs [[buffer(2)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]])
{
    if (state->converged != 0) return;
    float sum = 0.0f;
    for (uint i = tid; i < num_tgs; i += 256) {
        sum += partialSums[i];
    }
    sharedData[tid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) sharedData[tid] += sharedData[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) {
        float pAp = sharedData[0];
        state->pAp = pAp;
        if (abs(pAp) <= 1e-20f || !isfinite(pAp)) {
            state->converged = 1;
        } else {
            state->alpha = state->rho / pAp;
        }
    }
}

// CG Step 2: x += alpha * p, r -= alpha * Ap
kernel void cg_update_x_r(
    device float* x [[buffer(0)]],
    device float* r [[buffer(1)]],
    device const float* p [[buffer(2)]],
    device const float* Ap [[buffer(3)]],
    device const CGState* state [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= n) return;
    float alpha = state->alpha;
    x[id] += alpha * p[id];
    r[id] -= alpha * Ap[id];
}

// CG Step 3: Compute beta and check convergence on GPU
kernel void cg_compute_beta(
    device const float* partialSums [[buffer(0)]],
    device CGState* state [[buffer(1)]],
    constant uint& num_tgs [[buffer(2)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]])
{
    if (state->converged != 0) return;
    float sum = 0.0f;
    for (uint i = tid; i < num_tgs; i += 256) {
        sum += partialSums[i];
    }
    sharedData[tid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) sharedData[tid] += sharedData[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) {
        float rho_new = sharedData[0];
        state->rho_new = rho_new;
        state->iter++;
        if (rho_new / state->initial_rho < state->tolSq || state->iter >= state->maxIters) {
            state->converged = 1;
        } else {
            state->beta = rho_new / state->rho;
            state->rho = rho_new;
        }
    }
}

// CG Step 4: p = r + beta * p
kernel void cg_update_p(
    device float* p [[buffer(0)]],
    device const float* r [[buffer(1)]],
    device const CGState* state [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= n) return;
    float beta = state->beta;
    p[id] = r[id] + beta * p[id];
}

// BB Step 1: nextX = x - alpha * gradient, deltaX = -alpha * gradient
kernel void bb_update_x(
    device const float* x [[buffer(0)]],
    device const float* gradient [[buffer(1)]],
    device float* nextX [[buffer(2)]],
    device float* deltaX [[buffer(3)]],
    device const BBState* state [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= n) return;
    float alpha = state->alpha;
    float dx = -alpha * gradient[id];
    deltaX[id] = dx;
    nextX[id] = x[id] + dx;
}

// BB Step 2: SIMDgroup-Cooperative gradient computation: nextG = A * nextX - b
kernel void bb_calc_gradient_coop8(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const float* x [[buffer(3)]],
    device const float* b [[buffer(4)]],
    device float* gradient [[buffer(5)]],
    device const BBState* state [[buffer(6)]],
    constant uint& n [[buffer(7)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= n) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float sum = 0.0f;
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * x[col_indices[k]];
    }
    sum += simd_shuffle_down(sum, 4);
    sum += simd_shuffle_down(sum, 2);
    sum += simd_shuffle_down(sum, 1);
    if (lane == 0) {
        gradient[row] = sum - b[row];
    }
}

// BB Step 3: deltaG = nextG - gradient
kernel void bb_calc_delta_g(
    device const float* nextG [[buffer(0)]],
    device const float* gradient [[buffer(1)]],
    device float* deltaG [[buffer(2)]],
    device const BBState* state [[buffer(3)]],
    constant uint& n [[buffer(4)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= n) return;
    deltaG[id] = nextG[id] - gradient[id];
}

// BB Step 4: Fused quad dot reduction: (sTy, yTy, sTs, curGNormSq)
kernel void bb_quad_dot(
    device const float* deltaX [[buffer(0)]],
    device const float* deltaG [[buffer(1)]],
    device const float* nextG [[buffer(2)]],
    device float4* partialSums [[buffer(3)]],
    device const BBState* state [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    threadgroup float4* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state->converged != 0) return;
    float4 val = float4(0.0f);
    if (gid < n) {
        float dx = deltaX[gid];
        float dg = deltaG[gid];
        float ng = nextG[gid];
        val = float4(dx * dg, dg * dg, dx * dx, ng * ng);
    }
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

// BB Step 5: Compute alpha and check convergence on GPU
kernel void bb_compute_alpha(
    device const float4* partialSums [[buffer(0)]],
    device BBState* state [[buffer(1)]],
    constant uint& num_tgs [[buffer(2)]],
    threadgroup float4* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]])
{
    if (state->converged != 0) return;
    float4 sum = float4(0.0f);
    for (uint i = tid; i < num_tgs; i += 256) {
        sum += partialSums[i];
    }
    sharedData[tid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = 128; s > 0; s >>= 1) {
        if (tid < s) sharedData[tid] += sharedData[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) {
        float4 totals = sharedData[0];
        float sTy = totals.x;
        float yTy = totals.y;
        float sTs = totals.z;
        float curGNormSq = totals.w;
        state->curGNormSq = curGNormSq;
        state->iter++;

        if (curGNormSq / state->initGNormSq <= state->tolSq || state->iter >= state->maxIters) {
            state->converged = 1;
        } else {
            float alphaShort = (yTy > 1e-20f && isfinite(sTy / yTy)) ? (sTy / yTy) : state->learningRate;
            float alphaLong  = (abs(sTy) > 1e-20f && isfinite(sTs / sTy)) ? (sTs / sTy) : state->learningRate;
            float alpha = state->shortAlpha ? alphaShort : alphaLong;
            if (!isfinite(alpha)) alpha = state->learningRate;
            state->alpha = alpha;
        }
    }
}

// BB Solution synchronization
kernel void bb_sync_solution(
    device float* x0 [[buffer(0)]],
    device const float* x1 [[buffer(1)]],
    device const BBState* state [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= n) return;
    if ((state->iter % 2) == 1) {
        x0[gid] = x1[gid];
    }
}

// Center-of-Mass (Zero-Mode) Projection (3D)
kernel void com_reduce_3d(
    device const float* x [[buffer(0)]],
    device float3* partialCom [[buffer(1)]],
    constant uint& numNodes [[buffer(2)]],
    threadgroup float3* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    float3 val = (gid < numNodes) ? float3(x[3 * gid], x[3 * gid + 1], x[3 * gid + 2]) : float3(0.0f);
    sharedData[tid] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    for (uint s = tg_size / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sharedData[tid] += sharedData[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    
    if (tid == 0) {
        partialCom[tg_id] = sharedData[0];
    }
}

kernel void com_finalize_3d(
    device const float3* partialCom [[buffer(0)]],
    device float3* comResult [[buffer(1)]],
    constant uint& numPartialSums [[buffer(2)]],
    constant uint& numNodes [[buffer(3)]],
    threadgroup float3* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]])
{
    float3 sum = float3(0.0f);
    for (uint i = tid; i < numPartialSums; i += tg_size) {
        sum += partialCom[i];
    }
    sharedData[tid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    
    for (uint s = tg_size / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sharedData[tid] += sharedData[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    
    if (tid == 0) {
        comResult[0] = sharedData[0] / float(numNodes);
    }
}

kernel void com_project_3d(
    device float* x [[buffer(0)]],
    device const float3* comResult [[buffer(1)]],
    constant uint& numNodes [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= numNodes) return;
    float3 shift = comResult[0];
    x[3 * gid + 0] -= shift.x;
    x[3 * gid + 1] -= shift.y;
    x[3 * gid + 2] -= shift.z;
}

// Center-of-Mass Projection (1D Fallback)
kernel void com_reduce_1d(
    device const float* x [[buffer(0)]],
    device float* partialCom [[buffer(1)]],
    constant uint& n [[buffer(2)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    float val = (gid < n) ? x[gid] : 0.0f;
    sharedData[tid] = val;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = tg_size / 2; s > 0; s >>= 1) {
        if (tid < s) sharedData[tid] += sharedData[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) partialCom[tg_id] = sharedData[0];
}

kernel void com_finalize_1d(
    device const float* partialCom [[buffer(0)]],
    device float* comResult [[buffer(1)]],
    constant uint& numPartialSums [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint tg_size [[threads_per_threadgroup]])
{
    float sum = 0.0f;
    for (uint i = tid; i < numPartialSums; i += tg_size) sum += partialCom[i];
    sharedData[tid] = sum;
    threadgroup_barrier(mem_flags::mem_threadgroup);
    for (uint s = tg_size / 2; s > 0; s >>= 1) {
        if (tid < s) sharedData[tid] += sharedData[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) comResult[0] = sharedData[0] / float(n);
}

kernel void com_project_1d(
    device float* x [[buffer(0)]],
    device const float* comResult [[buffer(1)]],
    constant uint& n [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= n) return;
    x[gid] -= comResult[0];
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

  // Pipelines
  id<MTLComputePipelineState> getSpmvPipe() const { return pipeSpmv; }
  id<MTLComputePipelineState> getSpmvCoopPipe() const { return pipeSpmvCoop; }
  id<MTLComputePipelineState> getSpmvCoopStatelessPipe() const
  {
    return pipeSpmvCoopStateless;
  }
  id<MTLComputePipelineState> getDotPipe() const { return pipeDot; }

  // CG pipelines
  id<MTLComputePipelineState> getCgAlphaPipe() const { return pipeCgAlpha; }
  id<MTLComputePipelineState> getCgUpdateXRPipe() const
  {
    return pipeCgUpdateXR;
  }
  id<MTLComputePipelineState> getCgBetaPipe() const { return pipeCgBeta; }
  id<MTLComputePipelineState> getCgUpdatePPipe() const { return pipeCgUpdateP; }

  // BB pipelines
  id<MTLComputePipelineState> getBbUpdateXPipe() const { return pipeBbUpdateX; }
  id<MTLComputePipelineState> getBbCalcGradPipe() const
  {
    return pipeBbCalcGrad;
  }
  id<MTLComputePipelineState> getBbCalcDeltaGPipe() const
  {
    return pipeBbCalcDeltaG;
  }
  id<MTLComputePipelineState> getBbQuadDotPipe() const { return pipeBbQuadDot; }
  id<MTLComputePipelineState> getBbAlphaPipe() const { return pipeBbAlpha; }
  id<MTLComputePipelineState> getBbSyncSolPipe() const { return pipeBbSyncSol; }

  // COM projection pipelines
  id<MTLComputePipelineState> getComReduce3DPipe() const
  {
    return pipeComReduce3D;
  }
  id<MTLComputePipelineState> getComFinalize3DPipe() const
  {
    return pipeComFinalize3D;
  }
  id<MTLComputePipelineState> getComProject3DPipe() const
  {
    return pipeComProject3D;
  }
  id<MTLComputePipelineState> getComReduce1DPipe() const
  {
    return pipeComReduce1D;
  }
  id<MTLComputePipelineState> getComFinalize1DPipe() const
  {
    return pipeComFinalize1D;
  }
  id<MTLComputePipelineState> getComProject1DPipe() const
  {
    return pipeComProject1D;
  }

  // Persistent Buffer Pooling
  void ensureMatrixBuffers(size_t numRows, size_t nnz)
  {
    size_t reqRowPtr = numRows + 1;
    if (!bufRowPtr || capRowPtr < reqRowPtr) {
      size_t newCap = std::max(reqRowPtr, static_cast<size_t>(capRowPtr * 1.5));
      bufRowPtr = [device newBufferWithLength:newCap * sizeof(uint32_t)
                                      options:MTLResourceStorageModeShared];
      capRowPtr = newCap;
    }

    if (!bufColIdx || capColIdx < nnz) {
      size_t newCap = std::max(nnz, static_cast<size_t>(capColIdx * 1.5));
      bufColIdx = [device newBufferWithLength:newCap * sizeof(uint32_t)
                                      options:MTLResourceStorageModeShared];
      capColIdx = newCap;
    }

    if (!bufValues || capValues < nnz) {
      size_t newCap = std::max(nnz, static_cast<size_t>(capValues * 1.5));
      bufValues = [device newBufferWithLength:newCap * sizeof(float)
                                      options:MTLResourceStorageModeShared];
      capValues = newCap;
    }
  }

  void ensureCommonVectorBuffers(size_t N)
  {
    if (!bufB || capN < N) {
      size_t newCap = std::max(N, static_cast<size_t>(capN * 1.5));
      bufB = [device newBufferWithLength:newCap * sizeof(float)
                                 options:MTLResourceStorageModeShared];
      bufX = [device newBufferWithLength:newCap * sizeof(float)
                                 options:MTLResourceStorageModeShared];
      capN = newCap;
    }
  }

  void ensureCgBuffers(size_t N, size_t numTgs)
  {
    ensureCommonVectorBuffers(N);
    if (!bufR || capCgVectors < N) {
      size_t newCap = std::max(N, static_cast<size_t>(capCgVectors * 1.5));
      bufR = [device newBufferWithLength:newCap * sizeof(float)
                                 options:MTLResourceStorageModeShared];
      bufP = [device newBufferWithLength:newCap * sizeof(float)
                                 options:MTLResourceStorageModeShared];
      bufAp = [device newBufferWithLength:newCap * sizeof(float)
                                  options:MTLResourceStorageModeShared];
      capCgVectors = newCap;
    }

    if (!bufPartial || capPartial < numTgs) {
      size_t newCap = std::max(numTgs, static_cast<size_t>(capPartial * 1.5));
      bufPartial = [device newBufferWithLength:newCap * sizeof(float)
                                       options:MTLResourceStorageModeShared];
      capPartial = newCap;
    }

    if (!bufCGState) {
      bufCGState = [device newBufferWithLength:sizeof(CGState)
                                       options:MTLResourceStorageModeShared];
    }
  }

  void ensureBbBuffers(size_t N, size_t numTgs)
  {
    ensureCommonVectorBuffers(N);
    if (!bufNextX || capBbVectors < N) {
      size_t newCap = std::max(N, static_cast<size_t>(capBbVectors * 1.5));
      bufNextX = [device newBufferWithLength:newCap * sizeof(float)
                                     options:MTLResourceStorageModeShared];
      bufG = [device newBufferWithLength:newCap * sizeof(float)
                                 options:MTLResourceStorageModeShared];
      bufNextG = [device newBufferWithLength:newCap * sizeof(float)
                                     options:MTLResourceStorageModeShared];
      bufDeltaX = [device newBufferWithLength:newCap * sizeof(float)
                                      options:MTLResourceStorageModeShared];
      bufDeltaG = [device newBufferWithLength:newCap * sizeof(float)
                                      options:MTLResourceStorageModeShared];
      capBbVectors = newCap;
    }

    if (!bufQuadPartial || capQuadPartial < numTgs) {
      size_t newCap =
        std::max(numTgs, static_cast<size_t>(capQuadPartial * 1.5));
      bufQuadPartial =
        [device newBufferWithLength:newCap * sizeof(simd_float4)
                            options:MTLResourceStorageModeShared];
      capQuadPartial = newCap;
    }

    if (!bufPartial || capPartial < numTgs) {
      size_t newCap = std::max(numTgs, static_cast<size_t>(capPartial * 1.5));
      bufPartial = [device newBufferWithLength:newCap * sizeof(float)
                                       options:MTLResourceStorageModeShared];
      capPartial = newCap;
    }

    if (!bufBBState) {
      bufBBState = [device newBufferWithLength:sizeof(BBState)
                                       options:MTLResourceStorageModeShared];
    }
  }

  void ensureComBuffers(size_t numItems, size_t numTgs, bool is3D)
  {
    (void)numItems;
    (void)is3D;
    if (!bufComPartial || capComPartial < numTgs) {
      size_t newCap =
        std::max(numTgs, static_cast<size_t>(capComPartial * 1.5));
      bufComPartial = [device newBufferWithLength:newCap * sizeof(simd_float3)
                                          options:MTLResourceStorageModeShared];
      capComPartial = newCap;
    }
    if (!bufComResult) {
      bufComResult = [device newBufferWithLength:sizeof(simd_float3)
                                         options:MTLResourceStorageModeShared];
    }
  }

  // Persistent Buffer Getters
  id<MTLBuffer> getBufRowPtr() const { return bufRowPtr; }
  id<MTLBuffer> getBufColIdx() const { return bufColIdx; }
  id<MTLBuffer> getBufValues() const { return bufValues; }
  id<MTLBuffer> getBufB() const { return bufB; }
  id<MTLBuffer> getBufX() const { return bufX; }

  // CG buffers
  id<MTLBuffer> getBufR() const { return bufR; }
  id<MTLBuffer> getBufP() const { return bufP; }
  id<MTLBuffer> getBufAp() const { return bufAp; }
  id<MTLBuffer> getBufCGState() const { return bufCGState; }

  // BB buffers
  id<MTLBuffer> getBufNextX() const { return bufNextX; }
  id<MTLBuffer> getBufG() const { return bufG; }
  id<MTLBuffer> getBufNextG() const { return bufNextG; }
  id<MTLBuffer> getBufDeltaX() const { return bufDeltaX; }
  id<MTLBuffer> getBufDeltaG() const { return bufDeltaG; }
  id<MTLBuffer> getBufBBState() const { return bufBBState; }

  // Reductions
  id<MTLBuffer> getBufPartial() const { return bufPartial; }
  id<MTLBuffer> getBufQuadPartial() const { return bufQuadPartial; }
  id<MTLBuffer> getBufComPartial() const { return bufComPartial; }
  id<MTLBuffer> getBufComResult() const { return bufComResult; }

  // Encode Center-of-Mass projection onto encoder
  void encodeComProjection(id<MTLComputeCommandEncoder> enc,
                           id<MTLBuffer> bufTarget,
                           uint32_t N)
  {
    uint32_t tg_size = 256;
    if (N % 3 == 0) {
      uint32_t numNodes = N / 3;
      uint32_t num_tgs = (numNodes + tg_size - 1) / tg_size;
      ensureComBuffers(numNodes, num_tgs, true);

      [enc setComputePipelineState:pipeComReduce3D];
      [enc setBuffer:bufTarget offset:0 atIndex:0];
      [enc setBuffer:bufComPartial offset:0 atIndex:1];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:2];
      [enc setThreadgroupMemoryLength:sizeof(simd_float3) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pipeComFinalize3D];
      [enc setBuffer:bufComPartial offset:0 atIndex:0];
      [enc setBuffer:bufComResult offset:0 atIndex:1];
      [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:3];
      [enc setThreadgroupMemoryLength:sizeof(simd_float3) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pipeComProject3D];
      [enc setBuffer:bufTarget offset:0 atIndex:0];
      [enc setBuffer:bufComResult offset:0 atIndex:1];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:2];
      [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    } else {
      uint32_t num_tgs = (N + tg_size - 1) / tg_size;
      ensureComBuffers(N, num_tgs, false);

      [enc setComputePipelineState:pipeComReduce1D];
      [enc setBuffer:bufTarget offset:0 atIndex:0];
      [enc setBuffer:bufComPartial offset:0 atIndex:1];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:2];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pipeComFinalize1D];
      [enc setBuffer:bufComPartial offset:0 atIndex:0];
      [enc setBuffer:bufComResult offset:0 atIndex:1];
      [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:3];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      [enc setComputePipelineState:pipeComProject1D];
      [enc setBuffer:bufTarget offset:0 atIndex:0];
      [enc setBuffer:bufComResult offset:0 atIndex:1];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:2];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
  }

private:
  bool available = false;
  id<MTLDevice> device = nil;
  id<MTLCommandQueue> queue = nil;

  // Pipeline states
  id<MTLComputePipelineState> pipeSpmv = nil;
  id<MTLComputePipelineState> pipeSpmvCoop = nil;
  id<MTLComputePipelineState> pipeSpmvCoopStateless = nil;
  id<MTLComputePipelineState> pipeDot = nil;

  id<MTLComputePipelineState> pipeCgAlpha = nil;
  id<MTLComputePipelineState> pipeCgUpdateXR = nil;
  id<MTLComputePipelineState> pipeCgBeta = nil;
  id<MTLComputePipelineState> pipeCgUpdateP = nil;

  id<MTLComputePipelineState> pipeBbUpdateX = nil;
  id<MTLComputePipelineState> pipeBbCalcGrad = nil;
  id<MTLComputePipelineState> pipeBbCalcDeltaG = nil;
  id<MTLComputePipelineState> pipeBbQuadDot = nil;
  id<MTLComputePipelineState> pipeBbAlpha = nil;
  id<MTLComputePipelineState> pipeBbSyncSol = nil;

  id<MTLComputePipelineState> pipeComReduce3D = nil;
  id<MTLComputePipelineState> pipeComFinalize3D = nil;
  id<MTLComputePipelineState> pipeComProject3D = nil;
  id<MTLComputePipelineState> pipeComReduce1D = nil;
  id<MTLComputePipelineState> pipeComFinalize1D = nil;
  id<MTLComputePipelineState> pipeComProject1D = nil;

  // Persistent MTLBuffers
  id<MTLBuffer> bufRowPtr = nil;
  id<MTLBuffer> bufColIdx = nil;
  id<MTLBuffer> bufValues = nil;
  id<MTLBuffer> bufB = nil;
  id<MTLBuffer> bufX = nil;

  id<MTLBuffer> bufR = nil;
  id<MTLBuffer> bufP = nil;
  id<MTLBuffer> bufAp = nil;
  id<MTLBuffer> bufCGState = nil;

  id<MTLBuffer> bufNextX = nil;
  id<MTLBuffer> bufG = nil;
  id<MTLBuffer> bufNextG = nil;
  id<MTLBuffer> bufDeltaX = nil;
  id<MTLBuffer> bufDeltaG = nil;
  id<MTLBuffer> bufBBState = nil;

  id<MTLBuffer> bufPartial = nil;
  id<MTLBuffer> bufQuadPartial = nil;
  id<MTLBuffer> bufComPartial = nil;
  id<MTLBuffer> bufComResult = nil;

  // Tracked capacities
  size_t capRowPtr = 0;
  size_t capColIdx = 0;
  size_t capValues = 0;
  size_t capN = 0;
  size_t capCgVectors = 0;
  size_t capBbVectors = 0;
  size_t capPartial = 0;
  size_t capQuadPartial = 0;
  size_t capComPartial = 0;

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
      pipeSpmvCoop = [device newComputePipelineStateWithFunction:
                               [library newFunctionWithName:@"csr_spmv_coop8"]
                                                           error:&error];
      pipeSpmvCoopStateless =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"csr_spmv_coop8_stateless"]
                                              error:&error];
      pipeDot = [device newComputePipelineStateWithFunction:
                          [library newFunctionWithName:@"dot_product"]
                                                      error:&error];

      pipeCgAlpha = [device newComputePipelineStateWithFunction:
                              [library newFunctionWithName:@"cg_compute_alpha"]
                                                          error:&error];
      pipeCgUpdateXR = [device newComputePipelineStateWithFunction:
                                 [library newFunctionWithName:@"cg_update_x_r"]
                                                             error:&error];
      pipeCgBeta = [device newComputePipelineStateWithFunction:
                             [library newFunctionWithName:@"cg_compute_beta"]
                                                         error:&error];
      pipeCgUpdateP = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"cg_update_p"]
                                                            error:&error];

      pipeBbUpdateX = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"bb_update_x"]
                                                            error:&error];
      pipeBbCalcGrad =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_gradient_coop8"]
                                              error:&error];
      pipeBbCalcDeltaG =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_delta_g"]
                                              error:&error];
      pipeBbQuadDot = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"bb_quad_dot"]
                                                            error:&error];
      pipeBbAlpha = [device newComputePipelineStateWithFunction:
                              [library newFunctionWithName:@"bb_compute_alpha"]
                                                          error:&error];
      pipeBbSyncSol =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_sync_solution"]
                                              error:&error];

      pipeComReduce3D = [device newComputePipelineStateWithFunction:
                                  [library newFunctionWithName:@"com_reduce_3d"]
                                                              error:&error];
      pipeComFinalize3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"com_finalize_3d"]
                                              error:&error];
      pipeComProject3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"com_project_3d"]
                                              error:&error];

      pipeComReduce1D = [device newComputePipelineStateWithFunction:
                                  [library newFunctionWithName:@"com_reduce_1d"]
                                                              error:&error];
      pipeComFinalize1D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"com_finalize_1d"]
                                              error:&error];
      pipeComProject1D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"com_project_1d"]
                                              error:&error];

      if (pipeSpmv && pipeSpmvCoop && pipeSpmvCoopStateless && pipeDot &&
          pipeCgAlpha && pipeCgUpdateXR && pipeCgBeta && pipeCgUpdateP &&
          pipeBbUpdateX && pipeBbCalcGrad && pipeBbCalcDeltaG &&
          pipeBbQuadDot && pipeBbAlpha && pipeBbSyncSol && pipeComReduce3D &&
          pipeComFinalize3D && pipeComProject3D && pipeComReduce1D &&
          pipeComFinalize1D && pipeComProject1D) {
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
                             const Eigen::VectorXd& initialX,
                             bool projectZeroModes)
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
    id<MTLCommandQueue> queue = backend.getQueue();

    uint32_t N = static_cast<uint32_t>(n);
    uint32_t tg_size = 256;
    uint32_t num_tgs = (N + tg_size - 1) / tg_size;
    size_t nnz = static_cast<size_t>(A.nonZeros());

    // 1. Persistent Buffer Pooling: Matrix buffers
    backend.ensureMatrixBuffers(N, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();

    uint32_t* rowPtrsMem = static_cast<uint32_t*>([bufRowPtr contents]);
    for (uint32_t i = 0; i <= N; ++i) {
      rowPtrsMem[i] = static_cast<uint32_t>(A.outerIndexPtr()[i]);
    }
    uint32_t* colIdxMem = static_cast<uint32_t*>([bufColIdx contents]);
    float* valuesMem = static_cast<float*>([bufValues contents]);
    for (size_t i = 0; i < nnz; ++i) {
      colIdxMem[i] = static_cast<uint32_t>(A.innerIndexPtr()[i]);
      valuesMem[i] = static_cast<float>(A.valuePtr()[i]);
    }

    // 2. Persistent Buffer Pooling: BB vectors & state
    backend.ensureBbBuffers(N, num_tgs);
    id<MTLBuffer> bufB = backend.getBufB();
    id<MTLBuffer> bufX0 = backend.getBufX();
    id<MTLBuffer> bufX1 = backend.getBufNextX();
    id<MTLBuffer> bufG0 = backend.getBufG();
    id<MTLBuffer> bufG1 = backend.getBufNextG();
    id<MTLBuffer> bufDeltaX = backend.getBufDeltaX();
    id<MTLBuffer> bufDeltaG = backend.getBufDeltaG();
    id<MTLBuffer> bufQuadPartial = backend.getBufQuadPartial();
    id<MTLBuffer> bufPartial = backend.getBufPartial();
    id<MTLBuffer> bufState = backend.getBufBBState();

    float* bMem = static_cast<float*>([bufB contents]);
    for (uint32_t i = 0; i < N; ++i) {
      bMem[i] = static_cast<float>(b[i]);
    }

    float* x0Mem = static_cast<float*>([bufX0 contents]);
    if (initialX.size() == n) {
      for (uint32_t i = 0; i < N; ++i)
        x0Mem[i] = static_cast<float>(initialX[i]);
    } else {
      std::memset(x0Mem, 0, sizeof(float) * N);
    }

    // Compute initial gradient: g0 = A*x0 - b
    uint32_t totalSpmvThreads = N * 8;
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getBbCalcGradPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX0 offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufG0 offset:0 atIndex:5];
      [enc setBuffer:nil offset:0 atIndex:6]; // stateless
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:7];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      // Initial gNorm
      [enc setComputePipelineState:backend.getDotPipe()];
      [enc setBuffer:bufG0 offset:0 atIndex:0];
      [enc setBuffer:bufG0 offset:0 atIndex:1];
      [enc setBuffer:bufPartial offset:0 atIndex:2];
      [enc setBuffer:nil offset:0 atIndex:3];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    float* parts = static_cast<float*>([bufPartial contents]);
    float initGNormSq = 0.0f;
    for (uint32_t i = 0; i < num_tgs; ++i)
      initGNormSq += parts[i];

    float initGNorm = std::sqrt(initGNormSq);
    if (initGNorm < 1e-12f) {
      iterationsDone = 0;
      return initialX.size() == n ? initialX : Eigen::VectorXd::Zero(n);
    }

    BBState* statePtr = static_cast<BBState*>([bufState contents]);
    statePtr->initGNormSq = initGNormSq;
    statePtr->curGNormSq = initGNormSq;
    statePtr->alpha = static_cast<float>(learningRate);
    statePtr->learningRate = static_cast<float>(learningRate);
    statePtr->tolSq = static_cast<float>(tolerance * tolerance);
    statePtr->shortAlpha = shortAlpha ? 1 : 0;
    statePtr->converged = 0;
    statePtr->iter = 0;
    statePtr->maxIters = static_cast<uint32_t>(maxIterations);

    const uint32_t batchSize = 20;

    // 3. GPU-Resident Batched Iteration Loop
    while (statePtr->converged == 0 && statePtr->iter < statePtr->maxIters) {
      uint32_t remaining = statePtr->maxIters - statePtr->iter;
      uint32_t itersInBatch = std::min(batchSize, remaining);

      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

      for (uint32_t step = 0; step < itersInBatch; ++step) {
        uint32_t iterIdx = statePtr->iter + step;
        id<MTLBuffer> curX = (iterIdx % 2 == 0) ? bufX0 : bufX1;
        id<MTLBuffer> curG = (iterIdx % 2 == 0) ? bufG0 : bufG1;
        id<MTLBuffer> nextX = (iterIdx % 2 == 0) ? bufX1 : bufX0;
        id<MTLBuffer> nextG = (iterIdx % 2 == 0) ? bufG1 : bufG0;

        // Step 1: nextX = x - alpha * g, deltaX = -alpha * g
        [enc setComputePipelineState:backend.getBbUpdateXPipe()];
        [enc setBuffer:curX offset:0 atIndex:0];
        [enc setBuffer:curG offset:0 atIndex:1];
        [enc setBuffer:nextX offset:0 atIndex:2];
        [enc setBuffer:bufDeltaX offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 2: SIMDgroup-Cooperative nextG = A * nextX - b
        [enc setComputePipelineState:backend.getBbCalcGradPipe()];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:nextX offset:0 atIndex:3];
        [enc setBuffer:bufB offset:0 atIndex:4];
        [enc setBuffer:nextG offset:0 atIndex:5];
        [enc setBuffer:bufState offset:0 atIndex:6];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:7];
        [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 3: deltaG = nextG - curG
        [enc setComputePipelineState:backend.getBbCalcDeltaGPipe()];
        [enc setBuffer:nextG offset:0 atIndex:0];
        [enc setBuffer:curG offset:0 atIndex:1];
        [enc setBuffer:bufDeltaG offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 4: Fused quad dot reduction (sTy, yTy, sTs, curGNormSq)
        [enc setComputePipelineState:backend.getBbQuadDotPipe()];
        [enc setBuffer:bufDeltaX offset:0 atIndex:0];
        [enc setBuffer:bufDeltaG offset:0 atIndex:1];
        [enc setBuffer:nextG offset:0 atIndex:2];
        [enc setBuffer:bufQuadPartial offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
        [enc setThreadgroupMemoryLength:sizeof(simd_float4) * tg_size
                                atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 5: Compute alpha and check convergence on GPU
        [enc setComputePipelineState:backend.getBbAlphaPipe()];
        [enc setBuffer:bufQuadPartial offset:0 atIndex:0];
        [enc setBuffer:bufState offset:0 atIndex:1];
        [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
        [enc setThreadgroupMemoryLength:sizeof(simd_float4) * tg_size
                                atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Center-of-Mass Projection every 25 iterations (if enabled)
        if (projectZeroModes && ((iterIdx + 1) % 25 == 0)) {
          backend.encodeComProjection(enc, nextX, N);
        }
      }

      // Sync final solution to bufX0 if ended on odd step
      [enc setComputePipelineState:backend.getBbSyncSolPipe()];
      [enc setBuffer:bufX0 offset:0 atIndex:0];
      [enc setBuffer:bufX1 offset:0 atIndex:1];
      [enc setBuffer:bufState offset:0 atIndex:2];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:3];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      // Periodic gradient recomputation every batch to avoid single-precision
      // drift
      [enc setComputePipelineState:backend.getBbCalcGradPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX0 offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufG0 offset:0 atIndex:5];
      [enc setBuffer:bufState offset:0 atIndex:6];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:7];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    // Final COM projection if requested
    if (projectZeroModes) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      backend.encodeComProjection(enc, bufX0, N);
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    iterationsDone = static_cast<int>(statePtr->iter);
    float* finalXPtr = static_cast<float*>([bufX0 contents]);
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
                const Eigen::VectorXd& initialX,
                bool projectZeroModes)
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
    id<MTLCommandQueue> queue = backend.getQueue();

    uint32_t N = static_cast<uint32_t>(n);
    uint32_t tg_size = 256;
    uint32_t num_tgs = (N + tg_size - 1) / tg_size;
    size_t nnz = static_cast<size_t>(A.nonZeros());

    // 1. Persistent Buffer Pooling: Matrix buffers
    backend.ensureMatrixBuffers(N, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();

    uint32_t* rowPtrsMem = static_cast<uint32_t*>([bufRowPtr contents]);
    for (uint32_t i = 0; i <= N; ++i) {
      rowPtrsMem[i] = static_cast<uint32_t>(A.outerIndexPtr()[i]);
    }
    uint32_t* colIdxMem = static_cast<uint32_t*>([bufColIdx contents]);
    float* valuesMem = static_cast<float*>([bufValues contents]);
    for (size_t i = 0; i < nnz; ++i) {
      colIdxMem[i] = static_cast<uint32_t>(A.innerIndexPtr()[i]);
      valuesMem[i] = static_cast<float>(A.valuePtr()[i]);
    }

    // 2. Persistent Buffer Pooling: CG vectors & state
    backend.ensureCgBuffers(N, num_tgs);
    id<MTLBuffer> bufB = backend.getBufB();
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufR = backend.getBufR();
    id<MTLBuffer> bufP = backend.getBufP();
    id<MTLBuffer> bufAp = backend.getBufAp();
    id<MTLBuffer> bufPartial = backend.getBufPartial();
    id<MTLBuffer> bufState = backend.getBufCGState();

    float* bMem = static_cast<float*>([bufB contents]);
    for (uint32_t i = 0; i < N; ++i) {
      bMem[i] = static_cast<float>(b[i]);
    }

    float* xPtr = static_cast<float*>([bufX contents]);
    float* rPtr = static_cast<float*>([bufR contents]);
    float* pPtr = static_cast<float*>([bufP contents]);

    uint32_t totalSpmvThreads = N * 8;

    if (initialX.size() == n) {
      for (uint32_t i = 0; i < N; ++i)
        xPtr[i] = static_cast<float>(initialX[i]);

      // r = b - A*x
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getSpmvCoopStatelessPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufAp offset:0 atIndex:4];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];

      float* apPtr = static_cast<float*>([bufAp contents]);
      for (uint32_t i = 0; i < N; ++i) {
        float val = bMem[i] - apPtr[i];
        rPtr[i] = val;
        pPtr[i] = val;
      }
    } else {
      std::memset(xPtr, 0, sizeof(float) * N);
      std::memcpy(rPtr, bMem, sizeof(float) * N);
      std::memcpy(pPtr, bMem, sizeof(float) * N);
    }

    // Compute initial rho = dot(r, r)
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getDotPipe()];
      [enc setBuffer:bufR offset:0 atIndex:0];
      [enc setBuffer:bufR offset:0 atIndex:1];
      [enc setBuffer:bufPartial offset:0 atIndex:2];
      [enc setBuffer:nil offset:0 atIndex:3];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
      [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
      [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    float* parts = static_cast<float*>([bufPartial contents]);
    float initial_rho = 0.0f;
    for (uint32_t i = 0; i < num_tgs; ++i)
      initial_rho += parts[i];

    if (initial_rho < 1e-24f) {
      iterationsDone = 0;
      return initialX.size() == n ? initialX : Eigen::VectorXd::Zero(n);
    }

    CGState* statePtr = static_cast<CGState*>([bufState contents]);
    statePtr->iter = 0;
    statePtr->converged = 0;
    statePtr->maxIters = static_cast<uint32_t>(maxIterations);
    statePtr->tolSq = static_cast<float>(tolerance * tolerance);
    statePtr->rho = initial_rho;
    statePtr->initial_rho = initial_rho;

    const uint32_t batchSize = 25;

    // 3. GPU-Resident Batched Iteration Loop
    while (statePtr->converged == 0 && statePtr->iter < statePtr->maxIters) {
      uint32_t remaining = statePtr->maxIters - statePtr->iter;
      uint32_t itersInBatch = std::min(batchSize, remaining);

      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

      for (uint32_t step = 0; step < itersInBatch; ++step) {
        // Step 1: SIMDgroup-Cooperative SpMV Ap = A * p
        [enc setComputePipelineState:backend.getSpmvCoopPipe()];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:bufP offset:0 atIndex:3];
        [enc setBuffer:bufAp offset:0 atIndex:4];
        [enc setBuffer:bufState offset:0 atIndex:5];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:6];
        [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 2: Dot pAp = dot(p, Ap)
        [enc setComputePipelineState:backend.getDotPipe()];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:bufAp offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 3: Compute alpha on GPU
        [enc setComputePipelineState:backend.getCgAlphaPipe()];
        [enc setBuffer:bufPartial offset:0 atIndex:0];
        [enc setBuffer:bufState offset:0 atIndex:1];
        [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 4: Update x += alpha * p, r -= alpha * Ap
        [enc setComputePipelineState:backend.getCgUpdateXRPipe()];
        [enc setBuffer:bufX offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufP offset:0 atIndex:2];
        [enc setBuffer:bufAp offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 5: Dot rho_new = dot(r, r)
        [enc setComputePipelineState:backend.getDotPipe()];
        [enc setBuffer:bufR offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 6: Compute beta and check convergence on GPU
        [enc setComputePipelineState:backend.getCgBetaPipe()];
        [enc setBuffer:bufPartial offset:0 atIndex:0];
        [enc setBuffer:bufState offset:0 atIndex:1];
        [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 7: Update p = r + beta * p
        [enc setComputePipelineState:backend.getCgUpdatePPipe()];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufState offset:0 atIndex:2];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:3];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Center-of-Mass Projection every 25 iterations (if enabled)
        if (projectZeroModes && ((statePtr->iter + step + 1) % 25 == 0)) {
          backend.encodeComProjection(enc, bufX, N);
        }
      }

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    // Final COM projection if requested
    if (projectZeroModes) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      backend.encodeComProjection(enc, bufX, N);
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    iterationsDone = static_cast<int>(statePtr->iter);
    float* finalXPtr = static_cast<float*>([bufX contents]);
    Eigen::VectorXd res(n);
    for (uint32_t i = 0; i < N; ++i) {
      res[i] = static_cast<double>(finalXPtr[i]);
    }
    return res;
  }
}

} // namespace pylimer_tools::sim::mehp
