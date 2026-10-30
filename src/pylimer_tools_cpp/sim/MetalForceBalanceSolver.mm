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

// ============================================================================
// 1D SpMV Kernels
// ============================================================================

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

// ============================================================================
// Multi-RHS 3D Block-SpMV (3x Memory Bandwidth Reduction)
// ============================================================================

inline float3 simd_shuffle_down_float3(float3 v, ushort delta) {
    return float3(
        simd_shuffle_down(v.x, delta),
        simd_shuffle_down(v.y, delta),
        simd_shuffle_down(v.z, delta)
    );
}

kernel void csr_spmm_3rhs(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const packed_float3* x [[buffer(3)]],
    device packed_float3* y [[buffer(4)]],
    device const CGState* state [[buffer(5)]],
    constant uint& numNodes [[buffer(6)]],
    uint row [[thread_position_in_grid]])
{
    if (state != nullptr && state->converged != 0) return;
    if (row >= numNodes) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float3 sum = float3(0.0f);
    for (uint k = start; k < end; ++k) {
        sum += values[k] * float3(x[col_indices[k]]);
    }
    y[row] = packed_float3(sum);
}

kernel void csr_spmm_3rhs_coop8(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const packed_float3* x [[buffer(3)]],
    device packed_float3* y [[buffer(4)]],
    device const CGState* state [[buffer(5)]],
    constant uint& numNodes [[buffer(6)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= numNodes) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float3 sum = float3(0.0f);
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * float3(x[col_indices[k]]);
    }
    sum += simd_shuffle_down_float3(sum, 4);
    sum += simd_shuffle_down_float3(sum, 2);
    sum += simd_shuffle_down_float3(sum, 1);
    if (lane == 0) {
        y[row] = packed_float3(sum);
    }
}

kernel void csr_spmm_3rhs_stateless(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const packed_float3* x [[buffer(3)]],
    device packed_float3* y [[buffer(4)]],
    constant uint& numNodes [[buffer(5)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= numNodes) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float3 sum = float3(0.0f);
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * float3(x[col_indices[k]]);
    }
    sum += simd_shuffle_down_float3(sum, 4);
    sum += simd_shuffle_down_float3(sum, 2);
    sum += simd_shuffle_down_float3(sum, 1);
    if (lane == 0) {
        y[row] = packed_float3(sum);
    }
}

// ============================================================================
// Jacobi (Diagonal Scaling) Preconditioning Kernels
// ============================================================================

kernel void extract_inv_diagonal(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device float* invDiag [[buffer(3)]],
    constant uint& numRows [[buffer(4)]],
    uint row [[thread_position_in_grid]])
{
    if (row >= numRows) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float d = 0.0f;
    for (uint k = start; k < end; ++k) {
        if (col_indices[k] == row) {
            d = values[k];
            break;
        }
    }
    invDiag[row] = (abs(d) > 1e-12f) ? (1.0f / d) : 1.0f;
}

// 1D Jacobi Application and Dot Product Reduction: z = invDiag * r, sum(r * z)
kernel void pcg_apply_jacobi_and_dot_1d(
    device const float* invDiag [[buffer(0)]],
    device const float* r [[buffer(1)]],
    device float* z [[buffer(2)]],
    device float* partialSums [[buffer(3)]],
    device const CGState* state [[buffer(4)]],
    constant uint& n [[buffer(5)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    float sumVal = 0.0f;
    if (gid < n) {
        float invD = invDiag[gid];
        float r_val = r[gid];
        float z_val = invD * r_val;
        z[gid] = z_val;
        sumVal = r_val * z_val;
    }
    sharedData[tid] = sumVal;
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

// 3D Jacobi Application and Dot Product Reduction: z = invDiag * r, sum(dot(r, z))
kernel void pcg_apply_jacobi_and_dot_3d(
    device const float* invDiag [[buffer(0)]],
    device const packed_float3* r [[buffer(1)]],
    device packed_float3* z [[buffer(2)]],
    device float* partialSums [[buffer(3)]],
    device const CGState* state [[buffer(4)]],
    constant uint& numNodes [[buffer(5)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    float sumVal = 0.0f;
    if (gid < numNodes) {
        float invD = invDiag[gid];
        float3 r_val = float3(r[gid]);
        float3 z_val = invD * r_val;
        z[gid] = packed_float3(z_val);
        sumVal = dot(r_val, z_val);
    }
    sharedData[tid] = sumVal;
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

// ============================================================================
// Vector Reductions and Updates
// ============================================================================

// 1D Dot product with parallel reduction
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

// 3D Dot product with parallel reduction
kernel void dot_product_3d(
    device const packed_float3* a [[buffer(0)]],
    device const packed_float3* b [[buffer(1)]],
    device float* partialSums [[buffer(2)]],
    device const CGState* state [[buffer(3)]],
    constant uint& numNodes [[buffer(4)]],
    threadgroup float* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    float val = (gid < numNodes) ? dot(float3(a[gid]), float3(b[gid])) : 0.0f;
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

// CG Step 2: 1D x += alpha * p, r -= alpha * Ap
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

// CG Step 2: 3D x += alpha * p, r -= alpha * Ap
kernel void cg_update_x_r_3d(
    device packed_float3* x [[buffer(0)]],
    device packed_float3* r [[buffer(1)]],
    device const packed_float3* p [[buffer(2)]],
    device const packed_float3* Ap [[buffer(3)]],
    device const CGState* state [[buffer(4)]],
    constant uint& numNodes [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= numNodes) return;
    float alpha = state->alpha;
    x[id] = packed_float3(float3(x[id]) + alpha * float3(p[id]));
    r[id] = packed_float3(float3(r[id]) - alpha * float3(Ap[id]));
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
        if (abs(rho_new) / abs(state->initial_rho) < state->tolSq || state->iter >= state->maxIters) {
            state->converged = 1;
        } else {
            state->beta = rho_new / state->rho;
            state->rho = rho_new;
        }
    }
}

// CG Step 4: 1D p = r + beta * p
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

// CG Step 4: 1D PCG p = z + beta * p
kernel void pcg_update_p_1d(
    device float* p [[buffer(0)]],
    device const float* z [[buffer(1)]],
    device const CGState* state [[buffer(2)]],
    constant uint& n [[buffer(3)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= n) return;
    float beta = state->beta;
    p[id] = z[id] + beta * p[id];
}

// CG Step 4: 3D PCG p = z + beta * p
kernel void pcg_update_p_3d(
    device packed_float3* p [[buffer(0)]],
    device const packed_float3* z [[buffer(1)]],
    device const CGState* state [[buffer(2)]],
    constant uint& numNodes [[buffer(3)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= numNodes) return;
    float beta = state->beta;
    p[id] = packed_float3(float3(z[id]) + beta * float3(p[id]));
}

// ============================================================================
// Barzilai-Borwein Kernels
// ============================================================================

// BB Step 1: 1D nextX = x - alpha * gradient, deltaX = -alpha * gradient
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

// BB Step 1: 3D nextX = x - alpha * gradient, deltaX = -alpha * gradient
kernel void bb_update_x_3d(
    device const packed_float3* x [[buffer(0)]],
    device const packed_float3* gradient [[buffer(1)]],
    device packed_float3* nextX [[buffer(2)]],
    device packed_float3* deltaX [[buffer(3)]],
    device const BBState* state [[buffer(4)]],
    constant uint& numNodes [[buffer(5)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= numNodes) return;
    float alpha = state->alpha;
    float3 dx = -alpha * float3(gradient[id]);
    deltaX[id] = packed_float3(dx);
    nextX[id] = packed_float3(float3(x[id]) + dx);
}

// BB Step 2: 1D gradient computation: nextG = A * nextX - b
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

// BB Step 2: 3D gradient computation: nextG = A * nextX - b
kernel void bb_calc_gradient_3rhs_coop8(
    device const uint* row_ptrs [[buffer(0)]],
    device const uint* col_indices [[buffer(1)]],
    device const float* values [[buffer(2)]],
    device const packed_float3* x [[buffer(3)]],
    device const packed_float3* b [[buffer(4)]],
    device packed_float3* gradient [[buffer(5)]],
    device const BBState* state [[buffer(6)]],
    constant uint& numNodes [[buffer(7)]],
    uint gid [[thread_position_in_grid]],
    uint simd_lane [[thread_index_in_simdgroup]])
{
    if (state != nullptr && state->converged != 0) return;
    uint row = gid / 8;
    uint lane = simd_lane % 8;
    if (row >= numNodes) return;
    uint start = row_ptrs[row];
    uint end = row_ptrs[row + 1];
    float3 sum = float3(0.0f);
    for (uint k = start + lane; k < end; k += 8) {
        sum += values[k] * float3(x[col_indices[k]]);
    }
    sum += simd_shuffle_down_float3(sum, 4);
    sum += simd_shuffle_down_float3(sum, 2);
    sum += simd_shuffle_down_float3(sum, 1);
    if (lane == 0) {
        gradient[row] = packed_float3(sum - float3(b[row]));
    }
}

// BB Step 3: 1D deltaG = nextG - gradient
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

// BB Step 3: 3D deltaG = nextG - gradient
kernel void bb_calc_delta_g_3d(
    device const packed_float3* nextG [[buffer(0)]],
    device const packed_float3* gradient [[buffer(1)]],
    device packed_float3* deltaG [[buffer(2)]],
    device const BBState* state [[buffer(3)]],
    constant uint& numNodes [[buffer(4)]],
    uint id [[thread_position_in_grid]])
{
    if (state->converged != 0 || id >= numNodes) return;
    deltaG[id] = packed_float3(float3(nextG[id]) - float3(gradient[id]));
}

// BB Step 4: 1D Fused quad dot reduction: (sTy, yTy, sTs, curGNormSq)
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

// BB Step 4: 3D Fused quad dot reduction: (sTy, yTy, sTs, curGNormSq)
kernel void bb_quad_dot_3d(
    device const packed_float3* deltaX [[buffer(0)]],
    device const packed_float3* deltaG [[buffer(1)]],
    device const packed_float3* nextG [[buffer(2)]],
    device float4* partialSums [[buffer(3)]],
    device const BBState* state [[buffer(4)]],
    constant uint& numNodes [[buffer(5)]],
    threadgroup float4* sharedData [[threadgroup(0)]],
    uint tid [[thread_index_in_threadgroup]],
    uint gid [[thread_position_in_grid]],
    uint tg_id [[threadgroup_position_in_grid]],
    uint tg_size [[threads_per_threadgroup]])
{
    if (state->converged != 0) return;
    float4 val = float4(0.0f);
    if (gid < numNodes) {
        float3 dx = float3(deltaX[gid]);
        float3 dg = float3(deltaG[gid]);
        float3 ng = float3(nextG[gid]);
        val = float4(dot(dx, dg), dot(dg, dg), dot(dx, dx), dot(ng, ng));
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

// BB Solution synchronization (1D)
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

// BB Solution synchronization (3D)
kernel void bb_sync_solution_3d(
    device packed_float3* x0 [[buffer(0)]],
    device const packed_float3* x1 [[buffer(1)]],
    device const BBState* state [[buffer(2)]],
    constant uint& numNodes [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= numNodes) return;
    if ((state->iter % 2) == 1) {
        x0[gid] = x1[gid];
    }
}

// ============================================================================
// Center-of-Mass (Zero-Mode) Projection (3D)
// ============================================================================

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

// ============================================================================
// End-to-End GPU Relaxation (Spring Strain & Inactivity Evaluation Kernel)
// ============================================================================

kernel void evaluate_springs_strain_and_inactivity(
    device const packed_float3* coords [[buffer(0)]],
    device const uint* springIndexA [[buffer(1)]],
    device const uint* springIndexB [[buffer(2)]],
    device const packed_float3* springBoxOffset [[buffer(3)]],
    device const float* springContourLength [[buffer(4)]],
    device uint* springFlag [[buffer(5)]],
    device atomic_uint* flaggedCounters [[buffer(6)]],
    constant float& breakingThreshold [[buffer(7)]],
    constant float& inactiveCutoff [[buffer(8)]],
    constant uint& numSprings [[buffer(9)]],
    uint s [[thread_position_in_grid]])
{
    if (s >= numSprings) return;
    uint idxA = springIndexA[s];
    uint idxB = springIndexB[s];
    if (idxA == idxB) {
        springFlag[s] = 0;
        return;
    }

    float3 posA = float3(coords[idxA]);
    float3 posB = float3(coords[idxB]);
    float3 offset = float3(springBoxOffset[s]);
    float3 distVec = (posB - posA) + offset;
    float len = length(distVec);
    float contour = springContourLength[s];

    uint flag = 0;
    if (breakingThreshold > 0.0f) {
        float weightedLen = (contour > 0.0f) ? (len / contour) : 0.0f;
        if (weightedLen > breakingThreshold) {
            flag |= 1u;
            atomic_fetch_add_explicit(&flaggedCounters[0], 1, memory_order_relaxed);
        }
    }

    if (inactiveCutoff > 0.0f) {
        float threshold = inactiveCutoff * max(contour, 1.0f);
        if (len <= threshold) {
            flag |= 2u;
            atomic_fetch_add_explicit(&flaggedCounters[1], 1, memory_order_relaxed);
        }
    }

    springFlag[s] = flag;
    if (flag != 0) {
        atomic_fetch_add_explicit(&flaggedCounters[2], 1, memory_order_relaxed);
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

  // 1D SpMV Pipelines
  id<MTLComputePipelineState> getSpmvPipe() const { return pipeSpmv; }
  id<MTLComputePipelineState> getSpmvCoopPipe() const { return pipeSpmvCoop; }
  id<MTLComputePipelineState> getSpmvCoopStatelessPipe() const
  {
    return pipeSpmvCoopStateless;
  }
  id<MTLComputePipelineState> getDotPipe() const { return pipeDot; }

  // 3D Block-SpMV Pipelines
  id<MTLComputePipelineState> getSpmm3rhsPipe() const { return pipeSpmm3rhs; }
  id<MTLComputePipelineState> getSpmm3rhsCoopPipe() const
  {
    return pipeSpmm3rhsCoop;
  }
  id<MTLComputePipelineState> getSpmm3rhsStatelessPipe() const
  {
    return pipeSpmm3rhsStateless;
  }
  id<MTLComputePipelineState> getDot3DPipe() const { return pipeDot3D; }

  // Jacobi Preconditioning Pipelines
  id<MTLComputePipelineState> getExtractInvDiagPipe() const
  {
    return pipeExtractInvDiag;
  }
  id<MTLComputePipelineState> getPcgJacobiDot1DPipe() const
  {
    return pipePcgJacobiDot1D;
  }
  id<MTLComputePipelineState> getPcgJacobiDot3DPipe() const
  {
    return pipePcgJacobiDot3D;
  }
  id<MTLComputePipelineState> getPcgUpdateP1DPipe() const
  {
    return pipePcgUpdateP1D;
  }
  id<MTLComputePipelineState> getPcgUpdateP3DPipe() const
  {
    return pipePcgUpdateP3D;
  }

  // CG pipelines
  id<MTLComputePipelineState> getCgAlphaPipe() const { return pipeCgAlpha; }
  id<MTLComputePipelineState> getCgUpdateXRPipe() const
  {
    return pipeCgUpdateXR;
  }
  id<MTLComputePipelineState> getCgUpdateXR3DPipe() const
  {
    return pipeCgUpdateXR3D;
  }
  id<MTLComputePipelineState> getCgBetaPipe() const { return pipeCgBeta; }
  id<MTLComputePipelineState> getCgUpdatePPipe() const { return pipeCgUpdateP; }

  // BB pipelines
  id<MTLComputePipelineState> getBbUpdateXPipe() const { return pipeBbUpdateX; }
  id<MTLComputePipelineState> getBbUpdateX3DPipe() const
  {
    return pipeBbUpdateX3D;
  }
  id<MTLComputePipelineState> getBbCalcGradPipe() const
  {
    return pipeBbCalcGrad;
  }
  id<MTLComputePipelineState> getBbCalcGrad3rhsPipe() const
  {
    return pipeBbCalcGrad3rhs;
  }
  id<MTLComputePipelineState> getBbCalcDeltaGPipe() const
  {
    return pipeBbCalcDeltaG;
  }
  id<MTLComputePipelineState> getBbCalcDeltaG3DPipe() const
  {
    return pipeBbCalcDeltaG3D;
  }
  id<MTLComputePipelineState> getBbQuadDotPipe() const { return pipeBbQuadDot; }
  id<MTLComputePipelineState> getBbQuadDot3DPipe() const
  {
    return pipeBbQuadDot3D;
  }
  id<MTLComputePipelineState> getBbAlphaPipe() const { return pipeBbAlpha; }
  id<MTLComputePipelineState> getBbSyncSolPipe() const { return pipeBbSyncSol; }
  id<MTLComputePipelineState> getBbSyncSol3DPipe() const
  {
    return pipeBbSyncSol3D;
  }

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

  // End-to-end Relaxation Pipeline
  id<MTLComputePipelineState> getEvalSpringsPipe() const
  {
    return pipeEvalSprings;
  }

  // Persistent Buffer Pooling: Matrix buffers
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

    if (!bufInvDiag || capInvDiag < numRows) {
      size_t newCap = std::max(numRows, static_cast<size_t>(capInvDiag * 1.5));
      bufInvDiag = [device newBufferWithLength:newCap * sizeof(float)
                                       options:MTLResourceStorageModeShared];
      capInvDiag = newCap;
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
      bufZ = [device newBufferWithLength:newCap * sizeof(float)
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

  void ensureSpringBuffers(size_t numSprings)
  {
    if (!bufSpringIndexA || capSprings < numSprings) {
      size_t newCap =
        std::max(numSprings, static_cast<size_t>(capSprings * 1.5));
      bufSpringIndexA =
        [device newBufferWithLength:newCap * sizeof(uint32_t)
                            options:MTLResourceStorageModeShared];
      bufSpringIndexB =
        [device newBufferWithLength:newCap * sizeof(uint32_t)
                            options:MTLResourceStorageModeShared];
      bufSpringBoxOffset =
        [device newBufferWithLength:newCap * 3 * sizeof(float)
                            options:MTLResourceStorageModeShared];
      bufSpringContour =
        [device newBufferWithLength:newCap * sizeof(float)
                            options:MTLResourceStorageModeShared];
      bufSpringFlags =
        [device newBufferWithLength:newCap * sizeof(uint32_t)
                            options:MTLResourceStorageModeShared];
      capSprings = newCap;
    }
    if (!bufFlaggedCounters) {
      bufFlaggedCounters =
        [device newBufferWithLength:3 * sizeof(uint32_t)
                            options:MTLResourceStorageModeShared];
    }
  }

  // Persistent Buffer Getters
  id<MTLBuffer> getBufRowPtr() const { return bufRowPtr; }
  id<MTLBuffer> getBufColIdx() const { return bufColIdx; }
  id<MTLBuffer> getBufValues() const { return bufValues; }
  id<MTLBuffer> getBufInvDiag() const { return bufInvDiag; }
  id<MTLBuffer> getBufB() const { return bufB; }
  id<MTLBuffer> getBufX() const { return bufX; }

  // CG buffers
  id<MTLBuffer> getBufR() const { return bufR; }
  id<MTLBuffer> getBufP() const { return bufP; }
  id<MTLBuffer> getBufAp() const { return bufAp; }
  id<MTLBuffer> getBufZ() const { return bufZ; }
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

  // Relaxation buffers
  id<MTLBuffer> getBufSpringIndexA() const { return bufSpringIndexA; }
  id<MTLBuffer> getBufSpringIndexB() const { return bufSpringIndexB; }
  id<MTLBuffer> getBufSpringBoxOffset() const { return bufSpringBoxOffset; }
  id<MTLBuffer> getBufSpringContour() const { return bufSpringContour; }
  id<MTLBuffer> getBufSpringFlags() const { return bufSpringFlags; }
  id<MTLBuffer> getBufFlaggedCounters() const { return bufFlaggedCounters; }

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

  // Pipeline states (1D)
  id<MTLComputePipelineState> pipeSpmv = nil;
  id<MTLComputePipelineState> pipeSpmvCoop = nil;
  id<MTLComputePipelineState> pipeSpmvCoopStateless = nil;
  id<MTLComputePipelineState> pipeDot = nil;

  // Pipeline states (3D Block-SpMV)
  id<MTLComputePipelineState> pipeSpmm3rhs = nil;
  id<MTLComputePipelineState> pipeSpmm3rhsCoop = nil;
  id<MTLComputePipelineState> pipeSpmm3rhsStateless = nil;
  id<MTLComputePipelineState> pipeDot3D = nil;

  // Jacobi pipelines
  id<MTLComputePipelineState> pipeExtractInvDiag = nil;
  id<MTLComputePipelineState> pipePcgJacobiDot1D = nil;
  id<MTLComputePipelineState> pipePcgJacobiDot3D = nil;
  id<MTLComputePipelineState> pipePcgUpdateP1D = nil;
  id<MTLComputePipelineState> pipePcgUpdateP3D = nil;

  // CG pipelines
  id<MTLComputePipelineState> pipeCgAlpha = nil;
  id<MTLComputePipelineState> pipeCgUpdateXR = nil;
  id<MTLComputePipelineState> pipeCgUpdateXR3D = nil;
  id<MTLComputePipelineState> pipeCgBeta = nil;
  id<MTLComputePipelineState> pipeCgUpdateP = nil;

  // BB pipelines
  id<MTLComputePipelineState> pipeBbUpdateX = nil;
  id<MTLComputePipelineState> pipeBbUpdateX3D = nil;
  id<MTLComputePipelineState> pipeBbCalcGrad = nil;
  id<MTLComputePipelineState> pipeBbCalcGrad3rhs = nil;
  id<MTLComputePipelineState> pipeBbCalcDeltaG = nil;
  id<MTLComputePipelineState> pipeBbCalcDeltaG3D = nil;
  id<MTLComputePipelineState> pipeBbQuadDot = nil;
  id<MTLComputePipelineState> pipeBbQuadDot3D = nil;
  id<MTLComputePipelineState> pipeBbAlpha = nil;
  id<MTLComputePipelineState> pipeBbSyncSol = nil;
  id<MTLComputePipelineState> pipeBbSyncSol3D = nil;

  // COM projection pipelines
  id<MTLComputePipelineState> pipeComReduce3D = nil;
  id<MTLComputePipelineState> pipeComFinalize3D = nil;
  id<MTLComputePipelineState> pipeComProject3D = nil;
  id<MTLComputePipelineState> pipeComReduce1D = nil;
  id<MTLComputePipelineState> pipeComFinalize1D = nil;
  id<MTLComputePipelineState> pipeComProject1D = nil;

  // End-to-end Relaxation Pipeline
  id<MTLComputePipelineState> pipeEvalSprings = nil;

  // Persistent MTLBuffers
  id<MTLBuffer> bufRowPtr = nil;
  id<MTLBuffer> bufColIdx = nil;
  id<MTLBuffer> bufValues = nil;
  id<MTLBuffer> bufInvDiag = nil;
  id<MTLBuffer> bufB = nil;
  id<MTLBuffer> bufX = nil;

  id<MTLBuffer> bufR = nil;
  id<MTLBuffer> bufP = nil;
  id<MTLBuffer> bufAp = nil;
  id<MTLBuffer> bufZ = nil;
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

  id<MTLBuffer> bufSpringIndexA = nil;
  id<MTLBuffer> bufSpringIndexB = nil;
  id<MTLBuffer> bufSpringBoxOffset = nil;
  id<MTLBuffer> bufSpringContour = nil;
  id<MTLBuffer> bufSpringFlags = nil;
  id<MTLBuffer> bufFlaggedCounters = nil;

  // Tracked capacities
  size_t capRowPtr = 0;
  size_t capColIdx = 0;
  size_t capValues = 0;
  size_t capInvDiag = 0;
  size_t capN = 0;
  size_t capCgVectors = 0;
  size_t capBbVectors = 0;
  size_t capPartial = 0;
  size_t capQuadPartial = 0;
  size_t capComPartial = 0;
  size_t capSprings = 0;

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

      // 1D SpMV pipelines
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

      // 3D Block-SpMV pipelines
      pipeSpmm3rhs = [device
        newComputePipelineStateWithFunction:[library
                                              newFunctionWithName:@"csr_spmm_3rhs"]
                                       error:&error];
      pipeSpmm3rhsCoop =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"csr_spmm_3rhs_coop8"]
                                              error:&error];
      pipeSpmm3rhsStateless =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"csr_spmm_3rhs_stateless"]
                                              error:&error];
      pipeDot3D = [device newComputePipelineStateWithFunction:
                            [library newFunctionWithName:@"dot_product_3d"]
                                                        error:&error];

      // Jacobi Preconditioning
      pipeExtractInvDiag =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"extract_inv_diagonal"]
                                              error:&error];
      pipePcgJacobiDot1D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"pcg_apply_jacobi_and_dot_1d"]
                                              error:&error];
      pipePcgJacobiDot3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"pcg_apply_jacobi_and_dot_3d"]
                                              error:&error];
      pipePcgUpdateP1D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"pcg_update_p_1d"]
                                              error:&error];
      pipePcgUpdateP3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"pcg_update_p_3d"]
                                              error:&error];

      // CG pipelines
      pipeCgAlpha = [device newComputePipelineStateWithFunction:
                              [library newFunctionWithName:@"cg_compute_alpha"]
                                                          error:&error];
      pipeCgUpdateXR = [device newComputePipelineStateWithFunction:
                                 [library newFunctionWithName:@"cg_update_x_r"]
                                                              error:&error];
      pipeCgUpdateXR3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"cg_update_x_r_3d"]
                                              error:&error];
      pipeCgBeta = [device newComputePipelineStateWithFunction:
                             [library newFunctionWithName:@"cg_compute_beta"]
                                                         error:&error];
      pipeCgUpdateP = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"cg_update_p"]
                                                            error:&error];

      // BB pipelines
      pipeBbUpdateX = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"bb_update_x"]
                                                            error:&error];
      pipeBbUpdateX3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_update_x_3d"]
                                              error:&error];
      pipeBbCalcGrad =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_gradient_coop8"]
                                              error:&error];
      pipeBbCalcGrad3rhs =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_gradient_3rhs_coop8"]
                                              error:&error];
      pipeBbCalcDeltaG =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_delta_g"]
                                              error:&error];
      pipeBbCalcDeltaG3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_calc_delta_g_3d"]
                                              error:&error];
      pipeBbQuadDot = [device newComputePipelineStateWithFunction:
                                [library newFunctionWithName:@"bb_quad_dot"]
                                                            error:&error];
      pipeBbQuadDot3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_quad_dot_3d"]
                                              error:&error];
      pipeBbAlpha = [device newComputePipelineStateWithFunction:
                              [library newFunctionWithName:@"bb_compute_alpha"]
                                                          error:&error];
      pipeBbSyncSol =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_sync_solution"]
                                              error:&error];
      pipeBbSyncSol3D =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"bb_sync_solution_3d"]
                                              error:&error];

      // COM projection pipelines
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

      // End-to-end Relaxation Pipeline
      pipeEvalSprings =
        [device newComputePipelineStateWithFunction:
                  [library newFunctionWithName:@"evaluate_springs_strain_and_inactivity"]
                                              error:&error];

      if (pipeSpmv && pipeSpmvCoop && pipeSpmvCoopStateless && pipeDot &&
          pipeSpmm3rhs && pipeSpmm3rhsCoop && pipeSpmm3rhsStateless && pipeDot3D &&
          pipeExtractInvDiag && pipePcgJacobiDot1D && pipePcgJacobiDot3D &&
          pipePcgUpdateP1D && pipePcgUpdateP3D && pipeCgAlpha &&
          pipeCgUpdateXR && pipeCgUpdateXR3D && pipeCgBeta && pipeCgUpdateP &&
          pipeBbUpdateX && pipeBbUpdateX3D && pipeBbCalcGrad &&
          pipeBbCalcGrad3rhs && pipeBbCalcDeltaG && pipeBbCalcDeltaG3D &&
          pipeBbQuadDot && pipeBbQuadDot3D && pipeBbAlpha && pipeBbSyncSol &&
          pipeBbSyncSol3D && pipeComReduce3D && pipeComFinalize3D &&
          pipeComProject3D && pipeComReduce1D && pipeComFinalize1D &&
          pipeComProject1D && pipeEvalSprings) {
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

// ============================================================================
// Multi-RHS 3D Block Barzilai-Borwein Solver
// ============================================================================

Eigen::VectorXd
solveSLEMetalBarzilaiBorwein3D(const Eigen::SparseMatrix<double>& scalarA,
                              const Eigen::VectorXd& b,
                              double learningRate,
                              double tolerance,
                              int maxIterations,
                              int& iterationsDone,
                              bool shortAlpha,
                              const Eigen::VectorXd& initialX,
                              bool projectZeroModes)
{
  Eigen::Index numNodesIndex = scalarA.rows();
  uint32_t numNodes = static_cast<uint32_t>(numNodesIndex);
  uint32_t N = numNodes * 3;

  if (numNodes == 0) {
    iterationsDone = 0;
    return Eigen::VectorXd::Zero(0);
  }
  if (b.isZero(1e-9)) {
    iterationsDone = 1;
    return Eigen::VectorXd::Zero(N);
  }

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t tg_size = 256;
    uint32_t num_tgs = (numNodes + tg_size - 1) / tg_size;
    size_t nnz = static_cast<size_t>(scalarA.nonZeros());

    // 1. Matrix Buffers (N x N scalar Laplacian)
    backend.ensureMatrixBuffers(numNodes, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();

    uint32_t* rowPtrsMem = static_cast<uint32_t*>([bufRowPtr contents]);
    for (uint32_t i = 0; i <= numNodes; ++i) {
      rowPtrsMem[i] = static_cast<uint32_t>(scalarA.outerIndexPtr()[i]);
    }
    uint32_t* colIdxMem = static_cast<uint32_t*>([bufColIdx contents]);
    float* valuesMem = static_cast<float*>([bufValues contents]);
    for (size_t i = 0; i < nnz; ++i) {
      colIdxMem[i] = static_cast<uint32_t>(scalarA.innerIndexPtr()[i]);
      valuesMem[i] = static_cast<float>(scalarA.valuePtr()[i]);
    }

    // 2. Vector Buffers (3N floats = numNodes packed_float3)
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
    if (initialX.size() == N) {
      for (uint32_t i = 0; i < N; ++i)
        x0Mem[i] = static_cast<float>(initialX[i]);
    } else {
      std::memset(x0Mem, 0, sizeof(float) * N);
    }

    // Compute initial gradient: g0 = A*x0 - b using 3-RHS SpMM
    uint32_t totalSpmvThreads = numNodes * 8;
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getBbCalcGrad3rhsPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX0 offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufG0 offset:0 atIndex:5];
      [enc setBuffer:nil offset:0 atIndex:6]; // stateless
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:7];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      // Initial gNorm using 3D dot product
      [enc setComputePipelineState:backend.getDot3DPipe()];
      [enc setBuffer:bufG0 offset:0 atIndex:0];
      [enc setBuffer:bufG0 offset:0 atIndex:1];
      [enc setBuffer:bufPartial offset:0 atIndex:2];
      [enc setBuffer:nil offset:0 atIndex:3];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
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
      return initialX.size() == N ? initialX : Eigen::VectorXd::Zero(N);
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

        // Step 1: nextX = x - alpha * g, deltaX = -alpha * g (3D)
        [enc setComputePipelineState:backend.getBbUpdateX3DPipe()];
        [enc setBuffer:curX offset:0 atIndex:0];
        [enc setBuffer:curG offset:0 atIndex:1];
        [enc setBuffer:nextX offset:0 atIndex:2];
        [enc setBuffer:bufDeltaX offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 2: 3-RHS SIMDgroup-Cooperative nextG = A * nextX - b
        [enc setComputePipelineState:backend.getBbCalcGrad3rhsPipe()];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:nextX offset:0 atIndex:3];
        [enc setBuffer:bufB offset:0 atIndex:4];
        [enc setBuffer:nextG offset:0 atIndex:5];
        [enc setBuffer:bufState offset:0 atIndex:6];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:7];
        [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 3: deltaG = nextG - curG (3D)
        [enc setComputePipelineState:backend.getBbCalcDeltaG3DPipe()];
        [enc setBuffer:nextG offset:0 atIndex:0];
        [enc setBuffer:curG offset:0 atIndex:1];
        [enc setBuffer:bufDeltaG offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
        [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 4: Fused 3D quad dot reduction (sTy, yTy, sTs, curGNormSq)
        [enc setComputePipelineState:backend.getBbQuadDot3DPipe()];
        [enc setBuffer:bufDeltaX offset:0 atIndex:0];
        [enc setBuffer:bufDeltaG offset:0 atIndex:1];
        [enc setBuffer:nextG offset:0 atIndex:2];
        [enc setBuffer:bufQuadPartial offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
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

        // COM Projection every 25 iterations
        if (projectZeroModes && ((iterIdx + 1) % 25 == 0)) {
          backend.encodeComProjection(enc, nextX, N);
        }
      }

      // Sync final solution to bufX0 if ended on odd step
      [enc setComputePipelineState:backend.getBbSyncSol3DPipe()];
      [enc setBuffer:bufX0 offset:0 atIndex:0];
      [enc setBuffer:bufX1 offset:0 atIndex:1];
      [enc setBuffer:bufState offset:0 atIndex:2];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:3];
      [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

      // Periodic gradient recomputation every batch
      [enc setComputePipelineState:backend.getBbCalcGrad3rhsPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX0 offset:0 atIndex:3];
      [enc setBuffer:bufB offset:0 atIndex:4];
      [enc setBuffer:bufG0 offset:0 atIndex:5];
      [enc setBuffer:bufState offset:0 atIndex:6];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:7];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

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
    Eigen::VectorXd res(N);
    for (uint32_t i = 0; i < N; ++i) {
      res[i] = static_cast<double>(finalXPtr[i]);
    }
    return res;
  }
}

// ============================================================================
// Multi-RHS 3D Block Conjugate Gradient Solver (Fused 3-RHS SpMM + Jacobi PCG)
// ============================================================================

Eigen::VectorXd
solveSLEMetalCG3D(const Eigen::SparseMatrix<double>& scalarA,
                  const Eigen::VectorXd& b,
                  double tolerance,
                  int maxIterations,
                  int& iterationsDone,
                  const Eigen::VectorXd& initialX,
                  bool projectZeroModes,
                  bool useJacobi)
{
  Eigen::Index numNodesIndex = scalarA.rows();
  uint32_t numNodes = static_cast<uint32_t>(numNodesIndex);
  uint32_t N = numNodes * 3;

  if (numNodes == 0) {
    iterationsDone = 0;
    return Eigen::VectorXd::Zero(0);
  }
  if (b.isZero(1e-9)) {
    iterationsDone = 1;
    return Eigen::VectorXd::Zero(N);
  }

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(),
                  "Metal GPU backend is not available on this device");

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t tg_size = 256;
    uint32_t num_tgs = (numNodes + tg_size - 1) / tg_size;
    size_t nnz = static_cast<size_t>(scalarA.nonZeros());

    // 1. Matrix Buffers (N x N scalar Laplacian)
    backend.ensureMatrixBuffers(numNodes, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();
    id<MTLBuffer> bufInvDiag = backend.getBufInvDiag();

    uint32_t* rowPtrsMem = static_cast<uint32_t*>([bufRowPtr contents]);
    for (uint32_t i = 0; i <= numNodes; ++i) {
      rowPtrsMem[i] = static_cast<uint32_t>(scalarA.outerIndexPtr()[i]);
    }
    uint32_t* colIdxMem = static_cast<uint32_t*>([bufColIdx contents]);
    float* valuesMem = static_cast<float*>([bufValues contents]);
    for (size_t i = 0; i < nnz; ++i) {
      colIdxMem[i] = static_cast<uint32_t>(scalarA.innerIndexPtr()[i]);
      valuesMem[i] = static_cast<float>(scalarA.valuePtr()[i]);
    }

    // 2. Vector Buffers (3N floats = numNodes packed_float3)
    backend.ensureCgBuffers(N, num_tgs);
    id<MTLBuffer> bufB = backend.getBufB();
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufR = backend.getBufR();
    id<MTLBuffer> bufP = backend.getBufP();
    id<MTLBuffer> bufAp = backend.getBufAp();
    id<MTLBuffer> bufZ = backend.getBufZ();
    id<MTLBuffer> bufPartial = backend.getBufPartial();
    id<MTLBuffer> bufState = backend.getBufCGState();

    float* bMem = static_cast<float*>([bufB contents]);
    for (uint32_t i = 0; i < N; ++i) {
      bMem[i] = static_cast<float>(b[i]);
    }

    float* xPtr = static_cast<float*>([bufX contents]);
    float* rPtr = static_cast<float*>([bufR contents]);
    float* pPtr = static_cast<float*>([bufP contents]);

    uint32_t totalSpmvThreads = numNodes * 8;

    // Extract inverse diagonal on GPU if Jacobi preconditioning is enabled
    if (useJacobi) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getExtractInvDiagPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufInvDiag offset:0 atIndex:3];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
      [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    } else {
      float* invDMem = static_cast<float*>([bufInvDiag contents]);
      std::fill_n(invDMem, numNodes, 1.0f);
    }

    if (initialX.size() == N) {
      for (uint32_t i = 0; i < N; ++i)
        xPtr[i] = static_cast<float>(initialX[i]);

      // r = b - A*x using 3-RHS stateless SpMM
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getSpmm3rhsStatelessPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufAp offset:0 atIndex:4];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
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

    // Initial preconditioner & rho calculation: z = invDiag * r, rho = dot(r, z)
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      if (useJacobi) {
        [enc setComputePipelineState:backend.getPcgJacobiDot3DPipe()];
        [enc setBuffer:bufInvDiag offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufZ offset:0 atIndex:2];
        [enc setBuffer:bufPartial offset:0 atIndex:3];
        [enc setBuffer:nil offset:0 atIndex:4];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      } else {
        [enc setComputePipelineState:backend.getDot3DPipe()];
        [enc setBuffer:bufR offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBuffer:nil offset:0 atIndex:3];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      }
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    float* parts = static_cast<float*>([bufPartial contents]);
    float initial_rho = 0.0f;
    for (uint32_t i = 0; i < num_tgs; ++i)
      initial_rho += parts[i];

    if (std::abs(initial_rho) < 1e-24f) {
      iterationsDone = 0;
      return initialX.size() == N ? initialX : Eigen::VectorXd::Zero(N);
    }

    // Copy z to p for initial search direction if Jacobi is enabled
    if (useJacobi) {
      float* zPtr = static_cast<float*>([bufZ contents]);
      std::memcpy(pPtr, zPtr, sizeof(float) * N);
    }

    CGState* statePtr = static_cast<CGState*>([bufState contents]);
    statePtr->iter = 0;
    statePtr->converged = 0;
    statePtr->maxIters = static_cast<uint32_t>(maxIterations);
    statePtr->tolSq = static_cast<float>(tolerance * tolerance);
    statePtr->rho = initial_rho;
    statePtr->initial_rho = initial_rho;

    const uint32_t batchSize = 25;

    // Batched GPU iteration loop
    while (statePtr->converged == 0 && statePtr->iter < statePtr->maxIters) {
      uint32_t remaining = statePtr->maxIters - statePtr->iter;
      uint32_t itersInBatch = std::min(batchSize, remaining);

      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

      for (uint32_t step = 0; step < itersInBatch; ++step) {
        // Step 1: SIMDgroup-Cooperative 3-RHS SpMM Ap = A * p
        [enc setComputePipelineState:backend.getSpmm3rhsCoopPipe()];
        [enc setBuffer:bufRowPtr offset:0 atIndex:0];
        [enc setBuffer:bufColIdx offset:0 atIndex:1];
        [enc setBuffer:bufValues offset:0 atIndex:2];
        [enc setBuffer:bufP offset:0 atIndex:3];
        [enc setBuffer:bufAp offset:0 atIndex:4];
        [enc setBuffer:bufState offset:0 atIndex:5];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:6];
        [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 2: Dot pAp = dot(p, Ap) (3D)
        [enc setComputePipelineState:backend.getDot3DPipe()];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:bufAp offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
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

        // Step 4: Update x += alpha * p, r -= alpha * Ap (3D)
        [enc setComputePipelineState:backend.getCgUpdateXR3DPipe()];
        [enc setBuffer:bufX offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufP offset:0 atIndex:2];
        [enc setBuffer:bufAp offset:0 atIndex:3];
        [enc setBuffer:bufState offset:0 atIndex:4];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
        [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        // Step 5: Preconditioning & dot: z = invDiag * r, rho_new = dot(r, z)
        if (useJacobi) {
          [enc setComputePipelineState:backend.getPcgJacobiDot3DPipe()];
          [enc setBuffer:bufInvDiag offset:0 atIndex:0];
          [enc setBuffer:bufR offset:0 atIndex:1];
          [enc setBuffer:bufZ offset:0 atIndex:2];
          [enc setBuffer:bufPartial offset:0 atIndex:3];
          [enc setBuffer:bufState offset:0 atIndex:4];
          [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
          [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
          [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        } else {
          [enc setComputePipelineState:backend.getDot3DPipe()];
          [enc setBuffer:bufR offset:0 atIndex:0];
          [enc setBuffer:bufR offset:0 atIndex:1];
          [enc setBuffer:bufPartial offset:0 atIndex:2];
          [enc setBuffer:bufState offset:0 atIndex:3];
          [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:4];
          [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
          [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        }
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

        // Step 7: Update p = z + beta * p (3D)
        [enc setComputePipelineState:backend.getPcgUpdateP3DPipe()];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:useJacobi ? bufZ : bufR offset:0 atIndex:1];
        [enc setBuffer:bufState offset:0 atIndex:2];
        [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:3];
        [enc dispatchThreads:MTLSizeMake(numNodes, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(numNodes, 256u), 1, 1)];
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
    Eigen::VectorXd res(N);
    for (uint32_t i = 0; i < N; ++i) {
      res[i] = static_cast<double>(finalXPtr[i]);
    }
    return res;
  }
}

// ============================================================================
// General / 1D Solvers
// ============================================================================

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
  // Check if A is an N x N scalar Laplacian for a 3N vector b
  if (b.size() == A.rows() * 3 && A.rows() > 0) {
    return solveSLEMetalBarzilaiBorwein3D(A,
                                         b,
                                         learningRate,
                                         tolerance,
                                         maxIterations,
                                         iterationsDone,
                                         shortAlpha,
                                         initialX,
                                         projectZeroModes);
  }

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
      [enc setBuffer:nil offset:0 atIndex:6];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:7];
      [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];

      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

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

        [enc setComputePipelineState:backend.getBbCalcDeltaGPipe()];
        [enc setBuffer:nextG offset:0 atIndex:0];
        [enc setBuffer:curG offset:0 atIndex:1];
        [enc setBuffer:bufDeltaG offset:0 atIndex:2];
        [enc setBuffer:bufState offset:0 atIndex:3];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
        [enc dispatchThreads:MTLSizeMake(N, 1, 1)
          threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

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

        [enc setComputePipelineState:backend.getBbAlphaPipe()];
        [enc setBuffer:bufQuadPartial offset:0 atIndex:0];
        [enc setBuffer:bufState offset:0 atIndex:1];
        [enc setBytes:&num_tgs length:sizeof(uint32_t) atIndex:2];
        [enc setThreadgroupMemoryLength:sizeof(simd_float4) * tg_size
                                atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(1, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

        if (projectZeroModes && ((iterIdx + 1) % 25 == 0)) {
          backend.encodeComProjection(enc, nextX, N);
        }
      }

      [enc setComputePipelineState:backend.getBbSyncSolPipe()];
      [enc setBuffer:bufX0 offset:0 atIndex:0];
      [enc setBuffer:bufX1 offset:0 atIndex:1];
      [enc setBuffer:bufState offset:0 atIndex:2];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:3];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];

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
                bool projectZeroModes,
                bool useJacobi)
{
  // Check if A is an N x N scalar Laplacian for a 3N vector b
  if (b.size() == A.rows() * 3 && A.rows() > 0) {
    return solveSLEMetalCG3D(A,
                             b,
                             tolerance,
                             maxIterations,
                             iterationsDone,
                             initialX,
                             projectZeroModes,
                             useJacobi);
  }

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

    backend.ensureMatrixBuffers(N, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();
    id<MTLBuffer> bufInvDiag = backend.getBufInvDiag();

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

    backend.ensureCgBuffers(N, num_tgs);
    id<MTLBuffer> bufB = backend.getBufB();
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufR = backend.getBufR();
    id<MTLBuffer> bufP = backend.getBufP();
    id<MTLBuffer> bufAp = backend.getBufAp();
    id<MTLBuffer> bufZ = backend.getBufZ();
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

    // Extract inverse diagonal on GPU if Jacobi preconditioning is enabled
    if (useJacobi) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getExtractInvDiagPipe()];
      [enc setBuffer:bufRowPtr offset:0 atIndex:0];
      [enc setBuffer:bufColIdx offset:0 atIndex:1];
      [enc setBuffer:bufValues offset:0 atIndex:2];
      [enc setBuffer:bufInvDiag offset:0 atIndex:3];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
      [enc dispatchThreads:MTLSizeMake(N, 1, 1)
        threadsPerThreadgroup:MTLSizeMake(std::min(N, 256u), 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    } else {
      float* invDMem = static_cast<float*>([bufInvDiag contents]);
      std::fill_n(invDMem, N, 1.0f);
    }

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

    // Initial preconditioner & rho calculation: z = invDiag * r, rho = dot(r, z)
    {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      if (useJacobi) {
        [enc setComputePipelineState:backend.getPcgJacobiDot1DPipe()];
        [enc setBuffer:bufInvDiag offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufZ offset:0 atIndex:2];
        [enc setBuffer:bufPartial offset:0 atIndex:3];
        [enc setBuffer:nil offset:0 atIndex:4];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      } else {
        [enc setComputePipelineState:backend.getDotPipe()];
        [enc setBuffer:bufR offset:0 atIndex:0];
        [enc setBuffer:bufR offset:0 atIndex:1];
        [enc setBuffer:bufPartial offset:0 atIndex:2];
        [enc setBuffer:nil offset:0 atIndex:3];
        [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
        [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
        [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
            threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      }
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    float* parts = static_cast<float*>([bufPartial contents]);
    float initial_rho = 0.0f;
    for (uint32_t i = 0; i < num_tgs; ++i)
      initial_rho += parts[i];

    if (std::abs(initial_rho) < 1e-24f) {
      iterationsDone = 0;
      return initialX.size() == n ? initialX : Eigen::VectorXd::Zero(n);
    }

    if (useJacobi) {
      float* zPtr = static_cast<float*>([bufZ contents]);
      std::memcpy(pPtr, zPtr, sizeof(float) * N);
    }

    CGState* statePtr = static_cast<CGState*>([bufState contents]);
    statePtr->iter = 0;
    statePtr->converged = 0;
    statePtr->maxIters = static_cast<uint32_t>(maxIterations);
    statePtr->tolSq = static_cast<float>(tolerance * tolerance);
    statePtr->rho = initial_rho;
    statePtr->initial_rho = initial_rho;

    const uint32_t batchSize = 25;

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

        // Step 5: Preconditioning & dot: z = invDiag * r, rho_new = dot(r, z)
        if (useJacobi) {
          [enc setComputePipelineState:backend.getPcgJacobiDot1DPipe()];
          [enc setBuffer:bufInvDiag offset:0 atIndex:0];
          [enc setBuffer:bufR offset:0 atIndex:1];
          [enc setBuffer:bufZ offset:0 atIndex:2];
          [enc setBuffer:bufPartial offset:0 atIndex:3];
          [enc setBuffer:bufState offset:0 atIndex:4];
          [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
          [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
          [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        } else {
          [enc setComputePipelineState:backend.getDotPipe()];
          [enc setBuffer:bufR offset:0 atIndex:0];
          [enc setBuffer:bufR offset:0 atIndex:1];
          [enc setBuffer:bufPartial offset:0 atIndex:2];
          [enc setBuffer:bufState offset:0 atIndex:3];
          [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];
          [enc setThreadgroupMemoryLength:sizeof(float) * tg_size atIndex:0];
          [enc dispatchThreadgroups:MTLSizeMake(num_tgs, 1, 1)
              threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
        }
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

        // Step 7: Update p = z + beta * p
        [enc setComputePipelineState:backend.getPcgUpdateP1DPipe()];
        [enc setBuffer:bufP offset:0 atIndex:0];
        [enc setBuffer:useJacobi ? bufZ : bufR offset:0 atIndex:1];
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

// ============================================================================
// End-to-End GPU Relaxation Evaluation
// ============================================================================

SpringEvaluationResult
evaluateSpringStrainAndInactivityMetal(size_t numNodes,
                                       const SpringRelaxationData& springData)
{
  SpringEvaluationResult result;
  result.numBroken = 0;
  result.numInactive = 0;
  result.totalFlagged = 0;
  result.springFlags.resize(springData.numSprings, 0);

  if (springData.numSprings == 0) {
    return result;
  }

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  if (!backend.isAvailable()) {
    return result;
  }

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t S = static_cast<uint32_t>(springData.numSprings);

    backend.ensureSpringBuffers(S);
    id<MTLBuffer> bufCoords = backend.getBufX();
    id<MTLBuffer> bufIdxA = backend.getBufSpringIndexA();
    id<MTLBuffer> bufIdxB = backend.getBufSpringIndexB();
    id<MTLBuffer> bufBoxOffset = backend.getBufSpringBoxOffset();
    id<MTLBuffer> bufContour = backend.getBufSpringContour();
    id<MTLBuffer> bufFlags = backend.getBufSpringFlags();
    id<MTLBuffer> bufCounters = backend.getBufFlaggedCounters();

    uint32_t* idxAMem = static_cast<uint32_t*>([bufIdxA contents]);
    uint32_t* idxBMem = static_cast<uint32_t*>([bufIdxB contents]);
    float* boxOffsetMem = static_cast<float*>([bufBoxOffset contents]);
    float* contourMem = static_cast<float*>([bufContour contents]);

    for (size_t s = 0; s < springData.numSprings; ++s) {
      idxAMem[s] = static_cast<uint32_t>(springData.springIndexA[s]);
      idxBMem[s] = static_cast<uint32_t>(springData.springIndexB[s]);
      boxOffsetMem[3 * s + 0] = static_cast<float>(springData.springBoxOffset[3 * s + 0]);
      boxOffsetMem[3 * s + 1] = static_cast<float>(springData.springBoxOffset[3 * s + 1]);
      boxOffsetMem[3 * s + 2] = static_cast<float>(springData.springBoxOffset[3 * s + 2]);
      contourMem[s] = static_cast<float>(springData.springContourLength[s]);
    }

    uint32_t* countersMem = static_cast<uint32_t*>([bufCounters contents]);
    countersMem[0] = 0;
    countersMem[1] = 0;
    countersMem[2] = 0;

    float breakingThreshold = static_cast<float>(springData.breakingThreshold);
    float inactiveCutoff = static_cast<float>(springData.inactiveCutoff);

    id<MTLCommandBuffer> cmd = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:backend.getEvalSpringsPipe()];
    [enc setBuffer:bufCoords offset:0 atIndex:0];
    [enc setBuffer:bufIdxA offset:0 atIndex:1];
    [enc setBuffer:bufIdxB offset:0 atIndex:2];
    [enc setBuffer:bufBoxOffset offset:0 atIndex:3];
    [enc setBuffer:bufContour offset:0 atIndex:4];
    [enc setBuffer:bufFlags offset:0 atIndex:5];
    [enc setBuffer:bufCounters offset:0 atIndex:6];
    [enc setBytes:&breakingThreshold length:sizeof(float) atIndex:7];
    [enc setBytes:&inactiveCutoff length:sizeof(float) atIndex:8];
    [enc setBytes:&S length:sizeof(uint32_t) atIndex:9];

    uint32_t tg_size = 256;
    [enc dispatchThreads:MTLSizeMake(S, 1, 1)
      threadsPerThreadgroup:MTLSizeMake(std::min(S, tg_size), 1, 1)];

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    result.numBroken = countersMem[0];
    result.numInactive = countersMem[1];
    result.totalFlagged = countersMem[2];

    if (result.totalFlagged > 0) {
      uint32_t* flagsMem = static_cast<uint32_t*>([bufFlags contents]);
      for (size_t s = 0; s < springData.numSprings; ++s) {
        result.springFlags[s] = flagsMem[s];
      }
    }

    return result;
  }
}

Eigen::VectorXd
getGpuResidentCoordinates(size_t numCoords)
{
  MetalFBSBackend& backend = MetalFBSBackend::instance();
  if (!backend.isAvailable()) {
    return Eigen::VectorXd::Zero(numCoords);
  }

  id<MTLBuffer> bufX = backend.getBufX();
  if (!bufX) {
    return Eigen::VectorXd::Zero(numCoords);
  }

  float* xPtr = static_cast<float*>([bufX contents]);
  Eigen::VectorXd res(numCoords);
  for (size_t i = 0; i < numCoords; ++i) {
    res[i] = static_cast<double>(xPtr[i]);
  }
  return res;
}

Eigen::VectorXd
spmvMetal3D(const Eigen::SparseMatrix<double>& scalarA, const Eigen::VectorXd& x)
{
  uint32_t numNodes = static_cast<uint32_t>(scalarA.rows());
  uint32_t N = numNodes * 3;
  RUNTIME_EXP_IFN(static_cast<size_t>(x.size()) == static_cast<size_t>(N),
                  "Vector x size must equal scalarA.rows() * 3 for 3D SpMV");

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(), "Metal backend unavailable");

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t tg_size = 256;
    size_t nnz = static_cast<size_t>(scalarA.nonZeros());

    backend.ensureMatrixBuffers(numNodes, nnz);
    id<MTLBuffer> bufRowPtr = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx = backend.getBufColIdx();
    id<MTLBuffer> bufValues = backend.getBufValues();

    uint32_t* rowPtrsMem = static_cast<uint32_t*>([bufRowPtr contents]);
    for (uint32_t i = 0; i <= numNodes; ++i) {
      rowPtrsMem[i] = static_cast<uint32_t>(scalarA.outerIndexPtr()[i]);
    }
    uint32_t* colIdxMem = static_cast<uint32_t*>([bufColIdx contents]);
    float* valuesMem = static_cast<float*>([bufValues contents]);
    for (size_t i = 0; i < nnz; ++i) {
      colIdxMem[i] = static_cast<uint32_t>(scalarA.innerIndexPtr()[i]);
      valuesMem[i] = static_cast<float>(scalarA.valuePtr()[i]);
    }

    backend.ensureCgBuffers(N, (numNodes + tg_size - 1) / tg_size);
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufAp = backend.getBufAp();

    float* xPtr = static_cast<float*>([bufX contents]);
    for (uint32_t i = 0; i < N; ++i) {
      xPtr[i] = static_cast<float>(x[i]);
    }

    uint32_t totalSpmvThreads = numNodes * 8;
    id<MTLCommandBuffer> cmd = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:backend.getSpmm3rhsStatelessPipe()];
    [enc setBuffer:bufRowPtr offset:0 atIndex:0];
    [enc setBuffer:bufColIdx offset:0 atIndex:1];
    [enc setBuffer:bufValues offset:0 atIndex:2];
    [enc setBuffer:bufX offset:0 atIndex:3];
    [enc setBuffer:bufAp offset:0 atIndex:4];
    [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
    [enc dispatchThreads:MTLSizeMake(totalSpmvThreads, 1, 1)
      threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    float* apPtr = static_cast<float*>([bufAp contents]);
    Eigen::VectorXd y(N);
    for (uint32_t i = 0; i < N; ++i) {
      y[i] = static_cast<double>(apPtr[i]);
    }
    return y;
  }
}

Eigen::VectorXd
spmvMetal(const Eigen::SparseMatrix<double>& A, const Eigen::VectorXd& x)
{
  uint32_t N = static_cast<uint32_t>(A.rows());
  RUNTIME_EXP_IFN(static_cast<size_t>(x.size()) == static_cast<size_t>(N),
                  "Vector x size must equal A.rows()");

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  RUNTIME_EXP_IFN(backend.isAvailable(), "Metal backend unavailable");

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t tg_size = 256;
    size_t nnz = static_cast<size_t>(A.nonZeros());

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

    backend.ensureCgBuffers(N, (N + tg_size - 1) / tg_size);
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufAp = backend.getBufAp();

    float* xPtr = static_cast<float*>([bufX contents]);
    for (uint32_t i = 0; i < N; ++i) {
      xPtr[i] = static_cast<float>(x[i]);
    }

    uint32_t totalSpmvThreads = N * 8;
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
    Eigen::VectorXd y(N);
    for (uint32_t i = 0; i < N; ++i) {
      y[i] = static_cast<double>(apPtr[i]);
    }
    return y;
  }
}

SpmvBenchmarkResult
benchmarkSpMVMetal(const Eigen::SparseMatrix<double>& scalarA,
                   const Eigen::SparseMatrix<double>& fullA,
                   int numRepetitions)
{
  SpmvBenchmarkResult result;
  uint32_t numNodes = static_cast<uint32_t>(scalarA.rows());
  uint32_t N = numNodes * 3;
  size_t nnz3D = static_cast<size_t>(scalarA.nonZeros());
  size_t nnzFull = static_cast<size_t>(fullA.nonZeros());

  result.memory3D_bytes = sizeof(uint32_t) * (numNodes + 1) + (sizeof(uint32_t) + sizeof(float)) * nnz3D;
  result.memoryFull_bytes = sizeof(uint32_t) * (N + 1) + (sizeof(uint32_t) + sizeof(float)) * nnzFull;
  result.bandwidthReductionRatio = static_cast<double>(result.memoryFull_bytes) / static_cast<double>(result.memory3D_bytes);

  MetalFBSBackend& backend = MetalFBSBackend::instance();
  if (!backend.isAvailable()) return result;

  @autoreleasepool {
    id<MTLCommandQueue> queue = backend.getQueue();
    uint32_t tg_size = 256;

    // 1. Benchmark 3D Block-SpMV
    backend.ensureMatrixBuffers(numNodes, nnz3D);
    id<MTLBuffer> bufRowPtr3D = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdx3D = backend.getBufColIdx();
    id<MTLBuffer> bufValues3D = backend.getBufValues();
    uint32_t* rMem3D = static_cast<uint32_t*>([bufRowPtr3D contents]);
    for (uint32_t i = 0; i <= numNodes; ++i) rMem3D[i] = static_cast<uint32_t>(scalarA.outerIndexPtr()[i]);
    uint32_t* cMem3D = static_cast<uint32_t*>([bufColIdx3D contents]);
    float* vMem3D = static_cast<float*>([bufValues3D contents]);
    for (size_t i = 0; i < nnz3D; ++i) {
      cMem3D[i] = static_cast<uint32_t>(scalarA.innerIndexPtr()[i]);
      vMem3D[i] = static_cast<float>(scalarA.valuePtr()[i]);
    }
    backend.ensureCgBuffers(N, (numNodes + tg_size - 1) / tg_size);
    id<MTLBuffer> bufX = backend.getBufX();
    id<MTLBuffer> bufAp = backend.getBufAp();
    float* xPtr = static_cast<float*>([bufX contents]);
    for (uint32_t i = 0; i < N; ++i) xPtr[i] = 1.0f;

    // Warmup 3D
    for (int w = 0; w < 5; ++w) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getSpmm3rhsStatelessPipe()];
      [enc setBuffer:bufRowPtr3D offset:0 atIndex:0];
      [enc setBuffer:bufColIdx3D offset:0 atIndex:1];
      [enc setBuffer:bufValues3D offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufAp offset:0 atIndex:4];
      [enc setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
      [enc dispatchThreads:MTLSizeMake(numNodes * 8, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    auto start3D = std::chrono::high_resolution_clock::now();
    id<MTLCommandBuffer> cmd3D = [queue commandBuffer];
    id<MTLComputeCommandEncoder> enc3D = [cmd3D computeCommandEncoder];
    [enc3D setComputePipelineState:backend.getSpmm3rhsStatelessPipe()];
    [enc3D setBuffer:bufRowPtr3D offset:0 atIndex:0];
    [enc3D setBuffer:bufColIdx3D offset:0 atIndex:1];
    [enc3D setBuffer:bufValues3D offset:0 atIndex:2];
    [enc3D setBytes:&numNodes length:sizeof(uint32_t) atIndex:5];
    for (int rep = 0; rep < numRepetitions; ++rep) {
      id<MTLBuffer> inBuf = (rep % 2 == 0) ? bufX : bufAp;
      id<MTLBuffer> outBuf = (rep % 2 == 0) ? bufAp : bufX;
      [enc3D setBuffer:inBuf offset:0 atIndex:3];
      [enc3D setBuffer:outBuf offset:0 atIndex:4];
      [enc3D dispatchThreads:MTLSizeMake(numNodes * 8, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc3D memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    [enc3D endEncoding];
    [cmd3D commit];
    [cmd3D waitUntilCompleted];
    auto end3D = std::chrono::high_resolution_clock::now();
    result.time3D_ms = std::chrono::duration<double, std::milli>(end3D - start3D).count();

    // 2. Benchmark Full 3Nx3N SpMV
    backend.ensureMatrixBuffers(N, nnzFull);
    id<MTLBuffer> bufRowPtrFull = backend.getBufRowPtr();
    id<MTLBuffer> bufColIdxFull = backend.getBufColIdx();
    id<MTLBuffer> bufValuesFull = backend.getBufValues();
    uint32_t* rMemFull = static_cast<uint32_t*>([bufRowPtrFull contents]);
    for (uint32_t i = 0; i <= N; ++i) rMemFull[i] = static_cast<uint32_t>(fullA.outerIndexPtr()[i]);
    uint32_t* cMemFull = static_cast<uint32_t*>([bufColIdxFull contents]);
    float* vMemFull = static_cast<float*>([bufValuesFull contents]);
    for (size_t i = 0; i < nnzFull; ++i) {
      cMemFull[i] = static_cast<uint32_t>(fullA.innerIndexPtr()[i]);
      vMemFull[i] = static_cast<float>(fullA.valuePtr()[i]);
    }
    backend.ensureCgBuffers(N, (N + tg_size - 1) / tg_size);

    // Warmup Full
    for (int w = 0; w < 5; ++w) {
      id<MTLCommandBuffer> cmd = [queue commandBuffer];
      id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
      [enc setComputePipelineState:backend.getSpmvCoopStatelessPipe()];
      [enc setBuffer:bufRowPtrFull offset:0 atIndex:0];
      [enc setBuffer:bufColIdxFull offset:0 atIndex:1];
      [enc setBuffer:bufValuesFull offset:0 atIndex:2];
      [enc setBuffer:bufX offset:0 atIndex:3];
      [enc setBuffer:bufAp offset:0 atIndex:4];
      [enc setBytes:&N length:sizeof(uint32_t) atIndex:5];
      [enc dispatchThreads:MTLSizeMake(N * 8, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [enc endEncoding];
      [cmd commit];
      [cmd waitUntilCompleted];
    }

    auto startFull = std::chrono::high_resolution_clock::now();
    id<MTLCommandBuffer> cmdFull = [queue commandBuffer];
    id<MTLComputeCommandEncoder> encFull = [cmdFull computeCommandEncoder];
    [encFull setComputePipelineState:backend.getSpmvCoopStatelessPipe()];
    [encFull setBuffer:bufRowPtrFull offset:0 atIndex:0];
    [encFull setBuffer:bufColIdxFull offset:0 atIndex:1];
    [encFull setBuffer:bufValuesFull offset:0 atIndex:2];
    [encFull setBytes:&N length:sizeof(uint32_t) atIndex:5];
    for (int rep = 0; rep < numRepetitions; ++rep) {
      id<MTLBuffer> inBuf = (rep % 2 == 0) ? bufX : bufAp;
      id<MTLBuffer> outBuf = (rep % 2 == 0) ? bufAp : bufX;
      [encFull setBuffer:inBuf offset:0 atIndex:3];
      [encFull setBuffer:outBuf offset:0 atIndex:4];
      [encFull dispatchThreads:MTLSizeMake(N * 8, 1, 1) threadsPerThreadgroup:MTLSizeMake(tg_size, 1, 1)];
      [encFull memoryBarrierWithScope:MTLBarrierScopeBuffers];
    }
    [encFull endEncoding];
    [cmdFull commit];
    [cmdFull waitUntilCompleted];
    auto endFull = std::chrono::high_resolution_clock::now();
    result.timeFull_ms = std::chrono::duration<double, std::milli>(endFull - startFull).count();

    result.speedup = result.timeFull_ms / result.time3D_ms;
    return result;
  }
}

} // namespace pylimer_tools::sim::mehp
