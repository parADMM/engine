//=============================================================================
//
//  paradmm-gpu.h
//
// == Authors ==
// Ning Hao
// Nate Derbinsky
//
//=============================================================================

#ifndef __PARADMM_GPU_H
#define __PARADMM_GPU_H

#include "paradmm.h"

MY_EXTERN_C __device__ double d_abs_d(double x);

MY_EXTERN_C __global__ void d_updateX(graph* G);
MY_EXTERN_C __global__ void d_updateM(graph* G);
MY_EXTERN_C __global__ void d_updateZ(graph* G);
MY_EXTERN_C __global__ void d_updateU(graph* G);
MY_EXTERN_C __global__ void d_updateU_per_edge(graph* G);
MY_EXTERN_C __global__ void d_updateN(graph* G);
MY_EXTERN_C __global__ void d_updateN_per_edge(graph* G);
MY_EXTERN_C __global__ void d_updateG_XMZUN_single_block(graph* G, int numiterations);

MY_EXTERN_C __global__ void printKernel(graph* G);

MY_EXTERN_C void copyGraphFromCPUtoGPU(graph c_graph, graph* h_d_graph);

//////////////////////////////////////////////////
//////////////////////////////////////////////////

#endif
