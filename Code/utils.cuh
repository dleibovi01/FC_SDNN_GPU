/**
 * Header file for CUDA-implemented utilities needed by the neural net
 * @author Aadyot Bhatngar
 * @date April 22, 2018
 */

#pragma once

#include <cuda_runtime.h>

template<typename T> void cudaMemsetType(T *dev_ptr, T val, int n_vals);

float CrossEntropyLoss(float *pred_Y, float *true_Y,
    int n, int c, int h, int w);

float SoftThresholdAccuracy(float *pred_Y, float *true_Y,
    int n, int c, int h, int w);

// void cudaCallImCmpProdScaleKernel(const unsigned int blocks,
//         const unsigned int threadsPerBlock,
//         const double *raw_data,
//         const cufftDoubleComplex *impulse_v,
//         cufftDoubleComplex *out_data,
//         const unsigned int length, 
//         const unsigned int scaling);

// void cudaCallRestrictionKernel(const unsigned int blocks,
//         const unsigned int threadsPerBlock,
//         const double * dev_buffer_cont_y,
//         double * dev_buffer_y_der,
//         const unsigned int Ny, 
//         const unsigned int Cy, 
//         const unsigned int N);

__global__ void CrossEntropyKernel(float *pred_Y, float *true_Y, float *loss,
    int n, int c, int h, int w);

__global__ void SoftThresholdAccKernel(float *pred_Y, float *true_Y, float *acc,
    int n, int c, int h, int w);


