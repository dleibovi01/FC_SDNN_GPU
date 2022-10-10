   
#include <cstdio>
#include "VectorOperations.cuh"
#include <cuda_runtime.h>
#include <cufft.h>
#include <device_launch_parameters.h>
#include <thrust/device_ptr.h>
#include <thrust/fill.h>
#include "helper_cuda.h"




__global__
void
cudaImCmpProdScaleKernel(const double *raw_data, 
    const cufftDoubleComplex *impulse_v, cufftDoubleComplex *out_data, 
    int length, int scaling) 
{
    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    while(thread_index < length)
    {
        out_data[thread_index].x = 
            - raw_data[thread_index] * impulse_v[thread_index].y;
        out_data[thread_index].x /= scaling;

        out_data[thread_index].y = 
            raw_data[thread_index] * impulse_v[thread_index].x;

        out_data[thread_index].y /= scaling;

        thread_index += blockDim.x * gridDim.x;
    }
}




__global__
void 
cudaRestrictionKernel(const cufftDoubleComplex * dev_buffer_cont_y,
    double * dev_buffer_y_der, const unsigned int Ny, const unsigned int Cy,
    const unsigned int N)
{
    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    uint i;
    uint j;
    uint ind;
    uint fourPts = Ny + Cy;
    while(thread_index < N*Ny)
    {
        i = thread_index / Ny;
        j = thread_index % Ny;
        ind = i*fourPts + j;
        dev_buffer_y_der[thread_index] = dev_buffer_cont_y[ind].x;

        thread_index += blockDim.x * gridDim.x;
    }
}

__global__
void 
cudaCopyDblToCmpKernel(cufftDoubleComplex * dev_buffer_cont_y,
        const double * dev_buffer_cont_y_T, const unsigned int N)
{
    uint thread_index = blockIdx.x * blockDim.x + threadIdx.x;
    while(thread_index < N)
    {
        dev_buffer_cont_y[thread_index].x = dev_buffer_cont_y_T[thread_index];

        thread_index += blockDim.x * gridDim.x;
    }   
}





void cudaCallRestrictionKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        const cufftDoubleComplex * dev_buffer_cont_y,
        double * dev_buffer_y_der,
        const unsigned int Ny, 
        const unsigned int Cy, 
        const unsigned int N)
{
    cudaRestrictionKernel<<<blocks, threadsPerBlock>>>(dev_buffer_cont_y,
        dev_buffer_y_der, Ny, Cy, N);    
}



void cudaCallImCmpProdScaleKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        const double *raw_data,
        const cufftDoubleComplex *impulse_v,
        cufftDoubleComplex *out_data,
        const unsigned int length,
        const unsigned int scaling)
{
    cudaImCmpProdScaleKernel<<<blocks, threadsPerBlock>>>(raw_data, impulse_v, 
        out_data, length, scaling);
}

void cudaCallCopyDblToCmpKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        cufftDoubleComplex * dev_buffer_cont_y,
        const double * dev_buffer_cont_y_T, const unsigned int N)
{
    cudaCopyDblToCmpKernel<<<blocks, threadsPerBlock>>>(dev_buffer_cont_y,
        dev_buffer_cont_y_T, N);    
}




