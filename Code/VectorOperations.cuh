

#ifndef CUDA_FFT_CONVOLVE_CUH
#define CUDA_FFT_CONVOLVE_CUH



#include <cufft.h>
#include "helper_cuda.h"

void cudaCallImCmpProdScaleKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        const double *raw_data,
        const cufftDoubleComplex *impulse_v,
        cufftDoubleComplex *out_data,
        const unsigned int length, 
        const unsigned int scaling);

// void cudaCallRestrictionKernel(const unsigned int blocks,
//         const unsigned int threadsPerBlock,
//         const double * dev_buffer_cont_y,
//         double * dev_buffer_y_der,
//         const unsigned int Ny, 
//         const unsigned int Cy, 
//         const unsigned int N);

void cudaCallRestrictionKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        const cufftDoubleComplex * dev_buffer_cont_y,
        double * dev_buffer_y_der,
        const unsigned int Ny, 
        const unsigned int Cy, 
        const unsigned int N);

void cudaCallCopyDblToCmpKernel(const unsigned int blocks,
        const unsigned int threadsPerBlock,
        cufftDoubleComplex * dev_buffer_cont_y,
        const double * dev_buffer_cont_y_T, const unsigned int N);


#endif
