/* 3D FC spatial differentiation schemes */

// #include <complex>
// #define MKL_Complex16 std::complex<double>
// #define MKL_Complex8 std::complex<float>
// #include "mkl.h"
#include "FC_3D.h"


FC_3D::FC_3D(const std::string xx, const std::string yy, const std::string zz, 
        int _Nx, int dx, int Cx, double hx, int _Ny, int dy, int Cy, double hy,
        int _Nz, int dz, int Cz, double hz, int alpha, int p) :
        Nx{_Nx}, Ny{_Ny}, Nz{_Nz}
{
    if(xx.compare("DD") == 0)
        FC_x = new FC_1D_DD(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("DN") == 0)
        FC_x = new FC_1D_DN(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("ND") == 0)
        FC_x = new FC_1D_ND(Nx, dx, Cx, hx, alpha, p);
    else if(xx.compare("NN") == 0)
        FC_x = new FC_1D_NN(Nx, dx, Cx, hx, alpha, p);

    if(yy.compare("DD") == 0)
        FC_y = new FC_1D_DD(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("DN") == 0)
        FC_y = new FC_1D_DN(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("ND") == 0)
        FC_y = new FC_1D_ND(Ny, dy, Cy, hy, alpha, p);
    else if(yy.compare("NN") == 0)
        FC_y = new FC_1D_NN(Ny, dy, Cy, hy, alpha, p);    

    if(zz.compare("DD") == 0)
        FC_z = new FC_1D_DD(Nz, dz, Cz, hz, alpha, p);
    else if(zz.compare("DN") == 0)
        FC_z = new FC_1D_DN(Nz, dz, Cz, hz, alpha, p);
    else if(zz.compare("ND") == 0)
        FC_z = new FC_1D_ND(Nz, dz, Cz, hz, alpha, p);
    else if(zz.compare("NN") == 0)
        FC_z = new FC_1D_NN(Nz, dz, Cz, hz, alpha, p); 

    // Setting buffers
    container = new double[_Nx * _Ny * _Nz];
    buffer_x = new double[_Nx * _Ny * _Nz];
    buffer_z = new double[_Nx * _Ny * _Nz];
    buffer_cont_y = new double[_Nx * (_Ny + Cy) * _Nz];
    buffer_cont_y_T = new double[_Nx * (_Ny + Cy) * _Nz];
    buffer_bd_y = new double[_Nx * (2 * dy) * _Nz];
    buffer_bd_y_T = new double[_Nx * (2 * dy) * _Nz];
    contMat_y = new double[2 * dy * Cy];
    buffer_fft_y = new std::complex<double>[_Nx * (_Ny + Cy) * _Nz];
    buffer_ifft_y = new std::complex<double>[_Nx * (_Ny + Cy) * _Nz];

    // Setting derCoeffs for 3D differentiation
    ones_y = new double[Nx*Nz];
    derCoeffs3D_y = new double[(Ny + Cy)*Nx*Nz];
    std::fill(ones_y, ones_y + Nx*Nz, 1.0);
    cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, Ny + Cy, Nx*Nz, 1,
        1.0, getDerCoeffsy(), Ny + Cy, ones_y, Nx*Nz, 0.0, derCoeffs3D_y, 
        Ny + Cy);    
    

    // Setting descriptor for batch differentiation
    int status;
    status = DftiCreateDescriptor(&desc_handle_3D_y,
        DFTI_DOUBLE, DFTI_COMPLEX, 1, Ny + Cy);
    status = DftiSetValue(desc_handle_3D_y, DFTI_NUMBER_OF_TRANSFORMS, Nx*Ny);
    status = DftiSetValue(desc_handle_3D_y, DFTI_INPUT_DISTANCE, Ny + Cy);
    status = DftiCommitDescriptor(desc_handle_3D_y);
    setContMat_y();


}

FC_3D::FC_3D(const std::string xx, const std::string yy, const std::string zz, 
        int _Nx, int dx, int Cx, double hx, int _Ny, int dy, int Cy, double hy,
        int _Nz, int dz, int Cz, double hz, double delta) :
        Nx{_Nx}, Ny{_Ny}, Nz{_Nz}
{
    if(xx.compare("DD") == 0)
        FC_x = new FC_1D_DD(Nx, dx, Cx, hx, delta);
    else if(xx.compare("DN") == 0)
        FC_x = new FC_1D_DN(Nx, dx, Cx, hx, delta);
    else if(xx.compare("ND") == 0)
        FC_x = new FC_1D_ND(Nx, dx, Cx, hx, delta);
    else if(xx.compare("NN") == 0)
        FC_x = new FC_1D_NN(Nx, dx, Cx, hx, delta);

    if(yy.compare("DD") == 0)
        FC_y = new FC_1D_DD(Ny, dy, Cy, hy, delta);
    else if(yy.compare("DN") == 0)
        FC_y = new FC_1D_DN(Ny, dy, Cy, hy, delta);
    else if(yy.compare("ND") == 0)
        FC_y = new FC_1D_ND(Ny, dy, Cy, hy, delta);
    else if(yy.compare("NN") == 0)
        FC_y = new FC_1D_NN(Ny, dy, Cy, hy, delta);   

    if(zz.compare("DD") == 0)
        FC_z = new FC_1D_DD(Nz, dz, Cz, hz, delta);
    else if(zz.compare("DN") == 0)
        FC_z = new FC_1D_DN(Nz, dz, Cz, hz, delta);
    else if(zz.compare("ND") == 0)
        FC_z = new FC_1D_ND(Nz, dz, Cz, hz, delta);
    else if(zz.compare("NN") == 0)
        FC_z = new FC_1D_NN(Nz, dz, Cz, hz, delta);  

    container = new double[_Nx * _Ny * _Nz]; 
    buffer_x = new double[_Nx * _Ny * _Nz]; 
    buffer_z = new double[_Nx * _Ny * _Nz]; 
    buffer_cont_y = new double[_Nx * (_Ny + Cy) * _Nz];
    buffer_cont_y_T = new double[_Nx * (_Ny + Cy) * _Nz];
    buffer_bd_y = new double[_Nx * (2 * dy) * _Nz];
    buffer_bd_y_T = new double[_Nx * (2 * dy) * _Nz];
    contMat_y = new double[2 * dy * Cy];
    buffer_fft_y = new std::complex<double>[_Nx * (_Ny + Cy) * _Nz];
    buffer_ifft_y = new std::complex<double>[_Nx * (_Ny + Cy) * _Nz];

    // Setting derCoeffs for 3D differentiation
    ones_y = new double[Nx*Nz];
    derCoeffs3D_y = new double[(Ny + Cy)*Nx*Nz];
    std::fill(ones_y, ones_y + Nx*Nz, 1.0);
    cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, Ny + Cy, Nx*Nz, 1,
        1.0, getDerCoeffsy(), Ny + Cy, ones_y, Nx*Nz, 0.0, derCoeffs3D_y, 
        Ny + Cy);  

    // Setting descriptor for batch differentiation
    int status;
    status = DftiCreateDescriptor(&desc_handle_3D_y,
        DFTI_DOUBLE, DFTI_COMPLEX, 1, Ny + Cy);
    status = DftiSetValue(desc_handle_3D_y, DFTI_NUMBER_OF_TRANSFORMS, Nx*Ny);
    status = DftiSetValue(desc_handle_3D_y, DFTI_INPUT_DISTANCE, Ny + Cy);
    status = DftiCommitDescriptor(desc_handle_3D_y);

    setContMat_y();

    // allocating the GPU buffers
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_y, Nx*Ny*Nz * sizeof(double)));  
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_y_der, Nx*Ny*Nz*sizeof(double)));   
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_bd_y, 
        2*dy*Nx*Nz*sizeof(double)));   
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_bd_y_T, 
        2*dy*Nx*Nz*sizeof(double)));  
    // CUDA_CALL(cudaMalloc((void **) &dev_buffer_cont_y, 
    //     Nx*(Ny + Cy)*Nz * sizeof(double)));  
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_cont_y, 
        Nx*(Ny + Cy)*Nz * sizeof(cufftDoubleComplex)));
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_cont_y_T, 
        Nx*(Ny + Cy)*Nz * sizeof(double)));
    CUDA_CALL(cudaMalloc((void **) &dev_contMat_y, 2*dy*Cy*sizeof(double))); 
    CUDA_CALL(cudaMalloc((void **) &dev_derCoeffs3D_y,
        (Ny + Cy)*Nx*Nz*sizeof(double)));  
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_fft_y,
        (Ny + Cy)*Nx*Nz*sizeof(cufftDoubleComplex))); 
    CUDA_CALL(cudaMalloc((void **) &dev_buffer_ifft_y,
        (Ny + Cy)*Nx*Nz*sizeof(cufftDoubleComplex)));         

    // cufftPlan1d(&dev_desc_handle_3D_y, Ny + Cy, CUFFT_D2Z, Nx*Nz);
    cufftPlan1d(&dev_desc_handle_3D_y, Ny + Cy, CUFFT_Z2Z, Nx*Nz);


    // cufftPlanMany(&dev_desc_handle_3D_y, Ny + Cy, CUFFT_D2Z, Nx*Nz);

    CUBLAS_CALL(cublasCreate(&cublasHandle));
    CUDNN_CALL( cudnnCreate(&cudnnHandle) );

    // setting some GPU buffers
    cudaMemcpy(dev_derCoeffs3D_y, derCoeffs3D_y, (Ny + Cy)*Nx*Nz*sizeof(double), 
        cudaMemcpyHostToDevice);
}

FC_3D::~FC_3D()
{
    delete FC_x;
    delete FC_y;
    delete FC_z;
    delete [] container;
    delete [] buffer_x;
    delete [] buffer_z;
    delete [] buffer_cont_y;
    delete [] buffer_cont_y_T;
    delete [] buffer_bd_y;
    delete [] buffer_bd_y_T;
    delete [] contMat_y;
    delete [] buffer_fft_y;
    delete [] buffer_ifft_y;
    delete [] derCoeffs3D_y;
    delete [] ones_y;

    DftiFreeDescriptor(&desc_handle_3D_y);

    // clearing GPU buffers
    cudaFree(dev_buffer_y);
    cudaFree(dev_buffer_y_der);
    cudaFree(dev_buffer_bd_y);
    cudaFree(dev_buffer_bd_y_T);
    cudaFree(dev_buffer_cont_y);
    cudaFree(dev_buffer_cont_y_T);
    cudaFree(dev_contMat_y);
    cudaFree(dev_derCoeffs3D_y);
    cudaFree(dev_buffer_fft_y);
    cudaFree(dev_buffer_ifft_y);
    cufftDestroy(dev_desc_handle_3D_y);
    CUBLAS_CALL( cublasDestroy(cublasHandle) );
    CUDNN_CALL( cudnnDestroy(cudnnHandle) );
}


void FC_3D::diff_y(const double* y, double *y_der, double h, const double* bc_d,
    const double* bc_u) const
{
    // This loop can be parallized
    for(int j = 0; j < Nx*Nz; j++)
    {
        FC_y->diff(y + j*Ny, y_der + j*Ny, h, bc_d[j], bc_u[j]);
        // produce the continuations first, then batch ffts
    }
}

void FC_3D::diff_x(const double* y, double *y_der, double h, const double* bc_d,
    const double* bc_u) const
{
    // The following transpositions and loop can be parallized
    mkl_domatcopy('C', 'T', Ny, Nx*Nz, 1.0, y, Ny, buffer_x, Nx*Nz);
    for(int j = 0; j < Ny*Nz; j++)
    {
        FC_x->diff(buffer_x + j*Nx, buffer_x + j*Nx, h, bc_d[j], bc_u[j]);
    }
    mkl_domatcopy('C', 'T', Nx*Nz, Ny, 1.0, buffer_x, Nx*Nz, y_der, Ny);

}

void FC_3D::diff_z(const double* y, double *y_der, double h, const double* bc_d,
    const double* bc_u) const
{
    // The following transpositions and loop can be parallized
    mkl_domatcopy('C', 'T', Nx*Ny, Nz, 1.0, y, Nx*Ny, buffer_z, Nz);
    for(int j = 0; j < Ny*Nx; j++)
    {
        FC_z->diff(buffer_z + j*Nz, buffer_z + j*Nz, h, bc_d[j], bc_u[j]);
    }
    mkl_domatcopy('C', 'T', Nz, Nx*Ny, 1.0, buffer_z, Nz, y_der, Nx*Ny);
}



void FC_3D::FcontGramBlend3D_y(const double* y)
{
    int N = Nx*Ny*Nz;
    int dy = getDy();
    int Cy = getCy();

    // Transpose y into buffer_cont_y_T)
    mkl_domatcopy('C', 'T', Ny, Nx*Nz, 1.0, y, Ny, buffer_cont_y_T, Nx*Nz);

    // // copy the boundary elements into the buffer_bd_y
    std::copy(buffer_cont_y_T, buffer_cont_y_T + dy*Nx*Nz, buffer_bd_y);
    std::copy(buffer_cont_y_T + N - dy*Nx*Nz, buffer_cont_y_T + N,
        buffer_bd_y + dy*Nx*Nz);

    // // Produce the continuations by multiplying  buffer_bd_y by contMat_y
    cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, Nx*Nz, Cy, 2*dy,
        1.0, buffer_bd_y, Nx*Nz, contMat_y, 2*dy, 0.0, buffer_cont_y_T + N, 
        Nx*Nz); 

    // // Transpose into buffer_cont_y
    mkl_domatcopy('C', 'T', Nx*Nz, Ny + Cy, 1.0, buffer_cont_y_T, Nx*Nz, 
        buffer_cont_y, Ny + Cy);
}


void FC_3D::FcontGramBlend3D_y_GPU(const double* y)
{
    int N = Nx*Ny*Nz;
    int dy = getDy();
    int Cy = getCy();

    // Transpose y into buffer_cont_y_T
    mkl_domatcopy('C', 'T', Ny, Nx*Nz, 1.0, y, Ny, buffer_cont_y_T, Nx*Nz);

    // copy the boundary elements into the buffer_bd_y
    std::copy(buffer_cont_y_T, buffer_cont_y_T + dy*Nx*Nz, buffer_bd_y);
    std::copy(buffer_cont_y_T + N - dy*Nx*Nz, buffer_cont_y_T + N,
        buffer_bd_y + dy*Nx*Nz);

    // Produce the continuations by multiplying  buffer_bd_y by contMat_y
    cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, Nx*Nz, Cy, 2*dy,
        1.0, buffer_bd_y, Nx*Nz, contMat_y, 2*dy, 0.0, buffer_cont_y_T + N, 
        Nx*Nz); 

    // Transpose into buffer_cont_y
    mkl_domatcopy('C', 'T', Nx*Nz, Ny + Cy, 1.0, buffer_cont_y_T, Nx*Nz, 
        buffer_cont_y, Ny + Cy);
}





void FC_3D::diff3D_y(const double* y, double *y_der)
{
    int Cy = getCy();
    int status;
    double fourPts = getFourPtsDbly();
    double * der_coeffs = getDerCoeffsy();
    // FcontGramBlend3D_y(y);

    std::copy(buffer_cont_y, buffer_cont_y + Nx*Nz*(Ny + Cy), buffer_fft_y); 

    status = DftiComputeForward(desc_handle_3D_y, buffer_fft_y);

    for(int j = 0; j < (Ny + Cy)*Nx*Nz; j++)
        buffer_fft_y[j] = buffer_fft_y[j] / fourPts;
    VectorMulImCmp((Ny + Cy)*Nx*Nz, derCoeffs3D_y, buffer_fft_y, buffer_ifft_y);
    
    status = DftiComputeBackward(desc_handle_3D_y, buffer_ifft_y);
    for(int i = 0; i < Nx*Nz; i++)
    {
        for(int j = 0; j < Ny; j++)
        {
            y_der[i*Ny + j] = buffer_ifft_y[i*(Ny + Cy) + j].real();
        }
    }
}

float FC_3D::diff3D_y_GPU(const double* y, double *y_der)
{
    // std::cout <<"In diff3D_y" << std::endl;

    float gpu_time_nc;
    // Use the CUDA machinery for recording time
    cudaEvent_t start_gpu_nc, stop_gpu_nc;
    cudaEventCreate(&start_gpu_nc);
    cudaEventCreate(&stop_gpu_nc);

    

    int status;
    int Cy = getCy();
    double fourPts = getFourPtsDbly();
    const uint threads_per_block = 1024;
    const uint max_block_count = 65535; 
    uint block_count = std::min(max_block_count, 
        ((Ny + Cy)*Nx*Nz / (threads_per_block)));

    // FcontGramBlend3D_y_GPU(y);



    // Copy data from host to device
    cudaMemcpy(dev_buffer_cont_y_T, buffer_cont_y, 
        Nx*Nz*(Ny + Cy)*sizeof(double), cudaMemcpyHostToDevice);

    // Start the timer
    cudaEventRecord(start_gpu_nc);

    cudaCallCopyDblToCmpKernel(block_count, threads_per_block, 
        dev_buffer_cont_y, dev_buffer_cont_y_T, Nx*Nz*(Ny + Cy));
    cufftExecZ2Z(dev_desc_handle_3D_y, dev_buffer_cont_y, dev_buffer_fft_y, 
        CUFFT_FORWARD);

    cudaCallImCmpProdScaleKernel(block_count, threads_per_block, 
        dev_derCoeffs3D_y, dev_buffer_fft_y, dev_buffer_ifft_y, (Ny + Cy)*Nx*Nz,
        Ny + Cy);

    cufftExecZ2Z(dev_desc_handle_3D_y, dev_buffer_ifft_y, dev_buffer_cont_y, 
        CUFFT_INVERSE);


    cudaCallRestrictionKernel(block_count, threads_per_block, dev_buffer_cont_y,
        dev_buffer_y_der, Ny, Cy, Nx*Nz);


    // Stop timer
    cudaEventRecord(stop_gpu_nc);
    cudaEventSynchronize(stop_gpu_nc);
    cudaEventElapsedTime(&gpu_time_nc, start_gpu_nc, stop_gpu_nc);


    cudaMemcpy(y_der, dev_buffer_y_der, Nx*Nz*Ny*sizeof(double), 
        cudaMemcpyDeviceToHost);


    return gpu_time_nc;
    
}





void FC_3D::setContMat_y()
{
    int dy = getDy();
    int Cy = getCy();
    std::copy(getFAQFy(), getFAQFy() + dy*Cy, contMat_y);
    std::copy(getAQy(), getAQy() + dy*Cy, contMat_y  + dy*Cy);
    mkl_dimatcopy('C', 'T', Cy, 2*dy, 1.0, contMat_y, Cy, 2*dy);
}
