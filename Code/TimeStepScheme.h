/* Time stepping scheme */

#ifndef TIMESTEPSCHEME_H
#define TIMESTEPSCHEME_H

#include "VectorField.h"
#include "Mesh.h"
#include "SpatDiffScheme.h"
#include "FC_1D.h"
#define MKL_Complex16 std::complex<double>
#define MKL_Complex8 std::complex<float>
#include "mkl.h"
#include <algorithm> 
#include "printing.h"
#include "PDE.h"
#include <array>

#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "helper_cuda.h"



class Fwd_Euler{

double dt;
int stages = 1;

public:

double getStages() {return stages;}
Fwd_Euler(double _dt) : dt{_dt} {}

template<typename Sp_Diff>
VectorField advance(const VectorField & v, const Sp_Diff &sp_diff)
{
    VectorField u = v - dt*sp_diff.diff(v);
    return u;
}

};


class SSPRK_4{

std::vector<int> stage0;
std::vector<int> stage1;
std::vector<int> stage2;
std::vector<int> stage3;
std::vector<int> stage4;
std::vector<std::vector<int> > extractions1;
std::vector<std::vector<int> > extractions2;
std::vector<std::vector<int> > extractions3;
std::vector<std::vector<int> > extractions4;
std::vector<std::vector<int> > extractions5;
std::vector<double> coeffs1;
std::vector<double> coeffs2;
std::vector<double> coeffs3;
std::vector<double> coeffs4;
std::vector<double> coeffs5;
std::vector<double> coeffs5_bis;
int stages;
int unknowns;
VectorField flux;
VectorField flux3;

double* device_data = nullptr;
double* device_flux = nullptr;

static constexpr double a11 = 0.391752226571890;

static constexpr double a21 = 0.444370493651235;
static constexpr double a22 = 0.555629506348765;
static constexpr double a23 = 0.368410593050371;

static constexpr double a31 = 0.620101851488403;
static constexpr double a32 = 0.379898148511597;
static constexpr double a33 = 0.251891774271694;

static constexpr double a41 = 0.178079954393132;
static constexpr double a42 = 0.821920045606868;
static constexpr double a43 = 0.544974750228521;

static constexpr double a51 = 0.517231671970585;
static constexpr double a52 = 0.096059710526147;
static constexpr double a53 = 0.063692468666290;
static constexpr double a54 = 0.386708617503269;
static constexpr double a55 = 0.226007483236906;

/** cuBLAS library context */
cublasHandle_t cublasHandle;

public:

    SSPRK_4(int _unknowns, int _N) : unknowns{_unknowns}, flux{_unknowns, _N}, 
        flux3{_unknowns, _N}
    {
        for(int i = 0; i < unknowns; i++)
        {
            stage0.push_back(i);
            stage1.push_back(i + unknowns);
            stage2.push_back(i + 2*unknowns);
            stage3.push_back(i + 3*unknowns);
            stage4.push_back(i + 4*unknowns);
        }   
        stages = 5;   
        extractions1.push_back(stage0); 
        extractions2.push_back(stage0);
        extractions2.push_back(stage1);
        extractions3.push_back(stage0);
        extractions3.push_back(stage2);
        extractions4.push_back(stage0);
        extractions4.push_back(stage3);
        extractions5.push_back(stage2);
        extractions5.push_back(stage3);
        extractions5.push_back(stage4);

        coeffs1 = {1.0, 0.0};
        coeffs2 = {a21, a22, 0.0};
        coeffs3 = {a31, a32, 0.0};
        coeffs4 = {a41, a42, 0.0};
        coeffs5 = {a51, a52, a54, 0.0};
        coeffs5_bis = {1, 0.0};

        CUDA_CALL(cudaMalloc((void **) &device_data, 
            stages * _unknowns * _N * sizeof(double)));
        CUDA_CALL(cudaMalloc((void **) &device_flux, 
            2 * _unknowns * _N * sizeof(double)));
        CUBLAS_CALL( cublasCreate(&cublasHandle) );
    }

    ~SSPRK_4()
    {
        // CUDA_CALL( cudaFree(device_data) );
        // CUDA_CALL( cudaFree(device_flux) );
    }
    int getStages() const {return stages;}


    ////////////////////////////////////////////////////////////////////////////
    // 3D routines
    ////////////////////////////////////////////////////////////////////////////

    template<typename Sp_Diff, typename PDE>
    void advance3D(Patch3D* patch, const Sp_Diff &sp_diff, const PDE &pde,
        const double dt, const double t,
        const std::vector<double *> &bcx_l, const std::vector<double *> &bcx_r,
        const std::vector<double *> &bcy_l, const std::vector<double *> &bcy_r,
        const std::vector<double *> &bcz_l, const std::vector<double *> &bcz_r)

    {

        auto v = patch->getFlowPtr();
        const int N = v->getLength();
        const double hx = patch->getHx();
        const double hy = patch->getHy();
        const double hz = patch->getHz();

        // 1st stage
        coeffs1[1] = -a11*dt;
        double coeffF = -a11*dt; 
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 0, hx, hy, hz,
            t, bcx_l, bcx_r, bcy_l, bcy_r, bcz_l, bcz_r);

        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + i*N, v->getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(device_flux + i*N, flux.getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }

        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF, device_flux,
            1, device_data, 1));

        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL( cudaMemcpy(v->getField(unknowns + i), device_data + i*N,
                N * sizeof(double), cudaMemcpyDeviceToHost) );
        }


        // linComb(v, extractions1, stage1, coeffs1, flux, 
        //     sp_diff->getContainer());




        // 2nd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 1, hx, hy, hz,
            t, bcx_l, bcx_r, bcy_l, bcy_r, bcz_l, bcz_r);
        coeffs2[2] = -a23*dt;
        
        // linComb(v, extractions2, stage2, coeffs2, flux, 
        //     sp_diff->getContainer());

        coeffF = - a23*dt;
        for(int i = 0; i < 2*unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + i*N, v->getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_flux + i*N, flux.getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        CUBLAS_CALL(cublasDscal(cublasHandle, N*unknowns, &a21, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &a22, 
            device_data + N*unknowns, 1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF, device_flux,
            1, device_data, 1));
            
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL( cudaMemcpy(v->getField(2*unknowns + i), device_data + i*N,
                N * sizeof(double), cudaMemcpyDeviceToHost) );
        }




        // 3rd stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 2, hx, hy, hz,
            t, bcx_l, bcx_r, bcy_l, bcy_r, bcz_l, bcz_r);
        coeffs3[2] = -a33*dt;
        // linComb(v, extractions3, stage3, coeffs3, flux, 
        //     sp_diff->getContainer());
        coeffF = - a33*dt;
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + i*N, v->getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + (unknowns + i)*N,
                v->getField(i + 2*unknowns), 
                N*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(device_flux + i*N, flux.getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        CUBLAS_CALL(cublasDscal(cublasHandle, N*unknowns, &a31, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &a32, 
            device_data + N*unknowns, 1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF, device_flux,
            1, device_data, 1));
            
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL( cudaMemcpy(v->getField(3*unknowns + i), device_data + i*N,
                N * sizeof(double), cudaMemcpyDeviceToHost) );
        }

        // 4th stage
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 3, hx, hy, hz,
            t, bcx_l, bcx_r, bcy_l, bcy_r, bcz_l, bcz_r);
        coeffs4[2] = -a43*dt;
        // linComb(v, extractions4, stage4, coeffs4, flux, 
        //     sp_diff->getContainer());
        coeffF = - a43*dt;
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + i*N, v->getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + (unknowns + i)*N,
                v->getField(i + 3*unknowns), 
                N*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(device_flux + i*N, flux.getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        CUBLAS_CALL(cublasDscal(cublasHandle, N*unknowns, &a41, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &a42, 
            device_data + N*unknowns, 1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF, device_flux,
            1, device_data, 1));
            
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL( cudaMemcpy(v->getField(4*unknowns + i), device_data + i*N,
                N * sizeof(double), cudaMemcpyDeviceToHost) );
        }

        // Stepping   
        flux3 = flux;
        pde.Cons_to_der_flux(*v, &flux, sp_diff, stages, 4, hx, hy, hz,
            t, bcx_l, bcx_r, bcy_l, bcy_r, bcz_l, bcz_r);
        std::vector<VectorField*> fluxes5 = {&flux3, &flux};
        coeffs5[3] = -a53*dt;
        coeffs5_bis[1] = -a55*dt;
        // linComb(v, extractions5, stage0, coeffs5, flux3, 
        //     sp_diff->getContainer());   
        // linComb(v, extractions1, stage0, coeffs5_bis, flux, 
        //     sp_diff->getContainer()); 
        coeffF = - a53*dt;
        double coeffF2 = -a55*dt;
        for(int i = 0; i < 3*unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_data + i*N, v->getField(2*unknowns + i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
        }
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL(cudaMemcpy(device_flux + i*N, flux3.getField(i), 
                N*sizeof(double), cudaMemcpyHostToDevice));
            CUDA_CALL(cudaMemcpy(device_flux + (unknowns + i)*N, 
                flux.getField(i), N*sizeof(double), cudaMemcpyHostToDevice));
        }
        CUBLAS_CALL(cublasDscal(cublasHandle, N*unknowns, &a51, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &a52, 
            device_data + N*unknowns, 1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF, device_flux,
            1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &a54, 
            device_data + 2*unknowns*N, 1, device_data, 1));
        CUBLAS_CALL(cublasDaxpy(cublasHandle, N*unknowns, &coeffF2, 
            device_flux + unknowns*N, 1, device_data, 1));
            
        for(int i = 0; i < unknowns; i++)
        {
            CUDA_CALL( cudaMemcpy(v->getField(i), device_data + i*N,
                N * sizeof(double), cudaMemcpyDeviceToHost) );
        }
    } 
    
private:

    // void linComb1()


};




#endif 