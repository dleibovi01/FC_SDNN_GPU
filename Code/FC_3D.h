/* 3D FC spatial differentiation schemes */

#ifndef FC_3D_H
#define FC_3D_H

#include "FC.h"
#include "FC_1D.h"
#include "Patch3D.h"
#include "SpatDiffScheme.h"
#include <string.h>
#include "printing.h"
#include <cmath>
#include <cufft.h>
#include <cudnn.h>
#include <cuda_runtime.h>
#include <cublas_v2.h>
#include "helper_cuda.h"
#include "VectorOperations.cuh"

class FC_3D : public SpatDiffScheme<Patch3D>{


FC_1D * FC_x;

FC_1D * FC_y;

FC_1D * FC_z;

int Nx;

int Ny;

int Nz;

double * contMat_y = nullptr;

double * ones_y = nullptr;

double * derCoeffs3D_y = nullptr;

double * container = nullptr;

double * buffer_x = nullptr;

double * buffer_z = nullptr;

double * buffer_cont_y = nullptr;

double * buffer_cont_y_T = nullptr;

double * buffer_bd_y = nullptr;

double * buffer_bd_y_T = nullptr;

std::complex<double> * buffer_fft_y = nullptr;

std::complex<double> * buffer_ifft_y = nullptr;

DFTI_DESCRIPTOR_HANDLE desc_handle_3D_y = NULL;

// GPU buffers and FFT plans 

double * dev_buffer_y = nullptr;
 
double * dev_buffer_y_der = nullptr;

double * dev_buffer_bd_y = nullptr;

double * dev_buffer_bd_y_T = nullptr;

// double * dev_buffer_cont_y = nullptr;

cufftDoubleComplex * dev_buffer_cont_y = nullptr;

double * dev_buffer_cont_y_T = nullptr;

double * dev_contMat_y = nullptr;

double * dev_derCoeffs3D_y = nullptr;

cufftDoubleComplex * dev_buffer_fft_y = nullptr;

cufftDoubleComplex * dev_buffer_ifft_y = nullptr;

cufftHandle dev_desc_handle_3D_y = NULL;

/** cuBLAS library context */
cublasHandle_t cublasHandle;

/** cuDNN library context */
cudnnHandle_t cudnnHandle;


public:


    FC_3D(const std::string xx, const std::string yy, const std::string zz, 
        int _Nx, int dx, int Cx, double hx, int _Ny, int dy, int Cy, double hy,
        int _Nz, int dz, int Cz, double hz, int alpha, int p);

    FC_3D(const std::string xx, const std::string yy, const std::string zz, 
        int _Nx, int dx, int Cx, double hx, int _Ny, int dy, int Cy, double hy,
        int _Nz, int dz, int Cz, double hz, double delta);

   ~FC_3D();

   FC_1D * getFCx() const {return FC_x;}

   FC_1D * getFCy() const {return FC_y;}

   FC_1D * getFCz() const {return FC_z;}

   double * getContainer() const {return container;}

   int getNx() const {return Nx;}

   int getNy() const {return Ny;}

   int getNz() const {return Nz;}

   int getCx() const {return FC_x->getC();}

   int getCy() const {return FC_y->getC();}

   int getCz() const {return FC_z->getC();}

   int getDx() const {return FC_x->getD();}

   int getDy() const {return FC_y->getD();}

   int getDz() const {return FC_z->getD();}

   double getFourPtsDblx() const {return FC_x->getFourPts_dbl();}

   double getFourPtsDbly() const {return FC_y->getFourPts_dbl();}

   double getFourPtsDblz() const {return FC_z->getFourPts_dbl();}

   double * getAQx() const {return FC_x->getAQ();}

   double * getFAQFx() const {return FC_x->getFAQF();} 

   double * getAQy() const {return FC_y->getAQ();}

   double * getFAQFy() const {return FC_y->getFAQF();} 

   double * getAQz() const {return FC_z->getAQ();}

   double * getFAQFz() const {return FC_z->getFAQF();} 

   std::complex<double>* getShiftCoeffsx() const{return FC_x->getShiftCoeffs();}

   std::complex<double>* getShiftCoeffsy() const{return FC_y->getShiftCoeffs();}

   std::complex<double>* getShiftCoeffsz() const{return FC_z->getShiftCoeffs();}

   double * getDerCoeffsx() const{return FC_x->getDerCoeffs();}

   double * getDerCoeffsy() const{return FC_y->getDerCoeffs();}

   double * getDerCoeffsz() const{return FC_z->getDerCoeffs();}

   DFTI_DESCRIPTOR_HANDLE getDescHandlex() const {return FC_x->getDescHandle();}

   DFTI_DESCRIPTOR_HANDLE getDescHandley() const {return FC_y->getDescHandle();}

   DFTI_DESCRIPTOR_HANDLE getDescHandlez() const {return FC_z->getDescHandle();}

   DFTI_DESCRIPTOR_HANDLE getDescHandleOopx() const 
      {return FC_x->getDescHandleOop();}

   DFTI_DESCRIPTOR_HANDLE getDescHandleOopy() const 
      {return FC_y->getDescHandleOop();}

   DFTI_DESCRIPTOR_HANDLE getDescHandleOopz() const 
      {return FC_z->getDescHandleOop();}

   void diff_y(const double* y, double *y_der, double h, const double* bc_d,
      const double* bc_u) const; 

   void diff_x(const double* y, double *y_der, double h, const double* bc_l,
      const double* bc_r) const;

   void diff_z(const double* y, double *y_der, double h, const double* bc_l,
      const double* bc_r) const;

   void FcontGramBlend3D_y(const double * y);

   void diff3D_y(const double* y, double *y_der);

   void FcontGramBlend3D_y_GPU(const double * y);

   float diff3D_y_GPU(const double* y, double *y_der);

private:
   
   void setContMat_y();

};

#endif