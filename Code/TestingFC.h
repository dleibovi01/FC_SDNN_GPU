/* Tests for FC and FC spatial differentiation schemes routines. */

#ifndef TESTINGFC_H
#define TESTINGFC_H

#include "Mesh2D.h"
#include "Mesh3D.h"
#include "FC_2D.h"
#include "FC_3D.h"
#include "IC.h"
#include <iostream>
#include <time.h>
#include <cmath>

#include <cuda_runtime.h>

void Testing2DFC()
{
    std::cout << "Testing 2D FC differentiation" << std::endl;

    // Setting up a 2D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int patchsize_x = 20;
    int patchsize_y = 20;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    Mesh2DUniform mesh{0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y, patchsize_x, 
        patchsize_y, overlap, fringe, unknowns, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "Burgers_Square";
    IC ic{problem, 1};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    Print_Mesh2DUniform(mesh);

    // Setting up an FC_2D spatial differentiation scheme
    int dx = 5;
    int dy = 5;
    int Cx = 27;
    int Cy = 27;

    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;    
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);
    std::vector<double *> bc_d;
    std::vector<double *> bc_u;
    std::vector<double *> bc_r;
    std::vector<double *> bc_l;
    bc_d.push_back(new double[patchsize_x]);
    bc_u.push_back(new double[patchsize_x]);
    bc_l.push_back(new double[patchsize_y]);
    bc_r.push_back(new double[patchsize_y]);
    // [patchsize_x];
    // double bc_u[patchsize_x];
    // double bc_l[patchsize_y];
    // double bc_r[patchsize_y];
    memset(bc_d[0], 0.0, patchsize_x);
    memset(bc_u[0], 0.0, patchsize_x);
    memset(bc_l[0], 0.0, patchsize_y);
    memset(bc_r[0], 0.0, patchsize_y);

    FC_2D fc_2d{"DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy, hy,
        alpha0, p};

    // Differentiating the data on the patch
    auto patches = mesh.getPatches();
    // // fc_2d.diff_y(patches[0]->getFlowPtr()->getField(0), 
    // //     patches[0]->getFlowPtr()->getField(0), hy, bc_d[0], bc_u[0]);
    // fc_2d.diff_x(patches[0]->getFlowPtr()->getField(0), 
    //     patches[0]->getFlowPtr()->getField(0), hx, bc_l[0], bc_r[0]);

    // testing filtering, der and 2nd der

    std::vector<std::complex<double> * > fft_data_x;
    std::vector<std::complex<double> * > fft_data_y;
    std::vector<std::vector<int> > fft_locs;
    std::vector<int> loc;
    for(int j = 0; j < unknowns; j++)
    {
        fft_data_x.push_back(
            new std::complex<double> [mesh.getPatches()[0]->getNnodes() + Cx]);
        fft_data_y.push_back(
            new std::complex<double> [mesh.getPatches()[0]->getNnodes() + Cy]);
        loc.push_back(j);
    }
    fft_locs.push_back(loc);
    loc.clear();
    std::vector<int> filt_unknowns;
    for(int i = 0; i < unknowns; i++)
    {
        filt_unknowns.push_back(i);
    }

    fc_2d.filter_y(patches[0]->getFlowPtr(), filt_unknowns, &fft_data_y,
        fft_locs[0], hy, bc_d, bc_u);


    // Setting the values of the nodes and printing
    patches[0]->VectorFieldToNodes();
    mesh.setPatches(patches);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "After y-filtering" << std::endl;
    Print_Mesh2DUniform(mesh);

    fc_2d.filter_x(patches[0]->getFlowPtr(), filt_unknowns, &fft_data_x,
        fft_locs[0], hx, bc_l, bc_r);

    // Setting the values of the nodes and printing
    patches[0]->VectorFieldToNodes();
    mesh.setPatches(patches);
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "After y and x filtering" << std::endl;
    Print_Mesh2DUniform(mesh);
        
        

}



void Testing3DFC_dery()
{
    std::cout << "Testing 3D FC differentiation" << std::endl;

    // Setting up a 3D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int npatches_z = 1;
    int patchsize_x = 101;
    int patchsize_y = 101;
    int patchsize_z = 101;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    Mesh3DUniform mesh{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y,
        npatches_z, patchsize_x, patchsize_y, patchsize_z, overlap, fringe,
        unknowns, true, true, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "LA_3D";
    IC ic{problem, unknowns};
    ic(&mesh);
    // std::cout << "Initial condition" << std::endl;
    // printMesh3DUniform(mesh); 

    // Setting up an FC_3D spatial differentiation scheme
    int dx = 5;
    int dy = 5;
    int dz = 5;
    int Cx = 27;
    int Cy = 27;
    int Cz = 27;

    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;    
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);
    double hz = 1.0/(double(patchsize_z) - 1.0);

    FC_3D fc_3d{"DD", "DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy,
        hy, patchsize_z, dz, Cz, hz, delta};

    // Differentiating the data on the patch
    auto patches = mesh.getPatches();

    // double ux[patchsize_x * patchsize_y * patchsize_z];
    double * uy_cpu = new double[patchsize_x * patchsize_y * patchsize_z];
    double * uy_gpu = new double[patchsize_x * patchsize_y * patchsize_z];


    // Producing periodic extensions (could be parallelized as well, but didn't
    // have enough time)
    fc_3d.FcontGramBlend3D_y(patches[0]->getFlowPtr()->getField(0));

    // Use the CUDA machinery for recording time
    cudaEvent_t start_cpu, stop_cpu, start_gpu, stop_gpu;
    cudaEventCreate(&start_cpu);
    cudaEventCreate(&stop_cpu);
    cudaEventCreate(&start_gpu);
    cudaEventCreate(&stop_gpu);

    ////////////////////////////////////////////////////////////////////////////
    // Running and measuring the runtime for differentiation for 
    // the CPU-only code
    cudaEventRecord(start_cpu);

    fc_3d.diff3D_y(patches[0]->getFlowPtr()->getField(0), uy_cpu);

    // Stop timer
    cudaEventRecord(stop_cpu);
    cudaEventSynchronize(stop_cpu);

    ////////////////////////////////////////////////////////////////////////////
    // Running and measuring the runtime for differentiation for the
    // CPU/GPU code
    cudaEventRecord(start_gpu);

    float gpu_time_no_copying = 
        fc_3d.diff3D_y_GPU(patches[0]->getFlowPtr()->getField(0), uy_gpu);

    // Stop timer
    cudaEventRecord(stop_gpu);
    cudaEventSynchronize(stop_gpu);

    ////////////////////////////////////////////////////////////////////////////

    // Checkig that the answers match
    double discrepancy = 0.0;
    for(int i = 0; i < patchsize_x * patchsize_y * patchsize_z; i++)
    {
        discrepancy = std::max(discrepancy, std::abs(uy_cpu[i] - uy_gpu[i]));
    }

    std::cout << "The maximum discrepancy (Linf norm) between the CPU and GPU output is " << discrepancy << std::endl;
    std::cout << std::endl;

    ////////////////////////////////////////////////////////////////////////////
    float cpu_time_milliseconds, gpu_time_milliseconds;
    cudaEventElapsedTime(&cpu_time_milliseconds, start_cpu, stop_cpu);
    cudaEventElapsedTime(&gpu_time_milliseconds, start_gpu, stop_gpu);

    std::cout << std::endl;
    std::cout << "CPU time: " << cpu_time_milliseconds << " milliseconds" << std::endl;
    std::cout << "GPU time: " << gpu_time_milliseconds << " milliseconds" << std::endl;
    std::cout << "GPU time without copying data: " << gpu_time_no_copying << " milliseconds" << std::endl;
    std::cout << std::endl << "Speedup factor: " <<
        cpu_time_milliseconds / gpu_time_milliseconds << std::endl << std::endl;
    std::cout << std::endl << "Speedup factor no copying: " <<
        cpu_time_milliseconds / gpu_time_no_copying << std::endl << std::endl;


    delete [] uy_cpu;
    delete [] uy_gpu;

    // std::cout << "CPU solution" << std::endl;
    // Print_Mat(uy_cpu, patchsize_y, patchsize_x * patchsize_y);  
    // std::cout << std::endl;

    // std::cout << "GPU solution" << std::endl;
    // Print_Mat(uy_gpu, patchsize_y, patchsize_x * patchsize_y);  
    // std::cout << std::endl;   
}



void Testing3DFC()
{
    std::cout << "Testing 3D FC differentiation" << std::endl;

    // Setting up a 3D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int npatches_z = 1;
    int patchsize_x = 10;
    int patchsize_y = 10;
    int patchsize_z = 10;
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    Mesh3DUniform mesh{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y,
        npatches_z, patchsize_x, patchsize_y, patchsize_z, overlap, fringe,
        unknowns, true, true, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "LA_3D";
    IC ic{problem, unknowns};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    printMesh3DUniform(mesh); 

    // Setting up an FC_3D spatial differentiation scheme
    int dx = 5;
    int dy = 5;
    int dz = 5;
    int Cx = 27;
    int Cy = 27;
    int Cz = 27;

    double alpha0 = 10.0;
    double p_0 = 2.0;
    double p = 14.0;    
    double delta = 0.1;
    double hx = 1.0/(double(patchsize_x) - 1.0);
    double hy = 1.0/(double(patchsize_y) - 1.0);
    double hz = 1.0/(double(patchsize_z) - 1.0);

    FC_3D fc_3d{"DD", "DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy,
        hy, patchsize_z, dz, Cz, hz, delta};

    // Differentiating the data on the patch
    auto patches = mesh.getPatches();

    double ux[patchsize_x * patchsize_y * patchsize_z];
    double uy[patchsize_x * patchsize_y * patchsize_z];
    double uz[patchsize_x * patchsize_y * patchsize_z];

    double bc[patchsize_x*patchsize_y];
    memset(bc, 0.0, patchsize_x*patchsize_y*sizeof(double));

    std::cout << "y-dimension" << std::endl;
    fc_3d.diff_y(patches[0]->getFlowPtr()->getField(0), uy, hy, bc, bc);

    std::cout << "x-dimension" << std::endl;
    fc_3d.diff_x(patches[0]->getFlowPtr()->getField(0), ux, hx, bc, bc);

    std::cout << "z-dimension" << std::endl;
    fc_3d.diff_z(patches[0]->getFlowPtr()->getField(0), uz, hz, bc, bc);

    std::cout << "printing uy" << std::endl;
    Print_Mat(uy + patchsize_x * patchsize_y, patchsize_y, patchsize_x);
    std::cout << std::endl;
    std::cout << std::endl;

    std::cout << "printing ux" << std::endl;    
    Print_Mat(ux + patchsize_x * patchsize_y, patchsize_y, patchsize_x);
    std::cout << std::endl;
    std::cout << std::endl;    

    std::cout << "printing uz" << std::endl;
    Print_Mat(uz + patchsize_x * patchsize_y, patchsize_y, patchsize_x);     
    std::cout << std::endl;
    std::cout << std::endl;       

}

#endif