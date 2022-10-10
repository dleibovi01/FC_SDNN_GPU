/* PDE solver */

#ifndef SOLVER_H
#define SOLVER_H

#include "SpatDiffScheme.h"
#include "TimeStepScheme.h"
#include "Mesh.h"
#include "Mesh1D.h"
#include "PDE.h"
#include "printing.h"
#include "SDNN.h"
#include <chrono>
#include <array>
#include <string>
#include <omp.h>
#include "FC_2D.h"
#include "FC_3D.h"

template<typename Meshtype, typename PDE, typename TS, typename Sp_diff, 
    typename Filter>
class Solver{
    Meshtype mesh;
    PDE pde;
    TS ts;
    const std::vector<Sp_diff> sp_diffs;
    const std::vector<Filter> filters;

    // ArtVisc artvisc;
    public:

    Solver(const Meshtype &_mesh, const PDE &_pde, TS &_ts, 
        const std::vector<Sp_diff> &_sp_diff, 
        const std::vector<Filter> &_filter) : mesh{_mesh}, pde{_pde}, ts{_ts},
        sp_diffs{_sp_diff}, filters{_filter} {}

    Solver(const Meshtype &_mesh, const PDE &_pde, TS &_ts, 
        const std::vector<Sp_diff> &_sp_diff) : mesh{_mesh}, pde{_pde}, ts{_ts},
        sp_diffs{_sp_diff} {}

    // const Meshtype& getMesh() const {return mesh;}
    const Meshtype& getMesh() const {return mesh;}

    


void solve3D(double dt)
    {
        std::cout << "Solver" << std::endl;
        double t = 0.0;
        double T = pde.getT();
        auto patches = mesh.getPatches();
        int npatches = sp_diffs.size();
        int phys_unknowns = pde.getPhysUnknowns();
        int stages = ts.getStages();
        pde.getIC()(&mesh);

        auto t1 = std::chrono::high_resolution_clock::now();
      
        
        double hx = 0.0;
        double hy = 0.0;
        double hz = 0.0;
        
        std::vector<std::vector<double* > > xbc_l(npatches);
        std::vector<std::vector<double* > > xbc_r(npatches);
        std::vector<std::vector<double* > > ybc_l(npatches);
        std::vector<std::vector<double* > > ybc_r(npatches);  
        std::vector<std::vector<double* > > zbc_l(npatches);
        std::vector<std::vector<double* > > zbc_r(npatches); 

        // std::cout << "Setting containers for bcs" << std::endl;    
        for(int i = 0; i < npatches; i++)
        {
            for(int j = 0; j < phys_unknowns; j++)
            {
                xbc_l[i].push_back(
                    new double[patches[i]->getNy() * patches[i]->getNz()]);
                xbc_r[i].push_back(
                    new double[patches[i]->getNy() * patches[i]->getNz()]);
                ybc_l[i].push_back(
                    new double[patches[i]->getNx() * patches[i]->getNz()]);
                ybc_r[i].push_back(
                    new double[patches[i]->getNx() * patches[i]->getNz()]);                                     
                zbc_l[i].push_back(
                    new double[patches[i]->getNy() * patches[i]->getNx()]);
                zbc_r[i].push_back(
                    new double[patches[i]->getNy() * patches[i]->getNx()]);   
            }
        } 
        // Time loop
        while (t < T)
        {
            // std::cout << "t = " << t << std::endl;
            patches = mesh.getPatches();      

            
            if(t + dt > T)
            {
                dt = T - t;
                t = T; 
            }
            else
            {
                t += dt;
            }           

            // Advancing
            for(int i = 0; i < npatches; i++)
            {
                pde.getBC().getBCX_L(&xbc_l[i], patches[i], t);
                pde.getBC().getBCX_R(&xbc_r[i], patches[i], t);
                pde.getBC().getBCY_L(&ybc_l[i], patches[i], t);
                pde.getBC().getBCY_R(&ybc_r[i], patches[i], t);
                pde.getBC().getBCZ_L(&zbc_l[i], patches[i], t);
                pde.getBC().getBCZ_R(&zbc_r[i], patches[i], t);

                ts.advance3D(patches[i], sp_diffs[i], pde, dt, t, xbc_l[i],
                    xbc_r[i], ybc_l[i], ybc_r[i], zbc_l[i], zbc_r[i]);  
            } 

            mesh.setPatches(patches);
            pde.getBC()(&mesh, t);
            
        }
        patches = mesh.getPatches();
        for(int i = 0; i < npatches; i++)
        {
            patches[i]->VectorFieldToNodes();
        }
        mesh.setPatches(patches);

        auto t2 = std::chrono::high_resolution_clock::now();

        // Getting number of milliseconds as an integer. //
        auto ms_int = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1); 
        // Getting number of milliseconds as a double. //
        std::chrono::duration<double, std::milli> ms_double = t2 - t1;

        std::cout << ms_int.count() << "ms\n" << std::endl;
        std::cout << ms_double.count() << "ms" << std::endl;  


    }  

    
};


void initialize1DFCSDNN(std::vector<FC_1D* > * diff_schemes, 
    std::vector<FC_1D* > * filters, const Mesh1DUniform &mesh, int npatches,
    std::string xx, int N, int d, int C, double delta, double alpha0,
    double p_0, double p);

void initialize2DFCSDNN(std::vector<FC_2D* > * diff_schemes, 
    std::vector<FC_2D* > * filters, const Mesh2DUniform &mesh, 
    int npatches_x, int npatches_y, std::string left, std::string right,
    std::string down, std::string up, int Nx, int Ny, int dx, int dy, int Cx,
    int Cy, double delta, double alpha0, double p_0, double p);

template<typename sp_diff>
void freeFCSDNN(std::vector<sp_diff* > * diff_schemes, 
    std::vector<sp_diff* > * filters)
{
    while(!diff_schemes->empty())
    {
        delete diff_schemes->back();
        diff_schemes->pop_back();
    }
    while(!filters->empty())
    {
        delete filters->back();
        filters->pop_back();
    }
}



#endif 