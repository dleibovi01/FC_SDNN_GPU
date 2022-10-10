/* Boundary conditions for 3D Linear advection equations*/

#ifndef BC_LA_3D_H
#define BC_LA_3D_H


#include "BC.h"
#include <vector>
#include <math.h>


class BC_LA_3D_wave : public BC
{
    std::vector<double> bcx_l;
    std::vector<double> bcx_r;
    std::vector<double> bcy_l;
    std::vector<double> bcy_r;
    std::vector<double> bcz_l;
    std::vector<double> bcz_r;    
    

public:

    BC_LA_3D_wave() : BC{1} {}

    std::vector<int> getEnforceableBdries(const Node* node,
        const double t) const 
    {
        std::vector<int> enforceable_bdries;
        std::vector<double> x = node->getPosition();
        if((x[0] == 0.0) || (x[1] == 0.0) || (x[2] == 0.0))
        {
            enforceable_bdries.push_back(0);
        }
        return enforceable_bdries;
    }

    std::vector<double> enforceBC(const Node* node, const double t) 
    {
        double pi = 3.1415926535897932384686;
        std::vector<double> bdry_values;
        std::vector<double> x = node->getPosition();
        // double value = std::cos(2*pi*(x[0] - t)) * std::cos(2*pi*(x[1] - t)) *
        //     std::cos(2*pi*(x[2] - t)); 
        double value = std::cos(2*pi*(x[0] - t)) * std::sin(2*pi*(x[1] - t)) *
            std::cosh(x[2] - t); 
        bdry_values.push_back(value);
        return bdry_values;
    }

    void getBCX_L(std::vector<double* > * BC_l, const Patch3D * patch,
        const double t) const 
    {
        int Nz = patch->getNz();
        int Ny = patch->getNy();
        memset((*BC_l)[0], 0.0, Nz*Ny*sizeof(double));
    }

    void getBCX_R(std::vector<double* > * BC_r, const Patch3D * patch,
        const double t) const 
    {
        int Nz = patch->getNz();
        int Ny = patch->getNy();
        memset((*BC_r)[0], 0.0, Ny*sizeof(double));
    }

    void getBCY_L(std::vector<double* > * BC_d, const Patch3D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        int Nz = patch->getNz();
        memset((*BC_d)[0], 0.0, Nz*Nx*sizeof(double));
    }

    void getBCY_R(std::vector<double* > * BC_u, const Patch3D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        int Nz = patch->getNz();
        memset((*BC_u)[0], 0.0, Nz*Nx*sizeof(double));
    }

    void getBCZ_L(std::vector<double* > * BC_d, const Patch3D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        int Ny = patch->getNy();
        memset((*BC_d)[0], 0.0, Ny*Nx*sizeof(double));
    }

    void getBCZ_R(std::vector<double* > * BC_u, const Patch3D * patch,
        const double t) const 
    {
        int Nx = patch->getNx();
        int Ny = patch->getNy();
        memset((*BC_u)[0], 0.0, Ny*Nx*sizeof(double));
    }


};




#endif