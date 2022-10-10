/* 3D Patch */

#include "Patch3D.h"
#include "printing.h"
#include "VectorField.h"

Patch3D & Patch3D::operator=(const Patch3D &patch)
{
    v = patch.getFlow();
    std::vector<Node*> n = patch.getNodes();
    Nnodes = n.size();    
    nodes.clear();
    for(int i = 0; i < Nnodes; i++)
    {
        nodes.push_back(n[i]);
    }
    phys_bdry_nodes = patch.getPhysBdryNodes();
    hx = patch.getHx();
    hy = patch.getHy();
    hz = patch.hz;
    Nx = patch.getNx();
    Ny = patch.getNy();
    Nz = patch.Nz;
    inner_nodes = patch.inner_nodes;
    return *this; 
}


Patch3D::Patch3D(int _Nx, int _Ny, int _Nz, int _unknowns, double a, double b,
    double c, double d, double e, double f, bool xlb, bool xrb, bool ybb, bool yfb,
    bool zdb, bool zub, int _fringe) : Nx{_Nx}, Ny{_Ny}, Nz{_Nz},
    fringe{_fringe}
{
    Nnodes = Nx * Ny * Nz;
    hx = (b - a) / double(Nx - 1);
    hy = (d - c) / double(Ny - 1);
    hz = (f - e) / double(Nz - 1);
    std::vector<double> position(3);
    for(int i = 0; i < Nz; i++)
    {
        for(int j = 0; j < Nx; j++)
        {
            for(int k = 0; k < Ny; k++)
            {
                nodes.push_back(new Node(_unknowns));
                nodes[Nx * Ny * i + Ny * j + k]->setIndex(Ny * j + k);
                position[0] = a + j*hx;
                position[1] = c + k*hy;
                position[2] = e + i*hz;
                nodes[Nx * Ny * i + Ny * j + k]->setPosition(position);
            }
        }
    }

    // If the boundary at x = a is a domain boundary 
    if(xlb)
    {
        for(int k = 0; k < Nz; k++)
        {
            for(int j = 0; j < Ny; j++)
            {
                phys_bdry_nodes.push_back(Nx * Ny * k + j);
            }
        }
    }

    // If the boundary at x = b is a domain boundary
    if(xrb)
    {
        for(int k = 0; k < Nz; k++)
        {
            for(int j = 0; j < Ny; j++)
            {
                phys_bdry_nodes.push_back(Nx * Ny * k + Ny * (Nx - 1) + j);
            }
        }    
    }

    // If the boundary at y = c is a domain boundary
    if(ybb)
    {
        for(int k = 0; k < Nz; k++)
        {
            for(int i = 0; i < Nx; i++)
            {
                phys_bdry_nodes.push_back(Nx * Ny * k + Ny * i);
            }
        }            
    }

    // If the boundary at y = d is a domain boundary
    if(yfb)
    {
        for(int k = 0; k < Nz; k++)
        {        
            for(int i = 0; i < Nx; i++)
            {
                phys_bdry_nodes.push_back(Nx * Ny * k + Ny * (i + 1) - 1);
            }  
        }        
    }

    // If the boundary at z = e is a domain boundary
    if(zdb)
    {
        for(int i = 0; i < Nx; i++)
        {
            for(int j = 0; j < Ny; j++)
            {
                phys_bdry_nodes.push_back(Ny * i + j);
            }
        }          
    }

    // If the boundary at z = f is a domain boundary
    if(zub)
    {
        for(int i = 0; i < Nx; i++)
        {
            for(int j = 0; j < Ny; j++)
            {
                phys_bdry_nodes.push_back(Nx * Ny * (Nz - 1) + Ny * i + j);
            }
        }          
    }    

    // Setting a vectorfield of length Nnodes for _unknowns unknowns.   
    v = VectorField{_unknowns, Nnodes};    
}