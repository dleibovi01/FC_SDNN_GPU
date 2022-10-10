/* 3D patch */

#ifndef PATCH3D_H
#define PATCH3D_H

#include "Node.h"
#include <iostream>
#include "Patch.h"
#include <algorithm>

// #include "printing.h"

class Patch3D : public Patch{

protected:

int Nx;

int Ny;

int Nz;

double hx;

double hy;

double hz;

int fringe;

public:

    Patch3D(){};

    Patch3D(const Patch3D &patch) : Patch(patch), Nx{patch.Nx}, Ny{patch.Ny}, 
        Nz{patch.Nz}, fringe{patch.fringe} {}

    Patch3D(int _Nx, int _Ny, int _Nz, int _unknowns, double a, double b,
        double c, double d, double e, double f, bool lb, bool rb, bool db,
        bool ub, bool bb, bool fb, int _fringe);

    Patch3D & operator=(const Patch3D &patch);

    int getNx() const {return Nx;}

    int getNy() const{return Ny;}

    int getNz() const{return Nz;}

    double getHx() const {return hx;}

    double getHy() const {return hy;}

    double getHz() const {return hz;}

    double getH() const {return std::min(std::min(hx, hy), hz);}

    int getFringe() const {return fringe;}

    // void setInnerBdry(Patch2D* patch, int position, int fringe, 
    //     int overlap, int unknowns, int stage);
 
};

#endif 