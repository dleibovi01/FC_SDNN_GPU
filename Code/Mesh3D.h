/* 3D Mesh */

#ifndef MESH_3D_H
#define MESH_3D_H

#include "Mesh.h"


class Mesh3DUniform : public Mesh<Patch3D>
{

protected: 

    /* number of patches in the x-direction. */
    int Npatches_x;

    /* number of patches in the y-direction. */
    int Npatches_y;

    /* number of patches in the y-direction. */
    int Npatches_z;

    /* overlap between the patches. */
    int overlap;

public:

    /* Constructor */   
    Mesh3DUniform(double a, double b, double c, double d, double e, double f,
        int npatches_x, int npatches_y, int npatches_z, int patchsize_x,
        int patchsize_y, int patchsize_z, int _overlap, int fringe,
        int unknowns, bool xlb, bool xrb, bool ybb, bool yfb, bool zb,
        bool zub);  

    int getOverlap() const {return overlap;}

    int getNpatches_x() const {return Npatches_x;}

    int getNpatches_y() const {return Npatches_y;}

    int getNpatches_z() const {return Npatches_y;}

    // void setIntraPatchBC(int unknowns, int stage);
};


#endif