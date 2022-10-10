/* 3D Mesh */

#include "Mesh3D.h"
#include "printing.h"

Mesh3DUniform::Mesh3DUniform(double a, double b, double c, double d, double e,
    double f, int npatches_x, int npatches_y, int npatches_z, int patchsize_x,
    int patchsize_y, int patchsize_z, int _overlap, int fringe, int unknowns,
    bool xlb, bool xrb, bool ybb, bool yfb, bool zdb, bool zub) : 
    overlap{_overlap}, Npatches_x{npatches_x}, Npatches_y{npatches_y},
    Npatches_z{npatches_z}  
{
    int patchsize = patchsize_x * patchsize_y * patchsize_z;
    int Nx = npatches_x * patchsize_x - (npatches_x - 1) * overlap;
    int Ny = npatches_y * patchsize_y - (npatches_y - 1) * overlap;
    int Nz = npatches_z * patchsize_z - (npatches_z - 1) * overlap;
    double hx = (b - a) / (double (Nx - 1));
    double hy = (d - c) / (double (Ny - 1));
    double hz = (f - e) / (double (Nz - 1));


    bool xlb_patch = false;
    bool xrb_patch = false;
    bool ybb_patch = false;
    bool yfb_patch = false;
    bool zdb_patch = false;
    bool zub_patch = false;

    // std::vector<int> inner_nodes;
    for(int k = 0; k < npatches_z; k++)
    {
        if(k == 0)
            zdb_patch = zdb;
        if(k == npatches_z - 1)
            zub_patch = zub;
        for(int i = 0; i < npatches_x; i++)
        {
            if(i == 0)
                xlb_patch = xlb;               
            if(i == npatches_x - 1)
                xrb_patch = xrb;
            for(int j = 0; j < npatches_y; j++)
            {
                if(j == 0)
                    ybb_patch = ybb;
                if(j == npatches_y - 1)
                    yfb_patch = yfb;
                    patches.push_back(new Patch3D(patchsize_x, patchsize_y, 
                        patchsize_z, unknowns, 
                        a + double(i*(patchsize_x - overlap))*hx, 
                        a + double(i*(patchsize_x - overlap) + 
                        (patchsize_x - 1))*hx,
                        c + double(j*(patchsize_y - overlap))*hy,
                        c + double(j*(patchsize_y - overlap) + 
                        (patchsize_y - 1))*hy,
                        e + double(k*(patchsize_z - overlap))*hz,
                        e + double(k*(patchsize_z - overlap) + 
                        (patchsize_z - 1))*hz,
                        xlb_patch, xrb_patch, ybb_patch, yfb_patch, zdb_patch,
                        zub_patch, fringe));            
                    ybb_patch = false;
                    yfb_patch = false;
            }
            xlb_patch = false;
            xrb_patch = false; 
        }
        zdb_patch = false;
        zub_patch = false; 
    }
}
