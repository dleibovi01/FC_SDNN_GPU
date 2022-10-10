
#include "TestingSolver.h"
#include "TimeStepScheme.h"
#include "Mesh2D.h"
#include "Mesh3D.h"
#include "FC_2D.h"
#include "FC_3D.h"
#include "Solver.h"
#include "IC.h"
#include <iostream>
#include "BC_Euler_2D.h"
#include "BC_LA_3D.h"
#include "SDNN.h"
#include "PDE.h"
#include "printing.h"


void Testing3DFCSolver()
{

    // Setting up a 3D mesh
    int npatches_x = 1;
    int npatches_y = 1;
    int npatches_z = 1;
    int patchsize_x = 101; //101 for the big test
    int patchsize_y = 101; //101 for the big test
    int patchsize_z = 101; //101 for the big test
    int overlap = 4;
    int fringe = 2;
    int unknowns = 1;
    int stages = 5;
    int mesh_unknowns = unknowns*stages + 3;
    Mesh3DUniform mesh{0.0, 1.0, 0.0, 1.0, 0.0, 1.0, npatches_x, npatches_y,
        npatches_z, patchsize_x, patchsize_y, patchsize_z, overlap, fringe,
        mesh_unknowns, true, true, true, true, true, true}; 

    // Setting up an initial condition
    std::string problem = "LA_3D";
    IC ic{problem, unknowns};
    ic(&mesh);
    std::cout << "Initial condition" << std::endl;
    // printMesh3DUniform(mesh, unknowns); 

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



    // Setting up a time differentiation scheme
    SSPRK_4 TS{unknowns, patchsize_x * patchsize_y * patchsize_z};

    // Creating a PDE
    BC_LA_3D_wave bc;
    double T = 0.1; // 0.1 for the big test, 1 for the small test
    double a1 = 1.0;
    double a2 = 1.0;
    double a3 = 1.0;
    LA3D<BC_LA_3D_wave> pde{ic, bc, T, a1, a2, a3};

    
    // Creating spatial diff schemes
    FC_3D fc_3d{"DD", "DD", "DD", patchsize_x, dx, Cx, hx, patchsize_y, dy, Cy,
        hy, patchsize_z, dz, Cz, hz, delta};
    std::vector<FC_3D *> diff_schemes;
    diff_schemes.push_back(&fc_3d);

    // Create a solver
    Solver<Mesh3DUniform, LA3D<BC_LA_3D_wave>, SSPRK_4, FC_3D*, FC_3D* >
        slv{mesh, pde, TS, diff_schemes};

    // Run the solver
    double dt = 0.001; // 0.001 for the big test, 0.05 for the small test
    slv.solve3D(dt);  

    // Printing solution
    std::cout << "Solution" << std::endl;
    // Print_Mesh2DUniform(mesh, unknowns);
    // printMesh3DUniform(mesh, unknowns);


    // std::string result_file = "result3D.txt";
    // printMesh3DUniform(mesh, result_file);
}