
#include "printing.h"
#include <algorithm>




void Print_Mat(const int * A, int nrows, int ncols)
{
    
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i] << " ";
        }
        std::cout << std::endl;
    }   
}


void Print_Mat(const double * A, int nrows, int ncols)
{
    std::cout.precision(5);
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i] << " ";
        }
        std::cout << std::endl;
    }   
}

void Print_Mat(const double * A, int nrows, int ncols, bool flag)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j + i*ncols] << " ";
        }
        std::cout << std::endl;
    }   
}

void Print_Mat(const std::complex<double> * A, int nrows, int ncols)
{
    for(int i = 0; i < nrows; i++)
    {
        for(int j = 0; j < ncols; j++)
        {
            std::cout << A[j*nrows + i].real() << " + " << A[j*nrows + i].imag() << " I ";
        }
        std::cout << std::endl;
    }   
}

void Print_VectorField(const VectorField &field)
{
    int N = field.getUnknowns();
    int length = field.getLength();
    std::vector<double*> flow = field.getFlow();
    double data [length*N] ; 
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < length; j++)
        {
            data[i*length + j] = (flow[i])[j];
        }    
    }
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << data[j*length + i] << "  ";
        }
        std::cout << std::endl; 
    }
}

void Print_VectorField(const VectorField &field, int precision)
{
    int N = field.getUnknowns();
    int length = field.getLength();
    std::cout.precision(precision);
    std::vector<double*> flow = field.getFlow();
    double data [length*N] ; 
    for (int i = 0; i < N; i++)
    {
        for(int j = 0; j < length; j++)
        {
            data[i*length + j] = (flow[i])[j];
        }    
    }
    for(int i = 0; i < length; i++)
    {
        for(int j = 0; j < N; j++)
        {
            std::cout << data[j*length + i] << "  ";
        }
        std::cout << std::endl; 
    }
}

void Print_VectorField(const VectorField &field, bool print)
{
    if(print)
    {
        VectorField u{field.getUnknowns() + 1, field.getLength()};
        double temp[field.getLength()];
        for(int i = 0; i < field.getLength(); i++)
        {
            temp[i] = double(i + 1);
        }
        u.setField(0, field.getLength(), temp);
        for(int i = 0; i < field.getUnknowns(); i++)
        {
            u.setField(i+1, field.getLength(), field.getField(i));
        }
        Print_VectorField(u);
    }      
    else
    {
        Print_VectorField(field);
    } 
}


void Print_VectorField(const VectorField &field, bool print, int precision)
{
    if(print)
    {
        VectorField u{field.getUnknowns() + 1, field.getLength()};
        double temp[field.getLength()];
        for(int i = 0; i < field.getLength(); i++)
        {
            temp[i] = double(i + 1);
        }
        u.setField(0, field.getLength(), temp);
        for(int i = 0; i < field.getUnknowns(); i++)
        {
            u.setField(i+1, field.getLength(), field.getField(i));
        }
        Print_VectorField(u, precision);
    }      
    else
    {
        Print_VectorField(field, precision);
    } 
}

void Print_Node(const Node &node)
{
    std::cout << "Node position: " << node.getPosition()[0] << std::endl;
    std::cout << "Node index: " << node.getIndex() << std::endl;
    std::cout << "Unknowns: " << node.getUnknowns() << std::endl;
    std::cout << "Stored values:" << std::endl;
    //Print_Mat(node.getValues(), node.getUnknowns(), 1);
    double a[node.getUnknowns()];
    for(int i = 0; i < node.getUnknowns(); i++)
    {
        a[i] = node.getValues()[i];
    } 
    Print_Mat(a, node.getUnknowns(), 1);
}

void Print_SparseMatrix_csr(int rows, int cols, int* rows_start, int* rows_end,
    int* col_indx, double* values)
{
    double entries[rows*cols];
    int current_col = 0;
    int val_index = 0;
    for(int i = 0; i < rows*cols; i++)
    {
        entries[i] = 0.0;
    }
    for(int i = 0; i < rows; i++)
    {
        current_col = rows_start[i];
        for(int j = rows_start[i]; j < rows_end[i]; j++)
        {
            entries[rows*col_indx[j] + i] = values[j];
        }
    }
    Print_Mat(entries, rows, cols);

}

void Print_Patch2D(const Patch2D &patch)
{
    int unknowns = patch.getNode(0)->getUnknowns();
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int i = 0; i < patchsize_y; i++)
        {
            for(int j = 0; j < patchsize_x; j++)
            {
                std::cout << patch.getNode(patchsize_y*j + i)->getValues()[u] 
                    << " ";
            }
            std::cout << std::endl;
        }
    }
}

void Print_Patch2D(const Patch2D &patch, int unknowns)
{
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int i = 0; i < patchsize_y; i++)
        {
            for(int j = 0; j < patchsize_x; j++)
            {
                std::cout << patch.getNode(patchsize_y*j + i)->getValues()[u] 
                    << " ";
            }
            std::cout << std::endl;
        }
    }
}


void Print_Mesh2DUniform(const Mesh2DUniform &mesh)
{
    int npatches_x = mesh.getNpatches_x();
    int npatches_y = mesh.getNpatches_y();
    int npatches = npatches_x * npatches_y;
    auto patches = mesh.getPatches();
    for(int j = 0; j < npatches_x; j++)
    {
        for(int i = 0; i < npatches_y; i++)
        {
            std::cout << "PATCH " << npatches_y*j + i << ". position : (" << 
                i << ", " << j << ")" << std::endl;
            Print_Patch2D(*patches[npatches_y*j + i]);
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }   
}

void Print_Mesh2DUniform(const Mesh2DUniform &mesh, int unknowns)
{
    int npatches_x = mesh.getNpatches_x();
    int npatches_y = mesh.getNpatches_y();
    int npatches = npatches_x * npatches_y;
    auto patches = mesh.getPatches();
    for(int j = 0; j < npatches_x; j++)
    {
        for(int i = 0; i < npatches_y; i++)
        {
            std::cout << "PATCH " << npatches_y*j + i << ". position : (" << 
                i << ", " << j << ")" << std::endl;
            Print_Patch2D(*patches[npatches_y*j + i], unknowns);
            std::cout << std::endl;
            std::cout << std::endl;
        }
    }   
}

void printPatch3D(const Patch3D &patch)
{
    int unknowns = patch.getNode(0)->getUnknowns();
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    int patchsize_z = patch.getNz();
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int k = 0; k < patchsize_z; k++)
        {
            std::cout << "z = " << 
                patch.getNode(k * patchsize_x * patchsize_y)->getPosition()[2]
                << std::endl;
            for(int i = 0; i < patchsize_y; i++)
            {
                for(int j = 0; j < patchsize_x; j++)
                {
                    std::cout << 
                        patch.getNode(k * patchsize_x * patchsize_y + 
                        patchsize_y*j + i)->getValues()[u] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}



void printPatch3D(const Patch3D &patch, int unknowns)
{
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    int patchsize_z = patch.getNz();
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int k = 0; k < patchsize_z; k++)
        {
            std::cout << "z = " << 
                patch.getNode(k * patchsize_x * patchsize_y)->getPosition()[2]
                << std::endl;
            for(int i = 0; i < patchsize_y; i++)
            {
                for(int j = 0; j < patchsize_x; j++)
                {
                    std::cout << 
                        patch.getNode(k * patchsize_x * patchsize_y + 
                        patchsize_y*j + i)->getValues()[u] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}


void printPatch3D(const Patch3D &patch, const int x0, const int x1,
    const int y0, const int y1, const int z0, const int z1)
{
    int unknowns = patch.getNode(0)->getUnknowns();
    int patchsize_x = patch.getNx();
    int patchsize_y = patch.getNy();
    int patchsize_z = patch.getNz();
    int j0 = std::max(0, x0);
    int i0 = std::max(0, y0);
    int k0 = std::max(0, z0);
    int j1 = std::min(x1, patchsize_x);
    int i1 = std::min(y1, patchsize_y);
    int k1 = std::min(z1, patchsize_z);
    for(int u = 0; u < unknowns; u++)
    {
        std::cout << "Unknown " << u << std::endl;
        for(int k = k0; k < k1; k++)
        {
            std::cout << "z = " << 
                patch.getNode(k * patchsize_x * patchsize_y)->getPosition()[2]
                << std::endl;
            for(int i = i0; i < i1; i++)
            {
                for(int j = j0; j < j1; j++)
                {
                    std::cout << 
                        patch.getNode(k * patchsize_x * patchsize_y + 
                        patchsize_y*j + i)->getValues()[u] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
        std::cout << std::endl;
    }
}


void printMesh3DUniform(const Mesh3DUniform &mesh)
{
    int npatches_x = mesh.getNpatches_x();
    int npatches_y = mesh.getNpatches_y();
    int npatches_z = mesh.getNpatches_z();
    int npatches = npatches_x * npatches_y * npatches_z;
    auto patches = mesh.getPatches();
    for(int k = 0; k < npatches_z; k++)
    {
        for(int j = 0; j < npatches_x; j++)
        {
            for(int i = 0; i < npatches_y; i++)
            {
                std::cout << "PATCH " << npatches_y*npatches_z*k +
                    npatches_y*j + i << ". position : (" << i << ", " << j <<
                    ", " << k << ")" << std::endl;
                printPatch3D(*patches[npatches_y*npatches_x*k + npatches_y*j + i]);
                std::cout << std::endl;
                std::cout << std::endl;
            }
        }  
    }
}

void printMesh3DUniform(const Mesh3DUniform &mesh, std::string filename) 
{
    auto patches = mesh.getPatches();
    int unknowns = patches[0]->getNode(0)->getUnknowns();
    int patchsize_x;
    int patchsize_y;
    int patchsize_z;
    int N;
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        myfile << patches.size();
        myfile << "\n";   
        myfile << unknowns;
        myfile << "\n";              
        for(int i = 0; i < patches.size(); i++)
        {
            auto patch = patches[i];
            auto v = patches[i]->getFlow();
            patchsize_x = patches[i]->getNx();
            patchsize_y = patches[i]->getNy();
            patchsize_z = patches[i]->getNz();
            myfile << patchsize_x;
            myfile << "\n";     
            myfile << patchsize_y;
            myfile << "\n"; 
            myfile << patchsize_z;
            myfile << "\n";   
            N = patchsize_x * patchsize_y * patchsize_z;                              
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << 
                    patch->getNode(j)->getPosition()[0];
                myfile << "\n";
            }
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << 
                    patch->getNode(j)->getPosition()[1];
                myfile << "\n";
            }
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << 
                    patch->getNode(j)->getPosition()[2];
                myfile << "\n";
            }                        
            for(int j = 0; j < unknowns; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    myfile << std::fixed << std::setprecision(17) <<
                        v.getFieldValue(j, k);
                    myfile << "\n";
                }
            }     
        }
        myfile.close();
    }
    else std::cout << "Unable to open file"; 
}

void printMesh3DUniform(const Mesh3DUniform &mesh, int unknowns)
{
    int npatches_x = mesh.getNpatches_x();
    int npatches_y = mesh.getNpatches_y();
    int npatches_z = mesh.getNpatches_z();
    int npatches = npatches_x * npatches_y * npatches_z;
    auto patches = mesh.getPatches();
    for(int k = 0; k < npatches_z; k++)
    {
        for(int j = 0; j < npatches_x; j++)
        {
            for(int i = 0; i < npatches_y; i++)
            {
                std::cout << "PATCH " << npatches_y*npatches_z*k +
                    npatches_y*j + i << ". position : (" << i << ", " << j <<
                    ", " << k << ")" << std::endl;
                printPatch3D(
                    *patches[npatches_y*npatches_x*k + npatches_y*j + i], 
                    unknowns);
                std::cout << std::endl;
                std::cout << std::endl;
            }
        }  
    }    
}
