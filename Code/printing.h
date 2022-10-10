
#ifndef PRINTING_H
#define PRINTING_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <complex>
#include <string>
#include "VectorField.h"
#include "Node.h"
#include "Patch.h"
#include "Patch1D.h"
#include "Patch2D.h"
#include "Patch3D.h"
#include "Mesh.h"
#include "Mesh2D.h"
#include "Mesh3D.h"

void Print_Mat(const double * A, int nrows, int ncols);
void Print_Mat(const double * A, int nrows, int ncols, bool flag);
void Print_Mat(const std::complex<double> * A, int nrows, int ncols);
void Print_Mat(const int* A, int nrows, int ncols);

void Print_VectorField(const VectorField &field);
void Print_VectorField(const VectorField &field, int precision);
void Print_VectorField(const VectorField &field, bool print, int precision);
void Print_VectorField(const VectorField &field, bool print);
void Print_Node(const Node &node);
void Print_SparseMatrix_csr(int rows, int cols, int* rows_start, int* rows_end,
    int* col_indx, double* values);

template<class Patch>
void Print_Patch1D(const Patch &patch)
{
    int unknowns;
    std::cout << patch.getNnodes() << " nodes" << std::endl;
    for(int i = 0; i < patch.getNnodes(); i++)
    {
        std::cout << "Position = " << patch.getNode(i)->getPosition()[0] << " ";
        unknowns = patch.getNode(i)->getUnknowns();
        for(int j = 0; j < unknowns; j++)
        {
            std::cout << "Value " << j << " = " << (patch.getNode(i)->getValues())[j] << "  ";
        }
        std::cout << std::endl;
    }
}

template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh)
{
    std::vector<Patch*> patches = mesh.getPatches();
    for(int i = 0; i < patches.size(); i++)
    {
        std::cout << "Patch " << i + 1 << std::endl;
        Print_Patch1D(*patches[i]);
        std::cout << std::endl;
        std::cout << std::endl;

    }
}

template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh, int unknowns, int intrb,
    std::string filename)
{
    int N;
    std::vector<Patch*> patches = mesh.getPatches();
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        myfile << patches.size();
        myfile << "\n";   
        myfile << unknowns;
        myfile << "\n";   
        myfile << intrb;
        myfile << "\n";            
        for(int i = 0; i < patches.size(); i++)
        {
            auto patch = patches[i];
            auto v = patches[i]->getFlow();
            N = v.getLength();
            myfile << N;
            myfile << "\n";             
            std::cout << "Patch " << i + 1 << std::endl;
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << patch->getNode(j)->getPosition()[0];
                myfile << "\n";
            }
            for(int j = 0; j < unknowns; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    // myfile << v.getFieldValue(j, k);
                    myfile << std::fixed << std::setprecision(17) << v.getFieldValue(j, k);
                    myfile << "\n";
                }
            }     
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";    
}



template<class Patch>
void Print_Mesh1D(const Mesh<Patch> &mesh, int unknowns, int intrb, int vector,
    std::string filename)
{
    int N;
    unknowns = 1;
    std::vector<Patch*> patches = mesh.getPatches();
    std::ofstream myfile (filename);
    if (myfile.is_open())
    {
        myfile << patches.size();
        myfile << "\n";   
        myfile << unknowns;
        myfile << "\n";   
        myfile << intrb;
        myfile << "\n";            
        for(int i = 0; i < patches.size(); i++)
        {
            auto patch = patches[i];
            auto v = patches[i]->getFlow();
            N = v.getLength();
            myfile << N;
            myfile << "\n";             
            // std::cout << "Patch " << i + 1 << std::endl;
            for(int j = 0; j < N; j++)
            {
                myfile << std::fixed << std::setprecision(17) << patch->getNode(j)->getPos();
                myfile << "\n";
            }
            for(int j = 0; j < unknowns; j++)
            {
                for(int k = 0; k < N; k++)
                {
                    // myfile << v.getFieldValue(j, k);
                    myfile << std::fixed << std::setprecision(17) << v.getFieldValue(vector + j, k);
                    myfile << "\n";
                }
            }     
        }
        myfile.close();
    }
    else std::cout << "Unable to open file";    
}



void Print_Patch2D(const Patch2D &patch);

void Print_Patch2D(const Patch2D &patch, int unknowns);

void Print_Mesh2DUniform(const Mesh2DUniform &mesh);

void Print_Mesh2DUniform(const Mesh2DUniform &mesh, int unknowns);

void printPatch3D(const Patch3D &patch);

// void printPatch3D(const Patch3D &patch, std::string filename);

void printPatch3D(const Patch3D &patch, int unknowns);

void printPatch3D(const Patch3D &patch, const int x0, const int x1,
    const int y0, const int y1, const int z0, const int z1);

void printMesh3DUniform(const Mesh3DUniform &mesh);

void printMesh3DUniform(const Mesh3DUniform &mesh, std::string filename);

void printMesh3DUniform(const Mesh3DUniform &mesh, int unknowns);


#endif 