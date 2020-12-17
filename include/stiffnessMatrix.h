#pragma once

#include "mesh.h"

#include <iostream>



using namespace std;


class stiffnessMatrix{
    
    private:

    matrix _globalStiff;
    unsigned int _size;
    unsigned int _freenodes;
    bool * _cons;
    myVector _solution;
    myVector _loads;

    public:

    stiffnessMatrix();

    stiffnessMatrix(int ndofs);

    stiffnessMatrix(int ndofs, int fixeddofs);

    // double *operator[](int i);

    friend ostream &operator<<(ostream &os, stiffnessMatrix &mat);

    void assembleStiffness(mesh &msh, material &mat);

    void addStiffness( matrix &elementStiff, const isoQuad4 &element );
    
    void bcFixedTemp(std::vector<unsigned int> & nodes, double value);

    void bcHeatGeneration(mesh & msh, material & mat, std::vector<unsigned int> & nodes, double q);

    void bcHeatGeneration(mesh & msh, material & mat, double q);

    void bcHeatFlux(mesh & msh, std::vector<unsigned int> & eles, double q);
    
    void addloads(myVector & elementloads, const isoQuad4 & element);

    void solveSystem();

    void printSolution();

    void writeSolution(mesh & msh);

    ~stiffnessMatrix();
};