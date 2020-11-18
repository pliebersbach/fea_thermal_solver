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
    
    void addfixeddof(std::vector<unsigned int> & nodes, double value);

    void solveSystem();

    ~stiffnessMatrix();
};