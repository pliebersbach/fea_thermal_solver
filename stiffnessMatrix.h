#pragma once

#include "mesh.h"

#include <iostream>



using namespace std;


class stiffnessMatrix{
    
    private:

    double* _globalStiff = nullptr;
    unsigned int _size;

    public:

    stiffnessMatrix();

    stiffnessMatrix(int ndofs);

    double *operator[](int i);

    friend ostream &operator<<(ostream &os, stiffnessMatrix &mat);

    // void assembleStiffness(mesh &msh, material &mat);

    // void addStiffness( double(*)[4][4] matrix, const isoQuad4 &element );

    ~stiffnessMatrix();
};