#include "elements.h"
#include "stiffnessMatrix.h"
#include "mesh.h"

#include <iostream>

using namespace std;

stiffnessMatrix::stiffnessMatrix(){
    _cons = nullptr;
    _size = 0;
    _freenodes = 0;
};

stiffnessMatrix::stiffnessMatrix(int ndofs){

    _globalStiff.allocate(ndofs);
    _size = ndofs;
    _freenodes = ndofs;
    _solution.allocate(ndofs);
    _cons = new bool[ndofs]();
    _loads.allocate(ndofs);

};

stiffnessMatrix::stiffnessMatrix(int ndofs, int fixeddofs){

    _globalStiff.allocate(ndofs);
    _size = ndofs;
    _freenodes = ndofs;
    _cons = new bool[ndofs]();
    _solution.allocate(ndofs);
    _loads.allocate(ndofs);

};

void stiffnessMatrix::addfixeddof(std::vector<unsigned int> & nodes, double value){
    
    for(unsigned int x : nodes){
        _cons[x] = true;
        _solution.setValue(value, x);
        _freenodes--;
    }

};

// double *stiffnessMatrix::operator[](int i){
//     return &_globalStiff[i*_size];
// };

stiffnessMatrix::~stiffnessMatrix(){
    std::cout << "Calling matrix destructor for stiffness matrix" << std::endl;
    _globalStiff.~matrix();
    _solution.~myVector();
    _loads.~myVector();
    delete[] _cons;
    _cons = nullptr;
};

ostream &operator<<(ostream &os, stiffnessMatrix &global){
    os << global._globalStiff;
    return os;
};

void stiffnessMatrix::assembleStiffness(mesh &msh, material &mat){

    std::vector<isoQuad4> elementlist = msh.getElementList();
    std::vector<node2d> nodelist = msh.getNodeList();
    unsigned int neles = msh.getNumEles();



    for(int i = 0; i < neles; i++){
        matrix elementstiff = matrix(4);
        elementlist[i].assembleLocalStiff(nodelist, mat, elementstiff);
        addStiffness(elementstiff, elementlist[i]);
    }

    // for(isoQuad4 element : msh.elementList){
    //     matrix = element.getStiffness();
    //     addStiffness(element, matrix);
    // }

};

void stiffnessMatrix::addStiffness(matrix & elementStiff, const isoQuad4 &element){

    std::cout << "From addStiffness matrix function" << std::endl;
    std::cout << "local stiffness matrix from element.assembleStiff(): " << std::endl;
    std::cout << elementStiff;

    for(int i = 0; i < 4; i++){
        for(int j = 0; j<4; j++){
            (*this)._globalStiff[element._nodes[i]][element._nodes[j]] += elementStiff[i][j];
        }
    }

};

void stiffnessMatrix::solveSystem(){
    
    matrix subMatrix = matrix(_freenodes);
    myVector subloads = myVector(_freenodes);
    myVector subsoln = myVector(_freenodes);
    
    int k = 0;
    int l = 0;
    double tempLoad = 0.0;

    for(int i = 0; i < _size; i++){

        if(!_cons[i]){
            tempLoad = _loads[i];
            for(int j = 0; j < _size; j++){
                if(!_cons[j]){
                    subMatrix[k][l] = _globalStiff[i][j];
                    l++;
                }
                else{
                    tempLoad -= _globalStiff[i][j]*_solution[j];
                }
            }
            subloads.setValue(tempLoad, k);
            k++;
            l = 0;
        }
    }

    std::cout << "Sub matrix for solver:" << std::endl;
    std::cout << subMatrix;
    std::cout << "Sub loads for solver:" << std::endl;
    std::cout << subloads;

    subMatrix.gaussElimination(subloads, subsoln);
    // subMatrix.gaussElimination(subloads.getVector(), subsoln.getVector());

    std::cout << "Sub Solution from solver" << std::endl;
    std::cout << subsoln;


};