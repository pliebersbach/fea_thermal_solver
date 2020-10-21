#include "elements.h"
#include "stiffnessMatrix.h"
#include "mesh.h"

#include <iostream>

using namespace std;

stiffnessMatrix::stiffnessMatrix(){
    _globalStiff = new double[2*2];
    _size = 2;
    _globalStiff[0] = M_PI;
    _globalStiff[1] = -1.0;
    _globalStiff[2] = -2.0;
    _globalStiff[3] = M_PI;
};

stiffnessMatrix::stiffnessMatrix(int ndofs){

    _globalStiff = new double[ndofs*ndofs]();
    _size = ndofs;

};

double *stiffnessMatrix::operator[](int i){
    return &_globalStiff[i*_size];
};

stiffnessMatrix::~stiffnessMatrix(){
    delete []_globalStiff;
    _globalStiff = nullptr;
    std::cout << "Stiff ness matrix deleted" << std::endl;
};

ostream &operator<<(ostream &os, stiffnessMatrix &mat){
    for(int i = 0; i < mat._size; i++){
        for(int j = 0; j < mat._size; j++){
            os << mat[i][j] << ", ";
        }
        os<<std::endl;
    }
    return os;
};

// void stiffnessMatrix::assembleStiffness(mesh &msh, material &mat){

//     std::vector<isoQuad4> elementlist = msh.getElementList();
//     std::vector<node2d> nodelist = msh.getNodeList();
//     unsigned int neles = msh.getNumEles();

//     for(int i = 0; i < neles; i++){
//         elementlist[i].assembleLocalStiff(nodelist, mat, *this);
//     }

//     // for(isoQuad4 element : msh.elementList){
//     //     matrix = element.getStiffness();
//     //     addStiffness(element, matrix);
//     // }

// }

// void addStiffness(double** matrix, const isoQuad4 &element){

//     std::cout << "From addStiffness matrix function" << std::endl;
//     std::cout << "local stiffness matrix from element.assembleStiff(): " << std::endl;
//     for(int i = 0; i < 4; i++){
//         for(int j = 0; j < 4; j++){
//             std::cout << matrix[i][j] << ", ";
//         }
//         std::cout << std::endl;
//     }

//     std::cout << "element nodes: " << std::endl;
//     std::cout << element._nodes[0] << ", " << element._nodes[1] << element._nodes[2] << element_nodes[3] << std::endl;

//     std::cout << "global stiffness matrix:" << std::endl;
//     std::cout << *this;

// }