#include "../include/elements.h"
#include "../include/mesh.h"
#include "../include/stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>

int main(){

    mesh mymesh = mesh(9, 4, 8);

    mymesh.addNode(0.0, 0.0);
    mymesh.addNode(0.0, 1.0);
    mymesh.addNode(0.0, 2.0);    
    mymesh.addNode(2.0, 2.0);
    mymesh.addNode(2.0, 1.0);
    mymesh.addNode(2.0, 0.5);
    mymesh.addNode(4.0, 0.0);
    mymesh.addNode(4.0, 1.0);
    mymesh.addNode(4.0, 2.0);

    mymesh.addElement(1, 0, 5, 4);
    mymesh.addElement(2, 1, 4, 3);
    mymesh.addElement(4, 5, 6, 7);
    mymesh.addElement(3, 4, 7, 8);

    mymesh.addBoundaryElement(2, 1);
    mymesh.addBoundaryElement(1, 0);
    mymesh.addBoundaryElement(0, 5);
    mymesh.addBoundaryElement(5, 6);
    mymesh.addBoundaryElement(6, 7);
    mymesh.addBoundaryElement(7, 8);
    mymesh.addBoundaryElement(8, 3);
    mymesh.addBoundaryElement(3, 2);

    // mesh mymesh = mesh(4,1);

    // mymesh.addNode(0.0, 1.0);
    // mymesh.addNode(0.0, 0.0);
    // mymesh.addNode(2.0, 0.5);
    // mymesh.addNode(2.0, 1.0);
    // mymesh.addElement(0, 1, 2, 3);

    std::cout <<  mymesh;

    material mat = material(1.0, 5.0);

    isoQuad4 element = mymesh.getElement(0);
    
    const std::vector<node2d> & list = mymesh.getNodeList();
    
    stiffnessMatrix global = stiffnessMatrix(9);

    matrix stiff = matrix(4);
    
    element.assembleLocalStiff(list, mat, stiff);

    std::cout << "Stiffness matrix object test" << std::endl;

    std::cout << stiff;

    std::cout << "Global Assembly test" << std::endl;

    global.assembleStiffness(mymesh, mat);

    std::vector<unsigned int> left = {0, 1, 2};
    std::vector<unsigned int> right = {6, 7, 8};
    global.bcFixedTemp(left, 100.0);
    global.bcFixedTemp(right, 0.0);

    global.solveSystem();

    // for(int i = 0; i < 2; i++){
    //     for(int j = 0; j < 2; j++){
    //         std::cout << matrix[i][j] << ", ";
    //     }
    //     std::cout << std::endl;
    // }

    std::cout << global;

    global.printSolution();
    global.writeSolution(mymesh);

    std::cout << "Done" << std::endl;


    return 0;
};