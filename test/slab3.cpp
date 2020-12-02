#include "../include/mesh.h"
#include "../include/elements.h"
#include "../include/stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>

/*
Patch test for heat conduction through wall slab.
Boundary conditions: Heat  Flux on left  side q''
                     Fixed Temp on right side T_R
                     Heat Generation q_dot
                     Slab of length L and material with conductivity k
Exact Solution:
                T(x) = (q_dot/(2k))*(L^2 - x^2) + (q''/k)*(L-x) + T_R
            ___________________
    |-->   |                   |
    |-->   |                   |
 q''|-->   |       q_dot       | T_R    
    |-->   |                   |
    |-->   |___________________|

           |-->x               |
           0                   L
*/

/*
T_{4}\left(x\right)=\frac{q_{dot}}{2k}\left(L^{2}-x^{2}\right)+\frac{q}{k}\left(L-x\right)+T_{R}
*/
int main(){

    mesh mymesh = mesh(20, 12, 14);

    mymesh.addNode(0.00, 0.00);
    mymesh.addNode(0.25, 0.00);
    mymesh.addNode(0.50, 0.00);
    mymesh.addNode(0.75, 0.00);
    mymesh.addNode(1.00, 0.00);
    mymesh.addNode(0.00, 0.25);
    mymesh.addNode(0.25, 0.25);
    mymesh.addNode(0.50, 0.25);
    mymesh.addNode(0.75, 0.25);
    mymesh.addNode(1.00, 0.25);
    mymesh.addNode(0.00, 0.50);
    mymesh.addNode(0.25, 0.50);
    mymesh.addNode(0.50, 0.50);
    mymesh.addNode(0.75, 0.50);
    mymesh.addNode(1.00, 0.50);
    mymesh.addNode(0.00, 0.75);
    mymesh.addNode(0.25, 0.75);
    mymesh.addNode(0.50, 0.75);
    mymesh.addNode(0.75, 0.75);
    mymesh.addNode(1.00, 0.75);

    mymesh.addElement(0, 1, 6, 5);
    mymesh.addElement(1, 2, 7, 6);
    mymesh.addElement(2, 3, 8, 7);
    mymesh.addElement(3, 4, 9, 8);
    mymesh.addElement(5, 6, 11, 10);
    mymesh.addElement(6, 7, 12, 11);
    mymesh.addElement(7, 8, 13, 12);
    mymesh.addElement(8, 9, 14, 13);
    mymesh.addElement(10, 11, 16, 15);
    mymesh.addElement(11, 12, 17, 16);
    mymesh.addElement(12, 13, 18, 17);
    mymesh.addElement(13, 14, 19, 18);

    
    mymesh.addBoundaryElement(0, 1);
    mymesh.addBoundaryElement(1, 2);
    mymesh.addBoundaryElement(2, 3);
    mymesh.addBoundaryElement(3, 4);
    mymesh.addBoundaryElement(4, 9);
    mymesh.addBoundaryElement(9, 14);
    mymesh.addBoundaryElement(14, 19);
    mymesh.addBoundaryElement(19, 18);
    mymesh.addBoundaryElement(18, 17);
    mymesh.addBoundaryElement(17, 16);
    mymesh.addBoundaryElement(16, 15);
    mymesh.addBoundaryElement(15, 10);
    mymesh.addBoundaryElement(10, 5);
    mymesh.addBoundaryElement(5, 0);

    
    stiffnessMatrix global = stiffnessMatrix(20);
    material mat = material(1, 100.0);
    std::vector<unsigned int> temp = {4, 9, 14, 19};
    std::vector<unsigned int> temp2 = {0, 5, 10, 15};

    std::vector<unsigned int> heat = {11, 12, 13};

    global.bcFixedTemp(temp, 300.0);
    // global.bcFixedTemp(temp2, 300.0);
    global.bcHeatFlux(mymesh, heat, -2500.0);
    global.bcHeatGeneration(mymesh, mat, 10000.0);
    global.assembleStiffness(mymesh, mat);
    global.solveSystem();
    global.printSolution();
    global.writeSolution(mymesh);

    return 0;

}