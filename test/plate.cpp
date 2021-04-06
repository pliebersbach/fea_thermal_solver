#include "../include/mesh.h"
#include "../include/elements.h"
#include "../include/stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>

int main(){
    mesh mymesh = mesh(49, 36, 24);

    mymesh.addNode(1.000, 0.000);
    mymesh.addNode(0.944, 0.330);
    mymesh.addNode(0.833, 0.545);
    mymesh.addNode(0.707, 0.707);
    mymesh.addNode(0.545, 0.838);
    mymesh.addNode(0.330, 0.944);
    mymesh.addNode(0.000, 1.000);
    mymesh.addNode(1.500, 0.000);
    mymesh.addNode(1.420, 0.496);
    mymesh.addNode(1.260, 0.817);
    mymesh.addNode(1.060, 1.060);
    mymesh.addNode(0.817, 1.260);
    mymesh.addNode(0.496, 1.420);
    mymesh.addNode(0.000, 1.500);
    mymesh.addNode(2.500, 0.000);
    mymesh.addNode(2.360, 0.826);
    mymesh.addNode(2.100, 1.360);
    mymesh.addNode(1.770, 1.770);
    mymesh.addNode(1.360, 2.100);
    mymesh.addNode(0.826, 2.360);
    mymesh.addNode(0.000, 2.500);
    mymesh.addNode(4.000, 0.000);
    mymesh.addNode(3.770, 1.320);
    mymesh.addNode(3.350, 2.180);
    mymesh.addNode(2.830, 2.830);
    mymesh.addNode(2.180, 3.350);
    mymesh.addNode(1.320, 3.770);
    mymesh.addNode(0.000, 4.000);
    mymesh.addNode(6.000, 0.000);
    mymesh.addNode(5.470, 1.920);
    mymesh.addNode(5.120, 3.320);
    mymesh.addNode(4.740, 4.740);
    mymesh.addNode(3.320, 5.120);
    mymesh.addNode(1.920, 5.470);
    mymesh.addNode(0.000, 6.000);
    mymesh.addNode(8.000, 0.000);
    mymesh.addNode(7.740, 2.710);
    mymesh.addNode(7.040, 4.580);
    mymesh.addNode(7.210, 7.210);
    mymesh.addNode(4.580, 7.040);
    mymesh.addNode(2.710, 7.740);
    mymesh.addNode(0.000, 8.000);
    mymesh.addNode(10.00, 0.000);
    mymesh.addNode(10.00, 3.500);
    mymesh.addNode(10.00, 6.500);
    mymesh.addNode(10.00, 10.00);
    mymesh.addNode(6.500, 10.00);
    mymesh.addNode(3.500, 10.00);
    mymesh.addNode(0.000, 10.00);

    mymesh.addElement(0, 7, 8, 1);
    mymesh.addElement(1, 8, 9, 2);
    mymesh.addElement(2, 9, 10, 3);
    mymesh.addElement(3, 10, 11, 4);
    mymesh.addElement(4, 11, 12, 5);
    mymesh.addElement(5, 12, 13, 6);
    mymesh.addElement(7, 14, 15, 8);
    mymesh.addElement(8, 15, 16, 9);
    mymesh.addElement(9, 16, 17, 10);
    mymesh.addElement(10, 17, 18, 11);
    mymesh.addElement(11, 18, 19, 12);
    mymesh.addElement(12, 19, 20, 13);
    mymesh.addElement(14, 21, 22, 15);
    mymesh.addElement(15, 22, 23, 16);
    mymesh.addElement(16, 23, 24, 17);
    mymesh.addElement(17, 24, 25, 18);
    mymesh.addElement(18, 25, 26, 19);
    mymesh.addElement(19, 26, 27, 20);
    mymesh.addElement(21, 28, 29, 22);
    mymesh.addElement(22, 29, 30, 23);
    mymesh.addElement(23, 30, 31, 24);
    mymesh.addElement(24, 31, 32, 25);
    mymesh.addElement(25, 32, 33, 26);
    mymesh.addElement(26, 33, 34, 27);
    mymesh.addElement(28, 35, 36, 29);
    mymesh.addElement(29, 36, 37, 30);
    mymesh.addElement(30, 37, 38, 31);
    mymesh.addElement(31, 38, 39, 32);
    mymesh.addElement(32, 39, 40, 33);
    mymesh.addElement(33, 40, 41, 34);
    mymesh.addElement(35, 42, 43, 36);
    mymesh.addElement(36, 43, 44, 37);
    mymesh.addElement(37, 44, 45, 38);
    mymesh.addElement(38, 45, 46, 39);
    mymesh.addElement(39, 46, 47, 40);
    mymesh.addElement(40, 47, 48, 41);
    
    mymesh.addBoundaryElement(0, 7);
    mymesh.addBoundaryElement(7, 14);
    mymesh.addBoundaryElement(14, 21);
    mymesh.addBoundaryElement(21, 28);
    mymesh.addBoundaryElement(28, 35);
    mymesh.addBoundaryElement(35, 42);
    mymesh.addBoundaryElement(42, 43);
    mymesh.addBoundaryElement(43, 44);
    mymesh.addBoundaryElement(44, 45);
    mymesh.addBoundaryElement(45, 46);
    mymesh.addBoundaryElement(46, 47);
    mymesh.addBoundaryElement(47, 48);
    mymesh.addBoundaryElement(48, 41);
    mymesh.addBoundaryElement(41, 34);
    mymesh.addBoundaryElement(34, 27);
    mymesh.addBoundaryElement(27, 20);
    mymesh.addBoundaryElement(20, 13);
    mymesh.addBoundaryElement(13, 6);
    mymesh.addBoundaryElement(6, 5);
    mymesh.addBoundaryElement(5, 4);
    mymesh.addBoundaryElement(4, 3);
    mymesh.addBoundaryElement(3, 2);
    mymesh.addBoundaryElement(2, 1);
    mymesh.addBoundaryElement(1, 0);
    
    heatXFer global = heatXFer(49);
    material mat = material(1, 100.0);
    std::vector<unsigned int> temp = {42, 43, 44, 45, 46, 47, 48};
    std::vector<unsigned int> temp2 = {0, 1, 2, 3, 4, 5, 6};

    std::vector<unsigned int> heat = {18, 19, 20, 21, 22, 23};

    global.bcFixedTemp(temp, 300.0);
    // global.bcFixedTemp(temp2, 500.0);
    global.bcHeatFlux(mymesh, heat, 5000.0);
    global.assembleStiffness(mymesh, mat);
    // global.solveSystem("gauss seidel");
    global.solveSystem("sor");
    global.printSolution();
    global.writeSolution(mymesh);
    

    return 0;

}