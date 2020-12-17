#include "include/mesh.h"
#include "include/elements.h"
#include "include/stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>

int main()
{
    mesh mymesh = mesh();

    mymesh.read_mesh("points.csv", 2);

    // mymesh.read_mesh("adsflh.csv", 2);

    std::cout << "returned to main" << endl;
    
    stiffnessMatrix global = stiffnessMatrix(63);
    material mat = material(1, 100.0);
    // std::vector<unsigned int> temp = {19, 22, 21, 20, 15, 18, 17, 16, 9};
    std::vector<unsigned int> temp = mymesh.getNodeSet(0);
    // std::vector<unsigned int> heat = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<unsigned int> heat = mymesh.getSideSet(2);
    global.bcFixedTemp(temp, 300.0);
    global.bcHeatFlux(mymesh, heat, 5000.0);
    global.assembleStiffness(mymesh, mat);
    global.solveSystem();
    global.printSolution();
    global.writeSolution(mymesh);
    // cout << mymesh;

    return 0;
}