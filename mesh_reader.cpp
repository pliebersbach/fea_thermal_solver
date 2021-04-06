#include "include/mesh.h"
#include "include/elements.h"
#include "include/stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <time.h>

int main(int argc, char* argv[])
{   
    string filename = argv[1];

    mesh mymesh = mesh();

    // mymesh.read_mesh("mesh_.pl", 2);
    mymesh.read_mesh(filename, 2);

    // mymesh.read_mesh("adsflh.csv", 2);

    // std::cout << "returned to main" << endl;
    
    heatXFer global = heatXFer(mymesh);
    material mat = material(1, 100.0);
    
    // std::vector<unsigned int> temp = {19, 22, 21, 20, 15, 18, 17, 16, 9};
    std::vector<unsigned int> temp = mymesh.getNodeSet(0);
    // std::vector<unsigned int> heat = {0, 1, 2, 3, 4, 5, 6, 7};
    std::vector<unsigned int> heat = mymesh.getSideSet(0);

    struct timespec start;
    struct timespec end;
    struct timespec middle;

    float time = 0;
    float time2 = 0;

    clock_gettime(CLOCK_MONOTONIC, &start);

    global.bcFixedTemp(temp, 300.0);
    global.bcHeatFlux(mymesh, heat, 5000.0);
    global.assembleStiffness(mymesh, mat);
    global.initializeSolution(315.0);

    clock_gettime(CLOCK_MONOTONIC, &middle);

    global.solveSystem("gauss elimination");
    // global.solveSystem("jacobi iteration");
    // global.solveSystem("gauss seidel");
    // global.solveSystem("sor");
    // global.solveSystem("conjugate gradient");


    clock_gettime(CLOCK_MONOTONIC, &end);

    size_t duration_usec = (middle.tv_sec - start.tv_sec)*1000*1000;
    duration_usec += (middle.tv_nsec - start.tv_nsec)/1000;
    time +=  duration_usec;

    size_t duration_usec2 = (end.tv_sec - start.tv_sec)*1000*1000;
    duration_usec2 += (end.tv_nsec - start.tv_nsec)/1000;
    time2 += duration_usec2;
    
    std::cout << "Assembly time: " << time << " microseconds" << std::endl;

    std::cout << "Total solution time: " << time2 << " microseconds" << std::endl;
    // global.writeStiffnessMatrix("cg_stiff");
    // global.writeLoads("cg_loads");
    // global.printSolution();
    global.writeSolution(mymesh);
    // cout << mymesh;

    return 0;
}