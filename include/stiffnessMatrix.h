#pragma once

#include "mesh.h"

#include <iostream>



using namespace std;


class heatXFer: public mesh{
    
    private:

    matrix myGlobalStiffnessMatrix;
    unsigned int mySizeStiffnessMatrix;
    unsigned int myFreeDofs;
    bool * myConstrainedDof;
    myVector mySolution;
    myVector myLoads;

    void addStiffness( matrix &elementStiff, const isoQuad4 &element );

    void addloads(myVector & elementloads, const isoQuad4 & element);

    public:

    heatXFer();

    heatXFer(int ndofs);

    heatXFer(mesh & msh);

    heatXFer(string name, int dim);

    // double *operator[](int i);

    void allocate();

    friend ostream &operator<<(ostream &os, heatXFer &mat);

    void assembleStiffness(mesh &msh, material &mat);

    void assembleStiffness(material &mat);
    
    void bcFixedTemp(std::vector<unsigned int> & nodes, double value);

    void bcHeatGeneration(mesh & msh, material & mat, std::vector<unsigned int> & nodes, double q);

    void bcHeatGeneration(material & mat, std::vector<unsigned int> & nodes, double q);

    void bcHeatGeneration(mesh & msh, material & mat, double q);

    void bcHeatGeneration(material & mat, double q);

    void bcHeatFlux(mesh & msh, std::vector<unsigned int> & eles, double q);
    
    void bcHeatFlux(std::vector<unsigned int> &eles, double q);

    void initializeSolution(double temp);

    void solveSystem(string solver);

    void printSolution();

    void writeSolution(mesh & msh);

    void writeSolution();

    void writeStiffnessMatrix(string name);

    void writeLoads(string name);

    ~heatXFer();
};