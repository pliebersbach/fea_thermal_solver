#include "../include/elements.h"
#include "../include/stiffnessMatrix.h"
#include "../include/mesh.h"

#include <iostream>
#include <fstream>

using namespace std;

heatXFer::heatXFer(){
    /* 
    Default Constructor for a heatXFer type object.
    Object contains information to conduct 2D thermal analysis
    Constructor sets all variables to zero/default and pointers to null
     */
    myConstrainedDof = nullptr;
    mySizeStiffnessMatrix = 0;
    myFreeDofs = 0;
};

heatXFer::heatXFer(int ndofs){

    /* 
    heatXFer object constructor
    Constructor based on number of degrees of freedom (dofs) in the model.  The number
    of dofs is known ahead of time and passed to constructor.
    Inputs:
        ndofs - The number of degrees of freedom in the model
    */

    myGlobalStiffnessMatrix.allocate(ndofs);
    mySizeStiffnessMatrix = ndofs;
    myFreeDofs = ndofs;
    mySolution.allocate(ndofs);
    myConstrainedDof = new bool[ndofs]();
    myLoads.allocate(ndofs);

};

heatXFer::heatXFer(mesh & msh){
    /* 
    heatXFer object constructor
    Constructor based on a mesh class object.  Mesh object is created ahead of time
    and passed as arguement.  Constructor uses infomation from passed mesh object to
    complete initialization.
    Inputs:
        msh - A mesh class object that holds mesh information
    */
    mySizeStiffnessMatrix = msh.getNumNodes();
    myGlobalStiffnessMatrix.allocate(mySizeStiffnessMatrix);
    myFreeDofs = msh.getNumNodes();
    mySolution.allocate(mySizeStiffnessMatrix);
    myConstrainedDof = new bool[mySizeStiffnessMatrix]();
    myLoads.allocate(mySizeStiffnessMatrix);
}

heatXFer::heatXFer(string name, int dim){
    /* 
    heatXFer object constructor
    Constructor based on a mesh file and spatial dimension.  Constructor reads the mesh
    file referenced by the passed name and uses the read mesh information to complete
    initialization.
    Inputs:
        name - string containing file name of mesh data
        dim - spatial dimension of analysis; 2 for 2D and 3 for 3D
    */
    read_mesh(name, dim);

    mySizeStiffnessMatrix = getNumNodes();
    myGlobalStiffnessMatrix.allocate(mySizeStiffnessMatrix);
    myFreeDofs = getNumNodes();
    mySolution.allocate(mySizeStiffnessMatrix);
    myConstrainedDof = new bool[mySizeStiffnessMatrix]();
    myLoads.allocate(mySizeStiffnessMatrix);

}

void heatXFer::allocate(){

    /* 
    Allocates memory and initializes member variables for a
    heatXFer object created with the default constructor.  Requires
    that mesh has been constructed/read into the object already.
    */
    
    mySizeStiffnessMatrix = getNumNodes();
    myGlobalStiffnessMatrix.allocate(mySizeStiffnessMatrix);
    myFreeDofs = getNumNodes();
    mySolution.allocate(mySizeStiffnessMatrix);
    myConstrainedDof = new bool[mySizeStiffnessMatrix]();
    myLoads.allocate(mySizeStiffnessMatrix);

};

void heatXFer::bcFixedTemp(std::vector<unsigned int> & nodes, double value){
    /* 
    Applies a Fixed Temperature boundary condition to the specfified nodes
    of the mesh.
    Inputs:
        nodes - reference to vector of node ids to apply the bc to
        value - the temperature to fix on the specified nodes
    */
    std::cout << "Applying fixed temp bc's" << std::endl;

    for(unsigned int x : nodes){
        myConstrainedDof[x] = true;
        // mySolution.setValue(value, x);
        mySolution.myData[x] = value;
        myFreeDofs--;    
    }

    std::cout << "Finished applying fixed temp bc's" << std::endl;

};

void heatXFer::initializeSolution(double temp){
    for(int i = 0; i < mySizeStiffnessMatrix; i++){
        if(!myConstrainedDof[i]){
            mySolution.myData[i] = temp;
        }
    }
}

void heatXFer::bcHeatGeneration(mesh & msh, material & mat, std::vector<unsigned int> & eles,  double q){
    /* 
    Applies a load to the model due to internal heat generation at specified elements
    Inputs:
        msh - a reference to the mesh object
        mat - a reference to the material object
        eles - a reference to a vector of element ids subject to heat generation
        q - the amount of heat generated in standard units, W/m^3
    */
    const std::vector<isoQuad4> &elementlist = msh.getElementList();
    const std::vector<node2d> &nodelist = msh.getNodeList();
    for(int i : eles){
        isoQuad4 element = msh.getElement(i);
        myVector elementloads = myVector(4);
        element.assembleHeatGeneration(nodelist, mat, elementloads, q);
        addloads(elementloads, element);
    }
};

void heatXFer::bcHeatGeneration(material & mat, std::vector<unsigned int> & eles,  double q){
    /* 
    Applies a load to the model due to internal heat generation at specified elements.
    Inputs:
        mat - a reference to the material object
        eles - a reference to a vector of element ids subject to heat generation
        q - the amount of heat generated in standard units, W/m^3
    */
    const std::vector<isoQuad4> &elementlist = getElementList();
    const std::vector<node2d> &nodelist = getNodeList();
    for(int i : eles){
        isoQuad4 element = getElement(i);
        myVector elementloads = myVector(4);
        element.assembleHeatGeneration(nodelist, mat, elementloads, q);
        addloads(elementloads, element);
    }
};

void heatXFer::bcHeatGeneration(mesh & msh, material & mat,  double q){
    
    /* 
    Applies a load to the model due to internal heat generation in all elements
    Inputs:
        msh - a reference to the mesh object
        mat - a reference to the material object
        q - the amount of heat generated in standard units, W/m^3
    */
    const std::vector<isoQuad4> &elementlist = msh.getElementList();
    const std::vector<node2d> &nodelist = msh.getNodeList();
    for(isoQuad4 element : elementlist){
        myVector elementloads = myVector(4);
        element.assembleHeatGeneration(nodelist, mat, elementloads, q);
        addloads(elementloads, element);
    }
};

void heatXFer::bcHeatGeneration(material & mat,  double q){
    /* 
    Applies a load to the model due to internal heat generation in all elements
    Inputs:
        mat - a reference to the material object
        q - the amount of heat generated in standard units, W/m^3
    */
    const std::vector<isoQuad4> &elementlist = getElementList();
    const std::vector<node2d> &nodelist = getNodeList();
    for(isoQuad4 element : elementlist){
        myVector elementloads = myVector(4);
        element.assembleHeatGeneration(nodelist, mat, elementloads, q);
        addloads(elementloads, element);
    }
};

void heatXFer::bcHeatFlux(mesh & msh, std::vector<unsigned int> & eles, double q){
    /* 
    Applies a loads to the model due to heat flux on a boundary.
    Inputs:
        msh - a reference to a mesh class object
        eles - a reference to a vector of line element ids to apply the bc to
        q - the value of the heat flux applied to the boundary
    */
    double length;
    double heat;
    const std::vector<node2d> & nodelist = msh.getNodeList();
    for(unsigned int i : eles){
        line2 element  = msh.getBoundaryElement(i);
        double length = element.getLength(nodelist);
        heat = (q*length)/2.0;
        // (*this).myLoads.addValue(heat, element._nodes[0]);
        // (*this).myLoads.addValue(heat, element._nodes[1]);
        myLoads.myData[element._nodes[0]] += heat;
        myLoads.myData[element._nodes[1]] += heat;
    }
};

void heatXFer::bcHeatFlux(std::vector<unsigned int> & eles, double q){
    /* 
    Applies a loads to the model due to heat flux on a boundary.
    Inputs:
        eles - a reference to a vector of line element ids to apply the bc to
        q - the value of the heat flux applied to the boundary
    */
    
    std::cout << "Assembling Loads from heat flux bc's" << std::endl;

    double length;
    double heat;
    const std::vector<node2d> & nodelist = getNodeList();
    for(unsigned int i : eles){
        line2 element  = getBoundaryElement(i);
        double length = element.getLength(nodelist);
        heat = (q*length)/2.0;
        // (*this).myLoads.addValue(heat, element._nodes[0]);
        // (*this).myLoads.addValue(heat, element._nodes[1]);
        myLoads.myData[element._nodes[0]] += heat;
        myLoads.myData[element._nodes[1]] += heat;
    }

    std::cout << "Finished assembling loads form heat flux bc's" << std::endl;

};

// double *heatXFer::operator[](int i){
//     return &_globalStiff[i*_size];
// };

heatXFer::~heatXFer(){
    /* 
    Destructor for heatXFer class.
    Frees allocated memeory and resets pointers to null
    */
    std::cout << "Calling matrix destructor for stiffness matrix" << std::endl;
    myGlobalStiffnessMatrix.~matrix();
    mySolution.~myVector();
    myLoads.~myVector();
    delete[] myConstrainedDof;
    myConstrainedDof = nullptr;
};

ostream &operator<<(ostream &os, heatXFer &global){
    /* 
    Overloaded out operator to print contents of the heatXFer object to terminal.
    Prints the globalstiffness matrix, and the load and solution vectors.
    */
    os << "Global Stiffness Matrix: " << endl;
    os << global.myGlobalStiffnessMatrix;
    os << "Load Vector: " << endl;
    os << global.myLoads;
    os << "Solution Vector: " << endl;
    os << global.mySolution;
    return os;
};

void heatXFer::assembleStiffness(mesh &msh, material &mat){

    /*
    Assembles the global stiffness matrix for the analysis.
    Inputs:
        msh - A reference to a mesh class object
        mat - A reference to a material class object.
    */

    std::vector<isoQuad4> elementlist = msh.getElementList();
    std::vector<node2d> nodelist = msh.getNodeList();
    unsigned int neles = msh.getNumEles();



    for(int i = 0; i < neles; i++){
        matrix elementstiff = matrix(4);
        myVector elementloads = myVector(4);
        elementlist[i].assembleLocalStiff(nodelist, mat, elementstiff);
        addStiffness(elementstiff, elementlist[i]);        
    }

    // for(isoQuad4 element : msh.elementList){
    //     matrix = element.getStiffness();
    //     addStiffness(element, matrix);
    // }

};

void heatXFer::assembleStiffness(material &mat){
   
    /*
    Assembles the global stiffness matrix for the analysis.
    Inputs:
        msh - A reference to a mesh class object
        mat - A reference to a material class object.
    */

    std::vector<isoQuad4> elementlist = getElementList();
    std::vector<node2d> nodelist = getNodeList();
    unsigned int neles = getNumEles();



    for(int i = 0; i < neles; i++){
        matrix elementstiff = matrix(4);
        myVector elementloads = myVector(4);
        elementlist[i].assembleLocalStiff(nodelist, mat, elementstiff);
        addStiffness(elementstiff, elementlist[i]);        
    }

    // for(isoQuad4 element : msh.elementList){
    //     matrix = element.getStiffness();
    //     addStiffness(element, matrix);
    // }

};

void heatXFer::addStiffness(matrix & elementStiff, const isoQuad4 &element){

    /*
    Adds an element's stiffness matrix to the global stiffness matrix.
    Called from assembleStiffness()
    Inputs:
        elementStiff - A reference to an element's stiffness matrix
        element - A reference to isoQuad4 class object, owner of the matrix being added to global one
    */
    // std::cout << "From addStiffness matrix function" << std::endl;
    // std::cout << "local stiffness matrix from element.assembleStiff(): " << std::endl;
    // std::cout << elementStiff;

    for(int i = 0; i < 4; i++){
        for(int j = 0; j<4; j++){
            myGlobalStiffnessMatrix[element._nodes[i]][element._nodes[j]] += elementStiff[i][j];
        }
    }

};

void heatXFer::addloads(myVector & elementloads, const isoQuad4 & element){
    
    /* 
    Adds loads from element load vector to the global load vector.
    Inputs:
        elementLoads - reference to myVector class object containing element's loads
        element - a reference to a line class object, the owner of elements loads being added to global
    */
    
    for(int i = 0; i < 4; i++){
        // myLoads.addValue(elementloads[i], element._nodes[i]);
        myLoads.myData[element._nodes[i]] += elementloads[i];

    }
};

void heatXFer::solveSystem(string solver){

    /*
    Solves the linear system of equations using the Gauss Elmination method.
    */
    
    // std::cout << "Load vector from top of solver: " << std::endl;
    // std::cout << _loads;

    matrix subMatrix = matrix(myFreeDofs);
    myVector subloads = myVector(myFreeDofs);
    double *b = subloads.myData;
    myVector subsoln = myVector(myFreeDofs);
    double *x = subsoln.myData;
    
    int k = 0;
    int l = 0;
    double tempLoad = 0.0;

    for(int i = 0; i < mySizeStiffnessMatrix; i++){

        if(!myConstrainedDof[i]){
            tempLoad = myLoads[i];
            for(int j = 0; j < mySizeStiffnessMatrix; j++){
                if(!myConstrainedDof[j]){
                    subMatrix[k][l] = myGlobalStiffnessMatrix[i][j];
                    l++;
                }
                else{
                    tempLoad -= myGlobalStiffnessMatrix[i][j]*mySolution[j];
                }
            }
            // subloads.setValue(tempLoad, k);
            b[k] = tempLoad;

            k++;
            l = 0;
        }
    }

    subMatrix.printMatrix("cg_stiff");
    subloads.printVector("cg_loads");
    // std::cout << "Sub matrix for solver:" << std::endl;
    // std::cout << subMatrix;
    // std::cout << "Sub loads for solver:" << std::endl;
    // std::cout << subloads;
    if(solver == "gauss elimination"){
        subMatrix.gaussElimination(subloads, subsoln);
    }
    else if(solver == "jacobi iteration"){
        subMatrix.jacobiIteration(subloads, subsoln, 500);
    }
    else if(solver ==  "gauss seidel"){
        subMatrix.gaussSeidel(subloads, subsoln, 1000);
    }
    else if(solver == "sor"){
        subMatrix.sor(subloads, subsoln, 1000);
    }
    else if(solver == "conjugate gradient"){
        subMatrix.conjugateGradient(subloads, subsoln, 1000000);
    }
    else{
        cerr << "No Solver specified. Exiting" << endl;
        exit(1);
    }
    // subMatrix.gaussElimination(subloads.getVector(), subsoln.getVector());

    // std::cout << "Sub Solution from solver" << std::endl;
    // std::cout << subsoln;
    int j = 0;

    for(int i = 0; i < mySizeStiffnessMatrix; i++){
    
        if(!myConstrainedDof[i]){
            // mySolution.myData[i] = x[j];
            mySolution.myData[i] = subsoln.myData[j];
            j++;
        }
    
    }


};

void heatXFer::printSolution(){
    
    /*
    Prints the nodal solution to the terminal.
    */

    std::cout << "Nodal Solution:" << std::endl;
    std::cout << mySolution;
};

void heatXFer::writeSolution(mesh & msh){

    /*Writes the nodal solution, and mesh information to .mat files for postprocessing
    with Matlab/Octave
    Inputs:
        msh - A reference to a mesh class object
    */
    
    std::cout << "Writing Solution File." << std::endl;
    std::cout << "..." << std::endl;

    ofstream soln;
    soln.open("solution.mat", ios::trunc);
    for(int i = 0; i < mySizeStiffnessMatrix; i++){
        soln << mySolution[i] << endl;
    }
    soln.close();

    const std::vector<node2d> &nodelist = msh.getNodeList();
    ofstream points;
    points.open("points.mat", ios::trunc);
    for(int i = 0; i < msh.getNumNodes(); i++){
        points << nodelist[i]._x << ", " << nodelist[i]._y << endl;
    }
    points.close();

    const std::vector<isoQuad4> & elementlist = msh.getElementList();
    ofstream meshout;
    meshout.open("mesh.mat", ios::trunc);
    for(int i = 0; i < msh.getNumEles(); i++){
        meshout << elementlist[i]._nodes[0] << ", ";
        meshout << elementlist[i]._nodes[1] << ", ";
        meshout << elementlist[i]._nodes[2] << ", ";
        meshout << elementlist[i]._nodes[3] << endl;
    }

    std::cout << "Done writing solution files." << std::endl;
    std::cout << std::endl;
};

void heatXFer::writeSolution(){

    /*Writes the nodal solution, and mesh information to .mat files for postprocessing
    with Matlab/Octave
    */
    
    std::cout << "Writing Solution File." << std::endl;
    std::cout << "..." << std::endl;
    
    ofstream soln;
    soln.open("solution.mat", ios::trunc);
    for(int i = 0; i < mySizeStiffnessMatrix; i++){
        soln << mySolution[i] << endl;
    }
    soln.close();

    const std::vector<node2d> &nodelist = getNodeList();
    ofstream points;
    points.open("points.mat", ios::trunc);
    for(int i = 0; i < getNumNodes(); i++){
        points << nodelist[i]._x << ", " << nodelist[i]._y << endl;
    }
    points.close();

    const std::vector<isoQuad4> & elementlist = getElementList();
    ofstream meshout;
    meshout.open("mesh.mat", ios::trunc);
    for(int i = 0; i < getNumEles(); i++){
        meshout << elementlist[i]._nodes[0] << ", ";
        meshout << elementlist[i]._nodes[1] << ", ";
        meshout << elementlist[i]._nodes[2] << ", ";
        meshout << elementlist[i]._nodes[3] << endl;
    }

    std::cout << "Done writing solution files." << std::endl;
    std::cout << std::endl;
};

void heatXFer::writeStiffnessMatrix(string name){
    myGlobalStiffnessMatrix.printMatrix(name);
};

void heatXFer::writeLoads(string name){
    myLoads.printVector(name);
}