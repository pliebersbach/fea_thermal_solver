#include "../include/elements.h"
#include "../include/mesh.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>


mesh::mesh(){
    cout << "Default Constructor" << endl;
    _numnodes = 0;
    _numelements = 0;
    _numBoundaryEles = 0;

};

mesh::mesh(unsigned int nnodes, unsigned int nelements){
    cout << "Mesh Constructor Called.  List Locations:" << endl;
    nodeList.reserve(nnodes);
    elementList.reserve(nelements);
    _numnodes=nnodes;
    _numelements=nelements;
    cout << "Element List Location: " << &elementList << endl;
    cout << "Node List Location:    " << &nodeList << endl;

};

mesh::mesh(unsigned int nnodes, unsigned int nelements, unsigned int nboundaryeles){
    nodeList.reserve(nnodes);
    elementList.reserve(nelements);
    boundaryEleList.reserve(nboundaryeles);
    _numnodes = nnodes;
    _numelements = nelements;
    _numBoundaryEles = nboundaryeles;
}

/*
mesh::mesh(unsigned int nnodes, unsigned int nelements, int dim){

    if(dim==3){
         nodeList.reserve(nnodes);
    }

 elementList.reserve(nelements);

};
*/

void mesh::read_mesh(string name, unsigned int dim){
    std::ifstream inFile;
    inFile.open(name);
    if(!inFile){
        cerr << "Unable to open file";
        exit(1);
    }
    else{
        cout << "File opened successfully" << endl;
    }

    //Read coordinates from file
    std::string line;

    while(std::getline(inFile, line)){
        
        if(line.find("Coordinates")==0){
            std::getline(inFile, line);
            unsigned int num_nodes = stoi(line);

            nodeList.reserve(num_nodes);
            _numnodes = num_nodes;

            double x, y, z;
            // cout << "Num Nodes: " << num_nodes << endl;

            for(int i = 0; i < num_nodes; i++){
                inFile >> x >> y >> z;
                if (dim == 2)
                {
                    nodeList.push_back(node2d(x, y));
                }
                // else if (dim == 3)
                // {
                //     nodeList.push_back(node3d(x, y, z));
                // }
                // std::cout <<  "i: " << i << " X: " << x << " Y: " << y << " Z: " << z << std::endl;
            }
        }

        if(line.find("Elements") == 0){
            std::getline(inFile, line);
            unsigned int num_eles = stoi(line);
            elementList.reserve(num_eles);
            _numelements = num_eles;
            unsigned int n0, n1, n2, n3;
            // cout << "Num Elements: " << num_eles << endl;

            for(int i = 0; i < num_eles; i++){
                inFile >> n0 >> n1 >> n2 >> n3;
                elementList.push_back(isoQuad4(n0, n1, n2, n3));
                // cout << "i: " << i << " n0: " << n0 << " n1: " << n1 << " n2: " << n2 << " n3: " << n3 << endl;
            }
        }

        if(line.find("Edges")==0){
            std::getline(inFile, line);
            unsigned int num_eles = stoi(line);
            boundaryEleList.reserve(num_eles);
            _numBoundaryEles = num_eles;
            unsigned int n0, n1;
            // cout << "Num Boundary Eles: " << num_eles << endl;
            for(int i = 0; i < num_eles; i++){
                inFile >> n0 >> n1;
                boundaryEleList.push_back(line2(n0, n1));
                // cout << "i: " << i << " n0: " << n0 << " n1: " << n1 << endl;
            }
        }
    }

    inFile.close();
}

void mesh::addElement(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4){
    elementList.push_back(isoQuad4(n1, n2, n3, n4));
}

void mesh::addNode(double x, double y){
    nodeList.push_back(node2d(x, y));
};

void mesh::addBoundaryElement(unsigned int n1, unsigned int n2){
    boundaryEleList.push_back(line2(n1, n2));
}

isoQuad4 mesh::getElement(unsigned int i){
    return elementList[i];
}

line2 mesh::getBoundaryElement(unsigned int i){
    return boundaryEleList[i];
}

const std::vector<isoQuad4> & mesh::getElementList(){
    return elementList;
}

const std::vector<node2d> & mesh::getNodeList(){
    return nodeList;
}

const std::vector<line2> & mesh::getBoundaryElementList(){
    return boundaryEleList;
}

unsigned int mesh::getNumEles(){
    return _numelements;
};

unsigned int mesh::getNumNodes(){
    return _numnodes;
}

// mesh::~mesh(){
//     std::cout << "Mesh destructor" << std::endl;
// };

ostream & operator<<(ostream & os, const mesh &msh){

    os << "Node List" << endl;

    for(node2d node : msh.nodeList){
        os << "x: " << node._x << " y: " << node._y << "Location: " << &node << endl;
    }

    os << "Element List" << endl;

    // for(isoQuad4 element : msh.elementList){
    //     os << element._nodes[0] << " " << element._nodes[1] << " " << element._nodes[2] << " " << element._nodes[3] << "element location: " << &element << endl;
    // }

    for(int i = 0; i < msh._numelements; i++){
        os << msh.elementList[i]._nodes[0] << " " << msh.elementList[i]._nodes[1] << " " << msh.elementList[i]._nodes[2] << " " << msh.elementList[i]._nodes[3] << "element location: " << &msh.elementList[i] << endl;
    }

    os << "Boundary Element List" << endl;

    for(int i = 0; i < msh._numBoundaryEles; i++){
        os << msh.boundaryEleList[i]._nodes[0] << " " << msh.boundaryEleList[i]._nodes[1] << endl;
    }

    return os;

}

double mesh::getTotalArea(){
    double totalArea = 0.0;
    for(isoQuad4 element : elementList){
        totalArea += element.getArea(nodeList);
    }
    return totalArea;
}