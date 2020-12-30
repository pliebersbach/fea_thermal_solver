#include "../include/elements.h"
#include "../include/mesh.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>


mesh::mesh(){
    cout << "Default Mesh Constructor" << endl;
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
        cerr << "Unable to open file" << std::endl;
        exit(1);
    }
    else{
        cout << "File opened successfully" << endl;
    }

    std::cout << "Begin Reading Mesh" << std::endl;

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

        if (line.find("Nodesets") == 0)
        {
            std::getline(inFile, line);
            unsigned int num_nodesets = stoi(line);
            unsigned int num_nodes;
            for (int i = 0; i < num_nodesets; i++)
            {
                inFile >> num_nodes;
                // std::cout << "nodeset #:" << i << std::endl;
                std::vector<unsigned int> nodeids;
                unsigned int node_id;
                for (int j = 0; j < num_nodes; j++)
                {
                    inFile >> node_id;
                    // std::cout << "node id:" << node_id << std::endl;
                    nodeids.push_back(node_id);
                }

                nodeset.push_back(nodeids);
            }
        }

        if (line.find("Sidesets") == 0)
        {
            std::getline(inFile, line);
            unsigned int num_sidesets = stoi(line);
            unsigned int num_edges;
            for (int i = 0; i < num_sidesets; i++)
            {
                inFile >> num_edges;
                // std::cout << "sideset #:" << i << std::endl;
                std::vector<unsigned int> edgeids;
                unsigned int edge_id;
                for (int j = 0; j < num_edges; j++)
                {
                    inFile >> edge_id;
                    // std::cout << "edge id:" << edge_id << std::endl;
                    edgeids.push_back(edge_id);
                }

                sideset.push_back(edgeids);
            }
        }

    }

    inFile.close();

    std::cout << "Finished Reading Mesh." << std::endl;
    std::cout << std::endl;
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

const std::vector<unsigned int> & mesh::getNodeSet(int i){
    return nodeset[i];
}

const std::vector<unsigned int> & mesh::getSideSet(int i){
    return sideset[i];
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

    os << "Nodeset boundary condition list:" << endl;
    int i = 0;
    for(std::vector<unsigned int> id_list : msh.nodeset){
        os << "List " << i << endl;
        for(unsigned int j : id_list){
            os << j << endl;
        }
        i++;
    }

    os << "Sideset boundary condition list:" << endl;
    i = 0;
    for(std::vector<unsigned int> id_list : msh.sideset){
        os << "List " << i << endl;
        for(unsigned int j : id_list){
            os << j << endl;
        }
        i++;
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