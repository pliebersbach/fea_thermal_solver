#pragma once

#include "elements.h"

#include <vector>
#include <iostream>
#include <cmath>



using namespace std;


class mesh{
    
    private:

    std::vector<node2d> nodeList;
    std::vector<isoQuad4> elementList;
    std::vector<line2> boundaryEleList;

    int dim;
    unsigned int _numnodes;
    unsigned int _numelements;
    unsigned int _numBoundaryEles;
    
    public:
        // std::vector<node2d> nodeList;
        // std::vector<isoQuad4> elementList;
        // int dim;
        // unsigned int _numnodes;
        // unsigned int _numelements;
        mesh();

        mesh(unsigned int nnodes, unsigned int nelements);

        mesh(unsigned int nnodes, unsigned int nelements, unsigned int nboundaryeles);


        //mesh(unsigned int nnodes, unsigned int nelements, int dim);
        void read_mesh(string name, unsigned int dim);

        void addElement(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4);

        void addNode(double x, double y);

        void addBoundaryElement(unsigned int n1, unsigned int n2);

        isoQuad4 getElement(unsigned int i);

        line2 getBoundaryElement(unsigned int i);

        const std::vector<isoQuad4> & getElementList();

        const std::vector<node2d> & getNodeList();

        const std::vector<line2> & getBoundaryElementList();

        unsigned int getNumEles();

        unsigned int getNumNodes();

        double getTotalArea();

        friend ostream &operator<<(ostream &os, const mesh &msh);

        // ~mesh();
};