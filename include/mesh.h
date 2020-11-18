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
    int dim;
    unsigned int _numnodes;
    unsigned int _numelements;
    
    public:
        // std::vector<node2d> nodeList;
        // std::vector<isoQuad4> elementList;
        // int dim;
        // unsigned int _numnodes;
        // unsigned int _numelements;

        mesh(unsigned int nnodes, unsigned int nelements);

        //mesh(unsigned int nnodes, unsigned int nelements, int dim);

        void addElement(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4);

        void addNode(double x, double y);

        isoQuad4 getElement(unsigned int i);

        const std::vector<isoQuad4> & getElementList();

        const std::vector<node2d> & getNodeList();

        unsigned int getNumEles();

        double getTotalArea();

        friend ostream &operator<<(ostream &os, const mesh &msh);

        // ~mesh();
};