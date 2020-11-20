#include "../include/elements.h"
#include "../include/mesh.h"
#include "../include/math_utilities.h"
// #include "stiffnessMatrix.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <cmath>

node2d::node2d(double x, double y) : _x(x), _y(y){};

node2d::node2d()
{
    std::cout << "Default Constructor Called" << std::endl;
};

// node2d::~node2d(){
//     std::cout << "node2d destructor called" << std::endl;
// };

node3d::node3d(double x, double y, double z) : _x(x), _y(y), _z(z){};

isoQuad4::isoQuad4()
{
    /*
    Iso Parametric 4 node quad element.
    Node numbering is counter clockwise starting from lower left corner as:
                           3 __________ 2
                            |          |
                            |          |
                            |__________|
                            0           1
    Default Constructor; all nodes initialized to 0
    */
    _nodes[0] = 0;
    _nodes[1] = 0;
    _nodes[2] = 0;
    _nodes[3] = 0;
};

isoQuad4::isoQuad4(unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3)
{
        /*
    Iso Parametric 4 node quad element.
    Node numbering is counter clockwise starting from lower left corner as:
                           3 __________ 2
                            |          |
                            |          |
                            |__________|
                            0           1
    Class constructor; nodes indeces set to values passed to constructor
    Inputs n0 - index of node 0 coordinate
           n1 - index of node 1 coordinate
           n2 - index of node 2 coordinate
           n3 - index of node 3 coordinate
    */
    _nodes[0] = n0;
    _nodes[1] = n1;
    _nodes[2] = n2;
    _nodes[3] = n3;
};

isoQuad4::isoQuad4(const isoQuad4 &quad)
{
    /*
    Iso Parametric 4 node quad element.

    Copy Constructor; copies nodal values to new quad object
    */
    cout << "copy constructor called" << endl;
    _nodes[0] = quad._nodes[0];
    _nodes[1] = quad._nodes[1];
    _nodes[2] = quad._nodes[2];
    _nodes[3] = quad._nodes[3];
};

void isoQuad4::assembleLocalStiff(const std::vector<node2d> &nodeList, material &mat, matrix & localStiff)
{
        /*
    Iso Parametric 4 node quad element.
    Method to assemble elemental stiffness matrix for heat conduction.  Uses Full Gauss Quadrature (2x2 points) for integration
    Inputs: nodelist - reference to vector of mesh coordinates
            mat - reference to element material object, holds material thickness and thermal conductivity
            localStiff - reference to matrix object that will store the computed elemental stiffness matrix
    */
    double gqPoints4[4][2] = {{-1 / sqrt(3), -1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)}, {1 / sqrt(3), -1 / sqrt(3)}, {1 / sqrt(3), 1 / sqrt(3)}};
    matrix coords(4,2);
    for (int i = 0; i < 4; i++)
    {
        coords[i][0] = nodeList[_nodes[i]]._x;
        coords[i][1] = nodeList[_nodes[i]]._y;
    }

    double thickness = mat.getThickness();
    double thecon = mat.getTheCon();

    matrix jacobian(2, 2);
    matrix g(2,2);
    matrix w(2,4);
    matrix B(2,4);
    double eta, xi;
    double jj;

    for (int i = 0; i < 4; i++)
    {

        eta = gqPoints4[i][1];
        xi = gqPoints4[i][0];

        //evaluate shape function derivatives matrix
        w[0][0] = -(1.0 - eta)/4.0;
        w[0][1] =  (1.0 - eta)/4.0;
        w[0][2] =  (1.0 + eta)/4.0;
        w[0][3] = -(1.0 + eta)/4.0;
        w[1][0] = -(1.0 - xi)/4.0;
        w[1][1] = -(1.0 + xi)/4.0;
        w[1][2] =  (1.0 + xi)/4.0;
        w[1][3] =  (1.0 - xi)/4.0;

        //Evaluate jacobian matrix
        jacobian.multiply(w, coords);

        //evaluate jacobian determinant
        jj = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];

        //evaluate jacobian inverse matrix

        g[0][0] = jacobian[1][1] / jj;
        g[0][1] = -jacobian[0][1] / jj;
        g[1][0] = -jacobian[1][0] / jj;
        g[1][1] = jacobian[0][0] / jj;

        //evaluate B matrix
        B.multiply(g, w);

        //evaluate B'k*B and add to element stiffness matrix
        for (int m = 0; m < 4; m++)
        {
            for (int k = m; k < 4; k++)
            {
                double sum = 0.0;
                for (int l = 0; l < 2; l++)
                {
                    sum = sum + B[l][m] * B[l][k];
                }
                localStiff[m][k] = localStiff[m][k] + sum * thickness * thecon*jj;
                localStiff[k][m] = localStiff[m][k];
            }
        }
    }

};

void isoQuad4::assembleHeatGeneration(const std::vector<node2d> & nodeList, material & mat, myVector & loads, double q){
    /*
    Iso Parametric 4 node quad element.
    Method to calculate nodal load vector contribution from internal heat generation.  Uses Full Gauss Quadrature
    (2x2 points) for integration.
    Inputs: nodelist - reference to vector of mesh coordinates
            mat - reference to element material object, holds material thickness and thermal conductivity
            loads - reference to vector object that will hold the computed load vector contribution
            q - the value of the internal heat generation
            
    */
    double gqPoints4[4][2] = {{-1 / sqrt(3), -1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)}, {1 / sqrt(3), -1 / sqrt(3)}, {1 / sqrt(3), 1 / sqrt(3)}};
    matrix coords(4,2);
    
    for (int i = 0; i < 4; i++)
    {
        coords[i][0] = nodeList[_nodes[i]]._x;
        coords[i][1] = nodeList[_nodes[i]]._y;
    }

    myVector N = myVector(4);
    matrix jacobian = matrix(2,2);
    matrix w = matrix(2,4);
    
    double localloads[4];
    localloads[0] = 0.0;
    localloads[1] = 0.0;
    localloads[2] = 0.0;
    localloads[3] = 0.0;

    double eta;
    double xi;
    double jj;

    for(int i = 0; i < 4; i++){

        eta = gqPoints4[i][1];
        xi = gqPoints4[i][0];

        N.setValue(0.25*(1-xi)*(1-eta), 0);
        N.setValue(0.25*(1+xi)*(1-eta), 1);
        N.setValue(0.25*(1+xi)*(1+eta), 2);
        N.setValue(0.25*(1-xi)*(1+eta), 3);

        //evaluate shape function derivatives matrix
        w[0][0] = -(1.0 - eta)/4.0;
        w[0][1] =  (1.0 - eta)/4.0;
        w[0][2] =  (1.0 + eta)/4.0;
        w[0][3] = -(1.0 + eta)/4.0;
        w[1][0] = -(1.0 - xi)/4.0;
        w[1][1] = -(1.0 + xi)/4.0;
        w[1][2] =  (1.0 + xi)/4.0;
        w[1][3] =  (1.0 - xi)/4.0;

        jacobian.multiply(w, coords);
        jj = jacobian[0][0]*jacobian[1][1] - jacobian[0][1]*jacobian[1][0];

        localloads[0] += q*N[0]*jj;
        localloads[1] += q*N[1]*jj;
        localloads[2] += q*N[2]*jj;
        localloads[3] += q*N[3]*jj;

    }

    loads.setValue(localloads[0], 0);
    loads.setValue(localloads[1], 1);
    loads.setValue(localloads[2], 2);
    loads.setValue(localloads[3], 3);

};

double isoQuad4::getArea(std::vector<node2d> &nodeList)
{
    /* 
        Iso Parametric 4 node quad element.
        Method to calculate the area of quadrilateral element from its nodal coordinate
        Inputs: nodelist - reference to list of nodal coordinates of the mesh

        Outputs: double - the area of the element

     */
    // std::cout << "Node List Location in getArea: " << &nodeList << std::endl;

    double area = (nodeList[_nodes[0]]._x * nodeList[_nodes[1]]._y - nodeList[_nodes[0]]._y * nodeList[_nodes[1]]._x) +
                  (nodeList[_nodes[1]]._x * nodeList[_nodes[2]]._y - nodeList[_nodes[1]]._y * nodeList[_nodes[2]]._x) +
                  (nodeList[_nodes[2]]._x * nodeList[_nodes[3]]._y - nodeList[_nodes[2]]._y * nodeList[_nodes[3]]._x) +
                  (nodeList[_nodes[3]]._x * nodeList[_nodes[0]]._y - nodeList[_nodes[3]]._y * nodeList[_nodes[0]]._x);

    return fabs(area / 2.0);
};

line2::line2(){
    /* 
        Iso Parametric 2 node line element.

                                0---------1
        Default constructor; initialize node indeces to 0

     */
    _nodes[0] = 0;
    _nodes[1] = 0;
};

line2::line2(unsigned int n1, unsigned int n2){
        /* 
        Iso Parametric 2 node line element.

                                0---------1
        Class constructor; set the node indeces to passed values
        Inputs: n1 - index of node 0 coordinate
                n2 - index of node 1 coordinate

     */
    _nodes[0] = n1;
    _nodes[1] = n2;
};

double line2::getLength(const std::vector<node2d> & nodelist){
        /* 
        Iso Parametric 2 node line element.
        Method to calculate length of line2 element.
        Inputs: nodelist - reference to vector of coordinates in mesh
        Outputs: double - the length of the line2 element

     */
    return sqrt( (nodelist[_nodes[1]]._x - nodelist[_nodes[0]]._x)*(nodelist[_nodes[1]]._x - nodelist[_nodes[0]]._x) + (nodelist[_nodes[1]]._y - nodelist[_nodes[0]]._y)*(nodelist[_nodes[1]]._y - nodelist[_nodes[0]]._y) );
}

material::material()
{
        /* 
        Material class. Holds material properties to be used in analysis.
        Default Constructor: Initialize the material thickness to unit thickness 1m,
                             Initialize the thermal conductivity to 100 W/mK

     */
    _thickness = 1.0;
    _kMat[0][0] = 100.0;
    _kMat[0][1] = 100.0;
    _kMat[1][0] = 100.0;
    _kMat[1][1] = 100.0;
};

material::material(double thickness, double thecon)
{
    /* 
        Material class. Holds material properties to be used in analysis.
        Class Constructor: Initialize the material Properties to passed values
        Inputs: thickness - the thickness of the material in metric units m
                thecon - the thermal conductivity in metrix units, W/mk

     */
    _thickness = thickness;
    _kMat[0][0] = thecon;
    _kMat[0][1] = thecon;
    _kMat[1][0] = thecon;
    _kMat[1][1] = thecon;
};

double material::getThickness()
{
    /* 
        Material class. Holds material properties to be used in analysis.
        Method to return thickness of material
        Output: double - the thickness of the material

     */
    return _thickness;
};

double material::getTheCon()
{   
    /* 
        Material class. Holds material properties to be used in analysis.
        Method to return the thermal conductivity of material.
        Assumes isotropic material and returns one value.
        For anisotropic material corresponds to returning k_xx
        Output: double - thermal conductivity

     */
    return _kMat[0][0];
}