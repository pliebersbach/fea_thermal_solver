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
    _nodes[0] = 0;
    _nodes[1] = 0;
    _nodes[2] = 0;
    _nodes[3] = 0;
};

isoQuad4::isoQuad4(unsigned int n0, unsigned int n1, unsigned int n2, unsigned int n3)
{
    _nodes[0] = n0;
    _nodes[1] = n1;
    _nodes[2] = n2;
    _nodes[3] = n3;
};

isoQuad4::isoQuad4(const isoQuad4 &quad)
{
    cout << "copy constructor called" << endl;
    _nodes[0] = quad._nodes[0];
    _nodes[1] = quad._nodes[1];
    _nodes[2] = quad._nodes[2];
    _nodes[3] = quad._nodes[3];
}
// isoQuad4::~isoQuad4(){
//     std::cout << "isoQuad4 destructor called" << std::endl;
// };
void isoQuad4::assembleLocalStiff(const std::vector<node2d> &nodeList, material &mat, matrix & localStiff)
{

    double gqPoints4[4][2] = {{-1 / sqrt(3), -1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)}, {1 / sqrt(3), -1 / sqrt(3)}, {1 / sqrt(3), 1 / sqrt(3)}};
    matrix coords(4,2);
    for (int i = 0; i < 4; i++)
    {
        coords[i][0] = nodeList[_nodes[i]]._x;
        coords[i][1] = nodeList[_nodes[i]]._y;
    }
    // cout << "element coordinates" << endl;
    // cout << coords;
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

        // cout << "w matrix" << endl;

        // cout << w;
        //Evaluate jacobian matrix
        jacobian.multiply(w, coords);

        // cout << "jacobian matrix: " << endl;
        // cout << jacobian;
        //evaluate jacobian determinant
        jj = jacobian[0][0] * jacobian[1][1] - jacobian[0][1] * jacobian[1][0];

        //evaluate jacobian inverse matrix

        g[0][0] = jacobian[1][1] / jj;
        g[0][1] = -jacobian[0][1] / jj;
        g[1][0] = -jacobian[1][0] / jj;
        g[1][1] = jacobian[0][0] / jj;
        // cout << "jacobian inverse" << endl;
        // cout << g;
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
    // cout << localStiff;
    // for (int i = 0; i < 4; i++)
    // {
    //     std::cout << localStiff[i][0] << ", " << localStiff[i][1] << ", " << localStiff[i][2] << ", " << localStiff[i][3] << std::endl;
    // };
    
    //assemble local stiffness matrix into global stiffness matrix
    // global.addStiffness(&localStiff, *this);
};

void isoQuad4::assembleHeatGeneration(const std::vector<node2d> & nodeList, material & mat, myVector & loads, double q){
    
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

    // double q = 6.0;
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

    std::cout << "Loads from internal heat generation" << std::endl;
    std::cout << loads;

};
// void isoQuad4::assembleLocalStiff(const std::vector<node2d> &nodeList, material &mat)
// {

//     double gqPoints4[4][2] = {{-1 / sqrt(3), -1 / sqrt(3)}, {-1 / sqrt(3), 1 / sqrt(3)}, {1 / sqrt(3), -1 / sqrt(3)}, {1 / sqrt(3), 1 / sqrt(3)}};
//     double coords[4][2];
//     for (int i = 0; i < 4; i++)
//     {
//         coords[i][0] = nodeList[_nodes[i]]._x;
//         coords[i][1] = nodeList[_nodes[i]]._y;
//     }

//     double thickness = mat.getThickness();
//     double thecon = mat.getTheCon();

//     double localStiff[4][4];
//     double j[2][2];
//     double g[2][2];
//     double w[2][4];
//     double B[2][4];
//     double eta, xi;
//     double jj;

//     for (int i = 0; i < 4; i++)
//     {

//         eta = gqPoints4[i][1];
//         xi = gqPoints4[i][0];

//         //evaluate shape function derivatives matrix
//         w[0][0] = -(1.0 - eta)/4.0;
//         w[0][1] =  (1.0 - eta)/4.0;
//         w[0][2] =  (1.0 + eta)/4.0;
//         w[0][3] = -(1.0 + eta)/4.0;
//         w[1][0] = -(1.0 - xi)/4.0;
//         w[1][1] = -(1.0 + xi)/4.0;
//         w[1][2] =  (1.0 + xi)/4.0;
//         w[1][3] =  (1.0 - xi)/4.0;

//         //Evaluate jacobian matrix

//         j[0][0] = w[0][0] * coords[0][0] + w[0][1] * coords[1][0] + w[0][2] * coords[2][0] + w[0][3] * coords[3][0];
//         j[0][1] = w[0][0] * coords[0][1] + w[0][1] * coords[1][1] + w[0][2] * coords[2][1] + w[0][3] * coords[3][1];
//         j[1][0] = w[1][0] * coords[0][0] + w[1][1] * coords[1][0] + w[1][2] * coords[2][0] + w[1][3] * coords[3][0];
//         j[1][1] = w[1][0] * coords[0][1] + w[1][1] * coords[1][1] + w[1][2] * coords[2][1] + w[1][3] * coords[3][1];

//         //evaluate jacobian determinant
//         jj = j[0][0] * j[1][1] - j[0][1] * j[1][0];

//         //evaluate jacobian inverse matrix

//         g[0][0] = j[1][1] / jj;
//         g[0][1] = -j[0][1] / jj;
//         g[1][0] = -j[1][0] / jj;
//         g[1][1] = j[0][0] / jj;

//         //evaluate B matrix

//         B[0][0] = g[0][0] * w[0][0] + g[0][1] * w[1][0];
//         B[0][1] = g[0][0] * w[0][1] + g[0][1] * w[1][1];
//         B[0][2] = g[0][0] * w[0][2] + g[0][1] * w[1][2];
//         B[0][3] = g[0][0] * w[0][3] + g[0][1] * w[1][3];
//         B[1][0] = g[1][0] * w[0][0] + g[1][1] * w[1][0];
//         B[1][1] = g[1][0] * w[0][1] + g[1][1] * w[1][1];
//         B[1][2] = g[1][0] * w[0][2] + g[1][1] * w[1][2];
//         B[1][3] = g[1][0] * w[0][3] + g[1][1] * w[1][3];

//         //evaluate B'k*B and add to element stiffness matrix
//         for (int m = 0; m < 4; m++)
//         {
//             for (int k = m; k < 4; k++)
//             {
//                 double sum = 0.0;
//                 for (int l = 0; l < 2; l++)
//                 {
//                     sum = sum + B[l][m] * B[l][k];
//                 }
//                 localStiff[m][k] = localStiff[m][k] + sum * thickness * thecon*jj;
//                 localStiff[k][m] = localStiff[m][k];
//             }
//         }
//     }

//     for (int i = 0; i < 4; i++)
//     {
//         std::cout << localStiff[i][0] << ", " << localStiff[i][1] << ", " << localStiff[i][2] << ", " << localStiff[i][3] << std::endl;
//     };
    
//     //assemble local stiffness matrix into global stiffness matrix
//     // global.addStiffness(&localStiff, *this);
// };

double isoQuad4::getArea(std::vector<node2d> &nodeList)
{

    std::cout << "Node List Location in getArea: " << &nodeList << std::endl;

    double area = (nodeList[_nodes[0]]._x * nodeList[_nodes[1]]._y - nodeList[_nodes[0]]._y * nodeList[_nodes[1]]._x) +
                  (nodeList[_nodes[1]]._x * nodeList[_nodes[2]]._y - nodeList[_nodes[1]]._y * nodeList[_nodes[2]]._x) +
                  (nodeList[_nodes[2]]._x * nodeList[_nodes[3]]._y - nodeList[_nodes[2]]._y * nodeList[_nodes[3]]._x) +
                  (nodeList[_nodes[3]]._x * nodeList[_nodes[0]]._y - nodeList[_nodes[3]]._y * nodeList[_nodes[0]]._x);

    return fabs(area / 2.0);
};

line2::line2(){
    _nodes[0] = 0;
    _nodes[1] = 0;
};

line2::line2(unsigned int n1, unsigned int n2){
    _nodes[0] = n1;
    _nodes[1] = n2;
};

double line2::getLength(const std::vector<node2d> & nodelist){
    return sqrt( (nodelist[_nodes[1]]._x - nodelist[_nodes[0]]._x)*(nodelist[_nodes[1]]._x - nodelist[_nodes[0]]._x) + (nodelist[_nodes[1]]._y - nodelist[_nodes[0]]._y)*(nodelist[_nodes[1]]._y - nodelist[_nodes[0]]._y) );
}

material::material()
{
    _thickness = 1.0;
    _kMat[0][0] = 100.0;
    _kMat[0][1] = 100.0;
    _kMat[1][0] = 100.0;
    _kMat[1][1] = 100.0;
};

material::material(double thickness, double thecon)
{
    _thickness = thickness;
    _kMat[0][0] = thecon;
    _kMat[0][1] = thecon;
    _kMat[1][0] = thecon;
    _kMat[1][1] = thecon;
};

double material::getThickness()
{
    return _thickness;
};

double material::getTheCon()
{
    return _kMat[0][0];
}