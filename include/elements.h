#pragma once
#include "math_utilities.h"

#include <vector>
#include <iostream>
#include <cmath>

class node2d{
    
    public:

    node2d();

    node2d(double x, double y);

    // ~node2d();
    
    double _x;
    double _y;
};

class node3d{
    public:

    node3d(double x, double y, double z);

    // ~node3d();

    double _x;
    double _y;
    double _z;
};

class material{
    private:
    double _thickness;
    double _kMat[2][2];

    public:
    material();
    material(double thickness, double thecon);
    double getThickness();
    double getTheCon();
};

class isoQuad4{

    private:
    
    public:

    unsigned int _nodes[4];

    isoQuad4();

    isoQuad4(unsigned int n1, unsigned int n2, unsigned int n3, unsigned int n4);

    isoQuad4(const isoQuad4 &quad);
    // ~isoQuad4();

    void assembleLocalStiff(const std::vector<node2d> & nodeList, material &mat, matrix & localStiff);

    void assembleHeatGeneration(const std::vector<node2d> & nodeList, material &mat, myVector& loads, double q);
    
    double getArea(std::vector<node2d> & nodeList);


};

class line2{
    
    private:

    public:

    unsigned int _nodes[2];

    line2();

    line2(unsigned int n1, unsigned int n2);

    double getLength(const std::vector<node2d> & nodelist);

    // ~line2();

    
};