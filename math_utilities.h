#pragma once

#include <iostream>

using namespace std;

class matrix{
    private:
    double* _matrix = nullptr;
    unsigned int _rows;
    unsigned int _cols;

    public:
    matrix();
    matrix(unsigned int size);
    matrix(unsigned int rows, unsigned int cols);
    ~matrix();

    double* operator[](int i);

    // matrix & operator*(matrix & mat);

    void multiply(matrix & a, matrix & b);

    friend ostream & operator<<(ostream &os, matrix &mat);
};