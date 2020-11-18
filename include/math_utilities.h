#pragma once

#include <iostream>
#include <vector>

using namespace std;

class myVector;

class matrix{
    private:
    double* _matrix;
    unsigned int _rows;
    unsigned int _cols;

    public:
    matrix();
    matrix(unsigned int size);
    matrix(unsigned int rows, unsigned int cols);
    void allocate(unsigned int size);
    void allocate(unsigned int rows, unsigned int cols);
    ~matrix();

    double* operator[](int i);

    unsigned int getRows();

    unsigned int getCols();

    // matrix & operator*(matrix & mat);

    void multiply(matrix & a, matrix & b);

    friend ostream & operator<<(ostream &os, matrix &mat);

    void gaussElimination(double * rhs, double * soln);

    void gaussElimination(myVector & rhs, myVector & soln);
};

class myVector{
    private:    
    double* _vector;
    unsigned int _size;

    public:    

    myVector();
    
    myVector(unsigned int size);
    
    void allocate(unsigned int size);
    
    ~myVector();
    
    double operator[](int i);
    
    void setValue(double value, unsigned int i);
    
    void addValue(double value, unsigned int i);

    void vector_multiply(matrix & a, myVector & b);
    
    void vector_subtract(myVector & a, myVector & b);
    
    double * getVector();
    
    friend ostream &operator<<(ostream &os, myVector &vec);
};
// void vector_multiply(matrix & a, std::vector<double> & b,  std::vector<double> &result);