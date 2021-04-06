#pragma once

#include <iostream>
#include <vector>

using namespace std;
class myVector;
class matrix;

double vectorDotProduct(double *vec1, double *vec2, int size);

void matrixVectorProduct(matrix &mat, double *vec, double *res);

double l2norm(double *vec, int size);

double l2norm(double *vec1, double *vec2, int size);


class matrix{

    public:

    double* myData;
    unsigned int myRows;
    unsigned int myCols;

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

    void jacobiIteration(myVector & rhs, myVector & soln, int max_iterations);

    void gaussSeidel(myVector & rhs, myVector & soln, int max_iterations);

    void sor(myVector & rhs, myVector & soln, int max_iterations);

    void conjugateGradient(myVector & rhs, myVector & soln, int max_iterations);

    void printMatrix(string fileName);
};

class myVector{

    public:

    double* myData;
    unsigned int mySize;
    
    myVector();
    
    myVector(unsigned int size);
    
    void allocate(unsigned int size);
    
    ~myVector();
    
    double operator[](int i);
    
    void setValue(double value, unsigned int i);
    
    void addValue(double value, unsigned int i);

    void vector_multiply(matrix & a, myVector & b);
    
    void vector_subtract(myVector & a, myVector & b);

    double l2Norm(myVector & a);

    double * getDataPointer();
    
    friend ostream &operator<<(ostream &os, myVector &vec);

    void printVector(string fileName);

};
// void vector_multiply(matrix & a, std::vector<double> & b,  std::vector<double> &result);