#include "math_utilities.h"

#include <iostream>
using namespace std;

matrix::matrix(){
    _matrix = new  double[4];
    _matrix[0] = 1.0;
    _matrix[1] = 0.0;
    _matrix[2] = 0.0;
    _matrix[3] = 1.0;
    _rows = 2;
    _cols = 2;
};

matrix::matrix(unsigned int size){
    _matrix = new double[size*size]();
    _rows = size;
    _cols = size;
};

matrix::matrix(unsigned int rows, unsigned int cols){
    _matrix = new double[rows*cols]();
    _rows = rows;
    _cols = cols;
}

matrix::~matrix(){
    std::cout << "Matrix destroyed" << std::endl;
    delete[] _matrix;
    _matrix = nullptr;
};

double* matrix::operator[](int i){
    return &_matrix[i*_rows];
};

ostream & operator<<(ostream &os, matrix &mat){
    os << "Print matrix:" << std::endl;
    for(int i = 0; i < mat._rows; i++){
        for(int j = 0; j < mat._cols; j++){
            os << mat[i][j] << ", ";
        }
        os << std::endl;
    }
    return os;
};

void matrix::multiply(matrix & a, matrix & b){
    if(a._cols != b._rows){
        std::cout << "Error: Dimension Mismatch" << std::endl;
        return;
    };

    for(int i = 0; i < _rows; i++){
        for(int j = 0; j < _cols; j++){
            double sum = 0.0;
            for(int k = 0; k < b._cols; k++){
                sum += a[i][k] * b[k][j];
            }
            (*this)[i][j] = sum;
        }
    }
};

// matrix & matrix::operator*(matrix & mat){
//     if(_cols != mat._rows){
//         std::cout << "Error: Dimension Mismatch" << std::endl;
//     };

//     matrix result = matrix(_rows, mat._cols);

//     for(int i = 0; i < _rows; i++){
//         for(int j = 0; j < mat._cols; j++){
//             double sum = 0.0;
//             for(int k = 0; k < _cols; k++){
//                 sum += (*this)[i][k] * mat[k][j];
//             }
//             result[i][j] = sum;
//         }
//     }
//     return result;
// }