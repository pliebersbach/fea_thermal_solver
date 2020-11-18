#include "math_utilities.h"

#include <iostream>
#include <vector>
using namespace std;

matrix::matrix(){
    _matrix = nullptr;
    _rows = 0;
    _cols = 0;
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
};

void matrix::allocate(unsigned int size){
    _matrix = new double[size*size]();
    _rows = size;
    _cols = size;
};

void matrix::allocate(unsigned int rows, unsigned int cols){
    _matrix = new double[rows*cols]();
    _rows = rows;
    _cols = cols;
};

matrix::~matrix(){
    std::cout << "Matrix destroyed" << std::endl;
    delete[] _matrix;
    _matrix = nullptr;
};

double* matrix::operator[](int i){
    return &_matrix[i*_cols];
};

unsigned int matrix::getCols(){
    return _cols;
};

unsigned int matrix::getRows(){
    return _rows;
};

ostream & operator<<(ostream &os, matrix &mat){
    os << "Print matrix:" << std::endl;
    for(int i = 0; i < mat._rows; i++){
        for(int j = 0; j < mat._cols; j++){
            os << mat[i][j] << ", ";
        }
        os << std::endl;
    }
    // for(int i = 0; i < mat._rows*mat._cols; i++){
    //     os << mat._matrix[i] << std::endl;
    // }
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
            for(int k = 0; k < a._cols; k++){
                sum += a[i][k] * b[k][j];
            }
            (*this)[i][j] = sum;
        }
    }
};

void matrix::gaussElimination(double * rhs, double * x){
    double l = 0.0;

    //Forward Elimination
    for(int k = 0; k < _rows - 1; k++){

        for(int i = k + 1; i < _rows; i++){
            l = (*this)[i][k]/(*this)[k][k];

            for(int j = k; j < _rows; j++){
                (*this)[i][j] = (*this)[i][j] - l*(*this)[k][j];
            }

            rhs[i] = rhs[i] - l*rhs[k];
        }
    }

    //Back Substitution

    for(int k = _rows - 1; k>-1; k--){
        x[k] = rhs[k];

        for(int j = k+1; j < _rows; j++){
            x[k] = x[k] - (*this)[k][j]*x[j];
        }
        x[k] = x[k]/(*this)[k][k];
    }
};

void matrix::gaussElimination(myVector & rhs, myVector & x){
    double l = 0.0;

    //Forward Elimination
    for(int k = 0; k < _rows - 1; k++){

        for(int i = k + 1; i < _rows; i++){
            l = (*this)[i][k]/(*this)[k][k];

            for(int j = k; j < _rows; j++){
                (*this)[i][j] = (*this)[i][j] - l*(*this)[k][j];
            }

            // rhs[i] = rhs[i] - l*rhs[k];
            rhs.setValue( rhs[i] - l*rhs[k] , i);
        }
    }

    //Back Substitution

    for(int k = _rows - 1; k>-1; k--){
        // x[k] = rhs[k];
        x.setValue( rhs[k], k);

        for(int j = k+1; j < _rows; j++){
            // x[k] = x[k] - (*this)[k][j]*x[j];
            x.setValue( x[k] - (*this)[k][j]*x[j] , k);
        }
        // x[k] = x[k]/(*this)[k][k];
        x.setValue( x[k]/(*this)[k][k] , k);
    }
};

myVector::myVector(){
    _vector = nullptr;
    _size = 0;
};

myVector::myVector(unsigned int size){
    _vector = new double[size]();
    _size = size;
};

void myVector::allocate(unsigned int size){
    _vector = new double[size]();
    _size = size;
}

myVector::~myVector(){
    std::cout << "Vector deleted" << std::endl;
    delete[] _vector;
    _vector = nullptr;
}

double myVector::operator[](int i){
    return _vector[i];
};

void myVector::setValue(double value, unsigned int i){
    if(i>_size){
        std::cout << "Index out of bounds.  Operation aborted" << std::endl;
        return;
    }
    _vector[i] = value;
};

double * myVector::getVector(){
    return _vector;
}

ostream & operator<<(ostream &os, myVector & vec){
    os << "Print vector:" << std::endl;
    for(int i = 0; i < vec._size; i++){
            os << vec._vector[i] << std::endl;
    }
    return os;
};

void myVector::vector_multiply(matrix & a, myVector & b){
    if(a.getCols() != b._size){
        std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
        return;
    }

    if(_size != a.getRows()){
        std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
        return;        
    }

    for(int i = 0; i < _size; i++){
        
        for(int j = 0; j < b._size; j++){
            _vector[i] += a[i][j] * b[j];
        }
    }
};

void myVector::vector_subtract(myVector & a, myVector & b){
    if(a._size != b._size && a._size != _size){
        std::cout << "Error: Dimension Mismatch.  Subtract aborted" << std::endl;
        return;
    }
    for(int i = 0; i < _size; i++){
        _vector[i] = b[i] - a[i];

    }
};



// void vector_multiply(matrix & a, std::vector<double> & b, std::vector<double> & result){
//     if(a.getCols() != b.size()){
//         std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
//         return;
//     }

//     if(result.size() != a.getRows()){
//         std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
//         return;        
//     }

//     for(int i = 0; i < result.size(); i++){
        
//         for(int j = 0; j < b.size(); j++){
//             result[i] += a[i][j] * b[j];
//         }
//     }
// };

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