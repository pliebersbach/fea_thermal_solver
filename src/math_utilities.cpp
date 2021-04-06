#include "../include/math_utilities.h"

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
using namespace std;

double vectorDotProduct(double *vec1, double *vec2, int size){
    double product = 0.0;
    for(int i = 0; i < size; i++){
        product += vec1[i]*vec2[i];
    }
    return product;
};

void matrixVectorProduct(matrix &mat, double *vec, double *result){
    
    for(int i = 0; i < mat.myRows; i++){
        result[i] = 0.0;
        for(int j = 0; j < mat.myCols; j++){
            result[i] += mat[i][j]*vec[j];
        }
    }

};

double l2norm(double *vec, int size){
    double norm = 0.0;

    for(int i = 0; i < size; i++){
        norm += vec[i]*vec[i];
    }

    return sqrt(norm);

};

double l2norm(double *vec1, double *vec2, int size){
    double norm = 0.0;

    for(int i = 0; i < size; i++){
        norm += (vec1[i]-vec2[i])*(vec1[i]-vec2[i]);
    }

    return sqrt(norm);

};

matrix::matrix(){
    myData = nullptr;
    myRows = 0;
    myCols = 0;
};

matrix::matrix(unsigned int size){
    myData = new double[size*size]();
    myRows = size;
    myCols = size;
};

matrix::matrix(unsigned int rows, unsigned int cols){
    myData = new double[rows*cols]();
    myRows = rows;
    myCols = cols;
};

void matrix::allocate(unsigned int size){
    myData = new double[size*size]();
    myRows = size;
    myCols = size;
};

void matrix::allocate(unsigned int rows, unsigned int cols){
    myData = new double[rows*cols]();
    myRows = rows;
    myCols = cols;
};

matrix::~matrix(){
    // std::cout << "Matrix destroyed" << std::endl;
    delete[] myData;
    myData = nullptr;
};

double* matrix::operator[](int i){
    return &myData[i*myCols];
};

unsigned int matrix::getCols(){
    return myCols;
};

unsigned int matrix::getRows(){
    return myRows;
};

ostream & operator<<(ostream &os, matrix &mat){
    os << "Print matrix:" << std::endl;
    for(int i = 0; i < mat.myRows; i++){
        for(int j = 0; j < mat.myCols; j++){
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
    if(a.myCols != b.myRows){
        std::cout << "Error: Dimension Mismatch" << std::endl;
        return;
    };

    for(int i = 0; i < myRows; i++){
        for(int j = 0; j < myCols; j++){
            double sum = 0.0;
            for(int k = 0; k < a.myCols; k++){
                sum += a[i][k] * b[k][j];
            }
            (*this)[i][j] = sum;
        }
    }
};

void matrix::gaussElimination(double * rhs, double * x){
    double l = 0.0;

    //Forward Elimination
    for(int k = 0; k < myRows - 1; k++){

        for(int i = k + 1; i < myRows; i++){
            l = (*this)[i][k]/(*this)[k][k];

            for(int j = k; j < myRows; j++){
                (*this)[i][j] = (*this)[i][j] - l*(*this)[k][j];
            }

            rhs[i] = rhs[i] - l*rhs[k];
        }
    }

    //Back Substitution

    for(int k = myRows - 1; k>-1; k--){
        x[k] = rhs[k];

        for(int j = k+1; j < myRows; j++){
            x[k] = x[k] - (*this)[k][j]*x[j];
        }
        x[k] = x[k]/(*this)[k][k];
    }
};

void matrix::gaussElimination(myVector & rhs, myVector & x){
    if(myRows != rhs.mySize){
        cerr << "Matrix - Vector Dimensions Mismatched in Gauss Solver. Operation Aborted" << endl;
        exit(1);
    }
    double l = 0.0;
    double *xx  = x.myData; //use xx as pointer to vector data in x
    double *rhsvec = rhs.myData; //use rhsvec as pointer to vector data in rhs

    //Forward Elimination
    for(int k = 0; k < myRows - 1; k++){

        for(int i = k + 1; i < myRows; i++){
            l = (*this)[i][k]/(*this)[k][k];

            for(int j = k; j < myRows; j++){
                (*this)[i][j] = (*this)[i][j] - l*(*this)[k][j];
            }

            // rhs[i] = rhs[i] - l*rhs[k];
            // rhs.myData[i] = rhs[i] - l*rhs[k];
            rhsvec[i] = rhsvec[i] - l*rhsvec[k];

        }
    }

    //Back Substitution

    for(int k = myRows - 1; k>-1; k--){
        // x[k] = rhs[k];
        // x.myData[k] = rhs[k];
        xx[k] = rhsvec[k];

        for(int j = k+1; j < myRows; j++){
            // x[k] = x[k] - (*this)[k][j]*x[j];
            // x.myData[k] = x[k] - (*this)[k][j]*x[j];
            xx[k] = xx[k] - (*this)[k][j]*xx[j];

        }
        // x[k] = x[k]/(*this)[k][k];
        // x.myData[k] = x[k]/(*this)[k][k];
        xx[k] = xx[k]/(*this)[k][k];

    }
};

void matrix::jacobiIteration(myVector &rhs, myVector &soln, int max_iterations)
{
    cout << "Jacobi iteration solver" << endl;
    myVector intermediate_soln = myVector(rhs.mySize);
    double *xk1 = intermediate_soln.myData;
    double *xk = soln.myData;
    double *temp = nullptr;

    double tol = 1e-05;
    double error = 1.0;
    int k = 0;
    // int max_iterations = 25;

    // cout << "intermediate solution location: " << intermediate_soln.myData << " pointer to it: " << xk1 << endl;
    // cout << "solution location: " << soln.myData << " pointer to it: " << xk << endl;

    while (error > tol && k < max_iterations){    //iteration loop
        k++;
        for (int i = 0; i < myRows; i++){   //iterate rows

            double sum = 0.0;
            for (int j = 0; j < myCols; j++)    //iterate cols
            {
                if(j != i){
                    sum += (*this)[i][j] * xk[j];
                }
            }
            xk1[i] = (1 / (*this)[i][i]) * (rhs[i] - sum);
        }

        error = soln.l2Norm(intermediate_soln);

        // for(int i = 0; i < myRows; i++){
        //     xk[i] = xk1[i];
        // }

        temp = soln.myData;
        soln.myData = intermediate_soln.myData;
        intermediate_soln.myData = temp;
        xk = soln.myData;
        xk1 = intermediate_soln.myData;

        cout <<"Iteration # " << k << " Error: " << error << endl;
        // cout << soln;

        // // xk = intermediate_soln.myData;
        // // xk1 = soln.myData;
        // soln.myData = xk1;
        // cout << soln;
        // cout << intermediate_soln;

    }
    cout << "Final Error: " << error << endl;
    cout << "Total Iterations: " << k << endl;
    cout << "Final Solution: " << endl;
    // cout << soln;
};

void matrix::gaussSeidel(myVector &rhs, myVector &soln, int max_iterations)
{
    cout << "Gauss Seidel iteration solver" << endl;
    myVector intermediate_soln = myVector(rhs.mySize);
    double *xk1 = intermediate_soln.myData;
    double *xk = soln.myData;
    double *temp = nullptr;
    double sum1;
    double sum2;

    double tol = 1e-05;
    double error = 1.0;
    int k = 0;
    // int max_iterations = 4;

    // cout << "intermediate solution location: " << intermediate_soln.myData << " pointer to it: " << xk1 << endl;
    // cout << "solution location: " << soln.myData << " pointer to it: " << xk << endl;

    // while (error > tol && k < max_iterations){    //iteration loop
    while(error > tol){
        k++;
        for (int i = 0; i < myRows; i++){   //iterate rows

            sum1 = 0.0;
            sum2 = 0.0;

            for (int j = 0; j < i; j++)    //iterate cols
            {
                sum1 += (*this)[i][j] * xk1[j];
                
            }

            for(int j = i+1; j < myCols; j++){
                sum2 += (*this)[i][j] * xk[j];
            }

            xk1[i] = (1 / (*this)[i][i]) * (rhs[i] - sum1 - sum2);
        }

        error = soln.l2Norm(intermediate_soln);

        // for(int i = 0; i < myRows; i++){
        //     xk[i] = xk1[i];
        // }

        temp = soln.myData;
        soln.myData = intermediate_soln.myData;
        intermediate_soln.myData = temp;
        xk = soln.myData;
        xk1 = intermediate_soln.myData;

        cout <<"Iteration #, " << k << ", Error:, " << error << endl;
        // cout << soln;

        // // xk = intermediate_soln.myData;
        // // xk1 = soln.myData;
        // soln.myData = xk1;
        // cout << soln;
        // cout << intermediate_soln;

    }
    cout << "Final Error: " << error << endl;
    cout << "Total Iterations: " << k << endl;
    cout << "Final Solution: " << endl;
    // cout << soln;
};

void matrix::sor(myVector &rhs, myVector &soln, int max_iterations)
{
    cout << "SOR Iterative solver" << endl;
    myVector intermediate_soln = myVector(rhs.mySize);
    double *xk1 = intermediate_soln.myData;
    double *xk = soln.myData;
    double *temp = nullptr;
    double sum1;
    double sum2;
    double w = 1.75;

    double tol = 1e-05;
    double error = 1.0;
    int k = 0;
    // int max_iterations = 4;

    // cout << "intermediate solution location: " << intermediate_soln.myData << " pointer to it: " << xk1 << endl;
    // cout << "solution location: " << soln.myData << " pointer to it: " << xk << endl;

    // while (error > tol && k < max_iterations){    //iteration loop
    while (error > tol){    //iteration loop

        k++;
        for (int i = 0; i < myRows; i++){   //iterate rows

            sum1 = 0.0;
            sum2 = 0.0;

            for (int j = 0; j < i; j++)    //iterate cols
            {
                sum1 += (*this)[i][j] * xk1[j];
                
            }

            for(int j = i+1; j < myCols; j++){
                sum2 += (*this)[i][j] * xk[j];
            }

            xk1[i] = (1-w)*xk[i] + (w / (*this)[i][i]) * (rhs[i] - sum1 - sum2);
        }

        error = soln.l2Norm(intermediate_soln);

        // for(int i = 0; i < myRows; i++){
        //     xk[i] = xk1[i];
        // }

        temp = soln.myData;
        soln.myData = intermediate_soln.myData;
        intermediate_soln.myData = temp;
        xk = soln.myData;
        xk1 = intermediate_soln.myData;

        cout <<"Iteration #, " << k << ", Error:, " << error << endl;
        // cout << soln;

        // // xk = intermediate_soln.myData;
        // // xk1 = soln.myData;
        // soln.myData = xk1;
        // cout << soln;
        // cout << intermediate_soln;

    }
    cout << "Final Error: " << error << endl;
    cout << "Total Iterations: " << k << endl;
    cout << "Final Solution: " << endl;
    // cout << soln;
};

void matrix::conjugateGradient(myVector & rhs, myVector & soln, int max_iterations){
    cout << "Conjugate Gradient Solver" << endl;
    
    double *res = new double[myRows];
    double *p = new double[myRows];
    double *v = new double[myRows];
    double *resnew = new double[myRows];
    double *pnew = new double[myRows];
    double *vnew = new double[myRows];
    double *xnew = new double[myRows];

    double *temp = nullptr;
    
    double alpha, beta, error, residual, alpha1, alpha2, beta1, beta2;
    double tol = 1e-05;
    
    matrixVectorProduct((*this), soln.myData, v);

    for(int i = 0; i < myRows; i++){
        p[i] = rhs[i] - v[i];
        res[i] = p[i];
        v[i] = 0.0;
    }
    int k;
    for(k = 0; k < max_iterations; k++){

        // cout << "Iteration #: " << k << endl;

        // cout << "Residual: " << endl;

        // for(int j = 0; j < myRows; j++){
        //     cout << res[j] << endl;
        // }

        matrixVectorProduct((*this), p, v);

        // cout << "vector v: " << endl;
        // for(int j = 0; j < myRows; j++){
        //     cout << v[j] << endl;
        // }

        // cout << "vector p:" << endl;
        // for(int j = 0; j < myRows; j++){
        //     cout << p[j] << endl;
        // }

        alpha1 = vectorDotProduct(res, res, myRows);
        alpha2 = vectorDotProduct(p, v, myRows);

        // cout << "alpha1: " << alpha1 << endl;
        // cout << "alpha2: " << alpha2 << endl;
        alpha = alpha1/alpha2;
        // cout << "alpha: " << alpha << endl;

        // cout << "vector xnew | resnew: " << endl;
        for(int j = 0; j < soln.mySize; j++){
            xnew[j] = soln[j] + alpha*p[j];
            resnew[j] = res[j] - alpha*v[j];
            // cout << xnew[j] << " | " << resnew[j] << endl;
        }
        
        error = l2norm(xnew, soln.myData, myRows);

        residual = l2norm(resnew, myRows);

        cout <<"Iteration #, " << k << ", Error:, " << error << ", Residual:, " << residual << endl;
        if (error < tol)
        {
            temp = soln.myData;
            soln.myData = xnew;
            xnew = temp;
            break;
        }
        beta1 = vectorDotProduct(resnew, resnew, myRows);
        beta2 = vectorDotProduct(res, res, myRows);
        beta = beta1/beta2;
        // cout << "beta: " << beta << endl;
        for(int j = 0; j < myRows; j++){
            pnew[j] = resnew[j] + beta*p[j];
        }

        temp = res;
        res = resnew;
        resnew = temp;

        temp = p;
        p = pnew;
        pnew = temp;

        temp = v;
        v = vnew;
        vnew = temp;

        temp = soln.myData;
        soln.myData = xnew;
        xnew = temp;

    }

    delete[] res;
    delete[] p;
    delete[] v;
    delete[] resnew;
    delete[] pnew;
    delete[] vnew;
    delete[] xnew;

    res = nullptr;
    p = nullptr;
    v = nullptr;
    resnew = nullptr;
    pnew = nullptr;
    vnew = nullptr;
    xnew = nullptr;
    temp = nullptr;
    
    cout << "Final Error: " << error << endl;
    cout << "Total Iterations: " << k << endl;
    cout << "Final Solution: " << endl;
    // cout << soln;
};

void matrix::printMatrix(string fileName){
    std::cout << "Writing Matrix to File." << std::endl;
    std::cout << "..." << std::endl;

    ofstream soln;
    soln.open(fileName+".mat", ios::trunc);
    for(int i = 0; i < myRows; i++){
        for(int j = 0; j < myCols; j++){
            soln << (*this)[i][j] << ", ";
        }
        soln << endl;
    }
    soln.close();

    std::cout << "Done writing solution files." << std::endl;
    std::cout << std::endl;
};

myVector::myVector(){
    myData = nullptr;
    mySize = 0;
};

myVector::myVector(unsigned int size){
    myData = new double[size]();
    mySize = size;
};

void myVector::allocate(unsigned int size){
    myData = new double[size]();
    mySize = size;
}

myVector::~myVector(){
    // std::cout << "Vector deleted" << std::endl;
    delete[] myData;
    myData = nullptr;
}

double myVector::operator[](int i){
    return myData[i];
};

void myVector::setValue(double value, unsigned int i){
    if(i>mySize){
        std::cout << "Index out of bounds.  Operation aborted" << std::endl;
        return;
    }
    myData[i] = value;
};

void myVector::addValue(double value, unsigned int i){
    
    if (i > mySize)
    {
        std::cout << "Index out of bounds.  Operation aborted" << std::endl;
        return;
    }

    myData[i] += value;
};

double * myVector::getDataPointer(){
    return myData;
};

ostream & operator<<(ostream &os, myVector & vec){
    os << "Print vector:" << std::endl;
    for(int i = 0; i < vec.mySize; i++){
            os << vec.myData[i] << std::endl;
    }
    return os;
};

void myVector::vector_multiply(matrix & a, myVector & b){
    if(a.myCols != b.mySize){
        std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
        return;
    }

    if(mySize != a.myRows){
        std::cout << "Error: Dimension Mismatch.  Multiply aborted" << std::endl;
        return;        
    }

    for(int i = 0; i < mySize; i++){
        
        for(int j = 0; j < b.mySize; j++){
            myData[i] += a[i][j] * b[j];
        }
    }
};

void myVector::vector_subtract(myVector & a, myVector & b){
    if(a.mySize != b.mySize && a.mySize != mySize){
        std::cout << "Error: Dimension Mismatch.  Subtract aborted" << std::endl;
        return;
    }
    for(int i = 0; i < mySize; i++){
        myData[i] = b[i] - a[i];

    }
};

double myVector::l2Norm(myVector & a){
    if(a.mySize != mySize){
        cerr << "Error in L2 Norm: Vector dimension mismatch." << endl;
        exit(1);
    }
    double sum = 0.0;
    for(int i = 0; i < mySize; i++){
        sum += (a[i]-myData[i])*(a[i] - myData[i]);
    }
    return sqrt(sum);
}

void myVector::printVector(string fileName){
    std::cout << "Writing Matrix to File." << std::endl;
    std::cout << "..." << std::endl;

    ofstream soln;
    soln.open(fileName+".mat", ios::trunc);
    for(int i = 0; i < mySize; i++){
        soln << (*this)[i] << endl;
    }
    soln.close();

    std::cout << "Done writing loads to file." << std::endl;
    std::cout << std::endl;
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