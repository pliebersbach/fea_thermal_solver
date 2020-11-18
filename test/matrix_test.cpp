#include "../include/math_utilities.h"

#include <cstdlib>
#include <iostream>
#include <cmath>

int main(){

    // matrix a = matrix(4,2);
    // // matrix b = matrix(3,3);
    // // matrix c = matrix(3,3);

    // for(int i = 0; i < 4; i++){
    //     for(int j = 0; j < 2; j++){
    //         a[i][j] = i*2+j;
    //         // b[i][j] = (i*3+j)*2;
    //     }
    // }

    // cout << a;

    // myVector b = myVector(2);

    // b.setValue(2.0, 0);
    // b.setValue(3.0, 1);

    // myVector c = myVector(4);

    // c.vector_multiply(a,b);

    //  cout << b;

    // // c.multiply(a, b);

    //  cout << c;

    int n = 10;

    matrix A = matrix(n);
    myVector b; 
    myVector x;
    x.allocate(n);
    b.allocate(n);

    for (int i = 0; i < n; i++)
    {
        for(int j = 0; j < n; j++){
            A[i][j] = fmin(i+1, j+1);
        }
    }

    double dx = 1.0/(n-1);
    for(int i = 0; i < n; i++){
        
        b.setValue(1.0+sin(M_PI*i*dx), i);
    }

    // A.gaussElimination(b.getVector(), x.getVector());
    A.gaussElimination(b, x);

    std::cout << "Matrix: " << std::endl;
    std::cout << A;
    std::cout << "rhs vector:" << std::endl;
    std::cout << b;
    std::cout << "solution vector:" << std::endl;
    std::cout << x;


    return 0;
}