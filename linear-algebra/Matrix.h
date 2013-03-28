/**
 * @file Matrix.h
 * @brief Library for handling matrices
 * @author Ashish Ranjan (Jalan)
 * @version 0.1
 * @date 2013-03-28
 */

#ifndef MATRIX_H

#define MATRIX_H

#include <vector>

#define __ERROR_DOUBLE__ 0.0001

class Matrix
{
public:
    /********** Constructors **********/
    Matrix(int, int);
    Matrix(int, int, double);
    /********** Constructors **********/

    inline std::vector<double> &operator[](int r) //Matrix Element Access using [][]
    {
        return mat[r];
    }

    /********** Basic Matrix Operations **********/
    Matrix operator+(const Matrix &) const; // Matrix Addition        : Binary Operation
    Matrix operator-(const Matrix &) const; // Matrix Subtraction     : Binary Operation
    Matrix operator*(const Matrix &) const; // Matrix Multiplication  : Binary Operation
    Matrix operator*(double) const;         // Scalar Multiplication  : Binary Operation
    Matrix transpose() const;               // Transpose              : Unary Operation
    double det() const;                     // Determinant            : Unary Operation
    Matrix inv() const;                     // Inverse                : Unary Operation
    /********** Basic Matrix Operations **********/

    /********** Advanced Matrix Operations **********/
    Matrix operator^(int) const;            // Matrix Exponentiation          : Binary Operation
    Matrix operator/(const Matrix &) const; // Gaussian Elimination Solution  : Binary Operation
    /********** Advanced Matrix Operations **********/

    void display() const;
private:
    int rows, columns;
    std::vector< std::vector<double> > mat;
};

#endif /* end of include guard: MATRIX_H */

