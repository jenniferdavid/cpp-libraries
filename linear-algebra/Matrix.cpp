/**
 * @file Matrix.cpp
 * @brief Library for handling matrices
 * @author Ashish Ranjan (Jalan)
 * @version 0.1
 * @date 2013-03-28
 */

#include <iostream>
#include <iomanip>                       /* for setw() */
#include <cstdlib>                       /* for exit() */
#include <cmath>                         /* for pow()  */

#include "Matrix.h"

/**
 * @brief Constructor 1
 *
 * @param r
 * @param c
 */
Matrix::Matrix(int r, int c): mat(r, std::vector<double>(c))
{
    rows = r;
    columns = c;
}

/**
 * @brief Constructor 2
 *
 * @param r
 * @param c
 * @param init_val
 */
Matrix::Matrix(int r, int c, double init_val): mat(r, std::vector<double>(c, init_val))
{
    rows = r;
    columns = c;
}

/**
 * @brief Matrix Addition
 *
 * @param mtx
 *
 * @return
 */
Matrix Matrix::operator+(const Matrix &mtx) const
{
    //matrix addition is possible
    if (rows == mtx.rows && columns == mtx.columns)
    {
        Matrix m(rows, columns);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                m.mat[i][j] = mat[i][j] + mtx.mat[i][j];
        return m;
    }
    //matrix addition is not possible
    else
    {
        std::cout << std::endl << "Matrix Addition NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Matrix Subtraction
 *
 * @param mtx
 *
 * @return
 */
Matrix Matrix::operator-(const Matrix &mtx) const
{
    //matrix subtraction is possible
    if (rows == mtx.rows && columns == mtx.columns)
    {
        Matrix m(rows, columns);
        for (int i = 0; i < rows; i++)
            for (int j = 0; j < columns; j++)
                m.mat[i][j] = mat[i][j] - mtx.mat[i][j];
        return m;
    }
    //matrix subtraction is not possible
    else
    {
        std::cout << std::endl << "Matrix Subtraction NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Matrix Multiplication
 *
 * @param mtx
 *
 * @return
 */
Matrix Matrix::operator*(const Matrix &mtx) const
{
    //matrix multiplication is possible
    if (columns == mtx.rows)
    {
        Matrix m(rows, mtx.columns);
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < mtx.columns; j++)
                for (int k = 0, l = 0; k < columns && l < mtx.rows; k++, l++)
                    m.mat[i][j] += mat[i][k] * mtx.mat[l][j];
        }
        return m;
    }
    //matrix multiplication is not possible
    else
    {
        std::cout << std::endl << "Matrix Multiplication NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Scalar Multiplication
 *
 * @param num
 *
 * @return
 */
Matrix Matrix::operator*(double num) const
{
    Matrix m(rows, columns);
    for (int i = 0; i < rows; i++)
        for (int j = 0; j < columns; j++)
            m.mat[i][j] = mat[i][j] * num;
    return m;
}

/**
 * @brief Returns transpose of a matrix
 *
 * @return
 */
Matrix Matrix::transpose() const
{
    //a_ij transforms to a_ji
    Matrix m(columns, rows);
    for (int i = 0; i < m.rows; i++)
        for (int j = 0; j < m.columns; j++)
            m.mat[i][j] = mat[j][i];
    return m;
}

/**
 * @brief Returns determinant of a matrix
 *
 * @return
 */
double Matrix::det() const
{
    if (rows == columns)
    {
        //the matrix is a square matrix
        //the number of rows is equal to the number of columns but are used individually in diffenet places just for understanding
        if (rows == 2)
            return ((mat[0][0] * mat[1][1]) - (mat[0][1] * mat[1][0]));
        else if (rows == 1)
            return mat[0][0];
        else
        {
            double determinant = 0;
            for (int c = 0; c < columns; c++)
            {
                if (mat[0][c] == 0)
                    continue;
                //creating a new matrix with an order equal to the order of the parent matrix minus one
                Matrix m_new((rows - 1), (rows - 1));
                //initializing the new matrix with the values from the parent matrix minus one row and one column
                for (int i = 0; i < m_new.rows; i++)
                {
                    for (int j = 0; j < m_new.columns; j++)
                    {
                        if (j >= c)
                            m_new.mat[i][j] = mat[i + 1][j + 1];
                        else
                            m_new.mat[i][j] = mat[i + 1][j];
                    }
                }
                determinant += mat[0][c] * pow(-1, c) * m_new.det();
            }
            return determinant;
        }
    }
    else
    {
        std::cout << std::endl << "Matrix Determinant NOT possible!!!" << std::endl;
        exit(1);
    }

}

/**
 * @brief Returns inverse of a matrix
 *
 * @return
 */
Matrix Matrix::inv() const
{
    if (rows == columns)
    {
        //the matrix is a square matrix
        //the number of rows is equal to the number of columns but are used individually in diffenet places just for understanding
        Matrix m(rows, columns);
        double determinant = det();
        if (!(determinant > (-__ERROR_DOUBLE__) && determinant < (__ERROR_DOUBLE__)))       //if(determinant != 0)
        {
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < columns; j++)
                {
                    //creating a new matrix with an order equal to the order of the parent matrix minus one
                    Matrix m_new((m.rows - 1), (m.rows - 1));
                    //initializing the new matrix with the values from the parent matrix minus one row and one column
                    for (int k = 0; k < m_new.rows; k++)
                    {
                        for (int l = 0; l < m_new.columns; l++)
                        {
                            if (k >= i && l >= j)
                                m_new.mat[k][l] = mat[k + 1][l + 1];
                            else if (k >= i && l < j)
                                m_new.mat[k][l] = mat[k + 1][l];
                            else if (k < i && l >= j)
                                m_new.mat[k][l] = mat[k][l + 1];
                            else
                                m_new.mat[k][l] = mat[k][l];
                        }
                    }
                    m.mat[i][j] = pow(-1, (i + j)) * m_new.det();
                }
            }
            return ((m.transpose()) * (1 / determinant));
        }
        else
        {
            std::cout << std::endl << "Matrix Inverse NOT possible!!!" << std::endl;
            exit(1);
        }
    }
    else
    {
        std::cout << std::endl << "Matrix Inverse NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Matrix Exponentiation
 *
 * @param n
 *
 * @return
 */
Matrix Matrix::operator^(int n) const
{
    if (rows == columns && n == -1)
    {
        return this->inv();
    }
    else if (rows == columns && n >= 2)
    {
        Matrix m(rows, columns);
        m = (*this);
        for (int i = 2; i <= n; i++)
            m = (*this) * m;
        return m;
    }
    else
    {
        std::cout << std::endl << "Matrix Exponentiation NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Solution of Ax = b or x = b/A using Gaussian Elimination
 *
 * @param mtx
 *
 * @return
 */
Matrix Matrix::operator/(const Matrix &mtx) const
{
    //A = mtx
    //b = *this
    //if A is a square matrix AND if b is a vector AND if matrix multiplication of A and b is possible
    if (mtx.rows == mtx.columns && columns == 1 && mtx.columns == rows)
    {
        //solution matrix : x
        Matrix x(mtx.rows, 1);

        //augmented matrix : a
        Matrix a(mtx.rows, (mtx.columns + 1));
        for (int i = 0; i < a.rows; i++)
            for (int j = 0; j < (a.columns - 1); j++)
                a.mat[i][j] = mtx.mat[i][j];
        for (int i = 0; i < a.rows; i++)
            a.mat[i][(a.columns - 1)] = mat[i][0];

        /**********  Gaussian Elimination Step 1 : Forward Elimination **********/
        int counter = 0;
        for (int i = 0; i < (a.rows - 1); i++)
        {
            double pivot = a.mat[i][i];
            //pivots cannot be zero : By Definition
            if (pivot > (-__ERROR_DOUBLE__) && pivot < (__ERROR_DOUBLE__))              //if(pivot == 0)
            {
                counter += 1;
                //exchange rows : Row "i" and Row "i+counter"
                for (int c = 0; c < a.columns; c++)
                {
                    //swap elements of the two rows : Row "i" and Row "i+counter"
                    a.mat[i][c] = a.mat[i][c] + a.mat[i + counter][c];
                    a.mat[i + counter][c] = a.mat[i][c] - a.mat[i + counter][c];
                    a.mat[i][c] = a.mat[i][c] - a.mat[i + counter][c];
                }
                i--;
                continue;
            }
            else                                        //else
            {
                for (int j = i + 1; j < a.rows; j++)
                {
                    double factor = a.mat[j][i] / pivot;
                    a.mat[j][i] = 0;
                    for (int c = i + 1; c < a.columns; c++)
                        a.mat[j][c] -= factor * a.mat[i][c];
                }
            }
            counter = 0;
        }

        /**********  Gaussian Elimination Step 2 : Back Substitution **********/
        for (int i = (a.rows - 1); i >= 0; i--)
        {
            double num = a.mat[i][(a.columns - 1)];
            for (int j = (a.columns - 2); j > i; j--)
                num -= a.mat[i][j] * x.mat[j][0];
            x.mat[i][0] = num / a.mat[i][i];
        }

        return x;
    }
    else
    {
        std::cout << std::endl << "Gaussian Elimination NOT possible!!!" << std::endl;
        exit(1);
    }
}

/**
 * @brief Displays all the elements of a matrix
 */
void Matrix::display() const
{
    std::cout << std::endl;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            if (mat[i][j] > (-__ERROR_DOUBLE__) && mat[i][j] < (__ERROR_DOUBLE__))
                std::cout << std::setw(10) << 0;
            else
                std::cout << std::setw(10) << mat[i][j];
        }
        std::cout << std::endl;
    }
}

