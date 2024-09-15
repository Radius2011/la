#ifndef LA_H
#define LA_H

#include <cmath>
#include <iostream>

namespace la {
    class SMatrix {
    private:
        int dimensions_;
        double** cells_;
    public:
        SMatrix(int dimensions);
        ~SMatrix();
        double* operator[](int x);
        SMatrix(const SMatrix&) = delete;
        SMatrix& operator=(const SMatrix&) = delete;
        int getDimensions() const;
        void input();
        void print(char element_separator=' ', char row_separator='\n');
        void transpose();
        int trace();
        void scalarMult(double scalar);
        bool rowMultByNum(int rowIndex, double k);
        bool rowAdd(int rowIndexA, int rowIndexB, double k);
        bool rowSwap(int rowIndexA, int rowIndexB);
        bool toREF();
        void stabilizeZeros();
    };

    class Vector {
    private:
        int dimensions_;
        double* cells_;
    public:
        Vector(int dimensions);
        ~Vector();
        double& operator[](int y);
        Vector(const Vector&) = delete;
        Vector& operator=(const Vector&) = delete;
        int getDimensions() const;
        void input();
        void print(char element_separator=' ');
        void scalarMult(double scalar);
        double length();
        void normalize();
        void stabilizeZeros();
    };
    bool MatMult(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix);
    bool MatAdd(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix);
    bool LinTrans(SMatrix& matrix, Vector& vector, Vector& resultVector);
    bool det_Gauss(SMatrix& matrix, double& result);
    bool MatInv_Gauss_Jordan(SMatrix& matrix, SMatrix& resultMatrix);
    bool VectAdd(Vector& vector1, Vector& vector2, Vector& resultVector);
    bool dot_prod(Vector& vector1, Vector& vector2, double& result);
    bool cross_prod(Vector& vector1, Vector& vector2, Vector& resultVector);
    bool LinSys_Gauss(SMatrix& matrix, Vector& vector, Vector& resultVector);
    bool LinSys_LUdecomp(SMatrix& matrix, Vector& vector, Vector& resultVector);
    bool LinSys_Cram(SMatrix& matrix, Vector& vector, Vector& resultVector);
    bool Vect_IndepCheck(SMatrix& matrix);
    bool MatCopy(SMatrix& matrixA, SMatrix& matrixB);
}

#endif
