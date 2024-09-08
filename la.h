#ifndef LA_H
#define LA_H

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
    };
    void MatMult(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix);
    void MatAdd(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix);
    void LinTrans(SMatrix& matrix, Vector& vector, Vector& resultVector);
    void Transpose_matrix(SMatrix& matrix);
    double det_Gauss(SMatrix& matrix);
    void MatInv_Gauss_Jordan(SMatrix& matrix, SMatrix& resultMatrix);
    int MatTrace(SMatrix& matrix);
    void Mat_scalarMult(SMatrix& matrix, int scalar);
    void VectAdd(Vector& vector1, Vector& vector2, Vector& resultVector);
    void Vect_scalarMult(Vector& vector, int scalar);
    long long dot_prod(Vector& vector1, Vector& vector2);
    void cross_prod(Vector& vector1, Vector& vector2, Vector& resultVector);
    double VectLen(Vector& vector);
    void VectNorm(Vector& vector);
}

#endif