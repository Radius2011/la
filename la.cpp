#include "la.h"
#include <cmath>

namespace la {
    SMatrix::SMatrix(int dimensions) : dimensions_(dimensions) {
        cells_ = new double*[dimensions_];
        for (int i=0; i<dimensions_; i++) {
            cells_[i] = new double[dimensions_];
        }
    }
    SMatrix::~SMatrix(){
        for (int i=0; i<dimensions_; i++) {
            delete[] cells_[i];
        }
        delete[] cells_;
    }
    double* SMatrix::operator[](int x) {
        if (x<0||x>=dimensions_) {
            return nullptr;
        }
        return cells_[x];
    }
    int SMatrix::getDimensions() const { return dimensions_; }
    Vector::Vector(int dimensions) : dimensions_(dimensions) {
        cells_ = new double[dimensions_];
    }
    Vector::~Vector(){
        delete[] cells_;
    }
    double& Vector::operator[](int y) {
        static double errorValue = -1;
        if(y<0||y>=dimensions_){
            return errorValue;
        }
        return cells_[y];
    }
    int Vector::getDimensions() const { return dimensions_; }

    // Matrix multiplication function
    void MatMult(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix){
        double sum;
        int dimensions = matrix1.getDimensions();
        if(dimensions!=matrix2.getDimensions()||dimensions!=resultMatrix.getDimensions())  return;
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                sum=0;
                for(int k=0; k<dimensions; k++){
                    sum+=matrix1[k][i]*matrix2[j][k];
                }
                resultMatrix[j][i]=sum;
            }
        }
    }

    // Matrix addition function
    void MatAdd(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix){
        int dimensions = matrix1.getDimensions();
        if(dimensions!=matrix2.getDimensions()||dimensions!=resultMatrix.getDimensions())  return;
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                resultMatrix[i][j]=matrix1[i][j]+matrix2[i][j];
            }
        }
    }

    // Function to perform linear transformation
    void LinTrans(SMatrix& matrix, Vector& vector, Vector& resultVector){
        double sum;
        int dimensions = matrix.getDimensions();
        if(dimensions!=vector.getDimensions()||dimensions!=resultVector.getDimensions())  return;
        for(int i=0; i<dimensions; i++){
            sum=0;
            for(int k=0; k<dimensions; k++){
                sum+=matrix[k][i]*vector[k];
            }
            resultVector[i]=sum;
        }
    }

    // Function for matrix transpose
    void Transpose_matrix(SMatrix& matrix){
        int dimensions = matrix.getDimensions();
        SMatrix cmatrix(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                cmatrix[i][j]=matrix[i][j];
            }
        }
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                matrix[i][j]=cmatrix[j][i];
            }
        }
    }

    // Calculates a determinant of given matrix
    double det_Gauss(SMatrix& matrix) {
        double res = 1;
        int dimensions = matrix.getDimensions();
        SMatrix cmatrix(dimensions);
        for (int i=0; i<dimensions; i++) {
            for (int j=0; j<dimensions; j++) {
                cmatrix[i][j]=matrix[i][j];
            }
        }
        for (int i=0; i<(dimensions-1); i++) {
            if (cmatrix[i][i]==0) {
                int sw=-1;
                for (int u=i+1; u<dimensions; u++) {
                    if (cmatrix[i][u]!=0) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return 0;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = cmatrix[i][l];
                    cmatrix[i][l] = cmatrix[sw][l];
                    cmatrix[sw][l] = temp;
                }
                res = -res;
            }
            for (int j=i+1; j<dimensions; j++) {
                double k=cmatrix[i][j]/cmatrix[i][i];
                for (int l=i; l<dimensions; l++) {
                    cmatrix[l][j] -= k*cmatrix[l][i];
                }
            }
        }
        for (int i=0; i<dimensions; i++) {
            res*=cmatrix[i][i];
        }
        return res;
    }

    //Matrix inversion function, which works according to the Gauss-Jordan method
    void MatInv_Gauss_Jordan(SMatrix& matrix, SMatrix& resultMatrix){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultMatrix.getDimensions()) return;
        SMatrix cmatrix(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                cmatrix[i][j]=matrix[i][j];
                resultMatrix[i][j]=(i==j)?1:0;
            }
        }
        for (int i=0; i<dimensions; i++) {
            if (cmatrix[i][i]==0) {
                int sw=-1;
                for (int u=i+1; u<dimensions; u++) {
                    if (cmatrix[i][u]!=0) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = cmatrix[i][l];
                    cmatrix[i][l] = cmatrix[sw][l];
                    cmatrix[sw][l] = temp;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = resultMatrix[i][l];
                    resultMatrix[i][l] = resultMatrix[sw][l];
                    resultMatrix[sw][l] = temp;
                }
            }
            double k=cmatrix[i][i];
            for (int l=0; l<dimensions; l++) {
                cmatrix[l][i]/=k;
                resultMatrix[l][i]/=k;
            }
            for(int l=0; l<dimensions; l++){
                if(l!=i){
                    k=cmatrix[i][l];
                    for (int j=0; j<dimensions; j++) {
                        cmatrix[j][l]-=k*cmatrix[j][i];
                        resultMatrix[j][l]-=k*resultMatrix[j][i];
                    }
                }
            }
        }
    }

    //function which returns a trace of given matrix
    int MatTrace(SMatrix& matrix){
        int res=0;
        for(int i=0; i<matrix.getDimensions(); i++) res+=matrix[i][i];
        return res;
    }

    //Multiplies Matrix by scalar
    void Mat_scalarMult(SMatrix& matrix, int scalar){
        for(int i=0; i<matrix.getDimensions(); i++){
            for(int j=0; j<matrix.getDimensions(); j++){
                matrix[i][j]*=scalar;
            }
        }
    }

    //Vector Addition
    void VectAdd(Vector& vector1, Vector& vector2, Vector& resultVector){
        int dimensions=vector1.getDimensions();
        if(dimensions!=vector2.getDimensions()||dimensions!=resultVector.getDimensions()) return;
        for(int i=0; i<dimensions; i++){
            resultVector[i]=vector1[i]+vector2[i];
        }
    }

    //Multiplies Vector by scalar
    void Vect_scalarMult(Vector& vector, int scalar){
        for(int i=0; i<vector.getDimensions(); i++){
            vector[i]*=scalar;
        }
    }

    //Returns dot product of given Vectors
    long long dot_prod(Vector& vector1, Vector& vector2){
        long long sum=0;
        int dimensions=vector1.getDimensions();
        if(dimensions!=vector2.getDimensions()) return __INT_MAX__+1;
        for(int i=0; i<dimensions; i++){
            sum+=vector1[i]*vector2[i];
        }
        return sum;
    }

    //Returns cross product of given Vectors
    void cross_prod(Vector& vector1, Vector& vector2, Vector& resultVector){
        int dimensions = 3;
        if(dimensions!=vector1.getDimensions()||dimensions!=vector2.getDimensions()||dimensions!=resultVector.getDimensions()) return;
        resultVector[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
        resultVector[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];
        resultVector[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];
    }

    //Returns Vector length
    double VectLen(Vector& vector){
        long long sum=0;
        for(int i=0; i<vector.getDimensions(); i++) sum+=vector[i]*vector[i];
        return sqrt(sum);
    }

    //Normalizes Vector
    void VectNorm(Vector& vector){
        double len=VectLen(vector);
        for(int i=0; i<vector.getDimensions(); i++) vector[i]/=len;
    }

    void LinSys_Gauss(SMatrix& matrix, Vector& vector, Vector& resultVector){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultVector.getDimensions()||dimensions!=vector.getDimensions()) return;
        SMatrix cmatrix(dimensions);
        Vector cvector(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                cmatrix[i][j]=matrix[i][j];
            }
            cvector[i]=vector[i];
        }
        for (int i=0; i<dimensions; i++) {
            if (std::abs(cmatrix[i][i])<1e-9) {
                int sw=-1;
                for (int u=i+1; u<dimensions; u++) {
                    if (std::abs(cmatrix[i][u])>=1e-9) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = cmatrix[l][i];
                    cmatrix[l][i] = cmatrix[l][sw];
                    cmatrix[l][sw] = temp;
                }
                double temp = cvector[i];
                cvector[i] = cvector[sw];
                cvector[sw] = temp;
            }
            double k=cmatrix[i][i];
            for (int l=0; l<dimensions; l++) {
                cmatrix[l][i]/=k;
            }
            cvector[i]/=k;
            for(int l=0; l<dimensions; l++){
                if(l!=i){
                    k=cmatrix[i][l];
                    for (int j=0; j<dimensions; j++) {
                        cmatrix[j][l]-=k*cmatrix[j][i];
                    }
                    cvector[l]-=k*cvector[i];
                }
            }
        }
        for(int i=0; i<dimensions; i++) resultVector[i]=cvector[i];
    }
}
