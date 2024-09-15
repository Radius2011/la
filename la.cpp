#include "la.h"

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
    void SMatrix::input(){
        for(int i=0;  i<dimensions_; i++) {
            for(int j=0; j<dimensions_; j++){
                std::cin >> cells_[j][i];
            }
        }
    }
    void SMatrix::print(char element_separator, char row_separator){
        stabilizeZeros();
        for(int i=0;  i<dimensions_; i++) {
            for(int j=0; j<dimensions_; j++){
                std::cout << cells_[j][i] << element_separator;
            }
            std::cout << row_separator;
        }
    }
    void SMatrix::transpose(){
        SMatrix cmatrix(dimensions_);
        for(int i=0; i<dimensions_; i++) for(int j=0; j<dimensions_; j++) cmatrix[i][j]=cells_[i][j];
        for(int i=0; i<dimensions_; i++){
            for(int j=0; j<dimensions_; j++){
                cells_[i][j]=cmatrix[j][i];
            }
        }
    }
    int SMatrix::trace(){
        int res=0;
        for(int i=0; i<dimensions_; i++) res+=cells_[i][i];
        return res;
    }
    void SMatrix::scalarMult(double scalar){
        for(int i=0; i<dimensions_; i++){
            for(int j=0; j<dimensions_; j++){
                cells_[i][j]*=scalar;
            }
        }
        stabilizeZeros();
    }
    bool SMatrix::rowMultByNum(int rowIndex, double k){
        if(dimensions_<rowIndex||rowIndex<0) return false;
        for (int l=0; l<dimensions_; l++) cells_[l][rowIndex]*=k;
        stabilizeZeros();
        return true;
    }
    bool SMatrix::rowAdd(int rowIndexA, int rowIndexB, double k){
        if(dimensions_<rowIndexA||rowIndexA<0||dimensions_<rowIndexB||rowIndexB<0) return false;
        for (int l=0; l<dimensions_; l++) cells_[l][rowIndexB]+=k*cells_[l][rowIndexA];
        stabilizeZeros();
        return true;
    }
    bool SMatrix::rowSwap(int rowIndexA, int rowIndexB){
        if(dimensions_<rowIndexA||rowIndexA<0||dimensions_<rowIndexB||rowIndexB<0) return false;
        for (int l=0; l<dimensions_; l++) {
            double temp = cells_[l][rowIndexA];
            cells_[l][rowIndexA] = cells_[l][rowIndexB];
            cells_[l][rowIndexB] = temp;
        }
        return true;
    }
    bool SMatrix::toREF(){
        for (int i=0; i<dimensions_; i++) {
            if (cells_[i][i]==0) {
                int sw=-1;
                for (int u=i+1; u<dimensions_; u++) {
                    if (cells_[i][u]!=0) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return false;
                }
                rowSwap(i, sw);
            }
            double k=1/cells_[i][i];
            rowMultByNum(i,k);
            for(int l=i+1; l<dimensions_; l++){
                k=-1*cells_[i][l];
                rowAdd(i, l, k);
            }
        }
        stabilizeZeros();
        return true;
    }
    void SMatrix::stabilizeZeros(){ for(int i=0; i<dimensions_; i++) for(int j=0; j<dimensions_; j++) if(fabs(cells_[i][j])<1e-9) cells_[i][j]=0;}
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
    void Vector::input(){
        for(int i=0;  i<dimensions_; i++) std::cin >> cells_[i];
    }
    void Vector::print(char element_separator){
        stabilizeZeros();
        for(int i=0;  i<dimensions_; i++) {
            std::cout << cells_[i] << element_separator;
        }
        std::cout << std::endl;
    }
    void Vector::scalarMult(double scalar){
        for(int i=0; i<dimensions_; i++) cells_[i]*=scalar;
        stabilizeZeros();
    }
    double Vector::length(){
        double sum=0;
        for(int i=0; i<dimensions_; i++) sum+=cells_[i]*cells_[i];
        return sqrt(sum);
    }
    void Vector::normalize(){
        double len_1=1/length();
        scalarMult(len_1);
        stabilizeZeros();
    }
    void Vector::stabilizeZeros(){
        for(int i=0; i<dimensions_; i++){
            if(fabs(cells_[i])<1e-9) cells_[i]=0;
        }
    }

    // Matrix multiplication function
    bool MatMult(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix){
        double sum;
        int dimensions = matrix1.getDimensions();
        if(dimensions!=matrix2.getDimensions()||dimensions!=resultMatrix.getDimensions()) return false;
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                sum=0;
                for(int k=0; k<dimensions; k++){
                    sum+=matrix1[k][i]*matrix2[j][k];
                }
                resultMatrix[j][i]=sum;
            }
        }
        resultMatrix.stabilizeZeros();
        return true;
    }

    // Matrix addition function
    bool MatAdd(SMatrix& matrix1, SMatrix& matrix2, SMatrix& resultMatrix){
        int dimensions = matrix1.getDimensions();
        if(dimensions!=matrix2.getDimensions()||dimensions!=resultMatrix.getDimensions()) return false;
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                resultMatrix[i][j]=matrix1[i][j]+matrix2[i][j];
            }
        }
        resultMatrix.stabilizeZeros();
        return true;
    }

    // Function to perform linear transformation
    bool LinTrans(SMatrix& matrix, Vector& vector, Vector& resultVector){
        double sum;
        int dimensions = matrix.getDimensions();
        if(dimensions!=vector.getDimensions()||dimensions!=resultVector.getDimensions()) return false;
        for(int i=0; i<dimensions; i++){
            sum=0;
            for(int k=0; k<dimensions; k++){
                sum+=matrix[k][i]*vector[k];
            }
            resultVector[i]=sum;
        }
        resultVector.stabilizeZeros();
        return true;
    }


    // Calculates a determinant of given matrix
    bool det_Gauss(SMatrix& matrix, double& result) {
        result = 1;
        int dimensions = matrix.getDimensions();
        SMatrix cmatrix(dimensions);
        MatCopy(matrix, cmatrix);
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
                    return false;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = cmatrix[l][i];
                    cmatrix[l][i] = cmatrix[l][sw];
                    cmatrix[l][sw] = temp;
                }
                result = -result;
            }
            for (int j=i+1; j<dimensions; j++) {
                double k=-1*cmatrix[i][j]/cmatrix[i][i];
                cmatrix.rowAdd(i, j, k);
            }
        }
        for (int i=0; i<dimensions; i++) result*=cmatrix[i][i];
        return true;
    }

    //Matrix inversion function, which works according to the Gauss-Jordan method
    bool MatInv_Gauss_Jordan(SMatrix& matrix, SMatrix& resultMatrix){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultMatrix.getDimensions()) return false;
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
                    return false;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = cmatrix[l][i];
                    cmatrix[l][i] = cmatrix[l][sw];
                    cmatrix[l][sw] = temp;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = resultMatrix[l][i];
                    resultMatrix[l][i] = resultMatrix[l][sw];
                    resultMatrix[l][sw] = temp;
                }
            }
            double k=1/cmatrix[i][i];
            cmatrix.rowMultByNum(i, k);
            resultMatrix.rowMultByNum(i, k);
            for(int l=0; l<dimensions; l++){
                if(l!=i){
                    k=-1*cmatrix[i][l];
                    cmatrix.rowAdd(i, l, k);
                    resultMatrix.rowAdd(i, l, k);
                }
            }
        }
        resultMatrix.stabilizeZeros();
        return true;
    }

    //Vector Addition
    bool VectAdd(Vector& vector1, Vector& vector2, Vector& resultVector){
        int dimensions=vector1.getDimensions();
        if(dimensions!=vector2.getDimensions()||dimensions!=resultVector.getDimensions()) return false;
        for(int i=0; i<dimensions; i++){
            resultVector[i]=vector1[i]+vector2[i];
        }
        resultVector.stabilizeZeros();
        return true;
    }

    //Returns dot product of given Vectors
    bool dot_prod(Vector& vector1, Vector& vector2, double& result){
        result=0;
        int dimensions=vector1.getDimensions();
        if(dimensions!=vector2.getDimensions()) return false;
        for(int i=0; i<dimensions; i++){
            result+=vector1[i]*vector2[i];
        }
        return true;
    }

    //Returns cross product of given Vectors
    bool cross_prod(Vector& vector1, Vector& vector2, Vector& resultVector){
        int dimensions = 3;
        if(dimensions!=vector1.getDimensions()||dimensions!=vector2.getDimensions()||dimensions!=resultVector.getDimensions()) return false;
        resultVector[0]=vector1[1]*vector2[2]-vector1[2]*vector2[1];
        resultVector[1]=vector1[2]*vector2[0]-vector1[0]*vector2[2];
        resultVector[2]=vector1[0]*vector2[1]-vector1[1]*vector2[0];
        resultVector.stabilizeZeros();
        return true;
    }

    //Function that solves a system of linear equations using the Gauss's method
    bool LinSys_Gauss(SMatrix& matrix, Vector& vector, Vector& resultVector){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultVector.getDimensions()||dimensions!=vector.getDimensions()) return false;
        SMatrix cmatrix(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                cmatrix[i][j]=matrix[i][j];
            }
            resultVector[i]=vector[i];
        }
        for (int i=0; i<dimensions; i++) {
            if (cmatrix[i][i]==0) {
                int sw=-1;
                for (int u=i+1; u<dimensions; u++) {
                    if (cmatrix[i][u]==0) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return false;
                }
                cmatrix.rowSwap(i, sw);
                double temp = resultVector[i];
                resultVector[i] = resultVector[sw];
                resultVector[sw] = temp;
            }
            double k=cmatrix[i][i];
            for (int l=0; l<dimensions; l++) {
                cmatrix[l][i]/=k;
            }
            resultVector[i]/=k;
            for(int l=0; l<dimensions; l++){
                if(l!=i){
                    k=cmatrix[i][l];
                    for (int j=0; j<dimensions; j++) {
                        cmatrix[j][l]-=k*cmatrix[j][i];
                    }
                    resultVector[l]-=k*resultVector[i];
                }
            }
        }
        resultVector.stabilizeZeros();
        return true;
    }

    //Function that solves a system of linear equations using the LU decomposition method
    bool LinSys_LUdecomp(SMatrix& matrix, Vector& vector, Vector& resultVector){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultVector.getDimensions()||dimensions!=vector.getDimensions()) return false;
        SMatrix U(dimensions);
        SMatrix L(dimensions);
        Vector P(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                U[i][j]=matrix[i][j];
                L[i][j]=(i==j)?1:0;
            }
            P[i]=vector[i];
        }
        for (int i=0; i<(dimensions-1); i++) {
            if (U[i][i]==0) {
                int sw=-1;
                for (int u=i+1; u<dimensions; u++) {
                    if (U[i][u]!=0) {
                        sw=u;
                        break;
                    }
                }
                if (sw==-1) {
                    return false;
                }
                for (int l=0; l<dimensions; l++) {
                    double temp = U[l][i];
                    U[l][i] = U[l][sw];
                    U[l][sw] = temp;
                    temp = P[i];
                    P[i] = P[sw];
                    P[sw] = temp;
                }
            }
            for (int j=i+1; j<dimensions; j++) {
                double k=U[i][j]/U[i][i];
                for (int l=i; l<dimensions; l++) {
                    U[l][j] -= k*U[l][i];
                }
                L[i][j]=k;
            }
        }
        Vector y(dimensions);
        for (int i=0; i<dimensions; i++) {
            double sum=0;
            for (int j=0; j<i; j++) {
                sum+=L[j][i]*y[j];
            }
            y[i]=P[i]-sum;
        }
        for (int i=dimensions-1; i>-1; i--) {
            double sum=0;
            for (int j=i+1; j<dimensions; j++) {
                sum+=U[j][i]*resultVector[j];
            }
            resultVector[i]=(y[i]-sum)/U[i][i];
        }
        resultVector.stabilizeZeros();
        return true;
    }

    //Function that solves a system of linear equations using the Cramer's rule
    bool LinSys_Cram(SMatrix& matrix, Vector& vector, Vector& resultVector){
        int dimensions = matrix.getDimensions();
        if(dimensions!=resultVector.getDimensions()||dimensions!=vector.getDimensions()) return false;
        SMatrix proc(dimensions);
        for(int i=0; i<dimensions; i++){
            for(int j=0; j<dimensions; j++){
                proc[i][j]=matrix[i][j];
            }
        }
        double a;
        det_Gauss(matrix,a);
        double db=1/a;
        for(int i=0; i<dimensions; i++){
            for(int l=0; l<dimensions; l++){
                proc[i][l]=vector[l];
            }
            double o;
            det_Gauss(proc,o);
            resultVector[i]=o*db;
            for(int l=0; l<dimensions; l++){
                proc[i][l]=matrix[i][l];
            }
        }
        resultVector.stabilizeZeros();
        return true;
    }

    //Checks linear independence of n n-dimensional vectors represented as matrix
    bool Vect_IndepCheck(SMatrix& matrix){
        double a;
        det_Gauss(matrix,a);
        return a!=0;
    }

    //copy matrix(A) data to matrix(B)
    bool MatCopy(SMatrix& matrixA, SMatrix& matrixB){
        int dimensions=matrixA.getDimensions();
        if(dimensions!=matrixB.getDimensions()) return false;
        for(int i=0; i<dimensions; i++) for(int j=0; j<dimensions; j++) matrixB[i][j]=matrixA[i][j];
        return true;
    }
}
