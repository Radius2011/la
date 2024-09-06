#include <iostream>

namespace la {
    
    // Define a class for a matrix with dynamic allocation
    class RegisterMatrix {
        private:
            int dimensions_;    // Number of dimensions of the matrix
            int** cells_;   // Pointer to a 2D array representing the matrix
        public:
            RegisterMatrix(int dimensions) : dimensions_(dimensions) {
                cells_ = new int*[dimensions_];
                for (int i=0; i<dimensions_; i++) {
                    cells_[i] = new int[dimensions_];
                }
            }
            ~RegisterMatrix(){
                for (int i=0; i<dimensions_; i++) {
                    delete[] cells_[i];
                }
                delete[] cells_;
            }
            int* operator[](int x) {
                if (x<0||x>=dimensions_) {
                    std::cout << "X index out of bounds!";
                    return nullptr;
                }
                return cells_[x];
            }
            RegisterMatrix(const RegisterMatrix&) = delete;
            RegisterMatrix& operator=(const RegisterMatrix&) = delete;
            int getDimensions() const { return dimensions_; }
    };

    // Function to perform matrix multiplication
    void MatrixMultiplication(RegisterMatrix& matrix1, RegisterMatrix& matrix2, RegisterMatrix& resultMatrix){
        int sum, dimensions = matrix1.getDimensions();
        if(dimensions!=matrix2.getDimensions()||dimensions!=resultMatrix.getDimensions()) {std::cout << "The number of matrix dimensions does not match!"; return;}
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
}