# la

Simple namespace for everything below

## Matrix Coordinates:

The coordinates of the matrix elements are implemented in a slightly unusual way:

![Матриця](RMD/Image.png)
# Classes:
### SMatrix:

When creating an object, it takes one integer value - the number of dimensions of the matrix: 
```
la::SMatrix a(2); //creating 2D matrix named "a"
```
The class is created for simple use of matrix data by a function.  
Constructor performs dynamic memory allocation.  
The destructor automatically releases the used memory.  
The class also has a getDimensions() function that returns an integer - the number of dimensions in this matrix.  
Access to matrix elements is implemented through the [ ] operator.
```
a[0][2]=4;
```
### Vector:

When creating an object, it takes one integer value - the number of dimensions of the vector: 
```
la::Vector a(2); //creating 2D vector named "a"
```
The class is created for simple use of vector data by a function.  
Constructor performs dynamic memory allocation.  
The destructor automatically releases the used memory.  
The class also has a getDimensions() function that returns an integer - the number of dimensions in this vector.  
Access to vector elements is implemented through the [ ] operator.
```
a[0]=0;
```
# Functions:

## Operations with matrices:
### MatMult:

The name of the function speaks for itself.  
MatMult is a simple function for multiplying matrices with a time complexity of O(n^3), where n is the number of dimensions of the matrix.

### MatAdd:

MatAdd is a function for adding matrices with a time complexity of O(n^2), where n is the number of dimensions of the matrix.

### transpose_matrix:

Transpose_matrix is a function performs the transposition of a square matrix.  
It hastime complexity of O(n^2), where n is the number of dimensions of the matrix.

### det_Gauss:

det_Gauss function computes the determinant of a square matrix using Gaussian elimination.  
It has time complexity of O(n^3), where n is the number of dimensions of the matrix.

### MatInv_Gauss_Jordan:

Matrix inversion function, which works according to the Gauss-Jordan method.  
It has time complexity of O(n^3), where n is the number of dimensions of the matrix.

### MatTrace:
MatTrace is a simple function, which calculates the trace of given Matrix, with time complexity of O(n), where n is the number of dimensions of the matrix.

### Mat_scalarMult:

Mat_scalarMult function multiplies each element of a square matrix by a scalar.  
Time complexity of the Mat_scalarMult function is O(n^2), where n is the number of dimensions of the matrix.

## Operations with vectors:
### VectAdd:

VectAdd function computes the element-wise sum of two vectors and stores the result in a third vector.  
Time complexity of this function is O(n), where n is the number of dimensions of the vector.

### LinTrans:

LinTrans is a function to perform linear transformation with a time complexity of O(n^2), where n is the number of dimensions of the matrix.

### Vect_scalarMult:

Vect_scalarMult function multiplies each element of a vector by a scalar.  
Time complexity of this function is O(n), where n is the number of dimensions of the vector.

### dot_prod:

Returns dot product of given Vectors.
Time complexity of this function is O(n), where n is the number of dimensions of the vector.

### cross_prod:

Returns cross product of given Vectors.
Time complexity of this function is O(1).

### VectLen:

Returns Vector length.  
Time complexity of this function is O(n), where n is the number of dimensions of the vector.

### VectNorm:

Normalizes Vector.  
Time complexity of this function is O(n), where n is the number of dimensions of the vector.
