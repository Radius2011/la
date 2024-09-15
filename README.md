# la

Namespace for everything below

## Matrix Coordinates:

The coordinates of the matrix elements are implemented in a slightly unusual way:

![Матриця](RMD/Image.png)

## Return type of functions:

All functions below that have return type - bool, excluding Vect_IndepCheck, return true if the function executes successfully, if the function returns false check if you pass the arguments correctly.

# Classes:
## SMatrix:

When creating an object, it takes one integer value - the number of dimensions of the matrix: 
```
la::SMatrix a(2); //creating 2D matrix named "a"
```
The class is created for simple use of matrix data by a function.  
Constructor performs dynamic memory allocation.  
The destructor automatically releases the used memory.  
Access to vector elements is implemented through the [ ] operator.
```
a[0][1]=3;
```
### Methods:
- getDimensions - returns the number of dimensions in this matrix. ```a.getDimensions();```  
- input - row by row fills the matrix with numbers from std::cin. ```a.input();```
- print - outputs the matrix to the console using std::cout. Takes two arguments:
```char elementSeparator``` and ```char rowSeparator```, the character after each outputed element and the character at the end of each line, respectively(default ' ' and '\n'). For example, if you call the function like ```a.print(',','\n');```, the console will show  
```
    1,2,3,  
    4,5,6,
    7,8,9,
```
- transpose - transposes the matrix.
- trace - returns a trace of the matrix.
- scalarMult - multiplies every element of the matrix by given argument ```double scalar```.
- rowMultByNum - takes two arguments ```int rowIndex``` and ```double k``` - index of the row which should be multiplied by ```k``` and k.
- rowAdd - adds row with index ```rowIndexA``` multiplied by ```k``` to row with index ```rowIndexB```.
- rowSwap - swaps the row with the index ```rowIndexA``` with the row with the index ```rowIndexB```.
- toREF - transforms the matrix to a REF(Row Echelon Form).
- stabilizeZeros - removes the remainder from the division of binary numbers from matrix elements close to zero, so that when printing elements of the type ```-0.0000000007089``` it does not show ```-0```(this function is executed at the end of almost every function in this file).

## Vector:

When creating an object, it takes one integer value - the number of dimensions of the vector: 
```
la::Vector a(2); //creating 2D vector named "a"
```
The class is created for simple use of vector data by a function.  
Constructor performs dynamic memory allocation.  
The destructor automatically releases the used memory.    
Access to vector elements is implemented through the [ ] operator.
```
a[0]=0;
```
### Methods:
- getDimensions - returns the number of dimensions in this vector. ```a.getDimensions();```  
- input - row by row fills the vector with numbers from std::cin. ```a.input();```
- print - outputs the matrix to the console using std::cout. Takes one arguments:
```char elementSeparator``` - the character after each outputed element(default ' '). For example, if you call the function like ```a.print(',');```, the console will show  
```
    1,2,3,
    
```
- length - returns a length of vector.
- normalize - normalizes vector(divides every element by vector length).
- stabilizeZeros - removes the remainder from the division of binary numbers from vector elements close to zero, so that when printing elements of the type ```-0.0000000007089``` it does not show ```-0```(this function is executed at the end of almost every function in this file).

# Functions:

## Operations with matrices:
### MatMult:

The name of the function speaks for itself.  
MatMult is a simple function for multiplying matrices.

### MatAdd:

MatAdd is a function for adding matrices.

### MatCopy:

Copy matrix(A) data to matrix(B).

### det_Gauss:

det_Gauss function computes the determinant of a square matrix using Gaussian elimination.

### MatInv_Gauss_Jordan:

Matrix inversion function, which works according to the Gauss-Jordan method.

## Operations with vectors:
### VectAdd:

VectAdd function computes the element-wise sum of two vectors and stores the result in a third vector.

### LinTrans:

LinTrans is a function to perform linear transformation.

### dot_prod:

Returns dot product of given Vectors.

### cross_prod:

Returns cross product of given Vectors.

### Vect_IndepCheck:

Checks linear independence of n n-dimensional vectors represented as matrix.

## Systems of linear equations:
### LinSys_Gauss:

It takes three arguments matrix, which represents coefficients near variables, vector (numbers after =) and resultVector where stores result values of all variables. Uses the Gaussian method to find the roots of the equations.

### LinSys_LUdecomp:

It takes three arguments matrix, which represents coefficients near variables, vector (numbers after =) and resultVector where stores result values of all variables.Uses the LU-decomposition method to find the roots of the equations.

### LinSys_Cram:

It takes three arguments matrix, which represents coefficients near variables, vector (numbers after =) and resultVector where stores result values of all variables.Uses the Cramar's rule to find the roots of the equations.
