# Matrix Inversion Using Cholesky Decomposition Method

Matrix inversion techniques based on Cholesky decomposition and the related LDL decomposition are efficient techniques widely used for inversion of positive definite/symmetrich matrices across multiple fiedls.
If A is a positive definite Hermitian matrix, Cholesky decomposition factorises it into a lower triangular matrix and its conjugate transpose. Which is useful for efficient numerical solutions.  

 ![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/be170024-c0ba-41a7-ab01-ad6259517693)

Or equivalently, using an upper triangular matrix ![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/75ba34d5-3916-439c-9a34-64d33b7563f8) as 

![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/bdd531c6-5d3c-40e5-bef0-bd2f51b1b162)


Every symmetric positive definite matrix A can be decomposed into a product of a unique lower triangular matrix L and its transpose.

In a software implementation the upper triangular matrix is preferred as operations are row-wise and compatible with C programming language. 
Cholesky decomposition is of order![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/f2ea0e7f-4f58-4c63-90ef-c5788ecb7823) and reguires ![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/5f3801ce-4573-4353-a511-e7260e333350)
 operations. Mtarix inversion based on Cholesky decomposition in numerically stable for well conditioned matrices. 

The elements ![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/0db4d662-79e9-48c2-b720-2c0c412befff) are given as follows. 

Diagonal elements: 

![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/ba6d197b-7285-4ab3-b507-402a407cc270)

Upper triangular elements i.e i < j

![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/bb37356a-9504-44f2-85a3-7efda8aba46e)


Mathematical Description of the Algorithm ; 
 Input data: a symmetric positive definite matrix A whose elements are denoted by aij).

 Output data: the lower triangular matrix L whose elements are denoted by lij).

The Cholesky algorithm can be represented in the form

![image](https://github.com/Mserhatarslan/M-x-N-size-matrix-transpose/assets/63358327/fc6ef71b-c503-4649-aff3-0ebd69f10322)

# Pseudo Code of Cholesky Decomposition : 


procedure CholeskyFactorization(dim: uint16_t, A: float[], L: float[])

    for j := 0 to dim-1 do
        for i := j to dim-1 do
            sum := A[i * dim + j]
            for k := 0 to j-1 do
                sum := sum - L[i * dim + k] * L[j * dim + k]
            end for

            if i == j then
                L[i * dim + j] := sqrt(sum)
            else
                L[i * dim + j] := sum / L[j * dim + j]
            end if
        end for
    end for
end procedure

