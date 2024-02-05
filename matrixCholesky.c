#include <stdio.h>
#include <math.h>
#include <stdint.h>
#include <assert.h>
#define NDEBUG 

inline float calculateSqrt(float x) {  // Optimization Flag 
    return 1.0f / sqrtf(x);
}

void track_matrixCholesky(uint16_t dim, const float *A, float *G) {
    float v[4] = {0.f};
    uint16_t i, j, k;

    // Symmetric ve Positive definite Control

#ifndef NDEBUG  // Check only in debug mode

    for (i = 0; i < dim; i++) {
        for (j = 0; j < dim; j++) {
            assert(A[i * dim + j] == A[j * dim + i]);
        }
        assert(A[i * dim + i] > 0);
    }
#endif

    for (i = 0; i < dim * dim; i++) {
        G[i] = 0;
    }

    for (j = 0; j < dim; j++) {
        for (i = j; i < dim; i++) {
            v[i] = A[i * dim + j];
        }

        for (k = 0; k < j; k++) {
            for (i = j; i < dim; i++) {
                v[i] -= G[j * dim + k] * G[i * dim + k];
            }
        }
        float temp = calculateSqrt(v[j]);

        for (i = j; i < dim; i++) {
                G[i * dim + j] = v[i] * temp;
        }
    }
}

int main(void) {

    float A[] = {
        25.0f,  15.0f,  -5.0f,
        15.0f,  18.0f,    .0f,
        -5.0f,  0.0f,   11.0f
    };

    uint16_t dim = 3;
    float G[dim * dim];

    track_matrixCholesky(dim, A, G);

    printf("Cholesky faktörizasyonu:\n");
    for (int i = 0; i < dim; i++) {
        for (int j = 0; j < dim; j++) {
            printf("%f\t", G[i * dim + j]);
        }
        printf("\n");
    }

    return 0;
}

//    gcc -o matrix matrixCholesky.c -lm -O5
//   ./matrix