#include "Matrix.h"
#include <assert.h>


static Matrix* strassen(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows);

Matrix* do_strassen(Matrix* A, Matrix* B) {
    assert(A != NULL);
    assert(B != NULL);

    if ((A->rows != A->cols) || (B->cols != B->rows) || (B->cols != A->cols)) {
        fprintf(stderr, "Failed to calculate: wrong matrix size\nMatrixA.rows=%d, MatrixA.columns=%d\n"
                        "MatrixB.rows=%d. MatrixB.columns=%d\n", A->rows, A->cols, B->rows, B->cols);

        return NULL;
    }

    Matrix* dstMatrix = (Matrix*)calloc(1, sizeof(Matrix));
    strassen(dstMatrix, A, B, A->rows);

    return dstMatrix;
}


static Matrix* strassen(Matrix* dest, const Matrix* srcA, const Matrix* srcB, uint32_t rows){

    if (rows == 2) {
        return matrix_ijk_matmul(dest, srcA, srcB);
    }

    uint32_t len = rows / 2;

    matrix_autofree Matrix *a11 = matrix_new(len, len), *a12 = matrix_new(len, len), *a21 = matrix_new(len, len), *a22 = matrix_new(len, len),
                           *b11 = matrix_new(len, len), *b12 = matrix_new(len, len), *b21 = matrix_new(len, len), *b22 = matrix_new(len, len), *c11 = matrix_new(len, len), *c12 = matrix_new(len, len), *c21 = matrix_new(len, len), *c22 = matrix_new(len, len), *m1 = matrix_new(len, len), *m2 = matrix_new(len, len), *m3 = matrix_new(len, len), *m4 = matrix_new(len, len), *m5 = matrix_new(len, len), *m6 = matrix_new(len, len), *m7 = matrix_new(len, len), *temp1 = matrix_new(len, len), *temp2 = matrix_new(len, len);

    /* Divide matrix into four parts */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(a11, i, j) = MXY(srcA, i, j);
            MXY(a12, i, j) = MXY(srcA, i, j + len);
            MXY(a21, i, j) = MXY(srcA, i + len, j);
            MXY(a22, i, j) = MXY(srcA, i + len, j + len);

            MXY(b11, i, j) = MXY(srcB, i, j);
            MXY(b12, i, j) = MXY(srcB, i, j + len);
            MXY(b21, i, j) = MXY(srcB, i + len, j);
            MXY(b22, i, j) = MXY(srcB, i + len, j + len);
        }
    }

    /* Calculate seven formulas of strassen Algorithm */
    strassen(m1, matrix_addition(temp1, a11, a22), matrix_addition(temp2, b11, b22), len);
    strassen(m2, matrix_addition(temp1, a21, a22), b11, len);
    strassen(m3, a11, matrix_subtraction(temp1, b12, b22), len);
    strassen(m4, a22, matrix_subtraction(temp1, b21, b11), len);
    strassen(m5, matrix_addition(temp1, a11, a12), b22, len);
    strassen(m6, matrix_subtraction(temp1, a21, a11), matrix_addition(temp2, b11, b12), len);
    strassen(m7, matrix_subtraction(temp1, a12, a22), matrix_addition(temp2, b21, b22), len);

     /* Merge the answer of matrix dest */
    /* c11 = m1 + m4 - m5 + m7 = m1 + m4 - (m5 - m7) */
    matrix_subtraction(c11, matrix_addition(temp1, m1, m4), matrix_subtraction(temp2, m5, m7));
    matrix_addition(c12, m3, m5);
    matrix_addition(c21, m2, m4);
    matrix_addition(c22, matrix_subtraction(temp1, m1, m2), matrix_addition(temp2, m3, m6));

    /* Store the answer of matrix multiplication */
    FOR(i, len)
    {
        FOR(j, len)
        {
            MXY(dest, i, j) = MXY(c11, i, j);
            MXY(dest, i, j + len) = MXY(c12, i, j);
            MXY(dest, i + len, j) = MXY(c21, i, j);
            MXY(dest, i + len, j + len) = MXY(c22, i, j);
        }
    }

    return dest;
}