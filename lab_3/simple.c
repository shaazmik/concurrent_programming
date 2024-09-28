#include "Matrix.h"

Matrix* do_strassen(Matrix* A, Matrix* B);

void initTestData(Matrix* mtx1, Matrix* mtx2) {
    assert(mtx1 != NULL);
    assert(mtx2 != NULL);


    for (int i = 0; i < mtx1->rows; i++) {
        for (int j = 0; j < mtx1->cols; j++) {
            mtx1->values[j + i*mtx1->cols] = j + i*mtx1->cols;
        }
    }
    for (int i = 0; i < mtx2->rows; i++) {
        for (int j = 0; j < mtx2->cols; j++) {
            mtx2->values[j + i*mtx2->cols] = mtx2->cols*mtx2->rows - j - i*mtx2->cols;
        }
    }
}

int main() {

    Matrix* mtx1 = matrix_new(2, 2);
    Matrix* mtx2 = matrix_new(2, 2);

    initTestData(mtx1, mtx2);

    matrix_print("Matrix1:", mtx1);
    matrix_print("Matrix2:", mtx2);


    Matrix* mtx3 = NULL;//do_strassen(mtx1, mtx2);
    
    fprintf(stdout, "%p\n", mtx3);
    matrix_free(&mtx1);
    matrix_free(&mtx2);
}