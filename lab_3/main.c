#include "Matrix.h"
#include "time.h"
#include "omp.h"


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

    Matrix* mtx1 = matrix_new(2048, 2048);
    Matrix* mtx2 = matrix_new(2048, 2048);

    initTestData(mtx1, mtx2);

   // matrix_print("Matrix1:", mtx1);
   // matrix_print("Matrix2:", mtx2);

    clock_t start_time = clock();

    Matrix* mtxClassic = matrix_new(2048, 2048);
    mtxClassic = matrix_ijk_matmul(mtxClassic, mtx1, mtx2);

    clock_t end_time = clock();
    double time_spent = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("Running time ijk_matmul algo: %f seconds\n", time_spent);
   
    //matrix_print("Answer Classic:", mtxClassic);

    double start_time_parallel = omp_get_wtime();  // Запуск таймера OpenMP

    Matrix* mtx3 = do_strassen(mtx1, mtx2);

    double end_time_parallel = omp_get_wtime();  // Остановка таймера OpenMP

    time_spent = end_time_parallel - start_time_parallel;  // Расчёт времени
    printf("Running time strassen algo:   %f seconds\n", time_spent);
  
    // matrix_print("Answer Strassen:", mtx3);

    matrix_free(&mtx1);
    matrix_free(&mtx2);
    matrix_free(&mtx3);
    matrix_free(&mtxClassic);
}