#include <stdio.h>
#include <mpi.h>

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    // Убедимся, что запущено хотя бы два процесса
    if (world_size != 2) {
        fprintf(stderr, "World size must be 2 for %s\n", argv[0]);
        MPI_Abort(MPI_COMM_WORLD, 1); 
    }

    const int MAX_NUMBER = 1000;
    int numbers[MAX_NUMBER];
    double start_time, end_time;

    if (world_rank == 0) {
        // Инициализация данных для отправки
        for(int i = 0; i < MAX_NUMBER; i++) {
            numbers[i] = i;
        }

        // Начало замера времени
        start_time = MPI_Wtime();

        // Отправляем массив чисел процессу 1
        MPI_Send(numbers, MAX_NUMBER, MPI_INT, 1, 0, MPI_COMM_WORLD);

        // Получаем массив чисел обратно от процесса 1
        MPI_Recv(numbers, MAX_NUMBER, MPI_INT, 1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Окончание замера времени
        end_time = MPI_Wtime();

        printf("Time of send and receive = %f / %d seconds\n", (end_time-start_time), 2 * MAX_NUMBER);
    } else if (world_rank == 1) {
        // Получаем массив чисел от процесса 0
        MPI_Recv(numbers, MAX_NUMBER, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        // Отправляем массив чисел обратно процессу 0
        MPI_Send(numbers, MAX_NUMBER, MPI_INT, 0, 0, MPI_COMM_WORLD);
    }

    MPI_Finalize();
    return 0;
}
