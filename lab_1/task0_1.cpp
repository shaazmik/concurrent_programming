#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <mpi.h>

/*
### Алгоритм Монте-Карло для вычисления числа π

Алгоритм Монте-Карло для вычисления числа π основан на принципе случайной выборки. Он использует 
статистический метод для аппроксимации значения π, который осуществляется следующим образом:

1. Подготовка: Визуализируйте единичный квадрат (стороны которого равны 1), внутри которого вписан круг 
радиусом 0.5 (площадь этого круга равна π/4, так как площадь круга = πr^2, а радиус r = 0.5).

2. Генерация случайных точек: Случайным образом генерируются точки в пределах квадрата. 
Каждая точка имеет координаты (x, y), где x и y – случайные числа в интервале 0, 1.

3. Определение положения точек: Для каждой сгенерированной точки проверяется, находится ли она внутри круга. 
Это осуществляется путем проверки условия: если x^2 + y^2 <= 1, точка находится внутри круга 
(расстояние от центра до точки не превышает радиус круга, равный 1).

4. Вычисление числа π: Отношение количества точек внутри круга к общему количеству сгенерированных точек 
приблизительно равно площади круга, деленной на площадь квадрата (π/4), так как вероятность попадания точки 
в круг равна отношению площадей круга и квадрата. Таким образом, умножив это отношение на 4, получаем 
приближенное значение π.

*/



int main(int argc, char* argv[])
{
    unsigned long long total_iterations = 1000000000; // Общее количество итераций
    unsigned long long i;
    int rank, size;
    double x, y, distance_squared, pi_estimate;
    int points_in_circle = 0, total_points_in_circle = 0;

    MPI_Init(&argc, &argv); // Инициализация MPI
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Получить ранг текущего процесса
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Получить общее количество процессов

    srand(time(NULL) + rank); // Инициализация генератора случайных чисел
    
    int iterations_per_process = total_iterations / size; // Итераций на процесс

    // Основной цикл
    for(i = 0; i < iterations_per_process; i++) {
        x = (double)rand() / RAND_MAX; // Генерируем случайную точку x
        y = (double)rand() / RAND_MAX; // Генерируем случайную точку y
        distance_squared = x * x + y * y; // Рассчитываем квадрат расстояния до центра
        if (distance_squared <= 1) {
            points_in_circle++; // Точка внутри круга
        }
    }

    // Собираем результаты от всех процессов
    MPI_Reduce(&points_in_circle, &total_points_in_circle, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) { // Если это главный процесс
        pi_estimate = 4.0 * total_points_in_circle / total_iterations; // Вычисляем приближение к π
        printf("Estimated pi value: %f\n", pi_estimate);
    }

    MPI_Finalize(); // Завершение работы с MPI

    return 0;
}
