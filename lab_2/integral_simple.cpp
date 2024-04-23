#include <iostream>
#include <pthread.h>
#include <cmath>
#include <vector>


struct ThreadData {
    double start, end, result;
    uint32_t steps;
};

double f(double x) {
    return sin(1/x);
}

void* worker(void* arg) {
    ThreadData* data = (ThreadData*)arg;
    double range = data->end - data->start;
    double step = range / data->steps;
    double integral = 0.0;
    for(int i = 0; i < data->steps; ++i) {
        double x = data->start + (i + 0.5) * step;
        integral += f(x);
    }
    data->result = integral * step;
    return 0;
}

int main(int argc, char** argv) {
    const int thread_count = 1;
    double a = 0.000001;
    double b = 1.0; 
    uint32_t steps = 100000000; 
    
    std::vector<ThreadData> data(thread_count);
    std::vector<pthread_t> threads(thread_count);
    double part = (b - a) / thread_count; 

    for(int i = 0; i < thread_count; ++i) {
        data[i].start = a + i * part;
        data[i].end = a + (i + 1) * part;
        data[i].steps = steps / thread_count;
        pthread_create(&threads[i], nullptr, &worker, &data[i]);
    }

    double total_integral = 0.0;
    for (int i = 0; i < thread_count; ++i) {
        pthread_join(threads[i], nullptr);
        total_integral += data[i].result;
    }

    std::cout.precision(17);
    std::cout << "Integral is " << total_integral << std::endl;

    return 0;
}
