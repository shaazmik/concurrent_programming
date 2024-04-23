#include <iostream>
#include <pthread.h>
#include <cmath>
#include <vector>
#include <stack>
#include <mutex>

struct Part {
    double a, b, f_a, f_b, s_ab;
};

struct ThreadData{
    double epsilon;
    int rank;
    int thread_count;
    int* global_nactive;
    double* global_result;
    std::vector<Part>* global_stk;
    std::mutex* global_result_mtx;
    std::mutex* global_stk_mtx;
    std::mutex* global_has_task_mtx;
};

double f(double x) {
    return sin(1/x);
}

void* worker(void* arg) {
    struct ThreadData* data = (ThreadData*)arg;

    std::vector<Part> stk;
    stk.reserve(1000);

    while ( 1 )
    {
        Part p;
        data->global_has_task_mtx->lock();

        {
            const std::lock_guard< std::mutex> lock( *data->global_stk_mtx);

            p = data->global_stk->back();
            data->global_stk->pop_back();

            if (data->global_stk->size() != 0 )
            {
                data->global_has_task_mtx->unlock();
            }

            if ( p.a <= p.b )
            {
                (*data->global_nactive)++;
            }
        }

        if ( p.a > p.b )
        {
            break; // stop 
        }

        double result = 0;
        while ( 1 )
        {
            double c = (p.a + p.b) / 2;
            double f_c = f(c);

            double s_ac = (p.f_a + f_c) * (c - p.a) / 2;
            double s_cb = (f_c + p.f_b) * (p.b - c) / 2;
            double s_acb = s_ac + s_cb;

            if ( std::abs( p.s_ab - s_acb) > data->epsilon * std::abs( s_acb) )
            {
                stk.push_back( Part{ p.a, c, p.f_a, f_c, s_ac});
                p.a = c;
                p.f_a = f_c;
                p.s_ab = s_cb;
            }
            else
            {
                result += s_acb;

                if ( stk.empty() )
                    break;

                p = stk.back();
                stk.pop_back();
            }

            const size_t kStackOffloadSize = 8;
            if ( stk.size() > kStackOffloadSize && data->global_stk->size() == 0 )
            {
                const std::lock_guard< std::mutex> lock( *data->global_stk_mtx);

                while ( stk.size() > 1 )
                {
                    data->global_stk->push_back( stk.back());
                    stk.pop_back();
                }

                data->global_has_task_mtx->unlock();
            }
        }

        {

            const std::lock_guard< std::mutex> lock( *data->global_result_mtx);
            *data->global_result += result;
        }

        {
            const std::lock_guard< std::mutex> lock( *data->global_stk_mtx);
            (*data->global_nactive)--;

            if ( *data->global_nactive == 0 && data->global_stk->size() == 0 )
            {
                for ( int i = 0; i < data->thread_count; i++ )
                    data->global_stk->push_back( Part{ 2.0, 1.0, 0.0, 0.0, 0.0});

                data->global_has_task_mtx->unlock();
            }
        }
    }
}

int main(int argc, char** argv) {
    const int thread_count = 4;
    double a = 0.001;
    double b = 1.0; 
    double epsilon = 1e-10;
    double f_a = f(a);
    double f_b = f(b);
    double s_ab = (f_a + f_b) * (b - a) / 2;
  
    int global_nactive = 0;
    double global_result = 0;
    std::vector<Part> global_stk;
    global_stk.reserve(1000);
    global_stk.push_back(Part{ a, b, f_a, f_b, s_ab});

    std::mutex global_result_mtx;
    std::mutex global_stk_mtx;
    std::mutex global_has_task_mtx;
    
    std::vector<pthread_t> threads(thread_count);

    for (int rank = 0; rank < thread_count; rank++)
    {
        pthread_create(&threads[rank], NULL, worker, new ThreadData{epsilon, rank, thread_count,
                                                                             &global_nactive, &global_result, &global_stk,
                                                                             &global_result_mtx, &global_stk_mtx, &global_has_task_mtx});
    }

    for (int rank = 0; rank < thread_count; rank++)
    {
        pthread_join(threads[rank], NULL);
    }

    std::cout.precision(10);
    std::cout << "Integral is " << global_result << std::endl;
    
    return 0;
}
