#include <iostream>
#include <vector>
#include <thread>
#include <mutex>
#include <functional>
#include <cmath>
#include <queue>

struct Task {
    double start;
    double end;
    uint32_t steps;
    Task(double s, double e, uint32_t st) : start(s), end(e), steps(st) {}
};

class IntegralSolver {
public:
    IntegralSolver(double a, double b, uint32_t totalSteps, uint32_t tasks) {
        double step = (b - a) / totalSteps;
        for (int i = 0; i < tasks; ++i) {
            double taskStart = a + (b - a) / tasks * i;
            double taskEnd = taskStart + (b - a) / tasks;
            queue.emplace(taskStart, taskEnd, (taskEnd - taskStart) / step);
        }
    }

    void start(int threadCount) {
        std::vector<std::thread> threads;
        for (int i = 0; i < threadCount; ++i) {
            threads.emplace_back([this]() { this->process(); });
        }

        for (auto &th : threads) th.join();
    }

    double getResult() const { return result; }

private:
    double result = 0.0;
    std::queue<Task> queue;
    std::mutex mutex;

    void process() {
        while (true) {
            Task task(0, 0, 0);
            // Безопасное извлечение задачи из очереди
            {
                std::lock_guard<std::mutex> lock(mutex); // RAII для блокировки
                if (queue.empty()) return;
                task = queue.front();
                queue.pop();
            }

            computeTask(task);
        }
    }

    void computeTask(const Task& task) {
        double stepSize = (task.end - task.start) / task.steps;
        double halfStepSum = 0.0, fullStepSum = 0.0;

        for (uint32_t i = 0; i < task.steps; ++i) {
            double xMiddle = task.start + (i + 0.5) * stepSize;
            fullStepSum += std::sin(1 / xMiddle);
        }
        fullStepSum *= stepSize;

        uint32_t halfSteps = task.steps * 2;
        double halfStepSize = (task.end - task.start) / halfSteps;
        for (uint32_t i = 0; i < halfSteps; ++i) {
            double xMiddle = task.start + (i + 0.5) * halfStepSize;
            halfStepSum += std::sin(1 / xMiddle);
        }
        halfStepSum *= halfStepSize;

        double resultTask = std::abs(halfStepSum - fullStepSum) < 0.00001 ? fullStepSum : halfStepSum;

        std::lock_guard<std::mutex> lock(mutex);
        result += resultTask;
    }
};

int main() {
    IntegralSolver solver(0.001, 1.0, 10000000, 100);
    solver.start(4);

    std::cout << "Integral result is: " << solver.getResult() << std::endl;
}

