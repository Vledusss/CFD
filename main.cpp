#include <iostream>
#include <iomanip>
#include <fstream>

#include "Godunov.hpp"
#include "Burgers.hpp"

template<indexType N>
std::array<double, N> func1(const double dx) {
    std::array<double, N> res;
    double x = 0;
    for (indexType i = 0; i < N; ++i) {
        res[i] = std::sin(4 * M_PI * x);
        x += dx;
    }
    return res;
}

template<typename T, indexType N>
void print(std::ofstream& file, const std::array<T, N>& arr) noexcept {
    for (const auto& elem : arr) { file << std::setprecision(15) << elem << ','; }
    file << std::endl;
}

int main() {
    const indexType M = 1000;  // по времени
    const indexType N = 3000;  // по пространству

    const double dt = 0.001;
    const double dx = 0.001;

    Burgers<N> eq(-1, 1);
    // Burgers<N> eq(func1<N>(dx));

    std::ofstream file1("res.csv");

    if (file1.is_open()) {
        for (indexType i = 0; i < M; ++i) {
            print(file1, eq.state);
            Godunov<Burgers<N>, N>::solve(eq, dx, dt, (i + 1) * dt);
        }
    }
    file1.close();
}