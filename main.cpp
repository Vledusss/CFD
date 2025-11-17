#include <iostream>
#include <iomanip>
#include <fstream>

#include "Godunov.hpp"
#include "Burgers/Burgers.hpp"
#include "Euler/EulerEquations.hpp"

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

void testBurgers() {
    const indexType M = 1000;  // по времени
    const indexType N = 3000;  // по пространству

    const double dt = 0.001;
    const double dx = 0.001;

    Burgers<N> eq(-1, 1);
    // Burgers<N> eq(func1<N>(dx));

    std::ofstream file("res.csv");

    if (file.is_open()) {
        for (indexType i = 0; i < M; ++i) {
            print(file, eq.state);
            Godunov<Burgers<N>, N>::solve(eq, dx, dt);
        }
    }
    file.close();
}

void testEuler() {
    const indexType N = 1000;         // по пространству
    Euler::EulerEquation<N> eq;       // gamma = 1.4
    
    // Инициализация условий Римана
    for (indexType i = 0; i < N/2; ++i) {
        // Левая часть: высокое давление
        eq.state[i] = Euler::EulerState(1.0, 2.0, 2.5);  // ρ=1, u=0, p=1.0
    }
    for (indexType i = N/2; i < N; ++i) {
        // Правая часть: низкое давление  
        eq.state[i] = Euler::EulerState(0.125, 0.0, 0.1); // ρ=0.125, u=0, p=0.1
    }
    
    const double dx = N / 10000.;
    const double dt = 0.1 * dx;  // CFL условие
    const double t_end = 100;

    double t = 0.0;
    
    std::ofstream file("res.csv");
    file << "t,x,rho,u,p" << std::endl;

    while (t < t_end) {
        Godunov<Euler::EulerEquation<N>, N>::solve(eq, dx, dt);
        for (indexType i = 0; i < N; ++i) {
            const double p = eq.getPressure(eq.state[i]);
            const double u = eq.getVelocity(eq.state[i]);
            const double rho = eq.state[i].rho;
        file << t << ',' << i * dx << ',' << rho << ',' << u << ',' << p << std::endl;
        }

        t += dt;
    }
    
    file.close();

}

int main() {
    testEuler();
}