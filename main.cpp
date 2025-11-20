#include <iostream>
#include <iomanip>
#include <fstream>

#include "Godunov.hpp"
#include "Burgers/Burgers.hpp"
#include "Euler/EulerEquations.hpp"

// template<indexType N>
// std::array<double, N> func1(const double dx) {
//     std::array<double, N> res;
//     double x = 0;
//     for (indexType i = 0; i < N; ++i) {
//         res[i] = std::sin(4 * M_PI * x);
//         x += dx;
//     }
//     return res;
// }

// template<typename T, indexType N>
// void print(std::ofstream& file, const std::array<T, N>& arr) noexcept {
//     for (const auto& elem : arr) { file << std::setprecision(15) << elem << ','; }
//     file << std::endl;
// }

// void testBurgers() {
//     const indexType M = 1000;  // по времени
//     const indexType N = 3000;  // по пространству

//     const double dt = 0.001;
//     const double dx = 0.001;

//     Burgers<N> eq(-1, 1);
//     // Burgers<N> eq(func1<N>(dx));

//     std::ofstream file("res.csv");

//     if (file.is_open()) {
//         for (indexType i = 0; i < M; ++i) {
//             print(file, eq.state);
//             Godunov<Burgers<N>, N>::solve(eq, dx, dt);
//         }
//     }
//     file.close();
// }

void testEuler() {
    const indexType N = 1000;                                   // по пространству
    std::vector<std::tuple<double, Euler::Equation<N>>> eq;     // gamma = 1.4
    Euler::Equation<N> initial;
    
    // НУ (ρ, u, p)
    for (indexType i = 0; i < N / 4; ++i) {
        initial.states[i] = Euler::State(1.0, 0.0, 1.0);
    }
    for (indexType i = N / 4; i < N; ++i) {
        initial.states[i] = Euler::State(0.2, 0.0, 0.5);
    }
    
    const double dx = 0.1;
    const double dt = 0.1 * dx;  // CFL условие
    const double startTime = 0;
    const double endTime = 100;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");
    file << "t,x,rho,u,p" << std::endl;

    Godunov<Euler::Equation<N>, N>::solve(eq, startTime, endTime, dx, dt);
    for (const auto& elem : eq) {
        for (indexType i = 0; i < N; ++i) {
            const auto U = std::get<1>(elem);
            const double rho = U.states[i].rho;
            const double u = U.getVelocity(U.states[i]);
            const double p = U.getPressure(U.states[i]);
            file << std::get<0>(elem) << ',' << i * dx << ',' << rho << ',' << u << ',' << p << std::endl;
        }
    }
    
    file.close();

}

int main() {
    testEuler();
}