#include <iostream>
#include <iomanip>
#include <fstream>

#include "Godunov.hpp"
#include "Burgers/Burgers.hpp"
#include "Euler/EulerEquations.hpp"

void testEuler() {
    const indexType N = 1000;                                   // по пространству
    std::vector<std::tuple<double, Euler::Equation<N>>> eq;     // gamma = 1.4
    
    Euler::Equation<N> initial;
    
    // НУ (ρ, u, p)
    for (indexType i = 0; i < N / 2; ++i) { initial.states[i] = Euler::State(1.0, 0.0, 2.0); }
    for (indexType i = N / 2; i < N; ++i) { initial.states[i] = Euler::State(0.8, 0.0, 0.5); }
    
    const double dx = 0.1;
    
    const double startTime = 0;
    const double endTime = 100;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");

    if (file.is_open()) {
        file << "t,x,rho,u,p" << std::endl;
        Godunov<Euler::Equation<N>, N>::solve(eq, startTime, endTime, dx, 0.4);
        // std::cout << "DONE" << std::endl;

        for (const auto& elem : eq) {
            for (indexType i = 0; i < N; ++i) {
                const auto U = std::get<1>(elem);
                const double rho = U.states[i].rho;
                const double u = U.getVelocity(U.states[i]);
                const double p = U.getPressure(U.states[i]);
                file << std::setprecision(4) << std::get<0>(elem) << ',' << i * dx << ',' << rho << ',' << u << ',' << p << std::endl;
            }
        }
    }
    
    file.close();
}

int main() {
    testEuler();
}