#include <iostream>
#include <iomanip>
#include <fstream>

#include "Core.hpp"

#include "Burgers/BurgersEquation.hpp"
#include "Euler/EulerEquations.hpp"

void testBurgers() {
    const indexType N = 1000; // по пространству
    std::vector<std::tuple<double, Burgers::Equation<N>>> eq;
    
    Burgers::Equation<N> initial(-10, 10);
    
    const double dx = 1;
    
    const double startTime = 0;
    const double endTime = 10;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");

    if (file.is_open()) {
        file << "t,x,u" << std::endl;
        Core<Burgers::Equation<N>, N>::solve(eq, startTime, endTime, dx);
        std::cout << "Recording..." << std::endl;

        for (const auto& elem : eq) {
            for (indexType i = 0; i < N; ++i) {
                const auto t = std::get<0>(elem);
                const auto u = std::get<1>(elem).states[i];

                file << std::setprecision(4) << t << ',' << i * dx << ',' << u << std::endl;
            }
        }
    }
    
    file.close();

    std::cout << "Recording complete" << std::endl;
}

void testEuler() {
    const indexType N = 1000;                                   // по пространству
    std::vector<std::tuple<double, Euler::Equation<N>>> eq;     // gamma = 1.4
    
    Euler::Equation<N> initial;
    
    for (indexType i = 0; i < N / 2; ++i) { initial.states[i] = Euler::State(1.0, 100.0, 1.0); }  // НУ (ρ, u, p)
    for (indexType i = N / 2; i < N; ++i) { initial.states[i] = Euler::State(0.8, 0.0, 0.5); }    // НУ (ρ, u, p)
    
    const double dx = 1;
    
    const double startTime = 0;
    const double endTime = 15;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");

    if (file.is_open()) {
        file << "t,x,rho,u,p" << std::endl;
        Core<Euler::Equation<N>, N>::solve(eq, startTime, endTime, dx, 0.4);
        std::cout << "Recording..." << std::endl;

        for (const auto& elem : eq) {
            for (indexType i = 0; i < N; ++i) {
                const auto t = std::get<0>(elem);
                const auto U = std::get<1>(elem);

                const auto rho = U.states[i].rho;
                const auto u = U.getVelocity(U.states[i]);
                const auto p = U.getPressure(U.states[i]);

                file << std::setprecision(4) << t << ',' << i * dx 
                     << ',' << rho << ',' << u << ',' << p << std::endl;
            }
        }
    }
    
    file.close();

    std::cout << "Recording complete" << std::endl;
}

int main() {
    // testEuler();
    testBurgers();
}