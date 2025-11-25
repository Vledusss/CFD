#include <iostream>
#include <iomanip>
#include <fstream>

#include "Core.hpp"

#include "Burgers/BurgersEquation.hpp"
#include "Euler/Simple/EulerEquations.hpp"
#include "Euler/RF/EulerEquations.hpp"

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

void testEulerSimple() {
    const indexType N = 1000;                                           // по пространству
    std::vector<std::tuple<double, Euler::Simple::Equation<N>>> eq;     // gamma = 1.4
    
    Euler::Simple::Equation<N> initial;
    
    for (indexType i = 0; i < N / 2; ++i) { initial.states[i] = Euler::Simple::State(1.0, 100.0, 1.0); }  // НУ (ρ, u, p)
    for (indexType i = N / 2; i < N; ++i) { initial.states[i] = Euler::Simple::State(0.8, 0.0, 0.5); }    // НУ (ρ, u, p)
    
    const double dx = 1;
    
    const double startTime = 0;
    const double endTime = 15;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");

    if (file.is_open()) {
        file << "t,x,rho,u,p" << std::endl;
        Core<Euler::Simple::Equation<N>, N>::solve(eq, startTime, endTime, dx, 0.4);
        std::cout << "Recording..." << std::endl;

        for (const auto& elem : eq) {
            for (indexType i = 0; i < N; ++i) {
                const auto t = std::get<0>(elem);
                const auto U = std::get<1>(elem);

                const auto rho = U.states[i].rho;
                const auto u = U.states[i].getVelocity();
                const auto p = U.states[i].getPressure();

                file << std::setprecision(4) << t << ',' << i * dx 
                     << ',' << rho << ',' << u << ',' << p << std::endl;
            }
        }
    }
    
    file.close();

    std::cout << "Recording complete" << std::endl;
}

void testEulerRF() {
    const indexType N = 1000;                                       // по пространству
    const indexType IGNITION_ZONE = N / 2;
    std::vector<std::tuple<double, Euler::RF::Equation<N>>> eq;     // gamma = 1.4
    
    Euler::RF::Equation<N> initial;
    
    for (indexType i = 0; i < IGNITION_ZONE; ++i) { initial.states[i] = Euler::RF::State(1.0, 0.0, 430500, 1, 1e5); }    // НУ (ρ, u, p, Yp, Q) -> HLLC
    for (indexType i = IGNITION_ZONE; i < N; ++i) { initial.states[i] = Euler::RF::State(0.8, 0.0, 101325, 0, 1e5); }    // НУ (ρ, u, p, Yp, Q) -> HLLC

    // for (indexType i = 0; i < IGNITION_ZONE; ++i) { initial.states[i] = Euler::RF::State(1.0, 0.0, 184000, 0, 1e6); }    // НУ (ρ, u, p, Yp, Q) -> HLL
    // for (indexType i = IGNITION_ZONE; i < N; ++i) { initial.states[i] = Euler::RF::State(0.8, 0.0, 86100, 0, 1e6); }     // НУ (ρ, u, p, Yp, Q) -> HLL

    // for (indexType i = 0; i < IGNITION_ZONE; ++i) { initial.states[i] = Euler::RF::State(1.37, 704.0, 430500, 1, 1e6); }    // НУ (ρ, u, p, Yp, Q) -> HLLC 
    // for (indexType i = IGNITION_ZONE; i < N; ++i) { initial.states[i] = Euler::RF::State(0.8, 0.0, 101325, 0, 1e6); }       // НУ (ρ, u, p, Yp, Q) -> HLLC
    
    const double dx = 1;
    
    const double startTime = 0;
    const double endTime = 0.8;

    eq.emplace_back(std::make_tuple(startTime, initial));
    
    std::ofstream file("res.csv");

    if (file.is_open()) {
        file << "t,x,rho,u,p,Yp,T" << std::endl;
        Core<Euler::RF::Equation<N>, N>::solve(eq, startTime, endTime, dx, 0.4);
        std::cout << "Recording..." << std::endl;

        for (const auto& elem : eq) {
            for (indexType i = 0; i < N; ++i) {
                if (i * dx >= 200 && i * dx <= 800) {
                    const auto t = std::get<0>(elem);
                    const auto U = std::get<1>(elem);

                    const auto rho = U.states[i].rho;
                    const auto u = U.states[i].getVelocity();
                    const auto p = U.states[i].getPressure();
                    const auto T = U.states[i].getTemperature();
                    const auto Yp = U.states[i].rho_Yp / rho;

                    file << std::setprecision(4) << t << ',' << i * dx - 200 
                        << ',' << rho << ',' << u << ',' << p << ',' << Yp << ',' << T << std::endl;
                }
            }
        }
    }
    
    file.close();

    std::cout << "Recording complete" << std::endl;
}

int main() {
    testEulerRF();

    // testEulerSimple();

    // testBurgers();
}