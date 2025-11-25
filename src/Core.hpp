#pragma once

#include <tuple>
#include <array>

#include "Burgers/BurgersEquation.hpp"
#include "Euler/Simple/EulerEquations.hpp"
#include "Euler/RF/EulerEquations.hpp"

#include "solvers/LaxFriedrichs.hpp"
#include "solvers/HLL.hpp"
#include "solvers/HLLC.hpp"
#include "solvers/ReactiveFlow.hpp"

using indexType = std::size_t;

template<typename Equation, typename Solver, indexType N>
struct Core {
    static void solve(std::vector<std::tuple<double, Equation>>& eq, 
                      const double startTime, const double endTime,
                      const double dx, const double CFL);
};

template<indexType N>
struct Core<Burgers::Equation<N>, Solvers::LaxFriedrichs<double, double, Burgers::Equation<N>>, N> {
    static void solve(std::vector<std::tuple<double, Burgers::Equation<N>>>& eq, 
                      const double startTime, const double endTime,
                      const double dx, const double CFL = 0.5) {
        assert(startTime < endTime);
        assert(dx > 0); 

        std::array<double, N + 1> F;
        std::array<double, N + 1> edgeStates;

        double t = startTime; 

        while (t < endTime) {
            Burgers::Equation<N> solution = std::get<1>(eq.back());;

            edgeStates.fill(0);
            F.fill(0);

            edgeStates.front() = solution.states.front(); // ГУ
            edgeStates.back() = solution.states.back();   // ГУ

            double maxVelocity = 0.0;

            for (const auto& u : solution.states) {
                maxVelocity = std::max(maxVelocity, u);
            }

            const double timeStep = CFL * dx / maxVelocity;
            std::cout << t << ' ' << maxVelocity << ' ' << timeStep << std::endl;
            assert(timeStep > 1e-8);

            for (indexType i = 0; i < N - 1; ++i) {
                edgeStates[i + 1] = solution.calcU(solution.states[i], solution.states[i + 1]);
            }

            using LFSolver = Solvers::LaxFriedrichs<double, double, Burgers::Equation<N>>;

            for (indexType i = 0; i < N + 1; ++i) {
                // F[i] = solution.calcFlux(edgeStates[i], i); // точное решение
                F[i] = LFSolver::solve(solution, i, dx, timeStep);
            }

            F[0] = solution.calcFlux(edgeStates[0], 0); // ГУ для LF
            F[N] = solution.calcFlux(edgeStates[N], N); // ГУ для LF

            for (indexType i = 0; i < N; ++i) {
                solution.states[i] -= timeStep / dx * (F[i + 1] - F[i]);
            }

            eq.emplace_back(std::make_tuple(t, solution));

            t += timeStep;
        }

        for (indexType i = 0; i < 26; ++i) {
            std::cout << '-';
        }
        std::cout << std::endl;
    }
};

template<typename Solver, indexType N>
struct Core<Euler::Simple::Equation<N>, Solver, N> {
    static void solve(std::vector<std::tuple<double, Euler::Simple::Equation<N>>>& eq, 
                      const double startTime, const double endTime,
                      const double dx, const double CFL = 0.5) {
        assert(startTime < endTime);
        assert(dx > 0);

        std::array<Euler::Simple::Flux, N + 1> F;    
        
        double t = startTime;  
        
        while (t < endTime) {
            Euler::Simple::Equation<N> solution = std::get<1>(eq.back());
            F.fill(Euler::Simple::Flux(0, 0, 0));

            F.front() = solution.calcFlux(solution.states.front()); // перенос из центра
            F.back() = solution.calcFlux(solution.states.back());   // перенос из центра

            double maxVelocity = 0.0;

            for (const auto& state : solution.states) {
                const double c = std::sqrt(state.getGamma() * state.getPressure() / state.rho);
                maxVelocity = std::max(maxVelocity, c + std::abs(state.getVelocity()));
            }

            const double timeStep = CFL * dx / maxVelocity;
            std::cout << t << ' ' << maxVelocity << ' ' << timeStep << std::endl;
            assert(timeStep > 1e-8);

            for (indexType i = 1; i < N; ++i) {
                F[i] = Solver::solve(solution, i);
            }

            for (indexType i = 0; i < N; ++i) {
                solution.states[i] -= (F[i + 1] - F[i]) * timeStep / dx;
            }

            eq.emplace_back(std::make_tuple(t, solution));

            // solution.smooth();

            t += timeStep;
        }

        for (indexType i = 0; i < 26; ++i) {
            std::cout << '-';
        }
        std::cout << std::endl;
    }
};


template<typename Solver, indexType N>
struct Core<Euler::RF::Equation<N>, Solver, N> {
    static void solve(std::vector<std::tuple<double, Euler::RF::Equation<N>>>& eqs, 
                      const double startTime, const double endTime,
                      const double dx, const double CFL = 0.5) {
        assert(startTime < endTime);
        assert(dx > 0);

        std::array<Euler::RF::Flux, N + 1> F;    
        
        double t = startTime;  
        
        while (t < endTime) {
            Euler::RF::Equation<N + 2> solution;

            const auto eq = std::get<1>(eqs.back());
            std::copy(eq.states.begin(), eq.states.end(), solution.states.begin() + 1);
            
            solution.states.front() = eq.states.front();
            solution.states.back() = eq.states.back();

            F.fill(Euler::RF::Flux(0, 0, 0));

            double maxVelocity = 0.0;

            for (const auto& state : solution.states) {
                const double c = std::sqrt(state.getGamma() * state.getPressure() / state.rho);
                maxVelocity = std::max(maxVelocity, c + std::abs(state.getVelocity()));
            }

            const double timeStep = CFL * dx / maxVelocity;
            std::cout << t << ' ' << maxVelocity << ' ' << timeStep << std::endl;
            assert(timeStep > 1e-8);

            using RFSolver = Solvers::ReactiveFlow<Euler::RF::State>;

            for (indexType i = 0; i < N + 1; ++i) {
                F[i] = Solver::solve(solution, i + 1);
            }

            for (indexType i = 0; i < N; ++i) {
                solution.states[i + 1] -= (F[i + 1] - F[i]) * timeStep / dx;
            }

            for (indexType i = 1; i < N + 1; ++i) {
                solution.states[i] = RFSolver::solve(solution.states[i], timeStep);
            }

            Euler::RF::Equation<N> result;
            std::copy(solution.states.begin() + 1, solution.states.end() - 1, result.states.begin());
            eqs.emplace_back(std::make_tuple(t, result));

            // solution.smooth();

            t += timeStep;
        }

        for (indexType i = 0; i < 26; ++i) {
            std::cout << '-';
        }
        std::cout << std::endl;
    }
};