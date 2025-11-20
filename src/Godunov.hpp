#pragma once

#include <tuple>
#include <array>

// #include "Core" // eigen
#include "Euler/EulerEquations.hpp"
#include "solvers/HLL.hpp"
#include "solvers/HLLC.hpp"

using indexType = std::size_t;

template<typename Equation, indexType N>
struct Godunov {
    static void solve(Equation& eq, const double dx, const double dt) {
        assert(dx > 0);
        assert(dt > 0);

        const double ratio = dt / dx;

        std::array<double, N + 1> tempState;
        std::array<double, N + 1> F;

        tempState.front() = eq.state.front(); // ГУ
        tempState.back() = eq.state.back();   // ГУ

        for (indexType j = 0; j < N - 1; ++j) {
            tempState[j + 1] = eq.calcU(eq.state[j], eq.state[j + 1]);
        }

        for (indexType j = 0; j < N + 1; ++j) {
            F[j] = eq.calcF(j, tempState[j]);
        }

        for (indexType j = 1; j <= N; ++j) {
            eq.state[j - 1] -= ratio * (F[j] - F[j - 1]);
        }
    }
};

template<indexType N>
struct Godunov<Euler::Equation<N>, N> {
    static void solve(std::vector<std::tuple<double, Euler::Equation<N>>>& eq, 
                      const double startTime, const double endTime,
                      const double dx, const double dt) {
        assert(startTime < endTime);
        assert(dx > 0);
        assert(dt > 0);

        double t = startTime;  
        std::array<Euler::Flux, N + 1> F;                 
        
        while (t < endTime)
        {
            Euler::Equation<N> solution = std::get<1>(eq.back());
            F.fill(Euler::Flux(0, 0, 0));

            solution.states[0] = solution.states[1];          // ГУ
            solution.states[N - 1] = solution.states[N - 2];  // ГУ

            F.front() = solution.calcFlux(solution.states.front()); // перенос из центра
            F.back() = solution.calcFlux(solution.states.back());   // перенос из центра

            double timeStep = dt;

            using HLLSolver = Solvers::HLL<Euler::State, Euler::Flux, Euler::Equation<N>>;
            using HLLCSolver = Solvers::HLLC<Euler::State, Euler::Flux, Euler::Equation<N>>;

            for (indexType i = 1; i < N; ++i) {
                F[i] = HLLSolver::solve(solution, solution.states[i - 1], solution.states[i], dx, timeStep);
            }

            for (indexType i = 0; i < N; ++i) {
                solution.states[i] -= (F[i + 1] - F[i]) * timeStep / dx;
                // if (solution.states[i].rho_E < 0) { std::cout << t << std::endl; }
            }

            // solution.smooth();

            eq.emplace_back(std::make_tuple(t, solution));

            t += timeStep;
        }
    }
};