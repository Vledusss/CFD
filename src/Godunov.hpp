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
                      const double dx, const double CFL = 0.2) {
        assert(startTime < endTime);
        assert(dx > 0);

        double t = startTime;  
        std::array<Euler::Flux, N + 1> F;                 
        
        while (t < endTime)
        {
            Euler::Equation<N> solution = std::get<1>(eq.back());
            F.fill(Euler::Flux(0, 0, 0));

            F.front() = solution.calcFlux(solution.states.front()); // перенос из центра
            F.back() = solution.calcFlux(solution.states.back());   // перенос из центра

            double maxVelocity = 0.0;

            for (const auto& state : solution.states) {
                const double c = std::sqrt(solution.getGamma() * solution.getPressure(state) / state.rho);
                maxVelocity = std::max(maxVelocity, c + std::abs(solution.getVelocity(state)));
            }

            const double timeStep = CFL * dx / maxVelocity;
            // std::cout << t << ' ' << maxVelocity << ' ' << timeStep << std::endl;

            using HLLSolver = Solvers::HLL<Euler::State, Euler::Flux, Euler::Equation<N>>;
            using HLLCSolver = Solvers::HLLC<Euler::State, Euler::Flux, Euler::Equation<N>>;

            for (indexType i = 1; i < N; ++i) {
                const double rhoL = solution.states[i - 1].rho;
                const double rhoR = solution.states[i].rho;

                const double pL = solution.getPressure(solution.states[i - 1]);
                const double pR = solution.getPressure(solution.states[i]);

                const bool densityCheck = std::max(rhoL / rhoR, rhoR / rhoL) > 5;
                const bool pressureCheck = std::max(pL / pR, rhoR / pL) > 3;

                F[i] = densityCheck || pressureCheck ?
                       HLLSolver::solve(solution, i) :
                       HLLCSolver::solve(solution, i);
            }

            for (indexType i = 0; i < N; ++i) {
                solution.states[i] -= (F[i + 1] - F[i]) * timeStep / dx;
            }

            eq.emplace_back(std::make_tuple(t, solution));

            // solution.smooth();

            t += timeStep;
        }
    }
};
