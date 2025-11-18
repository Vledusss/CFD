#pragma once

#include <tuple>
#include <array>
#include <vector>
#include <cassert>

#include "Core" // eigen
#include "Euler/EulerEquations.hpp"

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
struct Godunov<Euler::EulerEquation<N>, N> {
    static void solve(std::vector<std::tuple<double, Euler::EulerEquation<N>>>& eq, 
                      const double startTime, const double endTime,
                      const double dx, const double dt) {
        assert(startTime < endTime);
        assert(dx > 0);
        assert(dt > 0);

        double t = startTime;  
        std::array<Euler::EulerFlux, N + 1> F;                 
        
        while (t <= endTime)
        {
            Euler::EulerEquation<N> solution = std::get<1>(eq.back());
            F.fill(Euler::EulerFlux(0, 0, 0));

            F.front() = solution.calcF(0, solution.state.front()); // ГУ
            F.back() = solution.calcF(0, solution.state.back());   // ГУ

            double timeStep = dt;

            // Вычисление потоков на гранях
            for (indexType j = 1; j < N; ++j) {
                F[j] = solution.hllFlux(solution.state[j - 1], solution.state[j], dx, timeStep);
            }

            const double ratio = timeStep / dx;

            // Обновление консервативных переменных
            for (indexType j = 0; j < N; ++j) {
                solution.state[j] -= (F[j + 1] - F[j]) * ratio;
            }

            eq.emplace_back(std::make_tuple(t, solution));

            t += timeStep;
        }
    }
};