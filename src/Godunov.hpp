#pragma once

#include <array>
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
    static void solve(Euler::EulerEquation<N>& eq, const double dx, const double dt) {
        assert(dx > 0);
        assert(dt > 0);

        const double ratio = dt / dx;

        std::array<Euler::EulerFlux, N + 1> F;

        // Граничные условия (простейшие - фиксированные)
        F.front() = eq.calcF(0, eq.state.front());
        F.back() = eq.calcF(0, eq.state.back());

        // Вычисление потоков на гранях
        for (indexType j = 1; j < N; ++j) {
            F[j] = eq.hllFlux(eq.state[j - 1], eq.state[j]);
        }

        // Обновление консервативных переменных
        for (indexType j = 1; j <= N; ++j) {
            eq.state[j - 1] -= (F[j] - F[j - 1]) * ratio;
        }
    }
};