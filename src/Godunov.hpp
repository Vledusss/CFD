#pragma once

#include <array>
#include <cassert>

#include "Core" // eigen

using indexType = std::size_t;

template<typename Equation, indexType N>
struct Godunov {
    static void solve(Equation& eq, const double dx, const double dt, const double t) {
        assert(dx > 0);
        assert(dt > 0);

        const double ratio = dt / dx;

        std::array<double, N + 1> tempState;
        std::array<double, N + 1> F;

        tempState.front() = eq.state.front();
        tempState.back() = eq.state.back();

        for (indexType j = 0; j < N - 1; ++j) {
            tempState[j + 1] = eq.calcU(eq.state[j], eq.state[j + 1], (j + 1) * dx, t);
        }

        for (indexType j = 0; j < N + 1; ++j) {
            F[j] = eq.calcF(j, tempState[j]);
        }

        for (indexType j = 1; j <= N; ++j) {
            eq.state[j - 1] -= ratio * (F[j] - F[j - 1]);
        }
    }
};
