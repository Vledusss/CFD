#pragma once

#include <array>
#include <vector>
#include <utility>

using indexType = std::size_t;

template<indexType N>
struct Burgers
{
    std::array<double, N> state;

    explicit Burgers(const double lval = 1, const double rval = 0) noexcept {
        for (indexType i = 0; i < N / 2; ++i) { state[i] = lval; }
        for (indexType i = N / 2; i < N; ++i) { state[i] = rval; }
    }

    explicit Burgers(const std::array<double, N>& vals) noexcept : state(vals) {}

    Burgers(const std::vector<std::pair<indexType, double>>& states) {
        indexType start = 0;
        for (const auto& elem : states) {
            for (indexType i = start; i < elem.first + start; ++i) {
                state[i] = elem.second;
            }
            start += elem.first;
        }
    }

    double calcF(const indexType j, const double U) { 
        if (j == 0 || j == N) { return U * U / 2.; }
        if (state[j - 1] < state[j]) {
            if (state[j - 1] >= 0) { return state[j - 1] * state[j - 1] / 2.; }
            else if (state[j] <= 0) { return state[j] * state[j] / 2.; }
            else { return 0; }
        }
        return U * U / 2.; 
    }

    double calcU(const double UL, const double UR) {
        return UL + UR > 0 ? UL : UR;
    }
};