#pragma once

#include <array>
#include <vector>
#include <utility>

constexpr double PI = 3.14159265358979323846;
constexpr double S0 = 0.7420192566214;
constexpr double maxF = 0.8445788230378;

using indexType = std::size_t;

template<indexType N>
struct BuckleyLeverett
{
    std::array<double, N> state;

    explicit BuckleyLeverett(const double lval = 1, const double rval = 0) noexcept {
        for (indexType i = 0; i < N / 2; ++i) { state[i] = lval; }
        for (indexType i = N / 2; i < N; ++i) { state[i] = rval; }
    }

    explicit BuckleyLeverett(const std::array<double, N>& vals) noexcept : state(vals) {}

    BuckleyLeverett(const std::vector<std::pair<indexType, double>>& states) {
        indexType start = 0;
        for (const auto& elem : states) {
            for (indexType i = start; i < elem.first + start; ++i) {
                state[i] = elem.second;
            }
            start += elem.first;
        }
    }

    double calcB(const double S) const noexcept { return (1 - std::cos(PI * S)) / 2.; }

    double calcF(const indexType j, const double S) { 
        if (j == 0 || j == N) { return 0; }
        const double dS = state[j] - state[j - 1];
        const double W = (calcB(state[j]) - calcB(state[j - 1])) / (state[j]  - state[j - 1]);
        if (dS < 0) { return -dS < 1 - S0 ? maxF : calcB(S); }
        else { return dS < S0 ? maxF : calcB(S); }
    }

    double calcU(const double SL, const double SR, const double x, const double t) {
        return SL + SR > 0 ? SL : SR;
    }
};