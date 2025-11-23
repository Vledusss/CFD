#pragma once

#include <array>
#include <vector>
#include <utility>

using indexType = std::size_t;

namespace Burgers {

template<indexType N>
struct Equation
{
    std::array<double, N> states;

    explicit Equation(const double lval = 1, const double rval = 0) noexcept {
        for (indexType i = 0; i < N / 2; ++i) { states[i] = lval; }
        for (indexType i = N / 2; i < N; ++i) { states[i] = rval; }
    }

    explicit Equation(const std::array<double, N>& vals) noexcept : states(vals) {}

    Equation(const std::vector<std::pair<indexType, double>>& zones) {
        indexType start = 0;
        for (const auto& zone : zones) {
            for (indexType i = start; i < zone.first + start; ++i) {
                states[i] = zone.second;
            }
            start += zone.first;
        }
    }

    double calcFlux(const double U, const indexType j = 0) const { 
        if (j == 0 || j == N) { return U * U / 2.; }
        if (states[j - 1] < states[j]) {
            if (states[j - 1] >= 0) { return states[j - 1] * states[j - 1] / 2.; }
            else if (states[j] <= 0) { return states[j] * states[j] / 2.; }
            else { return 0; }
        }
        return U * U / 2.; 
    }

    static double calcU(const double UL, const double UR) {
        return UL + UR > 0 ? UL : UR;
    }
};

}  // namespace Burgers
