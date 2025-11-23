#pragma once

#include <cmath>
#include <vector>

namespace Solvers {

template<typename State, typename Flux, typename Equation>
struct HLL {
    static Flux solve(const Equation& eq, const indexType i) {
        const State left = eq.states[i - 1];
        const State right = eq.states[i];

        const double pL = left.getPressure();
        const double pR = right.getPressure();
        
        const double uL = left.getVelocity();
        const double uR = right.getVelocity();
        
        const double cL = std::sqrt(left.getGamma() * pL / left.rho);
        const double cR = std::sqrt(right.getGamma() * pR / right.rho);
        
        const double SL = std::min(uL - cL, uR - cR);
        const double SR = std::max(uL + cL, uR + cR);
        
        const Flux FL = eq.calcFlux(left);
        const Flux FR = eq.calcFlux(right);
        
        if (SL >= 0) { return FL; } 
        else if (SR <= 0) { return FR; } 
        else { return (FL * SR - FR * SL + (right - left) * SL * SR) / (SR - SL); }
    }
};

}  // namespace Solvers
