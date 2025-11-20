#pragma once

#include <cmath>
#include <vector>

namespace Solvers {

template<typename State, typename Flux, typename Equation>
struct HLL {
    static Flux solve(const Equation& eq, const indexType i, const double dx,
                      double& timeStep, const double CFL = 0.5) {
        const State left = eq.states[i - 1];
        const State right = eq.states[i];

        const double pL = eq.getPressure(left);
        const double pR = eq.getPressure(right);
        
        const double uL = eq.getVelocity(left);
        const double uR = eq.getVelocity(right);
        
        const double cL = std::sqrt(eq.getGamma() * pL / left.rho);
        const double cR = std::sqrt(eq.getGamma() * pR / right.rho);
        
        const double SL = std::min(uL - cL, uR - cR);
        const double SR = std::max(uL + cL, uR + cR);

        const double maxVelocity = std::max(std::abs(SR), std::abs(SL));

        if (maxVelocity * timeStep / dx > 1) {
            timeStep = std::max(CFL * dx / maxVelocity, 1e-3);
        }
        
        const Flux FL = eq.calcFlux(left);
        const Flux FR = eq.calcFlux(right);
        
        if (SL >= 0) { return FL; } 
        else if (SR <= 0) { return FR; } 
        else { return (FL * SR - FR * SL + (right - left) * SL * SR) / (SR - SL); }
    }
};

}  // namespace Solvers
