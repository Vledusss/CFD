#pragma once

#include <cmath>
#include <vector>

namespace Solvers {

template<typename State, typename Flux, typename Equation>
struct HLL
{
    static Flux solve(const Equation& eq, const State& left, const State& right, 
                      const double dx, double& timeStep, const double CFL = 0.8) {
        const double pL = eq.getPressure(left);
        const double pR = eq.getPressure(right);
        
        const double uL = eq.getVelocity(left);
        const double uR = eq.getVelocity(right);
        
        const double aL = std::sqrt(std::abs(eq.getGamma() * pL / left.rho));
        const double aR = std::sqrt(std::abs(eq.getGamma() * pR / right.rho));
        
        const double SL = std::min(uL - aL, uR - aR);
        const double SR = std::max(uL + aL, uR + aR);

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
