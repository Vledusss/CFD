#pragma once

#include <cmath>
#include <vector>
#include <cassert>
#include <algorithm>

namespace Solvers {

template<typename State, typename Flux, typename Equation>
struct HLLC {
    static Flux solve(const Equation& eq, const State& left, const State& right, 
                      const double dx, double& timeStep, const double CFL = 0.8) {
        const double pL = eq.getPressure(left);
        const double pR = eq.getPressure(right);
        
        const double uL = eq.getVelocity(left);
        const double uR = eq.getVelocity(right);

        const double rhoL = left.rho;
        const double rhoR = right.rho;
        
        const double cL = std::sqrt(eq.getGamma() * pL / left.rho);
        const double cR = std::sqrt(eq.getGamma() * pR / right.rho);
        
        const double SL = std::min(uL - cL, uR - cR);
        const double SR = std::max(uL + cL, uR + cR);

        const double denom = rhoL * (SL - uL) - rhoR * (SR - uR);
        const double SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom;

        const double pStarL = pL + rhoL * (SL - uL) * (SM - uL); 
        const double pStarR = pR + rhoR * (SR - uR) * (SM - uR); 
        assert(pStarL == pStarR);
        const double pStar = pStarL;

        const double maxVelocity = std::max({std::abs(SR), std::abs(SL), std::abs(SM)});

        if (maxVelocity * timeStep / dx > 1) {
            timeStep = std::max(CFL * dx / maxVelocity, 1e-3);
        }
        
        const Flux FL = eq.calcFlux(left);
        const Flux FR = eq.calcFlux(right);

        if (SL >= 0) { return FL; } 
        else if (SR <= 0) { return FR; } 
        else if (SL <= 0 && SM >= 0) {
            State UStarL;
            const double factor = rhoL * (SL - uL) / (SL - SM);
            UStarL.rho = factor;
            UStarL.rho_u = factor * SM;
            UStarL.rho_E = factor * (left.rho_E / rhoL + (SM - uL) * (SM + pL / rhoL / (SL - uL)));

            return FL + (UStarL - left) * SL; 
        }
        else {
            State UStarR;
            const double factor = rhoR * (SR - uR) / (SR - SM);
            UStarR.rho = factor;
            UStarR.rho_u = factor * SM;
            UStarR.rho_E = factor * (right.rho_E / rhoR + (SM - uR) * (SM + pR / rhoR / (SR - uR)));

            return FR + (UStarR - right) * SR; 
        }
    }
};

}  // namespace Solvers
