#pragma once

#include <cmath>
#include <vector>
#include <algorithm>

#include "HLL.hpp"

namespace Solvers {

template<typename State, typename Flux, typename Equation>
struct HLLC {
    static Flux solve(const Equation& eq, const indexType i) {
        const State left = eq.states[i - 1];
        const State right = eq.states[i];

        const double pL = eq.getPressure(left);
        const double pR = eq.getPressure(right);
        
        const double rhoL = left.rho;
        const double rhoR = right.rho;

        const bool highDensity = std::max(rhoL / rhoR, rhoR / rhoL) > 10;
        const bool highPressure = std::max(pL / pR, pR / pL) > 10;

        if (highDensity || highPressure) { 
            // std::cout << "HLL fallback!" << std::endl;
            return HLL<State, Flux, Equation>::solve(eq, i); 
        }
        
        const double uL = eq.getVelocity(left);
        const double uR = eq.getVelocity(right);
        
        const double cL = std::sqrt(eq.getGamma() * pL / left.rho);
        const double cR = std::sqrt(eq.getGamma() * pR / right.rho);
        
        const double SL = std::min(uL - cL, uR - cR);
        const double SR = std::max(uL + cL, uR + cR);

        const double denom = rhoL * (SL - uL) - rhoR * (SR - uR);
        const double SM = (pR - pL + rhoL * uL * (SL - uL) - rhoR * uR * (SR - uR)) / denom;

        const double pStarL = pL + rhoL * (SL - uL) * (SM - uL); 
        const double pStarR = pR + rhoR * (SR - uR) * (SM - uR); 
        const double pStar = (pStarL + pStarR) / 2.;
        
        const Flux FL = eq.calcFlux(left);
        const Flux FR = eq.calcFlux(right);

        if (SL >= 0) { return FL; } 
        else if (SR <= 0) { return FR; } 
        else if (SL <= 0 && SM >= 0) {
            State UStarL;
            const double factor = rhoL * (SL - uL) / (SL - SM);
            UStarL.rho = factor;
            UStarL.rho_u = factor * SM;
            UStarL.E = factor * (left.E / rhoL + (SM - uL) * (SM + pL / rhoL / (SL - uL)));

            return FL + (UStarL - left) * SL; 
        }
        else {
            State UStarR;
            const double factor = rhoR * (SR - uR) / (SR - SM);
            UStarR.rho = factor;
            UStarR.rho_u = factor * SM;
            UStarR.E = factor * (right.E / rhoR + (SM - uR) * (SM + pR / rhoR / (SR - uR)));

            return FR + (UStarR - right) * SR; 
        }
    }
};

}  // namespace Solvers
