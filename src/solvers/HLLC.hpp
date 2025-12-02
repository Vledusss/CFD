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

        const double pL = left.getPressure();
        const double pR = right.getPressure();
        
        const double rhoL = left.rho;
        const double rhoR = right.rho;

        const bool highDensity = std::max(rhoL / rhoR, rhoR / rhoL) > 10;
        const bool highPressure = std::max(pL / pR, pR / pL) > 10;

        if (highDensity || highPressure) { 
            // std::cout << "HLL fallback!" << std::endl;
            return HLL<State, Flux, Equation>::solve(eq, i); 
        }
        
        const double uL = left.getVelocity();
        const double uR = right.getVelocity();
        
        const double cL = std::sqrt(left.getGamma() * pL / rhoL);
        const double cR = std::sqrt(right.getGamma() * pR / rhoR);
        
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
            State UStarL = left;
            UStarL *= (SL - uL) / (SL - SM);
            UStarL.rho_u = UStarL.rho * SM;
            UStarL.E = UStarL.rho * (left.E / rhoL + (SM - uL) * (SM + pL / rhoL / (SL - uL)));

            return FL + (UStarL - left) * SL; 
        }
        else {
            State UStarR = right;
            UStarR *= (SR - uR) / (SR - SM);
            UStarR.rho_u = UStarR.rho * SM;
            UStarR.E = UStarR.rho * (right.E / rhoR + (SM - uR) * (SM + pR / rhoR / (SR - uR)));

            return FR + (UStarR - right) * SR; 
        }
    }
};

}  // namespace Solvers
