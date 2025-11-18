#pragma once

#include <vector>
#include <cassert>

#include "EulerState.hpp"
#include "EulerFlux.hpp"

using indexType = size_t;

namespace Euler {

template<indexType N>
class EulerEquation {
public:
    std::array<EulerState, N> state;
    double gamma;
    
    EulerEquation(double g = 1.4) : gamma(g) {
        for (auto& s : state) { s = EulerState(0, 0, 0); }
    }
    
    double getVelocity(const EulerState& U) const {
        return U.rho_u / U.rho;
    }

    double getPressure(const EulerState& U) const {
        const double vel = getVelocity(U);
        const double p = (gamma - 1.0) * (U.rho_E - 0.5 * U.rho * vel * vel);
        return std::max(p, 1e-10);
    }
    
    EulerFlux calcF(indexType, const EulerState& U) const {
        const double p = getPressure(U);
        
        return EulerFlux(
            U.rho_u,                             // ρu
            U.rho_u * U.rho_u + p,               // ρu² + p  
            U.rho_u * (U.rho_E + p)              // u(ρE + p)
        );
    }
    
EulerFlux hllFlux(const EulerState& left, const EulerState& right, 
                  const double dx, double& timeStep, 
                  const double CFL = 0.8) const {
    const double pL = getPressure(left);
    const double pR = getPressure(right);
    
    const double uL = getVelocity(left);
    const double uR = getVelocity(right);
    
    const double aL = std::sqrt(gamma * pL / left.rho);
    const double aR = std::sqrt(gamma * pR / right.rho);
    
    const double SL = std::min(uL - aL, uR - aR);
    const double SR = std::max(uL + aL, uR + aR);

    if (std::max(SR, -SL) * timeStep / dx > 1) {
        timeStep = CFL * dx / std::max(SR, -SL);
    }
    
    const EulerFlux FL = calcF(0, left);
    const EulerFlux FR = calcF(0, right);
    
    if (SL >= 0) {
        return FL;
    } else if (SR <= 0) {
        return FR;
    } else {
        const double factor = 1.0 / (SR - SL);
        return EulerFlux(
            (SR * FL.density - SL * FR.density + SL * SR * (right.rho - left.rho)) * factor,
            (SR * FL.momentum - SL * FR.momentum + SL * SR * (right.rho_u - left.rho_u)) * factor,
            (SR * FL.energy - SL * FR.energy + SL * SR * (right.rho_E - left.rho_E)) * factor
        );
    }
}
};

}  // namespace Euler
