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
        return std::isnan(p) || p <= 0 ? 1e-10 : p;
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
        
        const double aL = std::sqrt(std::abs(gamma * pL / left.rho));
        const double aR = std::sqrt(std::abs(gamma * pR / right.rho));
        
        const double SL = std::min(uL - aL, uR - aR);
        const double SR = std::max(uL + aL, uR + aR);

        if (std::max(SR, -SL) * timeStep / dx > 1) {
            timeStep = std::max(CFL * dx / std::max(SR, -SL), 1e-3);
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

    void smooth(const double alpha = 0.1) {
    std::array<EulerState, N> old_state = state;
    
    for (indexType i = 1; i < N - 1; ++i) {
        state[i].rho = (1.0 - alpha) * old_state[i].rho + 
                      0.5 * alpha * (old_state[i-1].rho + old_state[i+1].rho);
        
        state[i].rho_u = (1.0 - alpha) * old_state[i].rho_u + 
                        0.5 * alpha * (old_state[i-1].rho_u + old_state[i+1].rho_u);
        
        state[i].rho_E = (1.0 - alpha) * old_state[i].rho_E + 
                        0.5 * alpha * (old_state[i-1].rho_E + old_state[i+1].rho_E);
        }
    }
};

}  // namespace Euler
