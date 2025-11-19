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
    std::array<State, N> states;
    double gamma;
    
    EulerEquation(double g = 1.4) : gamma(g) {
        for (auto& state : states) { state = State(0, 0, 0); }
    }
    
    double getVelocity(const State& U) const {
        return U.rho_u / U.rho;
    }

    double getPressure(const State& state) const {
        const double vel = getVelocity(state);
        const double p = (gamma - 1) * (state.rho_E - 0.5 * state.rho * vel * vel);
        return std::isnan(p) || p <= 0 ? 1e-10 : p;
    }
    
    Flux calcF(indexType, const State& state) const {
        const double p = getPressure(state);
        
        return Flux(
            state.rho_u,                            // ρu
            state.rho_u * state.rho_u + p,          // ρu² + p  
            state.rho_u * (state.rho_E + p)         // u(ρE + p)
        );
    }
    
    Flux hllFlux(const State& left, const State& right, 
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
        
        const Flux FL = calcF(0, left);
        const Flux FR = calcF(0, right);
        
        if (SL >= 0) {
            return FL;
        } else if (SR <= 0) {
            return FR;
        } else {
            const double factor = 1.0 / (SR - SL);
            return (FL * SR - FR * SL + (right - left) * SL * SR) * factor;
        }
    }

    void smooth(const double alpha = 0.1) {
    std::array<State, N> old_state = states;
    
    for (indexType i = 1; i < N - 1; ++i) {
        states[i].rho = (1.0 - alpha) * old_state[i].rho + 
                      0.5 * alpha * (old_state[i-1].rho + old_state[i+1].rho);
        
        states[i].rho_u = (1.0 - alpha) * old_state[i].rho_u + 
                        0.5 * alpha * (old_state[i-1].rho_u + old_state[i+1].rho_u);
        
        states[i].rho_E = (1.0 - alpha) * old_state[i].rho_E + 
                        0.5 * alpha * (old_state[i-1].rho_E + old_state[i+1].rho_E);
        }
    }
};

}  // namespace Euler
