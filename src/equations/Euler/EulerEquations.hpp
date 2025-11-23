#pragma once

#include <cmath>
#include <vector>

#include "EulerState.hpp"

using indexType = size_t;

namespace Euler {

template<indexType N>
struct Equation {
    std::array<State, N> states;
    
    Equation(const double gamma = 1.4) {
        for (auto& state : states) { state = State(0, 0, 0, gamma); }
    }
    
    Flux calcFlux(const State& state) const {
        const double p = state.getPressure();
        const double u = state.getVelocity();
        
        return Flux(
            state.rho_u,              // ρu
            state.rho_u * u + p,      // ρu² + p  
            u * (state.E + p)         // u(E + p)
        );
    }

    void smooth(const double alpha = 0.1) {
    std::array<State, N> old_states = states;
    
    for (indexType i = 1; i < N - 1; ++i) {
        states[i].rho = (1.0 - alpha) * old_states[i].rho + 
                      0.5 * alpha * (old_states[i - 1].rho + old_states[i + 1].rho);
        
        states[i].rho_u = (1.0 - alpha) * old_states[i].rho_u + 
                        0.5 * alpha * (old_states[i - 1].rho_u + old_states[i + 1].rho_u);
        
        states[i].E = (1.0 - alpha) * old_states[i].E + 
                        0.5 * alpha * (old_states[i - 1].E + old_states[i + 1].E);
        }
    }
};

}  // namespace Euler
