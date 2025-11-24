#pragma once

#include <cmath>
#include <vector>

#include "Euler/RF/EulerState.hpp"

using indexType = size_t;

namespace Euler::RF {

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
            u * (state.E + p),        // u(E + p)
            u * state.rho_Yp          // ρuYp
        );
    }

    void smooth(const double alpha = 0.1) {
        std::array<State, N> oldStates = states;
        
        for (indexType i = 1; i < N - 1; ++i) {
            states[i].rho = (1.0 - alpha) * oldStates[i].rho + 
                0.5 * alpha * (oldStates[i - 1].rho + oldStates[i + 1].rho);
            
            states[i].rho_u = (1.0 - alpha) * oldStates[i].rho_u + 
                0.5 * alpha * (oldStates[i - 1].rho_u + oldStates[i + 1].rho_u);
            
            states[i].E = (1.0 - alpha) * oldStates[i].E + 
                0.5 * alpha * (oldStates[i - 1].E + oldStates[i + 1].E);

            states[i].rho_Yp = (1.0 - alpha) * oldStates[i].rho_Yp + 
                0.5 * alpha * (oldStates[i - 1].rho_Yp + oldStates[i + 1].rho_Yp);
            }
    }
};

}  // namespace Euler::RF
