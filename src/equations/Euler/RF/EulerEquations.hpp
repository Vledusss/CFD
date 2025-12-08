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
};

}  // namespace Euler::RF
