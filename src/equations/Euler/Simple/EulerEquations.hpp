#pragma once

#include <cmath>
#include <vector>

#include "Euler/Simple/EulerState.hpp"

using indexType = size_t;

namespace Euler::Simple {

template<indexType N>
struct Equation {
    std::array<State, N> states;
    
    Flux calcFlux(const State& state) const {
        const double p = state.getPressure();
        const double u = state.getVelocity();
        
        return Flux(
            state.rho_u,              // ρu
            state.rho_u * u + p,      // ρu² + p  
            u * (state.E + p)         // u(E + p)
        );
    }
};

}  // namespace Euler::Simple
