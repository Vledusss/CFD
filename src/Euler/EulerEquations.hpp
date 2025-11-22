#pragma once

#include <cmath>
#include <vector>

#include "EulerState.hpp"

using indexType = size_t;

namespace Euler {

template<indexType N>
class Equation {
private:
    double gamma;
public:
    std::array<State, N> states;
    
    Equation(double g = 1.4) : gamma(g) {
        for (auto& state : states) { state = State(0, 0, 0); }
    }
    
    static double getVelocity(const State& state) {
        // assert(state.rho > 0);
        return state.rho_u / std::max(state.rho, 1e-10);
    }

    double getPressure(const State& state) const {
        const double vel = getVelocity(state);
        const double p = (gamma - 1) * (state.E - 0.5 * state.rho * vel * vel);
        // assert(p > 0);
        return std::max(p, 1e-10);
    }

    double getGamma() const noexcept { return gamma; }
    
    Flux calcFlux(const State& state) const {
        const double p = getPressure(state);
        
        return Flux(
            state.rho_u,                            // ρu
            state.rho_u * state.rho_u + p,          // ρu² + p  
            state.rho_u * (state.E + p)         // u(ρE + p)
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
