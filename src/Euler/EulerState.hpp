#pragma once 

#include "EulerFlux.hpp"

namespace Euler {
struct State {
    double rho;             // ρ
    double rho_u;           // ρu  
    double rho_E;           // ρE
    
    State(const double density = 0, const double velocity = 0, 
          const double pressure = 0, const double gamma = 1.4) : rho(density) {
        rho_u = density * velocity;
        const double E = pressure / ((gamma - 1) * rho) + 0.5 * velocity * velocity;
        rho_E = rho * E;
    }

    State& operator-=(const Flux& flux) {
        rho -= flux.density;
        rho_u -= flux.momentum;
        rho_E -= flux.energy;
        return *this;
    }
    
    State operator*(const double scalar) const {
        return State(rho * scalar, rho_u * scalar, rho_E * scalar);
    }
};

}  // namespace Euler
