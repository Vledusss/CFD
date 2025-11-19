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
        rho_E = pressure / (gamma - 1) + 0.5 * rho * velocity * velocity;
    }

    State& operator-=(const Flux& flux) {
        rho -= flux.density;
        rho_u -= flux.momentum;
        rho_E -= flux.energy;
        return *this;
    }

    State operator+(const State& other) const {
        return State(
            rho + other.rho,
            rho_u + other.rho_u,
            rho_E + other.rho_E
        );
    }

    State operator-(const State& other) const {
        return State(
            rho - other.rho,
            rho_u - other.rho_u,
            rho_E - other.rho_E
        );
    }
    
    Flux operator*(const double scalar) const {
        return Flux(rho * scalar, rho_u * scalar, rho_E * scalar);
    }
};

}  // namespace Euler
