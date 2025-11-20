#pragma once 

#include "EulerFlux.hpp"

namespace Euler {

struct State {
    double rho;
    double rho_u;
    double E;
    
    State(const double density = 0, const double velocity = 0, 
          const double pressure = 0, const double gamma = 1.4) : rho(density) {
        rho_u = density * velocity;
        E = pressure / (gamma - 1) + 0.5 * rho * velocity * velocity;
    }

    State& operator-=(const Flux& flux) {
        rho -= flux.density;
        rho_u -= flux.momentum;
        E -= flux.energy;
        return *this;
    }

    State operator+(const State& other) const {
        return State(
            rho + other.rho,
            rho_u + other.rho_u,
            E + other.E
        );
    }

    State operator-(const State& other) const {
        return State(
            rho - other.rho,
            rho_u - other.rho_u,
            E - other.E
        );
    }
    
    Flux operator*(const double scalar) const {
        return Flux(rho * scalar, rho_u * scalar, E * scalar);
    }
};

}  // namespace Euler
