#pragma once 

#include "EulerFlux.hpp"

namespace Euler {
// Вектор состояния для Эйлера (ρ, ρu, ρE)
struct EulerState {
    double rho;     // плотность
    double rho_u;   // импульс  
    double rho_E;   // полная энергия
    
    EulerState(const double density = 0, 
               const double velocity = 0, 
               const double pressure = 0, 
               const double gamma = 1.4) : rho(density) {
        rho_u = density * velocity;
        const double E = pressure / ((gamma - 1) * rho) + 0.5 * velocity * velocity;
        rho_E = rho * E;
    }
    
    // Арифметические операции для схемы
    EulerState& operator-=(const EulerFlux& flux) {
        rho -= flux.f_rho;
        rho_u -= flux.f_momentum;
        rho_E -= flux.f_energy;
        return *this;
    }
    
    EulerState operator*(double scalar) const {
        return EulerState(rho * scalar, rho_u * scalar, rho_E * scalar);
    }
};

}  // namespace Euler
