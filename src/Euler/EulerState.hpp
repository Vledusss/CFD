#pragma once 

namespace Euler {
// Вектор состояния для Эйлера (ρ, ρu, ρE)
struct EulerState {
    double rho;     // плотность
    double rho_u;   // импульс  
    double rho_E;   // полная энергия
    
    EulerState(double r = 0, double ru = 0, double rE = 0) 
        : rho(r), rho_u(ru), rho_E(rE) {}
    
    // Арифметические операции для схемы
    EulerState& operator-=(const EulerState& other) {
        rho -= other.rho;
        rho_u -= other.rho_u;
        rho_E -= other.rho_E;
        return *this;
    }
    
    EulerState operator*(double scalar) const {
        return EulerState(rho * scalar, rho_u * scalar, rho_E * scalar);
    }
};

}  // namespace Euler
