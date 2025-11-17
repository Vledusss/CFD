#pragma once 

namespace Euler {
// Поток для Эйлера
struct EulerFlux {
    double f_rho;             // ρu
    double f_momentum;        // ρu² + p
    double f_energy;          // u(ρE + p)
    
    EulerFlux(double fr = 0, double fm = 0, double fe = 0)
        : f_rho(fr), f_momentum(fm), f_energy(fe) {}
    
    EulerFlux& operator-=(const EulerFlux& other) {
        f_rho -= other.f_rho;
        f_momentum -= other.f_momentum;
        f_energy -= other.f_energy;
        return *this;
    }

    EulerFlux operator-(const EulerFlux& other) const {
        return EulerFlux(
            f_rho - other.f_rho,
            f_momentum - other.f_momentum,
            f_energy - other.f_energy
        );
    }
    
    EulerFlux operator*(double scalar) const {
        return EulerFlux(f_rho * scalar, f_momentum * scalar, f_energy * scalar);
    }
};

}  // namespace Euler
