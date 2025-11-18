#pragma once 

namespace Euler {
struct EulerFlux {
    double density;               // ρu
    double momentum;          // ρu² + p
    double energy;            // u(ρE + p)
    
    EulerFlux(double fr = 0, double fm = 0, double fe = 0)
        : density(fr), momentum(fm), energy(fe) {}
    
    EulerFlux& operator-=(const EulerFlux& other) {
        density -= other.density;
        momentum -= other.momentum;
        energy -= other.energy;
        return *this;
    }

    EulerFlux operator-(const EulerFlux& other) const {
        return EulerFlux(
            density - other.density,
            momentum - other.momentum,
            energy - other.energy
        );
    }
    
    EulerFlux operator*(const double scalar) const {
        return {density * scalar, momentum * scalar, energy * scalar};
    }
};

}  // namespace Euler
