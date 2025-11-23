#pragma once 

#include <cassert>

namespace Euler::RF {

struct Flux {
    double density;     // ρu
    double momentum;    // ρu² + p
    double energy;      // u(E + p)
    double product;     // ρuYp
    
    Flux(const double fr = 0, const double fm = 0, 
         const double fe = 0, const double fp = 0)
        : density(fr), momentum(fm), energy(fe), product(fp) {}

    Flux operator+(const Flux& other) const {
        return Flux(
            density + other.density,
            momentum + other.momentum,
            energy + other.energy,
            product + other.product
        );
    }

    Flux operator-(const Flux& other) const {
        return Flux(
            density - other.density,
            momentum - other.momentum,
            energy - other.energy,
            product - other.product
        );
    }
    
    Flux operator*(const double scalar) const {
        return Flux(
            density * scalar, 
            momentum * scalar, 
            energy * scalar, 
            product * scalar
        );
    }

    Flux operator/(const double scalar) const {
        assert(scalar != 0.0);
        return Flux(
            density / scalar, 
            momentum / scalar, 
            energy / scalar, 
            product / scalar
        );
    }
};

}  // namespace Euler::RF
