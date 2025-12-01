#pragma once 

#include <cassert>

namespace Euler::RF {

struct Flux {
    double density;     // ρu
    double momentum;    // ρu² + p
    double energy;      // u(E + p)
    double product;     // ρuYp

    Flux operator+(const Flux& other) const {
        return {
            density + other.density,
            momentum + other.momentum,
            energy + other.energy,
            product + other.product
        };
    }

    Flux operator-(const Flux& other) const {
        return {
            density - other.density,
            momentum - other.momentum,
            energy - other.energy,
            product - other.product
        };
    }
    
    Flux operator*(const double scalar) const {
        return {
            density * scalar, 
            momentum * scalar, 
            energy * scalar, 
            product * scalar
        };
    }

    Flux operator/(const double scalar) const {
        assert(scalar != 0);
        return {
            density / scalar, 
            momentum / scalar, 
            energy / scalar, 
            product / scalar
        };
    }
};

}  // namespace Euler::RF
