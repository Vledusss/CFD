#pragma once 

#include <cassert>

namespace Euler::Simple {

struct Flux {
    double density;           // ρu
    double momentum;          // ρu² + p
    double energy;            // u(E + p)

    Flux operator+(const Flux& other) const {
        return {
            density + other.density,
            momentum + other.momentum,
            energy + other.energy
        };
    }

    Flux operator-(const Flux& other) const {
        return {
            density - other.density,
            momentum - other.momentum,
            energy - other.energy
        };
    }
    
    Flux operator*(const double scalar) const {
        return {
            density * scalar, 
            momentum * scalar, 
            energy * scalar 
        };
    }

    Flux operator/(const double scalar) const {
        assert(scalar != 0);
        return {
            density / scalar, 
            momentum / scalar, 
            energy / scalar 
        };
    }
};

}  // namespace Euler::Simple
