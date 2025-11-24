#pragma once 

#include <algorithm>

#include "Euler/RF/EulerFlux.hpp"

namespace Euler::RF {

struct State {
    double rho;
    double rho_u;
    double E;
    double rho_Yp;
    
    State(const double density = 0, const double velocity = 0, 
          const double pressure = 0, const double product = 0,
          const double reactionHeat = 1e6, 
          const double adiabaticIndex = 1.4, 
          const double Rspec = 287) 
          : rho(density), gamma(adiabaticIndex), R(Rspec), Q(reactionHeat) {
        rho_u = density * velocity;
        rho_Yp = density * product;
        E = pressure / (gamma - 1) + 0.5 * rho * velocity * velocity  + rho_Yp * Q;
    }

    
    double getGamma() const noexcept { return gamma; }

    double getR() const noexcept { return R; }

    double getQ() const noexcept { return Q; }
    
    double getVelocity() const {
        // assert(rho > 0);
        return rho_u / std::max(rho, 1e-10);
    }

    double getPressure() const {
        const double u = getVelocity();
        const double p = (gamma - 1) * (E - 0.5 * rho * u * u - rho_Yp * Q);
        // assert(p > 0);
        return std::max(p, 1e-10);
    }

    double getTemperature() const {
        return getPressure() / std::max(rho, 1e-10) / R;
    }

    State& operator-=(const Flux& flux) {
        const double eps = 1e-10; 

        const double temp_rho = rho - flux.density;
        const double temp_rho_u = rho_u - flux.momentum;
        const double temp_E = E - flux.energy;
        const double temp_rho_Yp = rho_Yp - flux.product;

        const double K = 0.5 * (temp_rho_u * temp_rho_u) / std::max(temp_rho, eps);
        const double e = temp_E - K; // ~ p
        E = e > eps / (gamma - 1) ? temp_E : eps / (gamma - 1) + K;

        rho = std::max(temp_rho, eps);
        rho_u = rho > eps ? temp_rho_u : 0.0;
        rho_Yp = rho > eps ? temp_rho_Yp : 0.0;
        rho_Yp = std::min(rho_Yp, rho);

        return *this;
    }

    State operator+(const State& other) const {
        return {
            rho + other.rho,
            rho_u + other.rho_u,
            E + other.E,
            rho_Yp + other.rho_Yp
        };
    }

    State operator-(const State& other) const {
        return {
            rho - other.rho,
            rho_u - other.rho_u,
            E - other.E,
            rho_Yp - other.rho_Yp
        };
    }
    
    Flux operator*(const double scalar) const {
        return Flux(
            rho * scalar, 
            rho_u * scalar, 
            E * scalar, 
            rho_Yp * scalar
        );
    }

    private:

    double gamma;
    double R;
    double Q;
};

}  // namespace Euler::RF
