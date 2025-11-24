#pragma once

#include "Eigen/Dense" 
#include "Eigen/LU"

namespace Solvers {

constexpr double A = 1e5;         // Фактор частоты
constexpr double Ea = 2.87e6;     // Энергия активации

template<typename State>
struct ReactiveFlow {
    static State solve(const State& state, const double dt, 
                       const double tol = 1e-12, 
                       const double eps = 1e-10) {
        const double rho = state.rho;
        const double gamma = state.getGamma();
        const double R = state.getR();
        const double Q = state.getQ();
        const double Ea_R = Ea / state.getR(); // ~ 10000К

        State currState = state;
        const Eigen::Vector2d U = {state.E, state.rho_Yp};

        for (indexType i = 0; i < 5; ++i) {
            const double T = currState.getTemperature();
            const double rho_Yp = currState.rho_Yp;
            
            const double rho_Yf = rho - rho_Yp; // (Yf = 1 - Yp)
            const double exp = std::exp(-Ea_R / T);
            const double R_chem = A * std::max(0.0, rho_Yf) * exp; // A * ρYf * exp(-Ea/RT)
             
            const double omega_E = R_chem * Q; 
            const double omega_Yp = R_chem;

            const Eigen::Vector2d Omega = {omega_E, omega_Yp};

            Eigen::Matrix2d J;
            
            J(1, 0) =  (gamma - 1.0) * R_chem * Ea_R / T / T / R / std::max(rho, eps); 
            J(1, 1) = -A * exp - Q * J(1, 0);
            J(0, 0) = Q * J(1, 0); 
            J(0, 1) = Q * J(1, 1);

            const Eigen::Matrix2d L = Eigen::Matrix2d::Identity() - dt * J; // L = I - dt * J
            const Eigen::Vector2d currU = {currState.E, currState.rho_Yp};
            
            const Eigen::Vector2d F = currU - U - dt * Omega;
            const Eigen::Vector2d deltaU = L.lu().solve(-F);
            
            currState.E += deltaU(0);
            currState.rho_Yp += deltaU(1);

            currState.rho_Yp = std::max(currState.rho_Yp, 0.0);
            currState.rho_Yp = std::min(currState.rho_Yp, rho);

            // std::cout << i << ' ' << deltaU.norm() << std::endl;
            
            if (deltaU.norm() < tol) { break; }
        }

        const double K = 0.5 * (currState.rho_u * currState.rho_u) / std::max(currState.rho, eps);
        const double e = currState.E - K; // ~ p
        currState.E = e > eps / (currState.getGamma() - 1) ? 
            currState.E : eps / (currState.getGamma() - 1) + K;

        currState.rho = std::max(currState.rho, eps);
        currState.rho_u = currState.rho > eps ? currState.rho_u : 0.0;
        
        return currState;
    }
};

}  // namespace Solvers