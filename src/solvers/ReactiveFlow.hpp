#pragma once

#include "Eigen/Dense" 
#include "Eigen/LU"

namespace Solvers {

template<typename State>
struct ReactiveFlow {
    static State solve(const State& state, const double dt, const double tol = 1e-6) {
        const double A = 1e5;      // Фактор частоты
        const double Ea_R = 1e4;   // E_a / R (Энергия активации / Газовая постоянная)
        const double Q = 1.0e6;    // Тепловыделение (Дж/кг)

        State currState = state;

        const Eigen::Vector2d U = {state.E, state.rho_Yp};

        for (indexType i = 0; i < 5; ++i) {
            const double T = currState.getTemperature();
            const double rho = currState.rho;
            const double rho_Yp = currState.rho_Yp;
            const double gamma = currState.getGamma();
            const double R_spec = currState.getR();
            
            const double rho_Yf = rho - rho_Yp; // (Yf = 1 - Yp)
            const double exp = std::exp(-Ea_R / T);
            const double R_chem = A * std::max(0.0, rho_Yf) * exp; // A * ρYf * exp(-Ea/RT)
            
            const double omega_P = R_chem; 
            const double omega_E = R_chem * Q; 

            const Eigen::Vector2d Omega = {omega_E, omega_P};
            
            const double dTdE = (gamma - 1.0) / (rho * R_spec); // ∂T/∂E
            const double Arr = Ea_R / (T * T);                  // Ea/(R*T^2)
            const double R_chem_div_rhoYf = R_chem / std::max(1e-10, rho_Yf);

            Eigen::Matrix2d J;
            
            J(0, 0) = Q * R_chem * Arr * dTdE; 
            J(0, 1) = -Q * A * exp;
            J(1, 0) = R_chem * Arr * dTdE; 
            J(1, 1) = -A * exp;

            const Eigen::Matrix2d L = Eigen::Matrix2d::Identity() - dt * J; // L = I - dt * J
            const Eigen::Vector2d currU = {currState.E, currState.rho_Yp};
            
            const Eigen::Vector2d F = currU - U - dt * Omega;
            const Eigen::Vector2d deltaU = L.lu().solve(-F);
            
            currState.E += deltaU(0);
            currState.rho_Yp += deltaU(1);

            currState.rho_Yp = std::max(currState.rho_Yp, 0.0);
            currState.rho_Yp = std::min(currState.rho_Yp, rho);

            if (deltaU.norm() < tol) { break; }
        }

        const double eps = 1e-10;

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