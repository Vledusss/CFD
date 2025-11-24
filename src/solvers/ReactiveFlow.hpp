#pragma once

#include "Eigen/Dense" 
#include "Eigen/LU"

namespace Solvers {

template<typename State>
State solveChemistry(const State& state, const double dt) {
    const double A = 1e8;      // Фактор частоты
    const double Ea_R = 10000; // E_a / R (Энергия активации / Газовая постоянная)
    const double Q = 1.0e6;    // Тепловыделение (Дж/кг)
    
    const double T = state.getTemperature();
    const double rho = state.rho;
    const double rho_Yp = state.rho_Yp;
    const double gamma = state.getGamma();
    const double R_spec = state.getR();
    
    const double rho_Yf = rho - rho_Yp; // (Yf = 1 - Yp)
    const double R_chem = A * std::max(0.0, rho_Yf) * std::exp(-Ea_R / T); // A * ρYf * exp(-Ea/RT)
    
    const double omega_P = R_chem; 
    const double omega_E = R_chem * Q; 

    Eigen::Vector2d Omega;
    Omega(0) = omega_E;
    Omega(1) = omega_P;
    
    const double dTdE = (gamma - 1.0) / (rho * R_spec); // ∂T/∂E
    const double Arr = Ea_R / (T * T);             // Ea/(R*T^2)
    const double R_chem_div_rhoYf = R_chem / std::max(1e-10, rho_Yf);

    Eigen::Matrix2d J;
    
    J(0, 0) = Q * R_chem * Arr * dTdE; 
    J(0, 1) = -Q * R_chem_div_rhoYf;
    J(1, 0) = R_chem * Arr * dTdE; 
    J(1, 1) = -R_chem_div_rhoYf;

    Eigen::Matrix2d L = Eigen::Matrix2d::Identity() - dt * J; // L = I - dt * J
    
    Eigen::Vector2d RHS; // RHS = U* + dt * Omega
    RHS(0) = state.E + dt * Omega(0);
    RHS(1) = state.rho_Yp + dt * Omega(1);

    Eigen::Vector2d U_chem_new = L.lu().solve(RHS);
    
    State newState = state;
    newState.E = U_chem_new(0);
    newState.rho_Yp = std::max(0.0, U_chem_new(1));

    newState.rho_Yp = std::min(newState.rho_Yp, newState.rho); 

    const double eps = 1e-10;

    const double K = 0.5 * (newState.rho_u * newState.rho_u) / std::max(newState.rho, eps);
    const double e = newState.E - K; // ~ p
    newState.E = e > eps / (gamma - 1) ? newState.E : eps / (gamma - 1) + K;

    newState.rho = std::max(newState.rho, eps);
    newState.rho_u = newState.rho > eps ? newState.rho_u : 0.0;
    newState.rho_Yp = newState.rho > eps ? newState.rho_Yp : 0.0;
    newState.rho_Yp = std::min(newState.rho_Yp, rho);
    
    return newState;
}

}  // namespace Solvers