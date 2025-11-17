#pragma once

#include <vector>

#include "EulerState.hpp"
#include "EulerFlux.hpp"

using indexType = size_t;

namespace Euler {

template<indexType N>
class EulerEquation {
public:
    // Массивы консервативных переменных
    std::array<EulerState, N> state;
    double gamma;  // показатель адиабаты
    
    EulerEquation(double g = 1.4) : gamma(g) {
        // Инициализация нулевым состоянием
        for (auto& s : state) {
            s = EulerState(0, 0, 0);
        }
    }
    
    // Восстановление примитивных переменных из консервативных
    double getVelocity(const EulerState& U) const {
        return U.rho_u / U.rho;  // ρu / ρ = u
    }

    double getPressure(const EulerState& U) const {
        const double vel = getVelocity(U);
        const double p = (gamma - 1.0) * (U.rho_E - 0.5 * U.rho * vel * vel);
        return std::max(p, 1e-10);  // избегаем отрицательного давления
    }
    
    double getEnergy(const EulerState& U) const {
        return U.rho_E / U.rho;  // ρE / ρ = E
    }
    
    // Вычисление потока F(u) для уравнений Эйлера
    EulerFlux calcF(indexType, const EulerState& U) const {
        const double p = getPressure(U);
        
        return EulerFlux(
            U.rho_u,                             // ρu
            U.rho_u * U.rho_u + p,               // ρu² + p  
            U.rho_u * (U.rho_E + p)              // u(ρE + p)
        );
    }
    
EulerFlux hllFlux(const EulerState& left, const EulerState& right) const {
    // Вычисляем скорости звука и характеристики
    const double pL = getPressure(left);
    const double pR = getPressure(right);
    
    const double uL = getVelocity(left);
    const double uR = getVelocity(right);
    
    const double aL = std::sqrt(gamma * pL / left.rho);
    const double aR = std::sqrt(gamma * pR / right.rho);
    
    // Оценки волновых скоростей (простой вариант)
    const double SL = std::min(uL - aL, uR - aR);
    const double SR = std::max(uL + aL, uR + aR);
    
    const EulerFlux FL = calcF(0, left);
    const EulerFlux FR = calcF(0, right);
    
    // HLL поток
    if (SL >= 0) {
        return FL;
    } else if (SR <= 0) {
        return FR;
    } else {
        // HLL центральная область
        const double factor = 1.0 / (SR - SL);
        return EulerFlux(
            (SR * FL.f_rho - SL * FR.f_rho + SL * SR * (right.rho - left.rho)) * factor,
            (SR * FL.f_momentum - SL * FR.f_momentum + SL * SR * (right.rho_u - left.rho_u)) * factor,
            (SR * FL.f_energy - SL * FR.f_energy + SL * SR * (right.rho_E - left.rho_E)) * factor
        );
    }
}
};

}  // namespace Euler
