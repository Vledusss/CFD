#pragma once

namespace Solvers {
    
template<typename State, typename Flux, typename Equation>
struct LaxFriedrichs {
    static Flux solve(const Equation& eq, const indexType i, 
                      const double dx, const double dt) {
        const State left = eq.states[i - 1];
        const State right = eq.states[i];

        return ((eq.calcFlux(left) + eq.calcFlux(right)) - (right - left) * dx / dt) / 2.;
    }
};

}  // namespace Solvers
