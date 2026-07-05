#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"

#include <cmath>
#include <stdexcept>

namespace bh {
RotationalEnergyResult rotational_energy(const double mass_kg, const double spin) {
    if (!std::isfinite(mass_kg) || mass_kg <= 0.0) {
        throw std::invalid_argument("mass must be finite and positive");
    }
    if (!std::isfinite(spin) || spin < 0.0 || spin >= 1.0) {
        throw std::invalid_argument("spin must satisfy 0 <= a_star < 1");
    }

    const double inner = std::sqrt(1.0 - spin * spin);
    const double irreducible_fraction = std::sqrt((1.0 + inner) / 2.0);
    const double mass_energy = mass_kg * speed_of_light_m_s * speed_of_light_m_s;
    return {mass_energy,
            mass_kg * irreducible_fraction,
            mass_energy * (1.0 - irreducible_fraction),
            1.0 - irreducible_fraction};
}
}
