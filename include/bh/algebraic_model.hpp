#pragma once

namespace bh {

// SI boundary model. Internally, the mass ratios are equivalent to G=c=1.
struct RotationalEnergyResult {
    double mass_energy_joules{};
    double irreducible_mass_kg{};
    double rotational_energy_joules{};
    double rotational_fraction{};
};

[[nodiscard]] RotationalEnergyResult rotational_energy(double mass_kg,
                                                        double dimensionless_spin);
}
