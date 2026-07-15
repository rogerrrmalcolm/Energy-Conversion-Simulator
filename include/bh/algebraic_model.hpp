#pragma once

namespace bh {

// SI boundary model. Internally, the mass ratios are equivalent to G=c=1.
struct MassRange {
    double lower_kg{};
    double central_kg{};
    double upper_kg{};
};

struct SpinRange {
    double lower{};
    double central{};
    double upper{};
};

struct RotationalEnergyInput {
    double mass_kg{};
    double dimensionless_spin{};
    SpinRange spin_uncertainty{};
};

struct RotationalEnergyResult {
    double mass_energy_joules{};
    double irreducible_mass_kg{};
    double irreducible_mass_fraction{};
    double rotational_energy_joules{};
    double rotational_fraction{};
    double rotational_energy_lower_joules{};
    double rotational_energy_upper_joules{};
    double d_rotational_energy_d_spin_joules{};
};

struct RotationalEnergyRangeInput {
    MassRange mass{};
    SpinRange spin{};
};

struct RotationalEnergyRangeResult {
    RotationalEnergyResult lower{};
    RotationalEnergyResult central{};
    RotationalEnergyResult upper{};
};

[[nodiscard]] RotationalEnergyResult rotational_energy(double mass_kg,
                                                        double dimensionless_spin);
[[nodiscard]] RotationalEnergyResult rotational_energy(const RotationalEnergyInput& input);
[[nodiscard]] RotationalEnergyRangeResult rotational_energy_range(
    const RotationalEnergyRangeInput& input);

[[nodiscard]] double irreducible_mass_fraction(double dimensionless_spin);
[[nodiscard]] double rotational_energy_fraction(double dimensionless_spin);
}
