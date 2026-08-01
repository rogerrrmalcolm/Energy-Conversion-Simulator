#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace bh {


namespace {
void validate_spin(const double spin) {
    if (!std::isfinite(spin) || spin < 0.0 || spin >= 1.0) {
        throw std::invalid_argument("spin must satisfy 0 <= a_star < 1");
    }
}

void validate_mass(const MassRange& mass) {
    if (!std::isfinite(mass.lower_kg) ||
        !std::isfinite(mass.central_kg) ||
        !std::isfinite(mass.upper_kg)) {
        throw std::invalid_argument("mass values must be finite");
    }

    if (mass.lower_kg <= 0.0 ||
        mass.central_kg <= 0.0 ||
        mass.upper_kg <= 0.0) {
        throw std::invalid_argument("mass values must be positive");
    }

    if (mass.lower_kg > mass.central_kg ||
        mass.central_kg > mass.upper_kg) {
        throw std::invalid_argument("mass range must satisfy lower <= central <= upper");
    }
}

void validate_mass(const double mass_kg) {
    if (!std::isfinite(mass_kg)) {
        throw std::invalid_argument("mass value must be finite");
    }
    if (mass_kg <= 0.0) {
        throw std::invalid_argument("mass value must be positive");
    }
}


void validate_spin_range(const SpinRange& range) {
    validate_spin(range.lower);
    validate_spin(range.central);
    validate_spin(range.upper);
    if (range.lower > range.central || range.central > range.upper) {
        throw std::invalid_argument("spin uncertainty must satisfy lower <= central <= upper");
    }
}

void validate_result(const RotationalEnergyResult& result) {
    if (!std::isfinite(result.mass_energy_joules) ||
        !std::isfinite(result.irreducible_mass_kg) ||
        !std::isfinite(result.irreducible_mass_fraction) ||
        !std::isfinite(result.rotational_energy_joules) ||
        !std::isfinite(result.rotational_fraction) ||
        !std::isfinite(result.rotational_energy_lower_joules) ||
        !std::isfinite(result.rotational_energy_upper_joules) ||
        !std::isfinite(result.d_rotational_energy_d_spin_joules)) {
        throw std::overflow_error("rotational-energy calculation overflowed");
    }
}

bool nearly_equal(const double left, const double right) {
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <=
           16.0 * std::numeric_limits<double>::epsilon() * scale;
}

double one_minus_spin_squared(const double spin) {
    return (1.0 - spin) * (1.0 + spin);
}

double rotational_sensitivity_fraction_per_spin(const double spin) {
    if (spin == 0.0) {
        return 0.0;
    }

    const double root = std::sqrt(one_minus_spin_squared(spin));
    return spin / (4.0 * root * irreducible_mass_fraction(spin));
}
}

double irreducible_mass_fraction(const double spin) {
    validate_spin(spin);
    const double inner = std::sqrt(one_minus_spin_squared(spin));
    return std::sqrt((1.0 + inner) / 2.0);
}

double rotational_energy_fraction(const double spin) {
    return 1.0 - irreducible_mass_fraction(spin);
}

RotationalEnergyResult rotational_energy(const double mass_kg, const double spin) {
    return rotational_energy({mass_kg, spin, {spin, spin, spin}});
}

RotationalEnergyResult rotational_energy(const RotationalEnergyInput& input) {
    validate_mass(input.mass_kg);
    validate_spin(input.dimensionless_spin);
    validate_spin_range(input.spin_uncertainty);
    if (!nearly_equal(input.spin_uncertainty.central, input.dimensionless_spin)) {
        throw std::invalid_argument("spin uncertainty central value must match dimensionless_spin");
    }

    const double irreducible_fraction = irreducible_mass_fraction(input.dimensionless_spin);
    const double mass_energy = input.mass_kg * speed_of_light_m_s * speed_of_light_m_s;
    const double rotational_fraction = 1.0 - irreducible_fraction;

    const RotationalEnergyResult result{
        mass_energy,
        input.mass_kg * irreducible_fraction,
        irreducible_fraction,
        mass_energy * rotational_fraction,
        rotational_fraction,
        mass_energy * rotational_energy_fraction(input.spin_uncertainty.lower),
        mass_energy * rotational_energy_fraction(input.spin_uncertainty.upper),
        mass_energy * rotational_sensitivity_fraction_per_spin(input.dimensionless_spin)};
    validate_result(result);
    return result;
}
RotationalEnergyRangeResult rotational_energy_range(
    const RotationalEnergyRangeInput& input) {
    validate_mass(input.mass);
    validate_spin_range(input.spin);

    const auto lower = rotational_energy(
        input.mass.lower_kg, input.spin.lower);
    const auto central = rotational_energy(
        input.mass.central_kg, input.spin.central);
    const auto upper = rotational_energy(
        input.mass.upper_kg, input.spin.upper);

    return {lower,
            central,
            upper,
            central.rotational_energy_joules - lower.rotational_energy_joules,
            upper.rotational_energy_joules - central.rotational_energy_joules};
}
}
