#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"

#include <cmath>
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
        mass.central_kg <= 0.0 || // the mass of a black hole must be greater than or equal to one another, if not there is no spin or error
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
    validate_spin(range.central); //what it meant was to take in those specific parameters and pass the range values to a different function that does the validating
    validate_spin(range.upper);
    if (range.lower > range.central || range.central > range.upper) { //if invalid clause
        throw std::invalid_argument("spin uncertainty must satisfy lower <= central <= upper");
    }
}

double rotational_sensitivity_fraction_per_spin(const double spin) {
    if (spin == 0.0) {
        return 0.0;
    }

    const double root = std::sqrt(1.0 - spin * spin);
    return spin / (4.0 * root * irreducible_mass_fraction(spin));
}
}

double irreducible_mass_fraction(const double spin) {
    validate_spin(spin);
    const double inner = std::sqrt(1.0 - spin * spin);
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
    if (input.spin_uncertainty.central != input.dimensionless_spin) {
        throw std::invalid_argument("spin uncertainty central value must match dimensionless_spin");
    }

    const double irreducible_fraction = irreducible_mass_fraction(input.dimensionless_spin);
    const double mass_energy = input.mass_kg * speed_of_light_m_s * speed_of_light_m_s;
    const double rotational_fraction = 1.0 - irreducible_fraction;

    return {mass_energy,
            input.mass_kg * irreducible_fraction,
            irreducible_fraction,
            mass_energy * rotational_fraction,
            rotational_fraction,
            mass_energy * rotational_energy_fraction(input.spin_uncertainty.lower),
            mass_energy * rotational_energy_fraction(input.spin_uncertainty.upper),
            mass_energy * rotational_sensitivity_fraction_per_spin(input.dimensionless_spin)};
}
RotationalEnergyRangeResult rotational_energy_range(
    const RotationalEnergyRangeInput& input) {
        
}
}
