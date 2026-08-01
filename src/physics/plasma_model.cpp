#include "bh/plasma_model.hpp"
#include "bh/constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bh {
PlasmaResult estimate_toy_plasma_transport(const PlasmaInput& in) {
    if (!std::isfinite(in.magnetic_field_tesla) ||
        !std::isfinite(in.mass_density_kg_m3) ||
        !std::isfinite(in.flow_area_m2) ||
        !std::isfinite(in.dimensionless_spin) ||
        !std::isfinite(in.duration_seconds) ||
        in.magnetic_field_tesla < 0.0 || in.mass_density_kg_m3 <= 0.0 ||
        in.flow_area_m2 < 0.0 || in.duration_seconds < 0.0 ||
        in.dimensionless_spin < 0.0 || in.dimensionless_spin >= 1.0) {
        throw std::invalid_argument("invalid plasma model input");
    }
    const double c2 = speed_of_light_m_s*speed_of_light_m_s;
    const double sigma = in.magnetic_field_tesla*in.magnetic_field_tesla /
                         (vacuum_permeability_si*in.mass_density_kg_m3*c2);
    const double alfven = speed_of_light_m_s*std::sqrt(sigma/(1.0+sigma));
    const double poynting = (in.magnetic_field_tesla*in.magnetic_field_tesla /
                            vacuum_permeability_si)*alfven*in.flow_area_m2;
    // Deliberately visible heuristic: spin availability times magnetic coupling.
    const double coupling = std::clamp(in.dimensionless_spin*in.dimensionless_spin *
                                       sigma/(1.0+sigma), 0.0, 1.0);
    const double outward_electromagnetic_power = poynting*coupling;
    const double outward_electromagnetic_energy =
        outward_electromagnetic_power*in.duration_seconds;
    if (!std::isfinite(sigma) || !std::isfinite(alfven) || !std::isfinite(poynting) ||
        !std::isfinite(coupling) || !std::isfinite(outward_electromagnetic_power) ||
        !std::isfinite(outward_electromagnetic_energy)) {
        throw std::overflow_error("plasma calculation overflowed");
    }
    return {sigma, alfven, poynting, coupling, outward_electromagnetic_power,
            outward_electromagnetic_energy};
}

PlasmaResult estimate_plasma_extraction(const PlasmaInput& input) {
    return estimate_toy_plasma_transport(input);
}
}
