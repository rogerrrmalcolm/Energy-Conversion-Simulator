#include "bh/plasma_model.hpp"
#include "bh/constants.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bh {
PlasmaResult estimate_plasma_extraction(const PlasmaInput& in) {
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
    return {sigma, alfven, poynting, coupling, poynting*coupling*in.duration_seconds};
}
}
