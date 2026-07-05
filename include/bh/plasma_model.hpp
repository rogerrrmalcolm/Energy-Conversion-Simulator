#pragma once

namespace bh {
// Transparent 0-D ideal-MHD scaling model; this is not a GRMHD simulation.
struct PlasmaInput {
    double magnetic_field_tesla{};
    double mass_density_kg_m3{};
    double flow_area_m2{};
    double dimensionless_spin{};
    double duration_seconds{};
};

struct PlasmaResult {
    double magnetization{};
    double alfven_speed_m_s{};
    double poynting_power_watts{};
    double spin_coupling_efficiency{};
    double idealized_extracted_energy_joules{};
};

[[nodiscard]] PlasmaResult estimate_plasma_extraction(const PlasmaInput& input);
}
