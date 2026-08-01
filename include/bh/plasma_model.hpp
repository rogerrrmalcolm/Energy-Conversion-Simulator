#pragma once

namespace bh {
// Transparent 0-D, ideal-MHD-inspired transport scaling; this does not solve MHD or GRMHD.
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
    // Raw electromagnetic flux scaling for a field perpendicular to the representative flow.
    double poynting_power_watts{};
    // A visible heuristic factor, not a measured physical conversion efficiency.
    double spin_coupling_efficiency{};
    double outward_electromagnetic_power_watts{};
    double outward_electromagnetic_energy_joules{};
};

[[nodiscard]] PlasmaResult estimate_toy_plasma_transport(const PlasmaInput& input);

// Compatibility name retained for callers built against the original API. The result is outward
// electromagnetic transport, not an astrophysical extraction or usable delivered energy estimate.
[[nodiscard]] PlasmaResult estimate_plasma_extraction(const PlasmaInput& input);
}
