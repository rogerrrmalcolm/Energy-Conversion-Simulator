#pragma once

#include "bh/kerr_geodesic.hpp"

#include <cstddef>
#include <string_view>

namespace bh {
// One fixed, equatorial split configuration. This is not a search or optimization state.
struct PenroseSplitParameters {
    double split_radius_over_m{};
    double incoming_lz_over_m_m{};
    double split_angle_rad{};
};

enum class PenroseEventStatus {
    outside_ergosphere,
    physics_invalid,
    captured_or_non_escaping,
    physically_feasible,
    integration_failed
};

struct EquatorialPenroseScenario {
    // Geometrized units: G = c = 1. The first baseline is a neutral two-body split.
    double black_hole_mass{1.0};
    double dimensionless_spin{0.999};
    double parent_rest_mass{1.0};
    double fragment_rest_mass{0.0};
    double incoming_specific_energy{1.0};
    double initial_radius_over_m{10.0};
    double escape_radius_over_m{20.0};
    double integration_step{0.002};
    std::size_t max_integration_steps{50'000};
    KerrIntegrationControl integration_control{};
    // Maximum dimensionless normalized conservation or initialization residual.
    double residual_tolerance{1.0e-8};
};

struct PenroseEventResult {
    PenroseEventStatus status{PenroseEventStatus::physics_invalid};
    PenroseSplitParameters split{};
    double horizon_radius{};
    double static_limit_radius{};
    double split_radius{};
    double input_energy{};
    double captured_energy{};
    double escaping_energy{};
    double eta_penrose{};
    double extracted_energy{};
    double four_momentum_residual{};
    double mass_shell_residual{};
    double energy_conservation_residual{};
    double angular_momentum_conservation_residual{};
    double geodesic_initialization_residual{};
    double maximum_normalized_residual{};
    // These are physical spacetime trajectories, unlike any future parameter search trace.
    Trajectory incoming_trajectory{};
    Trajectory captured_trajectory{};
    Trajectory escaping_trajectory{};
};

[[nodiscard]] constexpr double classical_penrose_efficiency_limit() {
    return 0.20710678118654752440;
}

[[nodiscard]] std::string_view penrose_event_status_name(PenroseEventStatus status);
[[nodiscard]] PenroseEventResult evaluate_equatorial_penrose_event(
    const EquatorialPenroseScenario& scenario, const PenroseSplitParameters& split);
}
