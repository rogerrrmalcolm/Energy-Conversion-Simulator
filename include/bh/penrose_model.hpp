#pragma once

#include "bh/kerr_geodesic.hpp"

#include <array>
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
    // Includes split, initialization, and trajectory radial-first-integral residuals.
    double maximum_normalized_residual{};
    // These are physical spacetime trajectories, unlike any future parameter search trace.
    Trajectory incoming_trajectory{};
    Trajectory captured_trajectory{};
    Trajectory escaping_trajectory{};
};

// Compact result used while scanning candidate parameter grids. It preserves
// classification data and trajectory outcomes without retaining trajectory points.
struct PenroseEventSummary {
    PenroseEventStatus status{PenroseEventStatus::physics_invalid};
    PenroseSplitParameters split{};
    double input_energy{};
    double captured_energy{};
    double escaping_energy{};
    double eta_penrose{};
    double extracted_energy{};
    double maximum_normalized_residual{};
    TrajectoryTermination incoming_termination{TrajectoryTermination::completed};
    TrajectoryTermination captured_termination{TrajectoryTermination::completed};
    TrajectoryTermination escaping_termination{TrajectoryTermination::completed};
};

struct PenroseEventBatch4Result {
    std::array<PenroseEventSummary, avx2_double_lanes> events{};
    // True only when this batch reached AVX2 arithmetic; early parent exits remain false.
    bool used_avx2{};
};

struct PenroseEnergyBatch4Input {
    DoubleBatch4 input_energies{};
    DoubleBatch4 escaping_energies{};
};

struct PenroseEnergyBatch4Result {
    DoubleBatch4 eta_penrose{};
    DoubleBatch4 extracted_energies{};
};

struct PenroseLocalMomentumBatch4 {
    DoubleBatch4 time{};
    DoubleBatch4 radial{};
    DoubleBatch4 azimuth{};
};

struct PenroseZamoGeometryBatch4 {
    DoubleBatch4 black_hole_masses{};
    DoubleBatch4 spin_lengths{};
    DoubleBatch4 radii{};
};

struct PenroseFragmentSplitBatch4Input {
    PenroseLocalMomentumBatch4 parent_unit_velocities{};
    PenroseLocalMomentumBatch4 radial_bases{};
    PenroseLocalMomentumBatch4 azimuth_bases{};
    DoubleBatch4 split_angles_rad{};
    DoubleBatch4 daughter_com_energies{};
    DoubleBatch4 daughter_com_momenta{};
};

struct PenroseFragmentSplitBatch4Result {
    PenroseLocalMomentumBatch4 first{};
    PenroseLocalMomentumBatch4 second{};
};

struct PenroseConservedConstantsBatch4 {
    DoubleBatch4 energies{};
    DoubleBatch4 angular_momenta{};
};

struct PenroseConservationBatch4Input {
    PenroseLocalMomentumBatch4 parent{};
    PenroseLocalMomentumBatch4 first{};
    PenroseLocalMomentumBatch4 second{};
    DoubleBatch4 fragment_rest_masses{};
    PenroseConservedConstantsBatch4 incoming_constants{};
    PenroseConservedConstantsBatch4 first_constants{};
    PenroseConservedConstantsBatch4 second_constants{};
};

struct PenroseConservationBatch4Result {
    DoubleBatch4 four_momentum_residuals{};
    DoubleBatch4 mass_shell_residuals{};
    DoubleBatch4 energy_residuals{};
    DoubleBatch4 angular_momentum_residuals{};
};

[[nodiscard]] constexpr double classical_penrose_efficiency_limit() {
    return 0.20710678118654752440;
}

[[nodiscard]] std::string_view penrose_event_status_name(PenroseEventStatus status);
// Computes only the four-lane energy ledger; geodesics still establish event validity.
[[nodiscard]] PenroseEnergyBatch4Result penrose_energy_extraction_batch4(
    const PenroseEnergyBatch4Input& input);
[[nodiscard]] PenroseLocalMomentumBatch4 coordinate_to_zamo_batch4(
    const PenroseZamoGeometryBatch4& geometry,
    const KerrFourMomentumBatch4& coordinate_momenta);
[[nodiscard]] KerrFourMomentumBatch4 zamo_to_coordinate_batch4(
    const PenroseZamoGeometryBatch4& geometry,
    const PenroseLocalMomentumBatch4& local_momenta);
[[nodiscard]] PenroseFragmentSplitBatch4Result split_penrose_fragments_batch4(
    const PenroseFragmentSplitBatch4Input& input);
[[nodiscard]] PenroseConservationBatch4Result penrose_conservation_residuals_batch4(
    const PenroseConservationBatch4Input& input);
[[nodiscard]] bool penrose_batch4_uses_avx2() noexcept;
// Evaluates four split angles that share one split radius and incoming Lz. The
// shared incoming geodesic is prepared once; selected candidates are still
// re-evaluated through the full scalar API below before being returned.
[[nodiscard]] PenroseEventBatch4Result evaluate_equatorial_penrose_angle_batch4(
    const EquatorialPenroseScenario& scenario,
    const std::array<PenroseSplitParameters, avx2_double_lanes>& splits);
[[nodiscard]] PenroseEventResult evaluate_equatorial_penrose_event(
    const EquatorialPenroseScenario& scenario, const PenroseSplitParameters& split);
}
