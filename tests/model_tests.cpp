#include "bh/algebraic_model.hpp"
#include "bh/dijkstra.hpp"
#include "bh/constants.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"
#include "bh/schwarzschild_geodesic.hpp"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>

namespace {
int failures = 0;

void check(bool condition, const char* message) {
    if (!condition) { std::cerr << "FAIL: " << message << '\n'; ++failures; }
}

void near(double actual, double expected, double tolerance, const char* message) {
    check(std::abs(actual-expected) <= tolerance, message);
}

void near_relative(double actual, double expected, double relative_tolerance,
                   const char* message) {
    const double scale = std::max({1.0, std::abs(actual), std::abs(expected)});
    check(std::abs(actual - expected) <= relative_tolerance * scale, message);
}
}

int main() {
    const bh::DijkstraGridKey grid_lower{0, 0, 0};
    const bh::DijkstraGridKey grid_upper{1, 1, 1};
    const bh::DijkstraGridKey grid_start{0, 0, 0};
    const bh::DijkstraGridKey grid_goal{1, 1, 1};
    const auto grid_path =
        bh::find_dijkstra_grid_path(grid_start, grid_goal, grid_lower, grid_upper);
    check(grid_path.found, "Dijkstra finds an in-bounds grid goal");
    check(grid_path.parameter_adjustment_path.size() == 4,
          "Dijkstra returns the three unit changes of a shortest path");
    check(grid_path.parameter_adjustment_path.front() == grid_start &&
              grid_path.parameter_adjustment_path.back() == grid_goal,
          "Dijkstra path starts and ends at the declared grid nodes");
    const auto repeated_grid_path =
        bh::find_dijkstra_grid_path(grid_start, grid_goal, grid_lower, grid_upper);
    check(repeated_grid_path.parameter_adjustment_path == grid_path.parameter_adjustment_path,
          "Dijkstra has deterministic tie breaking");
    const auto already_at_goal =
        bh::find_dijkstra_grid_path(grid_start, grid_start, grid_lower, grid_upper);
    check(already_at_goal.found && already_at_goal.parameter_adjustment_path.size() == 1,
          "Dijkstra returns a one-node path when start already meets the goal");
    check(!bh::find_dijkstra_grid_path(grid_start, {2, 0, 0}, grid_lower, grid_upper).found,
          "Dijkstra rejects an out-of-bounds goal");

    const auto nonrotating = bh::rotational_energy(bh::solar_mass_kg, 0.0);
    near(nonrotating.rotational_energy_joules, 0.0, 1.0, "Schwarzschild has no rotational reservoir");
    near(nonrotating.irreducible_mass_fraction, 1.0, 0.0, "Schwarzschild irreducible mass equals total mass");
    near(bh::rotational_energy_fraction(0.0), 0.0, 0.0, "zero-spin reservoir fraction is zero");

    const auto near_extremal = bh::rotational_energy(1.0, 0.999999999);
    near(near_extremal.rotational_fraction, 1.0-1.0/std::sqrt(2.0), 3.0e-5,
         "near-extremal reservoir approaches 29.29 percent");
    const auto closest_subextremal = bh::rotational_energy(1.0, std::nextafter(1.0, 0.0));
    check(std::isfinite(closest_subextremal.rotational_fraction) &&
              std::isfinite(closest_subextremal.d_rotational_energy_d_spin_joules),
          "closest representable sub-extremal spin has finite algebraic outputs");
    near(bh::classical_penrose_efficiency_limit(), (std::sqrt(2.0) - 1.0) / 2.0, 1e-15,
         "classical ideal Penrose efficiency limit is approximately 20.71 percent");

    const auto uncertain = bh::rotational_energy(
        {bh::solar_mass_kg, 0.9, {0.8, 0.9, 0.99}});
    check(uncertain.rotational_energy_lower_joules < uncertain.rotational_energy_joules,
          "lower spin bound produces lower rotational reservoir");
    check(uncertain.rotational_energy_joules < uncertain.rotational_energy_upper_joules,
          "upper spin bound produces upper rotational reservoir");
    check(uncertain.d_rotational_energy_d_spin_joules > 0.0,
          "rotational reservoir sensitivity to spin is positive");
    const double finite_difference_step = 1.0e-6;
    const double finite_difference = (
        bh::rotational_energy(bh::solar_mass_kg, 0.9 + finite_difference_step)
            .rotational_energy_joules -
        bh::rotational_energy(bh::solar_mass_kg, 0.9 - finite_difference_step)
            .rotational_energy_joules) /
        (2.0 * finite_difference_step);
    near_relative(uncertain.d_rotational_energy_d_spin_joules, finite_difference, 1.0e-8,
                  "analytic rotational-energy sensitivity matches a finite difference");

    bool roundoff_rejected = false;
    try {
        (void)bh::rotational_energy(
            {1.0, 0.9, {0.8, std::nextafter(0.9, 1.0), 0.99}});
    } catch (const std::invalid_argument&) { roundoff_rejected = true; }
    check(!roundoff_rejected,
          "spin uncertainty accepts a central value within floating-point roundoff");

    const auto ranged = bh::rotational_energy_range(
        {{8.0 * bh::solar_mass_kg, 10.0 * bh::solar_mass_kg, 12.0 * bh::solar_mass_kg},
         {0.8, 0.9, 0.99}});
    check(ranged.lower.rotational_energy_joules < ranged.central.rotational_energy_joules,
          "lower mass and spin give lower rotational reservoir");
    check(ranged.central.rotational_energy_joules < ranged.upper.rotational_energy_joules,
          "upper mass and spin give upper rotational reservoir");
    near(ranged.rotational_energy_uncertainty_minus_joules,
         ranged.central.rotational_energy_joules - ranged.lower.rotational_energy_joules,
         0.0,
         "range reports lower asymmetric uncertainty");
    near(ranged.rotational_energy_uncertainty_plus_joules,
         ranged.upper.rotational_energy_joules - ranged.central.rotational_energy_joules,
         0.0,
         "range reports upper asymmetric uncertainty");
    const auto doubled_mass = bh::rotational_energy(2.0 * bh::solar_mass_kg, 0.9);
    const auto single_mass = bh::rotational_energy(bh::solar_mass_kg, 0.9);
    near_relative(doubled_mass.rotational_energy_joules,
                  2.0 * single_mass.rotational_energy_joules, 1.0e-14,
                  "rotational reservoir is linear in mass at fixed spin");

    bool rejected = false;
    try { (void)bh::rotational_energy(1.0, 1.0); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "extremal spin is rejected by sub-extremal model");

    rejected = false;
    try {
        (void)bh::rotational_energy({1.0, 0.7, {0.6, 0.8, 0.9}});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "spin uncertainty central value must match the selected spin");

    rejected = false;
    try {
        (void)bh::rotational_energy_range({{2.0, 1.0, 3.0}, {0.2, 0.3, 0.4}});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "mass range rejects unordered bounds");

    const double r = 10.0;
    const double circular_l = std::sqrt(r*r/(r-3.0));
    const double circular_energy = (r - 2.0) / std::sqrt(r * (r - 3.0));
    const bh::SchwarzschildOrbit circular_orbit{1.0, circular_energy, circular_l};
    near(bh::schwarzschild_radial_potential(circular_orbit, r), 0.0, 1e-14,
         "Schwarzschild circular orbit satisfies the radial mass-shell relation");
    const auto circular = bh::integrate_schwarzschild(
        circular_orbit, {0.0, r, 0.0, 0.0, 0.0},
        1e-3, 1000, 20.0);
    near(circular.points.back().radius, r, 1e-10, "Schwarzschild circular orbit remains circular");
    check(circular.diagnostics.maximum_normalized_radial_residual <= 1.0e-12,
          "Schwarzschild circular orbit preserves the radial first integral");

    const bh::SchwarzschildOrbit radial_schwarzschild{1.0, 1.0, 0.0};
    const double radial_start = 10.0;
    const auto schwarzschild_infall = bh::integrate_schwarzschild(
        radial_schwarzschild,
        {0.0, radial_start, 0.0, 0.0, -std::sqrt(2.0 / radial_start)},
        1e-3, 30'000, 20.0);
    const double schwarzschild_horizon_event = 2.0 * (1.0 + 1.0e-6);
    const double expected_schwarzschild_infall =
        2.0 * (std::pow(radial_start, 1.5) -
               std::pow(schwarzschild_horizon_event, 1.5)) /
        (3.0 * std::sqrt(2.0));
    check(schwarzschild_infall.termination == bh::TrajectoryTermination::crossed_horizon,
          "Schwarzschild radial free fall reaches the horizon event");
    near_relative(schwarzschild_infall.points.back().affine_parameter,
                  expected_schwarzschild_infall, 1.0e-6,
                  "Schwarzschild RK4 free fall matches the analytic proper time");
    check(schwarzschild_infall.diagnostics.maximum_normalized_radial_residual <= 1.0e-8,
          "Schwarzschild infall keeps a small radial first-integral residual");

    rejected = false;
    try {
        (void)bh::integrate_schwarzschild(
            radial_schwarzschild, {0.0, radial_start, 0.0, 0.0, 0.0},
            1e-3, 1000, 20.0);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Schwarzschild integration rejects an inconsistent radial velocity");

    near(bh::kerr_outer_horizon(1.0, 0.0), 2.0, 1e-14,
         "Kerr horizon reduces to Schwarzschild horizon");
    near(bh::kerr_spin_length(2.0, 0.9), 1.8, 1e-14,
         "dimensionless spin converts to Kerr spin length");
    near(bh::kerr_static_limit_radius(1.0, 0.5, 1.57079632679489661923), 2.0, 1e-14,
         "equatorial Kerr static limit is 2M");
    near(bh::kerr_static_limit_radius(1.0, 0.5, 0.0), bh::kerr_outer_horizon(1.0, 0.5),
         1e-14, "Kerr static limit meets horizon at the pole");
    check(bh::kerr_is_within_equatorial_ergosphere(1.0, 0.5, 1.9),
          "equatorial point between horizon and static limit is in ergosphere");
    const bh::KerrOrbit outward{1.0, 0.5, 1.0, 0.0, 1.0, 1};
    check(bh::kerr_radial_potential(outward, 10.0) >= 0.0, "Kerr test orbit is admissible");
    const auto escaping = bh::integrate_kerr(outward, 10.0, 0.01, 10'000, 11.0);
    check(escaping.termination == bh::TrajectoryTermination::reached_escape_radius,
          "outward Kerr trajectory reaches escape radius");
    check(escaping.diagnostics.maximum_normalized_radial_residual <= 1.0e-12,
          "outward Kerr trajectory preserves its radial first integral");
    const auto target = bh::integrate_kerr_to_radius(outward, 10.0, 10.5, 0.01, 10'000);
    check(target.termination == bh::TrajectoryTermination::reached_target_radius,
          "outward Kerr trajectory reaches requested target radius");
    const bh::KerrOrbit inward{1.0, 0.5, 1.0, 0.0, 1.0, -1};
    const auto captured = bh::integrate_kerr(inward, 3.0, 0.001, 50'000, 20.0);
    check(captured.termination == bh::TrajectoryTermination::crossed_horizon,
          "inward Kerr trajectory crosses horizon");

    const double bound_specific_energy = 0.95;
    const double expected_outer_turning_radius = 2.0 / (1.0 -
                                                        bound_specific_energy *
                                                            bound_specific_energy);
    const bh::KerrOrbit bound_outward{1.0, 0.0, bound_specific_energy, 0.0, 1.0, 1};
    const auto turning = bh::integrate_kerr(
        bound_outward, 10.0, 0.05, 100'000, 30.0, {1.0e-9, 1.0e-9, 1.0e-6});
    check(turning.termination == bh::TrajectoryTermination::turning_point,
          "outward Kerr trajectory stops at a radial turning point before escape");
    near_relative(turning.points.back().radius, expected_outer_turning_radius, 1.0e-7,
                  "Kerr turning-point event localizes the analytic Schwarzschild root");
    near(turning.points.back().radial_velocity, 0.0, 0.0,
         "Kerr turning-point event has zero radial velocity");
    check(turning.diagnostics.rejected_steps > 0,
          "Kerr integrator refines a step before reporting a turning point");

    const bh::KerrOrbit azimuth_check{1.0, 0.0, 1.0, 2.0, 1.0, 1};
    const auto momentum = bh::kerr_equatorial_four_momentum(azimuth_check, 10.0);
    near(momentum.azimuth, 0.02, 1e-14,
         "Schwarzschild limit preserves conventional positive azimuthal momentum");

    rejected = false;
    try {
        (void)bh::kerr_equatorial_four_momentum({1.0, 0.0, -1.0, 0.0, 1.0, -1}, 10.0);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Kerr momentum rejects a past-directed Schwarzschild state");

    const bh::KerrOrbit radial_infall{1.0, 0.0, 1.0, 0.0, 1.0, -1};
    const bh::KerrIntegrationControl strict_control{1.0e-10, 1.0e-10, 0.0};
    const auto analytic_infall = bh::integrate_kerr_to_radius(
        radial_infall, 10.0, 8.0, 0.5, 100'000, strict_control);
    const double expected_infall_affine_parameter =
        2.0 * (std::pow(10.0, 1.5) - std::pow(8.0, 1.5)) / (3.0 * std::sqrt(2.0));
    check(analytic_infall.termination == bh::TrajectoryTermination::reached_target_radius,
          "adaptive radial Kerr infall reaches its target radius");
    near_relative(analytic_infall.points.back().affine_parameter,
                  expected_infall_affine_parameter, 1.0e-7,
                  "adaptive Kerr infall matches the analytic Schwarzschild proper time");
    check(analytic_infall.diagnostics.accepted_steps > 0 &&
              analytic_infall.diagnostics.maximum_normalized_error <= 1.0,
          "adaptive Kerr diagnostics report an accepted bounded-error solution");
    check(analytic_infall.diagnostics.maximum_normalized_radial_residual <= 1.0e-12,
          "adaptive Kerr infall preserves the radial first integral");
    const auto adaptive_probe = bh::integrate_kerr_to_radius(
        radial_infall, 10.0, 8.0, 1.0, 100'000, {1.0e-12, 1.0e-12, 0.0});
    check(adaptive_probe.termination == bh::TrajectoryTermination::reached_target_radius &&
              adaptive_probe.diagnostics.rejected_steps > 0 &&
              adaptive_probe.diagnostics.maximum_normalized_error <= 1.0,
          "adaptive Kerr integration rejects and refines an overly coarse trial step");

    rejected = false;
    try {
        (void)bh::integrate_kerr_to_radius(
            radial_infall, 10.0, 8.0, 0.5, 100, {0.0, 1.0e-9, 0.0});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Kerr integration rejects invalid error-control settings");

    const auto penrose = bh::evaluate_equatorial_penrose_event(
        {1.0, 0.999, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, {}, 1e-7},
        {1.10, 2.07, -2.0});
    if (penrose.status != bh::PenroseEventStatus::physically_feasible) {
        std::cerr << "Penrose evaluator status: "
                  << bh::penrose_event_status_name(penrose.status) << '\n';
    }
    check(penrose.status == bh::PenroseEventStatus::physically_feasible,
          "reference equatorial split is physically feasible");
    check(penrose.four_momentum_residual <= 1e-7,
          "reference split conserves local four-momentum");
    check(penrose.mass_shell_residual <= 1e-7,
          "reference split satisfies fragment mass shells");
    check(penrose.geodesic_initialization_residual <= 1e-7,
          "reference fragments reconstruct from their Kerr conserved quantities");
    check(penrose.maximum_normalized_residual <= 1e-7,
          "reference event satisfies its normalized residual tolerance");
    check(penrose.captured_trajectory.termination == bh::TrajectoryTermination::crossed_horizon,
          "negative-energy reference fragment crosses horizon");
    check(penrose.escaping_trajectory.termination ==
              bh::TrajectoryTermination::reached_escape_radius,
          "positive-energy reference fragment reaches configured escape radius");
    check(penrose.extracted_energy > 0.0 && penrose.eta_penrose > 0.0,
          "reference event reports positive net extracted energy");
    near_relative(penrose.extracted_energy,
                  penrose.escaping_energy - penrose.input_energy, 1.0e-14,
                  "Penrose extracted energy is the escaping energy minus input energy");
    near_relative(penrose.eta_penrose,
                  penrose.extracted_energy / penrose.input_energy, 1.0e-14,
                  "Penrose efficiency is the normalized net extracted energy");
    check(std::max({penrose.incoming_trajectory.diagnostics.maximum_normalized_radial_residual,
                    penrose.captured_trajectory.diagnostics.maximum_normalized_radial_residual,
                    penrose.escaping_trajectory.diagnostics.maximum_normalized_radial_residual}) <=
              1.0e-7,
          "Penrose trajectories preserve their radial first integrals");
    const auto scaled_penrose = bh::evaluate_equatorial_penrose_event(
        {2.0, 0.999, 2.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, {}, 1e-7},
        {1.10, 2.07, -2.0});
    check(scaled_penrose.status == bh::PenroseEventStatus::physically_feasible,
          "scaled normalized event remains physically feasible");
    near_relative(scaled_penrose.eta_penrose, penrose.eta_penrose, 1.0e-8,
                  "dimensionless Penrose efficiency is invariant under common scaling");
    check(scaled_penrose.maximum_normalized_residual <= 1e-7,
          "scaled event satisfies the normalized residual tolerance");
    const auto outside_ergosphere = bh::evaluate_equatorial_penrose_event(
        {1.0, 0.999, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, {}, 1e-7},
        {2.1, 2.07, -2.0});
    check(outside_ergosphere.status == bh::PenroseEventStatus::outside_ergosphere,
          "split outside the equatorial ergosphere is rejected");
    const auto nonrotating_penrose = bh::evaluate_equatorial_penrose_event(
        {1.0, 0.0, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, {}, 1e-7},
        {2.1, 2.07, -2.0});
    check(nonrotating_penrose.status == bh::PenroseEventStatus::outside_ergosphere,
          "nonrotating Kerr has no equatorial Penrose split domain");

    const auto loaded_event = bh::load_equatorial_penrose_event_input(
        std::filesystem::path(BH_SOURCE_DIR) / "scenarios" / "equatorial_penrose_reference.cfg");
    near(loaded_event.scenario.dimensionless_spin, 0.999, 0.0,
         "scenario loader preserves the declared spin");
    near(loaded_event.split.split_radius_over_m, 1.10, 0.0,
         "scenario loader preserves the declared split radius");
    near(loaded_event.scenario.integration_control.relative_tolerance, 1.0e-9, 0.0,
         "scenario loader preserves the declared integration tolerance");
    const auto loaded_penrose =
        bh::evaluate_equatorial_penrose_event(loaded_event.scenario, loaded_event.split);
    check(loaded_penrose.status == bh::PenroseEventStatus::physically_feasible,
          "versioned reference scenario evaluates as physically feasible");

    rejected = false;
    try {
        (void)bh::load_equatorial_penrose_event_input(
            std::filesystem::path(BH_SOURCE_DIR) / "tests" / "data" /
            "invalid_penrose_unknown_key.cfg");
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "scenario loader rejects unknown keys");

    rejected = false;
    try {
        (void)bh::load_equatorial_penrose_event_input(
            std::filesystem::path(BH_SOURCE_DIR) / "tests" / "data" /
            "invalid_penrose_duplicate_key.cfg");
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "scenario loader rejects duplicate keys");

    const auto weak = bh::estimate_toy_plasma_transport({0.0, 1.0, 10.0, 0.9, 2.0});
    near(weak.outward_electromagnetic_power_watts, 0.0, 0.0,
         "zero field produces zero toy outward electromagnetic power");
    near(weak.outward_electromagnetic_energy_joules, 0.0, 0.0,
         "zero field produces zero toy outward electromagnetic energy");
    const auto plasma = bh::estimate_toy_plasma_transport({1.0, 1e-8, 10.0, 0.9, 2.0});
    check(plasma.alfven_speed_m_s > 0.0 && plasma.alfven_speed_m_s < bh::speed_of_light_m_s,
          "relativistic Alfven speed is causal");
    check(plasma.spin_coupling_efficiency >= 0.0 && plasma.spin_coupling_efficiency <= 1.0,
          "toy spin coupling is bounded");
    near_relative(plasma.outward_electromagnetic_power_watts,
                  plasma.poynting_power_watts * plasma.spin_coupling_efficiency, 1.0e-12,
                  "toy outward electromagnetic power applies the visible coupling");
    near_relative(plasma.outward_electromagnetic_energy_joules,
                  plasma.outward_electromagnetic_power_watts * 2.0, 1.0e-12,
                  "toy outward electromagnetic energy scales with duration");
    const auto zero_spin_plasma = bh::estimate_toy_plasma_transport({1.0, 1e-8, 10.0, 0.0, 2.0});
    check(zero_spin_plasma.poynting_power_watts > 0.0 &&
              zero_spin_plasma.outward_electromagnetic_power_watts == 0.0,
          "toy spin factor removes outward coupled power for a nonrotating source");
    const auto compatibility_plasma = bh::estimate_plasma_extraction({1.0, 1e-8, 10.0, 0.9, 2.0});
    near_relative(compatibility_plasma.outward_electromagnetic_energy_joules,
                  plasma.outward_electromagnetic_energy_joules, 0.0,
                  "legacy plasma API preserves the toy transport calculation");

    rejected = false;
    try {
        (void)bh::estimate_toy_plasma_transport(
            {std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0, 0.5, 1.0});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "plasma model rejects non-finite public input");

    bool overflow_rejected = false;
    try {
        (void)bh::rotational_energy(std::numeric_limits<double>::max(), 0.5);
    } catch (const std::overflow_error&) { overflow_rejected = true; }
    check(overflow_rejected, "algebraic model rejects non-finite overflow results");

    overflow_rejected = false;
    try {
        (void)bh::estimate_toy_plasma_transport(
            {std::numeric_limits<double>::max(), 1.0, 1.0, 0.5, 1.0});
    } catch (const std::overflow_error&) { overflow_rejected = true; }
    check(overflow_rejected, "plasma model rejects non-finite overflow results");

    overflow_rejected = false;
    try {
        (void)bh::kerr_outer_horizon(
            std::numeric_limits<double>::max(),
            std::numeric_limits<double>::max() / 2.0);
    } catch (const std::overflow_error&) { overflow_rejected = true; }
    check(overflow_rejected, "Kerr geometry rejects non-finite overflow results");

    overflow_rejected = false;
    try { (void)bh::kerr_equatorial_sigma(std::numeric_limits<double>::max()); }
    catch (const std::overflow_error&) { overflow_rejected = true; }
    check(overflow_rejected, "Kerr sigma rejects non-finite overflow results");

    overflow_rejected = false;
    try {
        (void)bh::schwarzschild_radial_potential(
            {std::numeric_limits<double>::max(), 1.0, 0.0},
            std::numeric_limits<double>::max());
    } catch (const std::overflow_error&) { overflow_rejected = true; }
    check(overflow_rejected, "Schwarzschild geometry rejects non-finite horizon results");

    rejected = false;
    try {
        (void)bh::kerr_outer_horizon(1.0, std::numeric_limits<double>::quiet_NaN());
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Kerr model rejects non-finite public input");

    if (failures == 0) std::cout << "All model tests passed\n";
    return failures == 0 ? 0 : 1;
}
