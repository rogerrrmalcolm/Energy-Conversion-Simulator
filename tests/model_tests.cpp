#include "bh/algebraic_model.hpp"
#include "bh/dijkstra.hpp"
#include "bh/constants.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"
#include "bh/schwarzschild_geodesic.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <vector>

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

    bh::PenroseEnergyBatch4Input energy_batch;
    energy_batch.input_energies = {1.0, 2.0, 4.0, 8.0};
    energy_batch.escaping_energies = {1.15, 1.5, 4.0, 10.0};
    const auto energy_batch_result = bh::penrose_energy_extraction_batch4(energy_batch);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        const double difference =
            energy_batch.escaping_energies[lane] - energy_batch.input_energies[lane];
        near_relative(energy_batch_result.eta_penrose[lane],
                      difference / energy_batch.input_energies[lane], 1.0e-14,
                      "four-lane Penrose efficiency matches scalar arithmetic");
        near_relative(energy_batch_result.extracted_energies[lane],
                      std::max(0.0, difference), 1.0e-14,
                      "four-lane Penrose extraction matches scalar arithmetic");
    }

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
    bh::KerrRadialPotentialBatch4 radial_batch;
    radial_batch.black_hole_masses = {1.0, 1.0, 2.0, 1.0};
    radial_batch.spin_lengths = {0.5, 0.8, 1.0, 0.2};
    radial_batch.energies = {1.0, 0.95, 2.0, 1.1};
    radial_batch.angular_momenta = {0.0, 2.0, -1.0, 0.5};
    radial_batch.rest_masses = {1.0, 1.0, 0.0, 1.0};
    radial_batch.carter_constants = {0.0, 0.0, 0.0, 0.0};
    radial_batch.radii = {10.0, 5.0, 8.0, 3.0};
    const auto radial_batch_result = bh::kerr_radial_potential_batch4(radial_batch);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        const bh::KerrOrbit lane_orbit{
            radial_batch.black_hole_masses[lane], radial_batch.spin_lengths[lane],
            radial_batch.energies[lane], radial_batch.angular_momenta[lane],
            radial_batch.rest_masses[lane], 1, radial_batch.carter_constants[lane]};
        near_relative(radial_batch_result[lane],
                      bh::kerr_radial_potential(lane_orbit, radial_batch.radii[lane]),
                      1.0e-14,
                      "four-lane Kerr radial potential matches the scalar kernel");
    }

    bh::KerrFourMomentumBatch4Input momentum_batch;
    momentum_batch.states.black_hole_masses = {1.0, 1.0, 1.0, 1.0};
    momentum_batch.states.spin_lengths = {0.5, 0.5, 0.5, 0.5};
    momentum_batch.states.energies = {1.0, 1.0, 1.0, 1.0};
    momentum_batch.states.angular_momenta = {0.0, 0.1, -0.1, 0.2};
    momentum_batch.states.rest_masses = {1.0, 1.0, 1.0, 1.0};
    momentum_batch.states.carter_constants = {0.0, 0.0, 0.0, 0.0};
    momentum_batch.states.radii = {10.0, 11.0, 12.0, 13.0};
    momentum_batch.radial_directions = {1, -1, 1, -1};
    const auto momentum_batch_result =
        bh::kerr_equatorial_four_momentum_batch4(momentum_batch);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        const bh::KerrOrbit lane_orbit{
            momentum_batch.states.black_hole_masses[lane],
            momentum_batch.states.spin_lengths[lane], momentum_batch.states.energies[lane],
            momentum_batch.states.angular_momenta[lane],
            momentum_batch.states.rest_masses[lane], momentum_batch.radial_directions[lane],
            momentum_batch.states.carter_constants[lane]};
        const auto scalar_momentum = bh::kerr_equatorial_four_momentum(
            lane_orbit, momentum_batch.states.radii[lane]);
        near_relative(momentum_batch_result.coordinate_time[lane],
                      scalar_momentum.coordinate_time, 1.0e-14,
                      "four-lane Kerr coordinate time matches the scalar kernel");
        near_relative(momentum_batch_result.radial[lane], scalar_momentum.radial,
                      1.0e-14,
                      "four-lane Kerr radial momentum matches the scalar kernel");
        near_relative(momentum_batch_result.azimuth[lane], scalar_momentum.azimuth,
                      1.0e-14,
                      "four-lane Kerr azimuth momentum matches the scalar kernel");
    }

    const bh::PenroseZamoGeometryBatch4 zamo_geometry{
        momentum_batch.states.black_hole_masses, momentum_batch.states.spin_lengths,
        momentum_batch.states.radii};
    const auto local_momentum =
        bh::coordinate_to_zamo_batch4(zamo_geometry, momentum_batch_result);
    const auto coordinate_round_trip =
        bh::zamo_to_coordinate_batch4(zamo_geometry, local_momentum);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        near_relative(coordinate_round_trip.coordinate_time[lane],
                      momentum_batch_result.coordinate_time[lane], 1.0e-13,
                      "four-lane ZAMO coordinate-time round trip is stable");
        near_relative(coordinate_round_trip.radial[lane], momentum_batch_result.radial[lane],
                      1.0e-13, "four-lane ZAMO radial round trip is stable");
        near_relative(coordinate_round_trip.azimuth[lane],
                      momentum_batch_result.azimuth[lane], 1.0e-13,
                      "four-lane ZAMO azimuth round trip is stable");
    }

    bh::PenroseFragmentSplitBatch4Input split_batch;
    split_batch.parent_unit_velocities.time = {1.0, 1.0, 1.0, 1.0};
    split_batch.radial_bases.radial = {1.0, 1.0, 1.0, 1.0};
    split_batch.azimuth_bases.azimuth = {1.0, 1.0, 1.0, 1.0};
    split_batch.split_angles_rad = {0.0, 0.5 * std::acos(-1.0), std::acos(-1.0),
                                    -0.5 * std::acos(-1.0)};
    split_batch.daughter_com_energies = {1.0, 1.0, 1.0, 1.0};
    split_batch.daughter_com_momenta = {0.6, 0.6, 0.6, 0.6};
    const auto split_batch_result = bh::split_penrose_fragments_batch4(split_batch);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        near_relative(split_batch_result.first.time[lane], 1.0, 1.0e-14,
                      "four-lane split gives the first daughter its center energy");
        near_relative(split_batch_result.first.radial[lane],
                      0.6 * std::cos(split_batch.split_angles_rad[lane]), 1.0e-14,
                      "four-lane split radial component matches scalar arithmetic");
        near_relative(split_batch_result.first.azimuth[lane],
                      0.6 * std::sin(split_batch.split_angles_rad[lane]), 1.0e-14,
                      "four-lane split azimuth component matches scalar arithmetic");
        near_relative(split_batch_result.second.radial[lane],
                      -split_batch_result.first.radial[lane], 1.0e-14,
                      "four-lane split produces opposite daughter radial momenta");
        near_relative(split_batch_result.second.azimuth[lane],
                      -split_batch_result.first.azimuth[lane], 1.0e-14,
                      "four-lane split produces opposite daughter azimuth momenta");
    }

    bh::PenroseConservationBatch4Input conservation_batch;
    conservation_batch.parent.time = {2.0, 2.0, 2.0, 2.0};
    conservation_batch.first = split_batch_result.first;
    conservation_batch.second = split_batch_result.second;
    conservation_batch.fragment_rest_masses = {0.8, 0.8, 0.8, 0.8};
    conservation_batch.incoming_constants.energies = {2.0, 2.0, 2.0, 2.0};
    conservation_batch.first_constants.energies = {1.0, 1.1, 1.2, 1.3};
    conservation_batch.second_constants.energies = {1.0, 0.9, 0.7, 0.7};
    conservation_batch.incoming_constants.angular_momenta = {0.0, 0.0, 0.0, 0.0};
    conservation_batch.first_constants.angular_momenta = {0.2, 0.3, 0.4, 0.5};
    conservation_batch.second_constants.angular_momenta = {-0.2, -0.4, -0.4, -0.3};
    const auto conservation =
        bh::penrose_conservation_residuals_batch4(conservation_batch);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        near(conservation.four_momentum_residuals[lane], 0.0, 1.0e-14,
             "four-lane split conserves local four-momentum");
        near(conservation.mass_shell_residuals[lane], 0.0, 1.0e-14,
             "four-lane split preserves the daughter mass shell");
        near_relative(
            conservation.energy_residuals[lane],
            std::abs(conservation_batch.first_constants.energies[lane] +
                     conservation_batch.second_constants.energies[lane] -
                     conservation_batch.incoming_constants.energies[lane]),
            1.0e-14, "four-lane energy residual matches scalar arithmetic");
        near_relative(
            conservation.angular_momentum_residuals[lane],
            std::abs(conservation_batch.first_constants.angular_momenta[lane] +
                     conservation_batch.second_constants.angular_momenta[lane] -
                     conservation_batch.incoming_constants.angular_momenta[lane]),
            1.0e-14, "four-lane angular-momentum residual matches scalar arithmetic");
    }
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

    const bh::KerrOrbitBatch4 integration_orbits{{outward, inward, outward, inward}};
    const bh::DoubleBatch4 integration_initial_radii{{10.0, 3.0, 10.0, 3.0}};
    const bh::DoubleBatch4 integration_escape_radii{{10.2, 20.0, 10.3, 20.0}};
    const bh::KerrLaneMaskBatch4 integration_lanes{{true, true, true, true}};
    const auto integration_batch = bh::integrate_kerr_batch4(
        integration_orbits, integration_initial_radii, 0.005, 50'000,
        integration_escape_radii, integration_lanes);
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        const auto scalar_trajectory = bh::integrate_kerr(
            integration_orbits[lane], integration_initial_radii[lane], 0.005,
            50'000, integration_escape_radii[lane]);
        check(integration_batch[lane].termination == scalar_trajectory.termination,
              "four-lane adaptive Kerr integration matches scalar termination");
        near_relative(integration_batch[lane].points.back().radius,
                      scalar_trajectory.points.back().radius, 1.0e-11,
                      "four-lane adaptive Kerr integration matches scalar final radius");
        near_relative(integration_batch[lane].points.back().affine_parameter,
                      scalar_trajectory.points.back().affine_parameter, 1.0e-10,
                      "four-lane adaptive Kerr integration matches scalar affine time");
        check(integration_batch[lane].diagnostics.maximum_normalized_error <= 1.0 &&
                  integration_batch[lane]
                          .diagnostics.maximum_normalized_radial_residual <=
                      1.0e-7,
              "four-lane adaptive Kerr integration preserves error and radial bounds");
    }

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

    const bh::EquatorialPenroseScenario reference_penrose_scenario{
        1.0, 0.999, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, {}, 1e-7};
    const auto penrose = bh::evaluate_equatorial_penrose_event(
        reference_penrose_scenario, {1.10, 2.07, -2.0});
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

    const std::array<bh::PenroseSplitParameters, bh::avx2_double_lanes> angle_batch_splits{{
        {1.10, 2.07, -2.00},
        {1.10, 2.07, -1.99},
        {1.10, 2.07, -1.98},
        {1.10, 2.07, -1.97}}};
    const auto angle_batch = bh::evaluate_equatorial_penrose_angle_batch4(
        reference_penrose_scenario, angle_batch_splits);
    check(angle_batch.used_avx2 == bh::penrose_batch4_uses_avx2(),
          "four-lane event evaluator reports whether AVX2 executed");
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        const auto scalar_event = bh::evaluate_equatorial_penrose_event(
            reference_penrose_scenario, angle_batch_splits[lane]);
        const auto& batch_event = angle_batch.events[lane];
        check(batch_event.status == scalar_event.status &&
                  batch_event.incoming_termination ==
                      scalar_event.incoming_trajectory.termination &&
                  batch_event.captured_termination ==
                      scalar_event.captured_trajectory.termination &&
                  batch_event.escaping_termination ==
                      scalar_event.escaping_trajectory.termination,
              "four-lane event status and trajectory outcomes match scalar evaluation");
        near_relative(batch_event.input_energy, scalar_event.input_energy, 1.0e-12,
                      "four-lane input energy matches scalar evaluation");
        near_relative(batch_event.captured_energy, scalar_event.captured_energy, 1.0e-12,
                      "four-lane captured energy matches scalar evaluation");
        near_relative(batch_event.escaping_energy, scalar_event.escaping_energy, 1.0e-12,
                      "four-lane escaping energy matches scalar evaluation");
        near_relative(batch_event.eta_penrose, scalar_event.eta_penrose, 1.0e-12,
                      "four-lane Penrose efficiency matches scalar evaluation");
        near_relative(batch_event.extracted_energy, scalar_event.extracted_energy, 1.0e-12,
                      "four-lane extracted energy matches scalar evaluation");
        near_relative(batch_event.maximum_normalized_residual,
                      scalar_event.maximum_normalized_residual, 1.0e-12,
                      "four-lane maximum residual matches scalar evaluation");
    }
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

    const auto loaded_search = bh::load_equatorial_penrose_dijkstra_input(
        std::filesystem::path(BH_SOURCE_DIR) / "tests" / "fixtures" /
        "equatorial_penrose_dijkstra_goal.cfg");
    near(loaded_search.search.start.split_radius_over_m, 1.09, 0.0,
         "search scenario loader preserves the declared start radius");
    near(loaded_search.search.step.split_radius_over_m, 0.01, 0.0,
         "search scenario loader preserves the declared radius step");
    near(loaded_search.search.eta_target, 0.04, 0.0,
         "internal goal fixture preserves its low test threshold");
    check(loaded_search.search.edge_costs[0] == 1 &&
              loaded_search.search.edge_costs[1] == 1 &&
              loaded_search.search.edge_costs[2] == 1 &&
              loaded_search.search.max_evaluations == 125 &&
              loaded_search.search.time_budget.count() == 0 &&
              loaded_search.search.algorithm == bh::PenroseSearchAlgorithm::dijkstra_h_zero,
          "search scenario loader preserves explicit unit-cost Dijkstra controls");
    check(bh::max_penrose_search_nodes == 25'000,
          "Penrose search exposes the declared 25000-node limit");
    const bh::PenroseSearchWindowSummary fixture_window =
        bh::describe_penrose_search_window(loaded_search.scenario, loaded_search.search);
    check(fixture_window.dimension_sizes == std::array<std::size_t, 3>{5, 5, 5} &&
              fixture_window.candidates == 125,
          "search-window description reports the exact discrete grid shape");
    bh::require_complete_penrose_search_window(loaded_search.scenario, loaded_search.search);

    const auto public_target_search = bh::load_equatorial_penrose_dijkstra_input(
        std::filesystem::path(BH_SOURCE_DIR) / "scenarios" /
        "equatorial_penrose_dijkstra_15_percent.cfg");
    near(public_target_search.search.eta_target, 0.15, 0.0,
         "public bounded search retains the declared 15 percent target");
    const bh::PenroseSearchWindowSummary public_target_window =
        bh::describe_penrose_search_window(
            public_target_search.scenario, public_target_search.search);
    check(public_target_window.dimension_sizes ==
              std::array<std::size_t, 3>{100, 10, 25} &&
              public_target_window.candidates == 25'000,
          "public 15 percent search contains exactly 25000 candidate nodes");
    bh::require_complete_penrose_search_window(
        public_target_search.scenario, public_target_search.search);

    auto exact_limit_search = public_target_search.search;
    exact_limit_search.max_evaluations = bh::max_penrose_search_nodes;
    const bh::PenroseSearchWindowSummary exact_limit_window =
        bh::describe_penrose_search_window(loaded_search.scenario, exact_limit_search);
    check(exact_limit_window.dimension_sizes ==
              std::array<std::size_t, 3>{100, 10, 25} &&
              exact_limit_window.candidates == bh::max_penrose_search_nodes,
          "search-window validation accepts exactly 25000 candidate nodes");
    bh::require_complete_penrose_search_window(loaded_search.scenario, exact_limit_search);

    rejected = false;
    try {
        auto oversized_search = exact_limit_search;
        oversized_search.upper_bound.split_angle_rad = -1.87;
        (void)bh::describe_penrose_search_window(loaded_search.scenario, oversized_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose search rejects a candidate window above 25000 nodes");

    rejected = false;
    try {
        auto incomplete_search = loaded_search.search;
        incomplete_search.max_evaluations = fixture_window.candidates - 1;
        bh::require_complete_penrose_search_window(
            loaded_search.scenario, incomplete_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "complete CLI window rejects a node budget that could skip a grid node");

    rejected = false;
    try {
        auto timed_search = loaded_search.search;
        timed_search.time_budget = std::chrono::milliseconds(1);
        bh::require_complete_penrose_search_window(loaded_search.scenario, timed_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "complete CLI window rejects a time budget that could skip a grid node");

    const auto penrose_search =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, loaded_search.search);
    check(penrose_search.found && penrose_search.target_reached &&
              penrose_search.status == bh::PenroseDijkstraSearchStatus::found_goal,
          "Penrose Dijkstra search finds a physics-valid target candidate");
    check(penrose_search.parameter_adjustment_path.size() == 2,
          "Penrose Dijkstra search returns the minimum one-step adjustment trace");
    check(penrose_search.parameter_adjustment_path.front().g_cost == 0 &&
              penrose_search.parameter_adjustment_path.front().h_cost == 0 &&
              penrose_search.parameter_adjustment_path.front().f_cost == 0 &&
              penrose_search.parameter_adjustment_path.back().g_cost == 1 &&
              penrose_search.parameter_adjustment_path.back().h_cost == 0 &&
              penrose_search.parameter_adjustment_path.back().f_cost == 1,
          "Penrose Dijkstra trace records the cumulative unit adjustment cost");
    check(penrose_search.parameter_adjustment_path.front().status !=
              bh::PenroseDijkstraNodeStatus::goal_feasible &&
              penrose_search.parameter_adjustment_path.back().status ==
                  bh::PenroseDijkstraNodeStatus::goal_feasible,
          "Penrose Dijkstra accepts only a candidate satisfying the physical goal predicate");
    check(!penrose_search.parameter_adjustment_path.front().parent_key &&
              penrose_search.parameter_adjustment_path.back().parent_key &&
              *penrose_search.parameter_adjustment_path.back().parent_key ==
                  penrose_search.parameter_adjustment_path.front().key &&
              penrose_search.parameter_adjustment_path.back().local_change ==
                  bh::DijkstraGridKey{1, 0, 0},
          "Penrose Dijkstra path retains parent and one-coordinate local-change metadata");
    near(penrose_search.parameter_adjustment_path.back().split.split_radius_over_m, 1.10, 1e-14,
         "Penrose Dijkstra derives the selected split from its integer grid key");
    check(penrose_search.selected_event.status == bh::PenroseEventStatus::physically_feasible &&
              penrose_search.selected_event.eta_penrose >= loaded_search.search.eta_target &&
              penrose_search.selected_event.captured_trajectory.termination ==
                  bh::TrajectoryTermination::crossed_horizon &&
              penrose_search.selected_event.escaping_trajectory.termination ==
                  bh::TrajectoryTermination::reached_escape_radius &&
              penrose_search.diagnostics.final_verification_evaluations == 1 &&
              penrose_search.diagnostics.nodes_evaluated == 7 &&
              penrose_search.diagnostics.duplicate_nodes_skipped > 0,
          "Penrose Dijkstra freshly verifies a captured-and-escaping physical goal event");
    const auto repeated_penrose_search =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, loaded_search.search);
    check(repeated_penrose_search.status == penrose_search.status &&
              repeated_penrose_search.diagnostics.nodes_generated ==
                  penrose_search.diagnostics.nodes_generated &&
              repeated_penrose_search.diagnostics.nodes_evaluated ==
                  penrose_search.diagnostics.nodes_evaluated &&
              repeated_penrose_search.diagnostics.nodes_expanded ==
                  penrose_search.diagnostics.nodes_expanded &&
              repeated_penrose_search.diagnostics.duplicate_nodes_skipped ==
                  penrose_search.diagnostics.duplicate_nodes_skipped &&
              repeated_penrose_search.parameter_adjustment_path.size() ==
                  penrose_search.parameter_adjustment_path.size() &&
              repeated_penrose_search.parameter_adjustment_path.back().key ==
                  penrose_search.parameter_adjustment_path.back().key,
          "Penrose Dijkstra repeats the same deterministic search outcome and trace");
    std::vector<bh::PenroseSearchProgress> dijkstra_progress;
    const auto observed_penrose_search = bh::find_penrose_dijkstra_path(
        loaded_search.scenario, loaded_search.search, {},
        {1, [&](const bh::PenroseSearchProgress& progress) {
             dijkstra_progress.push_back(progress);
         }});
    check(!dijkstra_progress.empty() &&
              dijkstra_progress.back().nodes_evaluated ==
                  observed_penrose_search.diagnostics.nodes_evaluated &&
              dijkstra_progress.back().best_eta_penrose &&
              observed_penrose_search.status == penrose_search.status &&
              observed_penrose_search.parameter_adjustment_path.back().key ==
                  penrose_search.parameter_adjustment_path.back().key,
          "Dijkstra progress reports aggregate state without changing traversal outcome");

    auto start_goal_search = loaded_search.search;
    start_goal_search.start = {1.10, 2.07, -2.0};
    start_goal_search.lower_bound = start_goal_search.start;
    start_goal_search.upper_bound = start_goal_search.start;
    start_goal_search.max_evaluations = 1;
    const auto start_goal =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, start_goal_search);
    check(start_goal.status == bh::PenroseDijkstraSearchStatus::found_goal &&
              start_goal.parameter_adjustment_path.size() == 1 &&
              start_goal.parameter_adjustment_path.front().g_cost == 0,
          "Penrose Dijkstra returns a one-node trace when the start already meets the goal");

    auto fallback_search = loaded_search.search;
    fallback_search.start = {1.09, 2.07, -2.0};
    fallback_search.lower_bound = fallback_search.start;
    fallback_search.upper_bound = {1.10, 2.07, -2.0};
    fallback_search.eta_target = 0.05;
    fallback_search.max_evaluations = 2;
    const auto fallback =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, fallback_search);
    check(fallback.found && !fallback.target_reached &&
              fallback.status ==
                  bh::PenroseDijkstraSearchStatus::best_feasible_below_target &&
              fallback.diagnostics.nodes_evaluated == 2 &&
              fallback.diagnostics.final_verification_evaluations == 1 &&
              fallback.parameter_adjustment_path.size() == 2 &&
              fallback.parameter_adjustment_path.back().g_cost == 1 &&
              fallback.parameter_adjustment_path.back().status ==
                  bh::PenroseDijkstraNodeStatus::escaping_without_target &&
              fallback.selected_event.status == bh::PenroseEventStatus::physically_feasible &&
              fallback.selected_event.eta_penrose < fallback_search.eta_target,
          "Penrose Dijkstra returns the greatest validated below-target extraction");
    near(fallback.parameter_adjustment_path.back().split.split_radius_over_m, 1.10, 1e-14,
         "Penrose fallback returns the best candidate in the completed bounded grid");

    auto fallback_tie_search = fallback_search;
    fallback_tie_search.start = {1.10, 2.07, -2.0};
    fallback_tie_search.lower_bound = fallback_tie_search.start;
    fallback_tie_search.upper_bound = fallback_tie_search.start;
    fallback_tie_search.step = {0.01, 0.01, 2.0 * std::acos(-1.0)};
    fallback_tie_search.upper_bound.split_angle_rad +=
        fallback_tie_search.step.split_angle_rad;
    const auto fallback_tie =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, fallback_tie_search);
    check(fallback_tie.status ==
                  bh::PenroseDijkstraSearchStatus::best_feasible_below_target &&
              fallback_tie.diagnostics.nodes_evaluated == 2 &&
              fallback_tie.parameter_adjustment_path.size() == 1 &&
              fallback_tie.parameter_adjustment_path.back().g_cost == 0,
          "Penrose fallback resolves equivalent extraction by minimum adjustment cost");

    auto no_solution_search = start_goal_search;
    no_solution_search.start = {1.99, 2.07, -2.0};
    no_solution_search.lower_bound = no_solution_search.start;
    no_solution_search.upper_bound = no_solution_search.start;
    no_solution_search.eta_target = 0.15;
    const auto no_solution =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, no_solution_search);
    check(!no_solution.found && !no_solution.target_reached &&
              no_solution.status ==
                  bh::PenroseDijkstraSearchStatus::no_solution_within_bounds,
          "Penrose Dijkstra retains no-solution status when no validated extraction exists");

    auto budget_limited_search = loaded_search.search;
    budget_limited_search.max_evaluations = 1;
    const auto budget_limited =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, budget_limited_search);
    check(!budget_limited.found &&
              budget_limited.status == bh::PenroseDijkstraSearchStatus::node_budget_exhausted &&
              budget_limited.diagnostics.nodes_evaluated == 1,
          "Penrose Dijkstra does not return a bounded fallback from a partial node budget");

    auto time_limited_search = loaded_search.search;
    time_limited_search.time_budget = std::chrono::milliseconds(1);
    const auto time_limited =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, time_limited_search);
    check(!time_limited.found &&
              time_limited.status == bh::PenroseDijkstraSearchStatus::time_budget_exhausted,
          "Penrose Dijkstra does not return a bounded fallback after a wall-clock timeout");

    auto unattainable_search = start_goal_search;
    unattainable_search.eta_target = bh::classical_penrose_efficiency_limit();
    const auto unattainable =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, unattainable_search);
    check(unattainable.status ==
                  bh::PenroseDijkstraSearchStatus::target_unattainable_under_model &&
              unattainable.diagnostics.nodes_evaluated == 0,
          "Penrose Dijkstra rejects the unattainable extremal single-split target before expansion");

    std::stop_source cancellation_source;
    cancellation_source.request_stop();
    const auto cancelled = bh::find_penrose_dijkstra_path(
        loaded_search.scenario, loaded_search.search, cancellation_source.get_token());
    check(cancelled.status == bh::PenroseDijkstraSearchStatus::cancelled &&
              cancelled.diagnostics.nodes_evaluated == 0,
          "Penrose Dijkstra supports cooperative cancellation before evaluation");

    auto coarse_radius_search = loaded_search.search;
    coarse_radius_search.lower_bound = {1.09, 2.07, -2.0};
    coarse_radius_search.upper_bound = {1.10, 2.07, -2.0};
    coarse_radius_search.step = {0.01, 0.01, 0.01};
    coarse_radius_search.start = {1.09, 2.07, -2.0};
    coarse_radius_search.max_evaluations = 2;
    const auto coarse_radius =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, coarse_radius_search);
    auto fine_radius_search = coarse_radius_search;
    fine_radius_search.step.split_radius_over_m = 0.005;
    fine_radius_search.max_evaluations = 3;
    const auto fine_radius =
        bh::find_penrose_dijkstra_path(loaded_search.scenario, fine_radius_search);
    check(coarse_radius.found && fine_radius.found &&
              coarse_radius.selected_event.eta_penrose >= coarse_radius_search.eta_target &&
              fine_radius.selected_event.eta_penrose >= fine_radius_search.eta_target,
          "coarse and refined radius grids both reproduce a validated target-reaching result");

    const auto small_phase_map =
        bh::evaluate_penrose_phase_map(loaded_search.scenario, coarse_radius_search);
    check(small_phase_map.complete &&
              small_phase_map.status == bh::PenrosePhaseMapStatus::completed &&
              small_phase_map.candidates.size() == 2 &&
              small_phase_map.diagnostics.nodes_evaluated == 2 &&
              small_phase_map.best_validated_candidate &&
              small_phase_map.best_validated_candidate->key == bh::DijkstraGridKey{1, 0, 0} &&
              small_phase_map.best_event.status == bh::PenroseEventStatus::physically_feasible &&
              small_phase_map.diagnostics.final_verification_evaluations == 1,
          "bounded phase map returns a freshly verified best extraction within its full grid");
    check(small_phase_map.best_validated_candidate &&
              fallback.parameter_adjustment_path.back().key ==
                  small_phase_map.best_validated_candidate->key &&
              fallback.selected_event.eta_penrose == small_phase_map.best_event.eta_penrose,
          "completed Dijkstra fallback agrees with the exhaustive bounded phase map");

    auto four_angle_phase_map_search = coarse_radius_search;
    four_angle_phase_map_search.lower_bound = {1.10, 2.07, -2.00};
    four_angle_phase_map_search.upper_bound = {1.10, 2.07, -1.97};
    four_angle_phase_map_search.step = {0.01, 0.01, 0.01};
    four_angle_phase_map_search.start = four_angle_phase_map_search.lower_bound;
    four_angle_phase_map_search.max_evaluations = 4;
    const auto four_angle_phase_map = bh::evaluate_penrose_phase_map(
        loaded_search.scenario, four_angle_phase_map_search);
    check(four_angle_phase_map.complete &&
              four_angle_phase_map.candidates.size() == 4 &&
              four_angle_phase_map.diagnostics.nodes_evaluated == 4 &&
              four_angle_phase_map.diagnostics.four_lane_batches == 1 &&
              four_angle_phase_map.diagnostics.four_lane_nodes == 4 &&
              four_angle_phase_map.diagnostics.scalar_nodes == 0 &&
              four_angle_phase_map.diagnostics.avx2_four_lane_batches ==
                  (bh::penrose_batch4_uses_avx2() ? 1U : 0U),
          "phase map dispatches four adjacent angle states through one batch");
    std::vector<bh::PenroseSearchProgress> phase_map_progress;
    const auto observed_four_angle_phase_map = bh::evaluate_penrose_phase_map(
        loaded_search.scenario, four_angle_phase_map_search, {},
        {1, [&](const bh::PenroseSearchProgress& progress) {
             phase_map_progress.push_back(progress);
         }});
    check(!phase_map_progress.empty() &&
              phase_map_progress.back().nodes_evaluated == 4 &&
              phase_map_progress.back().candidates_in_domain == 4 &&
              observed_four_angle_phase_map.best_validated_candidate &&
              four_angle_phase_map.best_validated_candidate &&
              observed_four_angle_phase_map.best_validated_candidate->key ==
                  four_angle_phase_map.best_validated_candidate->key,
          "phase-map progress reports aggregate state without changing selection");
    for (std::size_t lane = 0; lane < bh::avx2_double_lanes; ++lane) {
        check(four_angle_phase_map.candidates[lane].key ==
                  bh::DijkstraGridKey{0, 0, static_cast<int>(lane)},
              "four-lane phase map preserves deterministic angle-key order");
    }
    auto limited_phase_map_search = coarse_radius_search;
    limited_phase_map_search.max_evaluations = 1;
    const auto limited_phase_map =
        bh::evaluate_penrose_phase_map(loaded_search.scenario, limited_phase_map_search);
    check(!limited_phase_map.complete &&
              limited_phase_map.status == bh::PenrosePhaseMapStatus::node_budget_exhausted &&
              limited_phase_map.candidates.size() == 1,
          "bounded phase map never labels a partial scan as a full-grid maximum");

    rejected = false;
    try {
        auto invalid_search = loaded_search.search;
        invalid_search.step.split_angle_rad = 0.0;
        (void)bh::find_penrose_dijkstra_path(loaded_search.scenario, invalid_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose Dijkstra rejects a non-positive candidate-grid step");

    rejected = false;
    try {
        auto invalid_search = loaded_search.search;
        invalid_search.edge_costs[0] = 2;
        (void)bh::find_penrose_dijkstra_path(loaded_search.scenario, invalid_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose Dijkstra preserves the declared unit-cost baseline");

    rejected = false;
    try {
        auto invalid_search = loaded_search.search;
        invalid_search.start.split_radius_over_m = 1.095;
        (void)bh::find_penrose_dijkstra_path(loaded_search.scenario, invalid_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose Dijkstra rejects a start candidate off the canonical grid");

    rejected = false;
    try {
        auto invalid_search = loaded_search.search;
        invalid_search.lower_bound.split_radius_over_m = 1.04;
        (void)bh::find_penrose_dijkstra_path(loaded_search.scenario, invalid_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose Dijkstra rejects search radius bounds outside the ergosphere");

    rejected = false;
    try {
        auto invalid_search = loaded_search.search;
        invalid_search.step.split_angle_rad = std::numeric_limits<double>::quiet_NaN();
        (void)bh::find_penrose_dijkstra_path(loaded_search.scenario, invalid_search);
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Penrose Dijkstra rejects non-finite grid parameters");

    auto invalid_evaluator_scenario = loaded_search.scenario;
    invalid_evaluator_scenario.parent_rest_mass = 0.0;
    const auto evaluator_failure =
        bh::find_penrose_dijkstra_path(invalid_evaluator_scenario, start_goal_search);
    check(evaluator_failure.status == bh::PenroseDijkstraSearchStatus::evaluation_failure &&
              !evaluator_failure.failure_message.empty(),
          "Penrose Dijkstra preserves an evaluator exception as an explicit search failure");

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
