#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"
#include "bh/schwarzschild_geodesic.hpp"

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
}

int main() {
    const auto nonrotating = bh::rotational_energy(bh::solar_mass_kg, 0.0);
    near(nonrotating.rotational_energy_joules, 0.0, 1.0, "Schwarzschild has no rotational reservoir");
    near(nonrotating.irreducible_mass_fraction, 1.0, 0.0, "Schwarzschild irreducible mass equals total mass");
    near(bh::rotational_energy_fraction(0.0), 0.0, 0.0, "zero-spin reservoir fraction is zero");

    const auto near_extremal = bh::rotational_energy(1.0, 0.999999999);
    near(near_extremal.rotational_fraction, 1.0-1.0/std::sqrt(2.0), 3.0e-5,
         "near-extremal reservoir approaches 29.29 percent");
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
    const auto circular = bh::integrate_schwarzschild(
        {1.0, 0.956182887, circular_l}, {0.0, r, 0.0, 0.0, 0.0},
        1e-3, 1000, 20.0);
    near(circular.points.back().radius, r, 1e-10, "Schwarzschild circular orbit remains circular");

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
    const auto target = bh::integrate_kerr_to_radius(outward, 10.0, 10.5, 0.01, 10'000);
    check(target.termination == bh::TrajectoryTermination::reached_target_radius,
          "outward Kerr trajectory reaches requested target radius");
    const bh::KerrOrbit inward{1.0, 0.5, 1.0, 0.0, 1.0, -1};
    const auto captured = bh::integrate_kerr(inward, 3.0, 0.001, 50'000, 20.0);
    check(captured.termination == bh::TrajectoryTermination::crossed_horizon,
          "inward Kerr trajectory crosses horizon");
    const bh::KerrOrbit azimuth_check{1.0, 0.0, 1.0, 2.0, 1.0, 1};
    const auto momentum = bh::kerr_equatorial_four_momentum(azimuth_check, 10.0);
    near(momentum.azimuth, 0.02, 1e-14,
         "Schwarzschild limit preserves conventional positive azimuthal momentum");

    const auto penrose = bh::evaluate_equatorial_penrose_event(
        {1.0, 0.999, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, 1e-7},
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
    check(penrose.captured_trajectory.termination == bh::TrajectoryTermination::crossed_horizon,
          "negative-energy reference fragment crosses horizon");
    check(penrose.escaping_trajectory.termination ==
              bh::TrajectoryTermination::reached_escape_radius,
          "positive-energy reference fragment reaches configured escape radius");
    check(penrose.extracted_energy > 0.0 && penrose.eta_penrose > 0.0,
          "reference event reports positive net extracted energy");
    const auto outside_ergosphere = bh::evaluate_equatorial_penrose_event(
        {1.0, 0.999, 1.0, 0.0, 1.0, 10.0, 20.0, 0.002, 50'000, 1e-7},
        {2.1, 2.07, -2.0});
    check(outside_ergosphere.status == bh::PenroseEventStatus::outside_ergosphere,
          "split outside the equatorial ergosphere is rejected");

    const auto loaded_event = bh::load_equatorial_penrose_event_input(
        std::filesystem::path(BH_SOURCE_DIR) / "scenarios" / "equatorial_penrose_reference.cfg");
    near(loaded_event.scenario.dimensionless_spin, 0.999, 0.0,
         "scenario loader preserves the declared spin");
    near(loaded_event.split.split_radius_over_m, 1.10, 0.0,
         "scenario loader preserves the declared split radius");
    const auto loaded_penrose =
        bh::evaluate_equatorial_penrose_event(loaded_event.scenario, loaded_event.split);
    check(loaded_penrose.status == bh::PenroseEventStatus::physically_feasible,
          "versioned reference scenario evaluates as physically feasible");

    const auto weak = bh::estimate_plasma_extraction({0.0, 1.0, 10.0, 0.9, 2.0});
    near(weak.idealized_extracted_energy_joules, 0.0, 0.0, "zero field produces zero toy extraction");
    const auto plasma = bh::estimate_plasma_extraction({1.0, 1e-8, 10.0, 0.9, 2.0});
    check(plasma.alfven_speed_m_s > 0.0 && plasma.alfven_speed_m_s < bh::speed_of_light_m_s,
          "relativistic Alfven speed is causal");
    check(plasma.spin_coupling_efficiency >= 0.0 && plasma.spin_coupling_efficiency <= 1.0,
          "toy spin coupling is bounded");

    rejected = false;
    try {
        (void)bh::estimate_plasma_extraction(
            {std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0, 0.5, 1.0});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "plasma model rejects non-finite public input");

    rejected = false;
    try {
        (void)bh::kerr_outer_horizon(1.0, std::numeric_limits<double>::quiet_NaN());
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "Kerr model rejects non-finite public input");

    if (failures == 0) std::cout << "All model tests passed\n";
    return failures == 0 ? 0 : 1;
}
