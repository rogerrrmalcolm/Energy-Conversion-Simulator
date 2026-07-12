#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/plasma_model.hpp"
#include "bh/schwarzschild_geodesic.hpp"

#include <cmath>
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

    const auto uncertain = bh::rotational_energy(
        {bh::solar_mass_kg, 0.9, {0.8, 0.9, 0.99}});
    check(uncertain.rotational_energy_lower_joules < uncertain.rotational_energy_joules,
          "lower spin bound produces lower rotational reservoir");
    check(uncertain.rotational_energy_joules < uncertain.rotational_energy_upper_joules,
          "upper spin bound produces upper rotational reservoir");
    check(uncertain.d_rotational_energy_d_spin_joules > 0.0,
          "rotational reservoir sensitivity to spin is positive");

    bool rejected = false;
    try { (void)bh::rotational_energy(1.0, 1.0); } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "extremal spin is rejected by sub-extremal model");

    rejected = false;
    try {
        (void)bh::rotational_energy({1.0, 0.7, {0.6, 0.8, 0.9}});
    } catch (const std::invalid_argument&) { rejected = true; }
    check(rejected, "spin uncertainty central value must match the selected spin");

    const double r = 10.0;
    const double circular_l = std::sqrt(r*r/(r-3.0));
    const auto circular = bh::integrate_schwarzschild(
        {1.0, 0.956182887, circular_l}, {0.0, r, 0.0, 0.0, 0.0},
        1e-3, 1000, 20.0);
    near(circular.points.back().radius, r, 1e-10, "Schwarzschild circular orbit remains circular");

    near(bh::kerr_outer_horizon(1.0, 0.0), 2.0, 1e-14,
         "Kerr horizon reduces to Schwarzschild horizon");
    const bh::KerrOrbit outward{1.0, 0.5, 1.0, 0.0, 1.0, 1};
    check(bh::kerr_radial_potential(outward, 10.0) >= 0.0, "Kerr test orbit is admissible");
    const auto escaping = bh::integrate_kerr(outward, 10.0, 0.01, 10'000, 11.0);
    check(escaping.termination == bh::TrajectoryTermination::reached_escape_radius,
          "outward Kerr trajectory reaches escape radius");

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
