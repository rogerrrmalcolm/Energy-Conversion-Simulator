#include "bh/algebraic_model.hpp"
#include "bh/constants.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/plasma_model.hpp"
#include "bh/schwarzschild_geodesic.hpp"

#include <iomanip>
#include <iostream>

int main() {
    std::cout << std::scientific << std::setprecision(6);

    const auto energy = bh::rotational_energy(10.0 * bh::solar_mass_kg, 0.9);
    std::cout << "Algebraic model\n  mass energy: " << energy.mass_energy_joules
              << " J\n  rotational reservoir: " << energy.rotational_energy_joules
              << " J (" << 100.0*energy.rotational_fraction << "%)\n\n";

    bh::TrajectoryPoint initial{0.0, 10.0, 0.0, 0.0, -0.05};
    const auto schwarzschild = bh::integrate_schwarzschild(
        {1.0, 0.96, 3.8}, initial, 0.001, 1'000, 20.0);
    std::cout << "Schwarzschild geodesic\n  points: " << schwarzschild.points.size()
              << "\n  final radius: " << schwarzschild.points.back().radius << " M\n\n";

    const auto kerr = bh::integrate_kerr({1.0, 0.9, 1.0, 2.0, 1.0, -1},
                                         8.0, 0.001, 1'000, 20.0);
    std::cout << "Kerr geodesic\n  outer horizon: " << bh::kerr_outer_horizon(1.0, 0.9)
              << " M\n  final radius: " << kerr.points.back().radius << " M\n\n";

    const auto plasma = bh::estimate_plasma_extraction({1.0, 1e-8, 1e6, 0.9, 1.0});
    std::cout << "Toy plasma model\n  magnetization: " << plasma.magnetization
              << "\n  idealized power: " << plasma.poynting_power_watts
              << " W\n  modeled extracted energy: "
              << plasma.idealized_extracted_energy_joules << " J\n";
}
