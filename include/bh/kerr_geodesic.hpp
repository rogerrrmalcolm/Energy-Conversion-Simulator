#pragma once

#include "bh/trajectory.hpp"

#include <cstddef>

namespace bh {
// First-pass equatorial Kerr geodesic using Carter's separated equations.
struct KerrOrbit {
    double black_hole_mass{1.0};
    double spin_length{0.9};       // a, with |a| < M
    double energy{1.0};            // conserved energy per unit rest mass
    double angular_momentum{0.0};  // conserved axial angular momentum
    double rest_mass{1.0};         // 1 timelike, 0 null
    int radial_direction{-1};      // -1 inward, +1 outward; no turning-point switch yet
};

[[nodiscard]] double kerr_outer_horizon(double mass, double spin_length);
[[nodiscard]] double kerr_radial_potential(const KerrOrbit& orbit, double radius);
[[nodiscard]] Trajectory integrate_kerr(const KerrOrbit& orbit, double initial_radius,
                                        double step, std::size_t max_steps,
                                        double escape_radius);
}
