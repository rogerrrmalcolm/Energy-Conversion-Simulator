#pragma once

#include "bh/trajectory.hpp"

#include <cstddef>

namespace bh {
// Equatorial Kerr geodesic in geometrized units (G=c=1), with Q=0.
struct KerrOrbit {
    double black_hole_mass{1.0};
    double spin_length{0.9};       // a, with |a| < M
    double energy{1.0};            // Conserved energy at infinity.
    double angular_momentum{0.0};  // Conserved axial angular momentum L_z.
    double rest_mass{1.0};         // Mass-shell parameter: 1 timelike, 0 null.
    int radial_direction{-1};      // -1 inward, +1 outward.
    double carter_constant{0.0};   // Must remain zero in the equatorial baseline.
};

struct KerrFourMomentum {
    double coordinate_time{};  // p^t = dt/dlambda
    double radial{};           // p^r = dr/dlambda
    double azimuth{};          // p^phi = dphi/dlambda
};

struct KerrIntegrationControl {
    // Error scales are applied to r and phi. Boyer-Lindquist t is reported but
    // excluded because it is coordinate-singular at the horizon; dr/dlambda is
    // derived directly from the radial potential at each state.
    double absolute_tolerance{1.0e-10};
    double relative_tolerance{1.0e-9};
    // A value of zero derives a lower bound from the requested maximum step.
    double minimum_step{0.0};
};

[[nodiscard]] double kerr_spin_length(double mass, double dimensionless_spin);
[[nodiscard]] double kerr_inner_horizon(double mass, double spin_length);
[[nodiscard]] double kerr_outer_horizon(double mass, double spin_length);
[[nodiscard]] double kerr_static_limit_radius(double mass, double spin_length,
                                                double polar_angle_radians);
[[nodiscard]] bool kerr_is_within_equatorial_ergosphere(double mass, double spin_length,
                                                         double radius);
[[nodiscard]] double kerr_delta(double mass, double spin_length, double radius);
[[nodiscard]] double kerr_equatorial_sigma(double radius);
[[nodiscard]] double kerr_radial_potential(const KerrOrbit& orbit, double radius);
[[nodiscard]] KerrFourMomentum kerr_equatorial_four_momentum(const KerrOrbit& orbit,
                                                               double radius);
[[nodiscard]] Trajectory integrate_kerr(const KerrOrbit& orbit, double initial_radius,
                                        double step, std::size_t max_steps,
                                        double escape_radius,
                                        const KerrIntegrationControl& control = {});
[[nodiscard]] Trajectory integrate_kerr_to_radius(const KerrOrbit& orbit,
                                                   double initial_radius,
                                                   double target_radius,
                                                   double step,
                                                   std::size_t max_steps,
                                                   const KerrIntegrationControl& control = {});
}
