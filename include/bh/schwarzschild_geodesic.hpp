#pragma once

#include "bh/trajectory.hpp"

#include <cstddef>

namespace bh {
// Equatorial timelike Schwarzschild geodesic in geometrized units (G=c=1).
struct SchwarzschildOrbit {
    double black_hole_mass{1.0};
    double specific_energy{1.0};
    double specific_angular_momentum{0.0};
};

[[nodiscard]] Trajectory integrate_schwarzschild(
    const SchwarzschildOrbit& orbit, const TrajectoryPoint& initial,
    double step, std::size_t max_steps, double escape_radius);
}
