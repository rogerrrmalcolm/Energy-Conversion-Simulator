#include "bh/schwarzschild_geodesic.hpp"

#include <cmath>
#include <stdexcept>

namespace bh {
namespace {
struct Derivative { double r, phi, t, ur; };

Derivative derivative(const SchwarzschildOrbit& o, const TrajectoryPoint& s) {
    const double r = s.radius;
    const double m = o.black_hole_mass;
    const double l2 = o.specific_angular_momentum * o.specific_angular_momentum;
    return {s.radial_velocity,
            o.specific_angular_momentum / (r * r),
            o.specific_energy / (1.0 - 2.0 * m / r),
            -m / (r * r) + l2 / (r * r * r) - 3.0 * m * l2 / (r * r * r * r)};
}

TrajectoryPoint shifted(const TrajectoryPoint& s, const Derivative& d, double h) {
    return {s.affine_parameter + h, s.radius + h*d.r, s.phi + h*d.phi,
            s.coordinate_time + h*d.t, s.radial_velocity + h*d.ur};
}
}

Trajectory integrate_schwarzschild(const SchwarzschildOrbit& o, const TrajectoryPoint& initial,
                                   const double h, const std::size_t max_steps,
                                   const double escape_radius) {
    if (!std::isfinite(o.black_hole_mass) || !std::isfinite(o.specific_energy) ||
        !std::isfinite(o.specific_angular_momentum) || !std::isfinite(initial.radius) ||
        !std::isfinite(initial.radial_velocity) || !std::isfinite(h) ||
        !std::isfinite(escape_radius) || o.black_hole_mass <= 0.0 || h <= 0.0 ||
        initial.radius <= 2.0*o.black_hole_mass || escape_radius <= 2.0*o.black_hole_mass) {
        throw std::invalid_argument("invalid Schwarzschild integration parameters");
    }
    Trajectory out;
    out.points.reserve(max_steps + 1);
    out.points.push_back(initial);
    for (std::size_t i = 0; i < max_steps; ++i) {
        const auto& s = out.points.back();
        const auto k1 = derivative(o, s);
        const auto k2 = derivative(o, shifted(s, k1, h/2.0));
        const auto k3 = derivative(o, shifted(s, k2, h/2.0));
        const auto k4 = derivative(o, shifted(s, k3, h));
        Derivative k{(k1.r+2*k2.r+2*k3.r+k4.r)/6.0,
                     (k1.phi+2*k2.phi+2*k3.phi+k4.phi)/6.0,
                     (k1.t+2*k2.t+2*k3.t+k4.t)/6.0,
                     (k1.ur+2*k2.ur+2*k3.ur+k4.ur)/6.0};
        auto next = shifted(s, k, h);
        if (!std::isfinite(next.radius) || !std::isfinite(next.coordinate_time)) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        out.points.push_back(next);
        if (next.radius <= 2.0*o.black_hole_mass * (1.0 + 1e-6)) {
            out.termination = TrajectoryTermination::crossed_horizon;
            break;
        }
        if (next.radius >= escape_radius) {
            out.termination = TrajectoryTermination::reached_escape_radius;
            break;
        }
    }
    return out;
}
}
