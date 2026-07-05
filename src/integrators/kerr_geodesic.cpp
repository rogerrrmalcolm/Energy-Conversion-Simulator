#include "bh/kerr_geodesic.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace bh {
double kerr_outer_horizon(const double m, const double a) {
    if (m <= 0.0 || std::abs(a) >= m) throw std::invalid_argument("Kerr requires M>0 and |a|<M");
    return m + std::sqrt(m*m - a*a);
}

double kerr_radial_potential(const KerrOrbit& o, const double r) {
    const double delta = r*r - 2.0*o.black_hole_mass*r + o.spin_length*o.spin_length;
    const double p = o.energy*(r*r + o.spin_length*o.spin_length)
                   - o.spin_length*o.angular_momentum;
    const double offset = o.angular_momentum - o.spin_length*o.energy;
    return p*p - delta*(o.rest_mass*o.rest_mass*r*r + offset*offset); // Q=0 equatorial
}

Trajectory integrate_kerr(const KerrOrbit& o, const double initial_radius, const double h,
                          const std::size_t max_steps, const double escape_radius) {
    const double horizon = kerr_outer_horizon(o.black_hole_mass, o.spin_length);
    if (initial_radius <= horizon || h <= 0.0 || (o.radial_direction != -1 && o.radial_direction != 1)) {
        throw std::invalid_argument("invalid Kerr integration parameters");
    }
    Trajectory out;
    out.points.reserve(max_steps + 1);
    out.points.push_back({0.0, initial_radius, 0.0, 0.0, 0.0});
    for (std::size_t i = 0; i < max_steps; ++i) {
        const auto s = out.points.back();
        const double r = s.radius;
        const double delta = r*r - 2.0*o.black_hole_mass*r + o.spin_length*o.spin_length;
        const double sigma = r*r;
        const double potential = kerr_radial_potential(o, r);
        if (potential < -1e-12 || delta <= 0.0) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        const double p = o.energy*(r*r + o.spin_length*o.spin_length)
                       - o.spin_length*o.angular_momentum;
        const double offset = o.angular_momentum - o.spin_length*o.energy;
        const double dr = o.radial_direction*std::sqrt(std::max(0.0, potential))/sigma;
        const double dphi = (-offset + o.spin_length*p/delta)/sigma;
        const double dt = (-o.spin_length*offset + (r*r+o.spin_length*o.spin_length)*p/delta)/sigma;
        TrajectoryPoint next{s.affine_parameter+h, r+h*dr, s.phi+h*dphi,
                             s.coordinate_time+h*dt, dr};
        if (!std::isfinite(next.radius) || !std::isfinite(next.coordinate_time)) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        out.points.push_back(next);
        if (next.radius <= horizon*(1.0+1e-6)) {
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
