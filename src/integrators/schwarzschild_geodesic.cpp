#include "bh/schwarzschild_geodesic.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace bh {
namespace {
constexpr double radial_potential_tolerance = 1.0e-10;

double horizon_radius(const double mass) {
    const double horizon = 2.0 * mass;
    if (!std::isfinite(horizon)) {
        throw std::overflow_error("Schwarzschild horizon radius overflowed");
    }
    return horizon;
}

struct Derivative {
    double r;
    double phi;
    double t;
    double ur;
};

void validate_orbit(const SchwarzschildOrbit& orbit) {
    if (!std::isfinite(orbit.black_hole_mass) || !std::isfinite(orbit.specific_energy) ||
        !std::isfinite(orbit.specific_angular_momentum) || orbit.black_hole_mass <= 0.0 ||
        orbit.specific_energy <= 0.0) {
        throw std::invalid_argument("invalid Schwarzschild orbit parameters");
    }
}

Derivative derivative(const SchwarzschildOrbit& orbit, const TrajectoryPoint& state) {
    const double radius = state.radius;
    const double mass = orbit.black_hole_mass;
    const double angular_momentum_squared =
        orbit.specific_angular_momentum * orbit.specific_angular_momentum;
    return {state.radial_velocity,
            orbit.specific_angular_momentum / (radius * radius),
            orbit.specific_energy / (1.0 - 2.0 * mass / radius),
            -mass / (radius * radius) + angular_momentum_squared / (radius * radius * radius) -
                3.0 * mass * angular_momentum_squared /
                    (radius * radius * radius * radius)};
}

TrajectoryPoint shifted(const TrajectoryPoint& state, const Derivative& derivative,
                        const double step) {
    return {state.affine_parameter + step,
            state.radius + step * derivative.r,
            state.phi + step * derivative.phi,
            state.coordinate_time + step * derivative.t,
            state.radial_velocity + step * derivative.ur};
}

TrajectoryPoint interpolate_event(const TrajectoryPoint& start, const TrajectoryPoint& end,
                                  const double boundary) {
    const double denominator = end.radius - start.radius;
    const double fraction = denominator == 0.0
                                ? 0.0
                                : std::clamp((boundary - start.radius) / denominator, 0.0, 1.0);
    return {start.affine_parameter + fraction * (end.affine_parameter - start.affine_parameter),
            boundary,
            start.phi + fraction * (end.phi - start.phi),
            start.coordinate_time + fraction * (end.coordinate_time - start.coordinate_time),
            start.radial_velocity + fraction * (end.radial_velocity - start.radial_velocity)};
}

double normalized_radial_residual(const SchwarzschildOrbit& orbit,
                                  const TrajectoryPoint& point) {
    const double potential = schwarzschild_radial_potential(orbit, point.radius);
    const double velocity_squared = point.radial_velocity * point.radial_velocity;
    if (!std::isfinite(potential) || !std::isfinite(velocity_squared)) {
        return std::numeric_limits<double>::infinity();
    }
    const double scale = std::max({1.0, std::abs(potential), std::abs(velocity_squared)});
    return std::abs(velocity_squared - potential) / scale;
}

void record_radial_residual(Trajectory* trajectory, const SchwarzschildOrbit& orbit,
                            const TrajectoryPoint& point) {
    trajectory->diagnostics.maximum_normalized_radial_residual = std::max(
        trajectory->diagnostics.maximum_normalized_radial_residual,
        normalized_radial_residual(orbit, point));
}

double physical_radial_velocity(const SchwarzschildOrbit& orbit, const double radius,
                                const int direction) {
    const double potential = schwarzschild_radial_potential(orbit, radius);
    return static_cast<double>(direction) * std::sqrt(std::max(0.0, potential));
}
}  // namespace

double schwarzschild_radial_potential(const SchwarzschildOrbit& orbit, const double radius) {
    validate_orbit(orbit);
    if (!std::isfinite(radius) || radius <= horizon_radius(orbit.black_hole_mass)) {
        throw std::invalid_argument("Schwarzschild radial potential requires radius outside horizon");
    }
    const double lapse = 1.0 - 2.0 * orbit.black_hole_mass / radius;
    const double angular_term =
        orbit.specific_angular_momentum * orbit.specific_angular_momentum / (radius * radius);
    const double potential = orbit.specific_energy * orbit.specific_energy -
                             lapse * (1.0 + angular_term);
    if (!std::isfinite(potential)) {
        throw std::overflow_error("Schwarzschild radial potential overflowed");
    }
    return potential;
}

Trajectory integrate_schwarzschild(const SchwarzschildOrbit& orbit,
                                   const TrajectoryPoint& initial, const double step,
                                   const std::size_t max_steps, const double escape_radius) {
    validate_orbit(orbit);
    if (!std::isfinite(initial.affine_parameter) || !std::isfinite(initial.radius) ||
        !std::isfinite(initial.phi) || !std::isfinite(initial.coordinate_time) ||
        !std::isfinite(initial.radial_velocity) || !std::isfinite(step) ||
        !std::isfinite(escape_radius) || step <= 0.0 || max_steps == 0 ||
        initial.radius <= horizon_radius(orbit.black_hole_mass) ||
        escape_radius <= initial.radius) {
        throw std::invalid_argument("invalid Schwarzschild integration parameters");
    }

    const double initial_potential = schwarzschild_radial_potential(orbit, initial.radius);
    const double radial_velocity_squared = initial.radial_velocity * initial.radial_velocity;
    const double normalization_scale =
        std::max({1.0, std::abs(initial_potential), radial_velocity_squared});
    if (initial_potential < -radial_potential_tolerance ||
        std::abs(radial_velocity_squared - initial_potential) >
            radial_potential_tolerance * normalization_scale) {
        throw std::invalid_argument(
            "initial Schwarzschild radial velocity violates the timelike mass-shell relation");
    }

    Trajectory out;
    out.points.reserve(max_steps + 1);
    out.points.push_back(initial);
    record_radial_residual(&out, orbit, initial);
    out.diagnostics.final_step = step;
    const double horizon_event = horizon_radius(orbit.black_hole_mass) * (1.0 + 1.0e-6);
    if (!std::isfinite(horizon_event)) {
        throw std::overflow_error("Schwarzschild horizon event radius overflowed");
    }

    for (std::size_t i = 0; i < max_steps; ++i) {
        const TrajectoryPoint state = out.points.back();
        const Derivative k1 = derivative(orbit, state);
        const Derivative k2 = derivative(orbit, shifted(state, k1, step / 2.0));
        const Derivative k3 = derivative(orbit, shifted(state, k2, step / 2.0));
        const Derivative k4 = derivative(orbit, shifted(state, k3, step));
        const Derivative average{
            (k1.r + 2.0 * k2.r + 2.0 * k3.r + k4.r) / 6.0,
            (k1.phi + 2.0 * k2.phi + 2.0 * k3.phi + k4.phi) / 6.0,
            (k1.t + 2.0 * k2.t + 2.0 * k3.t + k4.t) / 6.0,
            (k1.ur + 2.0 * k2.ur + 2.0 * k3.ur + k4.ur) / 6.0};
        const TrajectoryPoint next = shifted(state, average, step);
        if (!std::isfinite(next.radius) || !std::isfinite(next.coordinate_time)) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        if (next.radius <= horizon_event) {
            TrajectoryPoint event = interpolate_event(state, next, horizon_event);
            event.radial_velocity = physical_radial_velocity(orbit, event.radius, -1);
            out.points.push_back(event);
            record_radial_residual(&out, orbit, event);
            ++out.diagnostics.accepted_steps;
            out.termination = TrajectoryTermination::crossed_horizon;
            break;
        }
        if (next.radius >= escape_radius) {
            TrajectoryPoint event = interpolate_event(state, next, escape_radius);
            event.radial_velocity = physical_radial_velocity(orbit, event.radius, 1);
            out.points.push_back(event);
            record_radial_residual(&out, orbit, event);
            ++out.diagnostics.accepted_steps;
            out.termination = TrajectoryTermination::reached_escape_radius;
            break;
        }
        out.points.push_back(next);
        record_radial_residual(&out, orbit, next);
        ++out.diagnostics.accepted_steps;
    }
    return out;
}
}  // namespace bh
