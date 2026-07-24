#include "bh/kerr_geodesic.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <optional>
#include <stdexcept>

namespace bh {
namespace {
constexpr double potential_tolerance = 1.0e-12;
constexpr double horizon_event_relative_tolerance = 1.0e-6;

void validate_geometry(const double mass, const double spin_length) {
    if (!std::isfinite(mass) || !std::isfinite(spin_length) || mass <= 0.0 ||
        std::abs(spin_length) >= mass) {
        throw std::invalid_argument("Kerr requires finite M > 0 and finite |a| < M");
    }
}

void validate_orbit(const KerrOrbit& orbit) {
    validate_geometry(orbit.black_hole_mass, orbit.spin_length);
    if (!std::isfinite(orbit.energy) || !std::isfinite(orbit.angular_momentum) ||
        !std::isfinite(orbit.rest_mass) || !std::isfinite(orbit.carter_constant) ||
        orbit.rest_mass < 0.0 || std::abs(orbit.carter_constant) > potential_tolerance ||
        (orbit.radial_direction != -1 && orbit.radial_direction != 1)) {
        throw std::invalid_argument("invalid equatorial Kerr orbit parameters");
    }
}

double radial_p(const KerrOrbit& orbit, const double radius) {
    return orbit.energy * (radius * radius + orbit.spin_length * orbit.spin_length) -
           orbit.spin_length * orbit.angular_momentum;
}

double angular_offset(const KerrOrbit& orbit) {
    return orbit.angular_momentum - orbit.spin_length * orbit.energy;
}

struct Derivative {
    double radial{};
    double azimuth{};
    double coordinate_time{};
};

bool crosses(const double start, const double end, const double boundary, const int direction) {
    return direction > 0 ? start < boundary && end >= boundary
                         : start > boundary && end <= boundary;
}

TrajectoryPoint event_point(const TrajectoryPoint& state, const Derivative& derivative,
                            const double step, const double boundary) {
    const double projected_radius = state.radius + step * derivative.radial;
    const double denominator = projected_radius - state.radius;
    const double fraction = denominator == 0.0
                                ? 0.0
                                : std::clamp((boundary - state.radius) / denominator, 0.0, 1.0);
    const double event_step = fraction * step;
    return {state.affine_parameter + event_step,
            boundary,
            state.phi + event_step * derivative.azimuth,
            state.coordinate_time + event_step * derivative.coordinate_time,
            derivative.radial};
}

std::optional<Derivative> derivative(const KerrOrbit& orbit, const TrajectoryPoint& state) {
    const double horizon = kerr_outer_horizon(orbit.black_hole_mass, orbit.spin_length);
    const double horizon_event = horizon * (1.0 + horizon_event_relative_tolerance);
    if (!std::isfinite(state.radius) || state.radius <= horizon_event) {
        return std::nullopt;
    }

    const double potential = kerr_radial_potential(orbit, state.radius);
    if (potential < -potential_tolerance) {
        return std::nullopt;
    }
    if (potential <= potential_tolerance) {
        return Derivative{};
    }

    const auto momentum = kerr_equatorial_four_momentum(orbit, state.radius);
    return Derivative{momentum.radial, momentum.azimuth, momentum.coordinate_time};
}

Trajectory integrate_kerr_impl(const KerrOrbit& orbit, const double initial_radius,
                               const double step, const std::size_t max_steps,
                               const std::optional<double> escape_radius,
                               const std::optional<double> target_radius) {
    validate_orbit(orbit);
    const double horizon = kerr_outer_horizon(orbit.black_hole_mass, orbit.spin_length);
    const double horizon_event = horizon * (1.0 + horizon_event_relative_tolerance);
    if (!std::isfinite(initial_radius) || !std::isfinite(step) || initial_radius <= horizon_event ||
        step <= 0.0) {
        throw std::invalid_argument("invalid Kerr integration parameters");
    }
    if (escape_radius.has_value() &&
        (!std::isfinite(*escape_radius) || *escape_radius <= horizon_event)) {
        throw std::invalid_argument("escape radius must be finite and outside the horizon");
    }
    if (target_radius.has_value()) {
        if (!std::isfinite(*target_radius) || *target_radius <= horizon_event ||
            (*target_radius - initial_radius) * orbit.radial_direction <= 0.0) {
            throw std::invalid_argument("target radius must lie in the selected radial direction");
        }
    }

    Trajectory out;
    out.points.reserve(max_steps + 1);
    out.points.push_back({0.0, initial_radius, 0.0, 0.0, 0.0});

    for (std::size_t i = 0; i < max_steps; ++i) {
        const TrajectoryPoint state = out.points.back();
        const auto k1 = derivative(orbit, state);
        if (!k1.has_value()) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        if (k1->radial == 0.0) {
            out.termination = TrajectoryTermination::turning_point;
            break;
        }

        const double predicted_radius = state.radius + step * k1->radial;
        if (orbit.radial_direction < 0 && crosses(state.radius, predicted_radius, horizon_event,
                                                  orbit.radial_direction)) {
            out.points.push_back(event_point(state, *k1, step, horizon_event));
            out.termination = TrajectoryTermination::crossed_horizon;
            break;
        }
        if (target_radius.has_value() &&
            crosses(state.radius, predicted_radius, *target_radius, orbit.radial_direction)) {
            out.points.push_back(event_point(state, *k1, step, *target_radius));
            out.termination = TrajectoryTermination::reached_target_radius;
            break;
        }
        if (escape_radius.has_value() && orbit.radial_direction > 0 &&
            crosses(state.radius, predicted_radius, *escape_radius, orbit.radial_direction)) {
            out.points.push_back(event_point(state, *k1, step, *escape_radius));
            out.termination = TrajectoryTermination::reached_escape_radius;
            break;
        }

        const auto midpoint_one = TrajectoryPoint{state.affine_parameter + step / 2.0,
                                                  state.radius + step * k1->radial / 2.0,
                                                  state.phi + step * k1->azimuth / 2.0,
                                                  state.coordinate_time + step * k1->coordinate_time / 2.0,
                                                  k1->radial};
        const auto k2 = derivative(orbit, midpoint_one);
        if (!k2.has_value()) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        const auto midpoint_two = TrajectoryPoint{state.affine_parameter + step / 2.0,
                                                  state.radius + step * k2->radial / 2.0,
                                                  state.phi + step * k2->azimuth / 2.0,
                                                  state.coordinate_time + step * k2->coordinate_time / 2.0,
                                                  k2->radial};
        const auto k3 = derivative(orbit, midpoint_two);
        if (!k3.has_value()) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        const auto endpoint = TrajectoryPoint{state.affine_parameter + step,
                                              state.radius + step * k3->radial,
                                              state.phi + step * k3->azimuth,
                                              state.coordinate_time + step * k3->coordinate_time,
                                              k3->radial};
        const auto k4 = derivative(orbit, endpoint);
        if (!k4.has_value()) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }

        const Derivative average{
            (k1->radial + 2.0 * k2->radial + 2.0 * k3->radial + k4->radial) / 6.0,
            (k1->azimuth + 2.0 * k2->azimuth + 2.0 * k3->azimuth + k4->azimuth) / 6.0,
            (k1->coordinate_time + 2.0 * k2->coordinate_time + 2.0 * k3->coordinate_time +
             k4->coordinate_time) /
                6.0};
        TrajectoryPoint next{state.affine_parameter + step,
                             state.radius + step * average.radial,
                             state.phi + step * average.azimuth,
                             state.coordinate_time + step * average.coordinate_time,
                             average.radial};
        if (!std::isfinite(next.radius) || !std::isfinite(next.phi) ||
            !std::isfinite(next.coordinate_time)) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        out.points.push_back(next);

        if (orbit.radial_direction < 0 && next.radius <= horizon_event) {
            out.points.back().radius = horizon_event;
            out.termination = TrajectoryTermination::crossed_horizon;
            break;
        }
        if (target_radius.has_value() &&
            crosses(state.radius, next.radius, *target_radius, orbit.radial_direction)) {
            out.points.back().radius = *target_radius;
            out.termination = TrajectoryTermination::reached_target_radius;
            break;
        }
        if (escape_radius.has_value() && orbit.radial_direction > 0 &&
            next.radius >= *escape_radius) {
            out.points.back().radius = *escape_radius;
            out.termination = TrajectoryTermination::reached_escape_radius;
            break;
        }
    }
    return out;
}
}  // namespace

double kerr_spin_length(const double mass, const double dimensionless_spin) {
    if (!std::isfinite(mass) || !std::isfinite(dimensionless_spin) || mass <= 0.0 ||
        dimensionless_spin < 0.0 || dimensionless_spin >= 1.0) {
        throw std::invalid_argument("dimensionless Kerr spin must satisfy 0 <= a_star < 1");
    }
    return mass * dimensionless_spin;
}

double kerr_inner_horizon(const double mass, const double spin_length) {
    validate_geometry(mass, spin_length);
    return mass - std::sqrt(mass * mass - spin_length * spin_length);
}

double kerr_outer_horizon(const double mass, const double spin_length) {
    validate_geometry(mass, spin_length);
    return mass + std::sqrt(mass * mass - spin_length * spin_length);
}

double kerr_static_limit_radius(const double mass, const double spin_length,
                                const double polar_angle_radians) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(polar_angle_radians)) {
        throw std::invalid_argument("polar angle must be finite");
    }
    const double cosine = std::cos(polar_angle_radians);
    return mass + std::sqrt(mass * mass - spin_length * spin_length * cosine * cosine);
}

bool kerr_is_within_equatorial_ergosphere(const double mass, const double spin_length,
                                          const double radius) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(radius)) {
        throw std::invalid_argument("radius must be finite");
    }
    return radius > kerr_outer_horizon(mass, spin_length) && radius < 2.0 * mass;
}

double kerr_delta(const double mass, const double spin_length, const double radius) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(radius)) {
        throw std::invalid_argument("radius must be finite");
    }
    return radius * radius - 2.0 * mass * radius + spin_length * spin_length;
}

double kerr_equatorial_sigma(const double radius) {
    if (!std::isfinite(radius)) {
        throw std::invalid_argument("radius must be finite");
    }
    return radius * radius;
}

double kerr_radial_potential(const KerrOrbit& orbit, const double radius) {
    validate_orbit(orbit);
    if (!std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument("radius must be finite and positive");
    }
    const double delta = kerr_delta(orbit.black_hole_mass, orbit.spin_length, radius);
    const double p = radial_p(orbit, radius);
    const double offset = angular_offset(orbit);
    return p * p - delta * (orbit.rest_mass * orbit.rest_mass * radius * radius +
                            offset * offset + orbit.carter_constant);
}

KerrFourMomentum kerr_equatorial_four_momentum(const KerrOrbit& orbit,
                                                const double radius) {
    validate_orbit(orbit);
    const double horizon = kerr_outer_horizon(orbit.black_hole_mass, orbit.spin_length);
    if (!std::isfinite(radius) || radius <= horizon) {
        throw std::invalid_argument("equatorial Kerr momentum requires radius outside horizon");
    }
    const double delta = kerr_delta(orbit.black_hole_mass, orbit.spin_length, radius);
    const double sigma = kerr_equatorial_sigma(radius);
    const double potential = kerr_radial_potential(orbit, radius);
    if (potential < -potential_tolerance) {
        throw std::invalid_argument("Kerr radial potential is forbidden at this radius");
    }
    const double p = radial_p(orbit, radius);
    const double offset = angular_offset(orbit);
    return {(orbit.spin_length * offset + (radius * radius + orbit.spin_length * orbit.spin_length) *
                                          p / delta) /
                sigma,
            static_cast<double>(orbit.radial_direction) *
                std::sqrt(std::max(0.0, potential)) / sigma,
            (offset + orbit.spin_length * p / delta) / sigma};
}

Trajectory integrate_kerr(const KerrOrbit& orbit, const double initial_radius,
                          const double step, const std::size_t max_steps,
                          const double escape_radius) {
    return integrate_kerr_impl(orbit, initial_radius, step, max_steps, escape_radius,
                               std::nullopt);
}

Trajectory integrate_kerr_to_radius(const KerrOrbit& orbit, const double initial_radius,
                                    const double target_radius, const double step,
                                    const std::size_t max_steps) {
    return integrate_kerr_impl(orbit, initial_radius, step, max_steps, std::nullopt,
                               target_radius);
}
}  // namespace bh
