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
constexpr double adaptation_safety = 0.9;

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

void validate_integration_control(const KerrIntegrationControl& control) {
    if (!std::isfinite(control.absolute_tolerance) ||
        !std::isfinite(control.relative_tolerance) ||
        !std::isfinite(control.minimum_step) || control.absolute_tolerance <= 0.0 ||
        control.relative_tolerance <= 0.0 || control.minimum_step < 0.0) {
        throw std::invalid_argument("invalid Kerr integration control");
    }
}

double finite_or_throw(const double value, const char* message) {
    if (!std::isfinite(value)) {
        throw std::overflow_error(message);
    }
    return value;
}

double horizon_root_fraction(const double mass, const double spin_length) {
    const double ratio = std::abs(spin_length) / mass;
    return std::sqrt(std::max(0.0, 1.0 - ratio * ratio));
}

double horizon_event_radius(const double horizon) {
    return finite_or_throw(horizon * (1.0 + horizon_event_relative_tolerance),
                           "Kerr horizon event radius overflowed");
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

enum class EvaluationStatus {
    valid,
    turning_point,
    forbidden_radial_region,
    invalid
};

struct DerivativeEvaluation {
    EvaluationStatus status{EvaluationStatus::invalid};
    Derivative value{};
};

struct StepResult {
    EvaluationStatus status{EvaluationStatus::invalid};
    TrajectoryPoint next{};
};

bool is_radial_boundary(const EvaluationStatus status) {
    return status == EvaluationStatus::turning_point ||
           status == EvaluationStatus::forbidden_radial_region;
}

bool crosses(const double start, const double end, const double boundary, const int direction) {
    return direction > 0 ? start < boundary && end >= boundary
                         : start > boundary && end <= boundary;
}

TrajectoryPoint advance(const TrajectoryPoint& state, const Derivative& derivative,
                        const double step) {
    return {state.affine_parameter + step,
            state.radius + step * derivative.radial,
            state.phi + step * derivative.azimuth,
            state.coordinate_time + step * derivative.coordinate_time,
            derivative.radial};
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

TrajectoryPoint horizon_event_point(const TrajectoryPoint& state, const Derivative& derivative,
                                    const double step, const double horizon_event) {
    const TrajectoryPoint projected = advance(state, derivative, step);
    return interpolate_event(state, projected, horizon_event);
}

std::optional<double> locate_turning_radius(const KerrOrbit& orbit,
                                            const double start_radius,
                                            const double end_radius) {
    if (!std::isfinite(start_radius) || !std::isfinite(end_radius) ||
        start_radius <= 0.0 || end_radius <= 0.0 || start_radius == end_radius) {
        return std::nullopt;
    }

    if (kerr_radial_potential(orbit, start_radius) <= potential_tolerance) {
        return start_radius;
    }

    // Sampling first catches a forbidden band even when R becomes positive again
    // before the horizon. Bisection then localizes the first R(r) = 0 boundary.
    constexpr int sample_count = 64;
    double allowed_radius = start_radius;
    for (int sample = 1; sample <= sample_count; ++sample) {
        const double fraction = static_cast<double>(sample) / sample_count;
        const double candidate_radius =
            start_radius + fraction * (end_radius - start_radius);
        if (kerr_radial_potential(orbit, candidate_radius) > potential_tolerance) {
            allowed_radius = candidate_radius;
            continue;
        }

        double blocked_radius = candidate_radius;
        for (int iteration = 0; iteration < 64; ++iteration) {
            const double midpoint_radius = (allowed_radius + blocked_radius) / 2.0;
            if (kerr_radial_potential(orbit, midpoint_radius) > potential_tolerance) {
                allowed_radius = midpoint_radius;
            } else {
                blocked_radius = midpoint_radius;
            }
        }
        return (allowed_radius + blocked_radius) / 2.0;
    }
    return std::nullopt;
}

std::optional<TrajectoryPoint> locate_turning_point(const KerrOrbit& orbit,
                                                     const TrajectoryPoint& state,
                                                     const TrajectoryPoint& probe) {
    const auto radius = locate_turning_radius(orbit, state.radius, probe.radius);
    if (!radius.has_value()) {
        return std::nullopt;
    }

    TrajectoryPoint event = interpolate_event(state, probe, *radius);
    event.radial_velocity = 0.0;
    return event;
}

DerivativeEvaluation derivative(const KerrOrbit& orbit, const TrajectoryPoint& state) {
    const double horizon = kerr_outer_horizon(orbit.black_hole_mass, orbit.spin_length);
    const double horizon_event = horizon_event_radius(horizon);
    if (!std::isfinite(state.affine_parameter) || !std::isfinite(state.radius) ||
        !std::isfinite(state.phi) || !std::isfinite(state.coordinate_time) ||
        state.radius <= horizon_event) {
        return {};
    }

    const double potential = kerr_radial_potential(orbit, state.radius);
    if (potential < -potential_tolerance) {
        return {EvaluationStatus::forbidden_radial_region, {}};
    }
    if (potential <= potential_tolerance) {
        return {EvaluationStatus::turning_point, {}};
    }

    const KerrFourMomentum momentum = kerr_equatorial_four_momentum(orbit, state.radius);
    if (!std::isfinite(momentum.radial) || !std::isfinite(momentum.azimuth) ||
        !std::isfinite(momentum.coordinate_time)) {
        return {};
    }
    return {EvaluationStatus::valid,
            {momentum.radial, momentum.azimuth, momentum.coordinate_time}};
}

StepResult rk4_step(const KerrOrbit& orbit, const TrajectoryPoint& state, const double step) {
    const DerivativeEvaluation k1 = derivative(orbit, state);
    if (k1.status != EvaluationStatus::valid) {
        return {k1.status, {}};
    }
    const DerivativeEvaluation k2 = derivative(orbit, advance(state, k1.value, step / 2.0));
    if (k2.status != EvaluationStatus::valid) {
        return {k2.status, {}};
    }
    const DerivativeEvaluation k3 = derivative(orbit, advance(state, k2.value, step / 2.0));
    if (k3.status != EvaluationStatus::valid) {
        return {k3.status, {}};
    }
    const DerivativeEvaluation k4 = derivative(orbit, advance(state, k3.value, step));
    if (k4.status != EvaluationStatus::valid) {
        return {k4.status, {}};
    }

    const Derivative average{
        (k1.value.radial + 2.0 * k2.value.radial + 2.0 * k3.value.radial +
         k4.value.radial) /
            6.0,
        (k1.value.azimuth + 2.0 * k2.value.azimuth + 2.0 * k3.value.azimuth +
         k4.value.azimuth) /
            6.0,
        (k1.value.coordinate_time + 2.0 * k2.value.coordinate_time +
         2.0 * k3.value.coordinate_time + k4.value.coordinate_time) /
            6.0};
    const TrajectoryPoint next = advance(state, average, step);
    if (!std::isfinite(next.radius) || !std::isfinite(next.phi) ||
        !std::isfinite(next.coordinate_time)) {
        return {};
    }
    return {EvaluationStatus::valid, next};
}

std::optional<TrajectoryPoint> locate_radius_event(const KerrOrbit& orbit,
                                                    const TrajectoryPoint& state,
                                                    const double step,
                                                    const double boundary,
                                                    const int direction) {
    double lower_step = 0.0;
    double upper_step = step;
    TrajectoryPoint event;
    bool found = false;

    for (int iteration = 0; iteration < 48; ++iteration) {
        const double midpoint_step = (lower_step + upper_step) / 2.0;
        const StepResult midpoint = rk4_step(orbit, state, midpoint_step);
        if (midpoint.status != EvaluationStatus::valid) {
            return std::nullopt;
        }
        if (crosses(state.radius, midpoint.next.radius, boundary, direction)) {
            upper_step = midpoint_step;
            event = midpoint.next;
            found = true;
        } else {
            lower_step = midpoint_step;
        }
    }

    if (!found) {
        return std::nullopt;
    }
    event.radius = boundary;
    const DerivativeEvaluation event_derivative = derivative(orbit, event);
    if (event_derivative.status == EvaluationStatus::valid) {
        event.radial_velocity = event_derivative.value.radial;
    }
    return event;
}

double normalized_component_error(const double left, const double right,
                                  const KerrIntegrationControl& control) {
    const double scale = control.absolute_tolerance +
                         control.relative_tolerance * std::max(std::abs(left), std::abs(right));
    return std::abs(left - right) / scale;
}

double normalized_rk4_error(const TrajectoryPoint& full_step,
                            const TrajectoryPoint& two_half_steps,
                            const KerrIntegrationControl& control) {
    // Step doubling estimates an RK4 local error; the 15 accounts for 2^4 - 1.
    return std::max(
        {normalized_component_error(full_step.radius, two_half_steps.radius, control) / 15.0,
         normalized_component_error(full_step.phi, two_half_steps.phi, control) / 15.0});
}

double normalized_radial_residual(const KerrOrbit& orbit, const TrajectoryPoint& point) {
    const double sigma = kerr_equatorial_sigma(point.radius);
    const double radial_term = sigma * point.radial_velocity;
    const double radial_term_squared = radial_term * radial_term;
    const double potential = kerr_radial_potential(orbit, point.radius);
    if (!std::isfinite(radial_term_squared) || !std::isfinite(potential)) {
        return std::numeric_limits<double>::infinity();
    }
    const double scale = std::max({1.0, std::abs(radial_term_squared), std::abs(potential)});
    return std::abs(radial_term_squared - potential) / scale;
}

void record_radial_residual(Trajectory* trajectory, const KerrOrbit& orbit,
                            const TrajectoryPoint& point) {
    trajectory->diagnostics.maximum_normalized_radial_residual = std::max(
        trajectory->diagnostics.maximum_normalized_radial_residual,
        normalized_radial_residual(orbit, point));
}

double physical_radial_velocity(const KerrOrbit& orbit, const double radius) {
    return kerr_equatorial_four_momentum(orbit, radius).radial;
}

double next_step_size(const double current_step, const double normalized_error,
                      const bool accepted, const double maximum_step) {
    if (normalized_error == 0.0) {
        return accepted ? maximum_step : current_step / 2.0;
    }
    const double raw_factor = adaptation_safety * std::pow(1.0 / normalized_error, 0.2);
    const double lower = accepted ? 0.2 : 0.1;
    const double upper = accepted ? 5.0 : 0.5;
    return std::min(maximum_step,
                    current_step * std::clamp(raw_factor, lower, upper));
}

double automatic_minimum_step(const double maximum_step, const double initial_radius) {
    const double scale = std::max({1.0, std::abs(maximum_step), std::abs(initial_radius)});
    return std::max(maximum_step * 1.0e-10,
                    32.0 * std::numeric_limits<double>::epsilon() * scale);
}

Trajectory integrate_kerr_impl(const KerrOrbit& orbit, const double initial_radius,
                               const double maximum_step, const std::size_t max_steps,
                               const std::optional<double> escape_radius,
                               const std::optional<double> target_radius,
                               const KerrIntegrationControl& control) {
    validate_orbit(orbit);
    validate_integration_control(control);

    const double horizon = kerr_outer_horizon(orbit.black_hole_mass, orbit.spin_length);
    const double horizon_event = horizon_event_radius(horizon);
    if (!std::isfinite(initial_radius) || !std::isfinite(maximum_step) ||
        initial_radius <= horizon_event || maximum_step <= 0.0 || max_steps == 0) {
        throw std::invalid_argument("invalid Kerr integration parameters");
    }
    const double minimum_step = control.minimum_step == 0.0
                                    ? automatic_minimum_step(maximum_step, initial_radius)
                                    : control.minimum_step;
    if (minimum_step > maximum_step) {
        throw std::invalid_argument("Kerr minimum step cannot exceed the maximum step");
    }
    if (escape_radius.has_value() &&
        (!std::isfinite(*escape_radius) || *escape_radius <= initial_radius)) {
        throw std::invalid_argument("escape radius must be finite and greater than the initial radius");
    }
    if (target_radius.has_value()) {
        if (!std::isfinite(*target_radius) || *target_radius <= horizon_event ||
            (*target_radius - initial_radius) * orbit.radial_direction <= 0.0) {
            throw std::invalid_argument("target radius must lie in the selected radial direction");
        }
    }

    Trajectory out;
    out.points.reserve(max_steps + 1);
    out.diagnostics.final_step = maximum_step;
    TrajectoryPoint initial{0.0, initial_radius, 0.0, 0.0, 0.0};
    const DerivativeEvaluation initial_derivative = derivative(orbit, initial);
    if (initial_derivative.status == EvaluationStatus::turning_point) {
        initial.radial_velocity = 0.0;
        out.points.push_back(initial);
        record_radial_residual(&out, orbit, initial);
        out.termination = TrajectoryTermination::turning_point;
        return out;
    }
    if (initial_derivative.status != EvaluationStatus::valid) {
        out.points.push_back(initial);
        out.termination = TrajectoryTermination::invalid_state;
        return out;
    }
    initial.radial_velocity = initial_derivative.value.radial;
    out.points.push_back(initial);
    record_radial_residual(&out, orbit, initial);

    double step = maximum_step;
    std::size_t attempts = 0;
    const std::size_t max_attempts =
        max_steps > std::numeric_limits<std::size_t>::max() / 16
            ? std::numeric_limits<std::size_t>::max()
            : max_steps * 16;

    while (out.diagnostics.accepted_steps < max_steps) {
        if (attempts++ >= max_attempts) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }

        const TrajectoryPoint state = out.points.back();
        const DerivativeEvaluation k1 = derivative(orbit, state);
        if (k1.status == EvaluationStatus::turning_point) {
            out.points.back().radial_velocity = 0.0;
            record_radial_residual(&out, orbit, out.points.back());
            out.termination = TrajectoryTermination::turning_point;
            break;
        }
        if (k1.status != EvaluationStatus::valid) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }

        const double projected_radius = state.radius + step * k1.value.radial;
        if (orbit.radial_direction < 0 &&
            crosses(state.radius, projected_radius, horizon_event, orbit.radial_direction)) {
            const TrajectoryPoint horizon_probe =
                interpolate_event(state, advance(state, k1.value, step), horizon_event);
            if (const auto turning = locate_turning_point(orbit, state, horizon_probe);
                turning.has_value()) {
                out.points.push_back(*turning);
                record_radial_residual(&out, orbit, *turning);
                ++out.diagnostics.accepted_steps;
                out.termination = TrajectoryTermination::turning_point;
                break;
            }
            TrajectoryPoint event = horizon_event_point(state, k1.value, step, horizon_event);
            event.radial_velocity = physical_radial_velocity(orbit, event.radius);
            out.points.push_back(event);
            record_radial_residual(&out, orbit, event);
            ++out.diagnostics.accepted_steps;
            out.termination = TrajectoryTermination::crossed_horizon;
            break;
        }

        const StepResult full_step = rk4_step(orbit, state, step);
        const StepResult first_half_step = rk4_step(orbit, state, step / 2.0);
        StepResult second_half_step;
        if (first_half_step.status == EvaluationStatus::valid) {
            second_half_step = rk4_step(orbit, first_half_step.next, step / 2.0);
        } else {
            second_half_step.status = first_half_step.status;
        }

        if (is_radial_boundary(full_step.status) || is_radial_boundary(first_half_step.status) ||
            is_radial_boundary(second_half_step.status)) {
            if (step > minimum_step) {
                step = std::max(minimum_step, step / 2.0);
                ++out.diagnostics.rejected_steps;
                continue;
            }
            const TrajectoryPoint probe = advance(state, k1.value, step);
            if (const auto turning = locate_turning_point(orbit, state, probe);
                turning.has_value()) {
                out.points.push_back(*turning);
                record_radial_residual(&out, orbit, *turning);
                ++out.diagnostics.accepted_steps;
                out.termination = TrajectoryTermination::turning_point;
                break;
            }
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        if (full_step.status != EvaluationStatus::valid ||
            first_half_step.status != EvaluationStatus::valid ||
            second_half_step.status != EvaluationStatus::valid) {
            if (step <= minimum_step) {
                out.termination = TrajectoryTermination::invalid_state;
                break;
            }
            step = std::max(minimum_step, step / 2.0);
            ++out.diagnostics.rejected_steps;
            continue;
        }

        const double error = normalized_rk4_error(full_step.next, second_half_step.next, control);
        if (!std::isfinite(error)) {
            out.termination = TrajectoryTermination::invalid_state;
            break;
        }
        if (error > 1.0) {
            if (step <= minimum_step) {
                out.termination = TrajectoryTermination::invalid_state;
                break;
            }
            step = std::max(minimum_step, next_step_size(step, error, false, maximum_step));
            ++out.diagnostics.rejected_steps;
            continue;
        }

        TrajectoryPoint next = second_half_step.next;
        next.radial_velocity = physical_radial_velocity(orbit, next.radius);
        out.diagnostics.maximum_normalized_error =
            std::max(out.diagnostics.maximum_normalized_error, error);
        ++out.diagnostics.accepted_steps;

        if (target_radius.has_value() &&
            crosses(state.radius, next.radius, *target_radius, orbit.radial_direction)) {
            const auto event = locate_radius_event(orbit, state, step, *target_radius,
                                                   orbit.radial_direction);
            TrajectoryPoint target_event = event.has_value()
                                               ? *event
                                               : interpolate_event(state, next, *target_radius);
            target_event.radial_velocity = physical_radial_velocity(orbit, target_event.radius);
            out.points.push_back(target_event);
            record_radial_residual(&out, orbit, target_event);
            out.termination = TrajectoryTermination::reached_target_radius;
            break;
        }
        if (escape_radius.has_value() && orbit.radial_direction > 0 &&
            crosses(state.radius, next.radius, *escape_radius, orbit.radial_direction)) {
            const auto event = locate_radius_event(orbit, state, step, *escape_radius,
                                                   orbit.radial_direction);
            TrajectoryPoint escape_event = event.has_value()
                                               ? *event
                                               : interpolate_event(state, next, *escape_radius);
            escape_event.radial_velocity = physical_radial_velocity(orbit, escape_event.radius);
            out.points.push_back(escape_event);
            record_radial_residual(&out, orbit, escape_event);
            out.termination = TrajectoryTermination::reached_escape_radius;
            break;
        }

        out.points.push_back(next);
        record_radial_residual(&out, orbit, next);
        step = std::max(minimum_step, next_step_size(step, error, true, maximum_step));
    }

    out.diagnostics.final_step = step;
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
    return finite_or_throw(mass * (1.0 - horizon_root_fraction(mass, spin_length)),
                           "Kerr inner horizon overflowed");
}

double kerr_outer_horizon(const double mass, const double spin_length) {
    validate_geometry(mass, spin_length);
    return finite_or_throw(mass * (1.0 + horizon_root_fraction(mass, spin_length)),
                           "Kerr outer horizon overflowed");
}

double kerr_static_limit_radius(const double mass, const double spin_length,
                                const double polar_angle_radians) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(polar_angle_radians)) {
        throw std::invalid_argument("polar angle must be finite");
    }
    const double ratio = (spin_length / mass) * std::cos(polar_angle_radians);
    const double static_limit = mass * (1.0 + std::sqrt(std::max(0.0, 1.0 - ratio * ratio)));
    return finite_or_throw(static_limit, "Kerr static limit overflowed");
}

bool kerr_is_within_equatorial_ergosphere(const double mass, const double spin_length,
                                           const double radius) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(radius)) {
        throw std::invalid_argument("radius must be finite");
    }
    const double static_limit = finite_or_throw(2.0 * mass, "Kerr equatorial static limit overflowed");
    return radius > kerr_outer_horizon(mass, spin_length) && radius < static_limit;
}

double kerr_delta(const double mass, const double spin_length, const double radius) {
    validate_geometry(mass, spin_length);
    if (!std::isfinite(radius)) {
        throw std::invalid_argument("radius must be finite");
    }
    return finite_or_throw(
        radius * radius - 2.0 * mass * radius + spin_length * spin_length,
        "Kerr delta overflowed");
}

double kerr_equatorial_sigma(const double radius) {
    if (!std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument("equatorial sigma requires a finite positive radius");
    }
    return finite_or_throw(radius * radius, "Kerr equatorial sigma overflowed");
}

double kerr_radial_potential(const KerrOrbit& orbit, const double radius) {
    validate_orbit(orbit);
    if (!std::isfinite(radius) || radius <= 0.0) {
        throw std::invalid_argument("radius must be finite and positive");
    }
    const double delta = kerr_delta(orbit.black_hole_mass, orbit.spin_length, radius);
    const double p = radial_p(orbit, radius);
    const double offset = angular_offset(orbit);
    return finite_or_throw(
        p * p - delta * (orbit.rest_mass * orbit.rest_mass * radius * radius +
                         offset * offset + orbit.carter_constant),
        "Kerr radial potential overflowed");
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
    const KerrFourMomentum momentum{
        (orbit.spin_length * offset +
         (radius * radius + orbit.spin_length * orbit.spin_length) * p / delta) /
            sigma,
        static_cast<double>(orbit.radial_direction) *
            std::sqrt(std::max(0.0, potential)) / sigma,
        (offset + orbit.spin_length * p / delta) / sigma};
    if (!std::isfinite(momentum.coordinate_time) || !std::isfinite(momentum.radial) ||
        !std::isfinite(momentum.azimuth)) {
        throw std::overflow_error("Kerr four-momentum overflowed");
    }
    if (momentum.coordinate_time <= 0.0) {
        throw std::invalid_argument("equatorial Kerr momentum must be future-directed");
    }
    return momentum;
}

Trajectory integrate_kerr(const KerrOrbit& orbit, const double initial_radius,
                          const double step, const std::size_t max_steps,
                          const double escape_radius,
                          const KerrIntegrationControl& control) {
    return integrate_kerr_impl(orbit, initial_radius, step, max_steps, escape_radius,
                               std::nullopt, control);
}

Trajectory integrate_kerr_to_radius(const KerrOrbit& orbit, const double initial_radius,
                                    const double target_radius, const double step,
                                    const std::size_t max_steps,
                                    const KerrIntegrationControl& control) {
    return integrate_kerr_impl(orbit, initial_radius, step, max_steps, std::nullopt,
                               target_radius, control);
}
}  // namespace bh
