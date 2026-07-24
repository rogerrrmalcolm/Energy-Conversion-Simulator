#include "bh/penrose_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace bh {
namespace {
struct LocalMomentum {
    double time{};
    double radial{};
    double azimuth{};
};

struct EquatorialMetric {
    double delta{};
    double g_tt{};
    double g_tphi{};
    double g_phiphi{};
    double g_rr{};
    double lapse{};
    double frame_dragging_omega{};
};

struct ConservedConstants {
    double energy{};
    double angular_momentum{};
};

double minkowski_dot(const LocalMomentum& left, const LocalMomentum& right) {
    return -left.time * right.time + left.radial * right.radial + left.azimuth * right.azimuth;
}

LocalMomentum add(const LocalMomentum& left, const LocalMomentum& right) {
    return {left.time + right.time, left.radial + right.radial, left.azimuth + right.azimuth};
}

LocalMomentum subtract(const LocalMomentum& left, const LocalMomentum& right) {
    return {left.time - right.time, left.radial - right.radial, left.azimuth - right.azimuth};
}

LocalMomentum scale(const LocalMomentum& value, const double factor) {
    return {value.time * factor, value.radial * factor, value.azimuth * factor};
}

double max_component_abs(const LocalMomentum& value) {
    return std::max({std::abs(value.time), std::abs(value.radial), std::abs(value.azimuth)});
}

double spacelike_norm(const LocalMomentum& value) {
    return std::sqrt(std::max(0.0, minkowski_dot(value, value)));
}

bool finite_split_parameters(const PenroseSplitParameters& split) {
    return std::isfinite(split.split_radius_over_m) &&
           std::isfinite(split.incoming_lz_over_m_m) &&
           std::isfinite(split.split_angle_rad);
}

void validate_scenario(const EquatorialPenroseScenario& scenario) {
    if (!std::isfinite(scenario.black_hole_mass) ||
        !std::isfinite(scenario.dimensionless_spin) ||
        !std::isfinite(scenario.parent_rest_mass) ||
        !std::isfinite(scenario.fragment_rest_mass) ||
        !std::isfinite(scenario.incoming_specific_energy) ||
        !std::isfinite(scenario.initial_radius_over_m) ||
        !std::isfinite(scenario.escape_radius_over_m) ||
        !std::isfinite(scenario.integration_step) || !std::isfinite(scenario.residual_tolerance) ||
        scenario.black_hole_mass <= 0.0 ||
        scenario.dimensionless_spin < 0.0 || scenario.dimensionless_spin >= 1.0 ||
        scenario.parent_rest_mass <= 0.0 || scenario.fragment_rest_mass < 0.0 ||
        2.0 * scenario.fragment_rest_mass > scenario.parent_rest_mass ||
        scenario.incoming_specific_energy <= 0.0 || scenario.initial_radius_over_m <= 0.0 ||
        scenario.escape_radius_over_m <= scenario.initial_radius_over_m ||
        scenario.integration_step <= 0.0 || scenario.max_integration_steps == 0 ||
        scenario.residual_tolerance <= 0.0) {
        throw std::invalid_argument("invalid equatorial Penrose scenario");
    }
}

EquatorialMetric equatorial_metric(const double mass, const double spin_length,
                                   const double radius) {
    const double delta = kerr_delta(mass, spin_length, radius);
    if (delta <= 0.0) {
        throw std::invalid_argument("equatorial metric requires a radius outside the horizon");
    }
    const double radius_squared = radius * radius;
    const double spin_squared = spin_length * spin_length;
    const double a_term = (radius_squared + spin_squared) * (radius_squared + spin_squared) -
                          spin_squared * delta;
    const double g_phiphi = a_term / radius_squared;
    return {delta,
            -(1.0 - 2.0 * mass / radius),
            -2.0 * mass * spin_length / radius,
            g_phiphi,
            radius_squared / delta,
            radius * std::sqrt(delta) / std::sqrt(a_term),
            2.0 * mass * spin_length * radius / a_term};
}

LocalMomentum coordinate_to_zamo(const double mass, const double spin_length,
                                 const double radius, const KerrFourMomentum& momentum) {
    const auto metric = equatorial_metric(mass, spin_length, radius);
    return {metric.lapse * momentum.coordinate_time,
            radius * momentum.radial / std::sqrt(metric.delta),
            std::sqrt(metric.g_phiphi) *
                (momentum.azimuth - metric.frame_dragging_omega * momentum.coordinate_time)};
}

KerrFourMomentum zamo_to_coordinate(const double mass, const double spin_length,
                                    const double radius, const LocalMomentum& momentum) {
    const auto metric = equatorial_metric(mass, spin_length, radius);
    const double coordinate_time = momentum.time / metric.lapse;
    return {coordinate_time,
            momentum.radial * std::sqrt(metric.delta) / radius,
            metric.frame_dragging_omega * coordinate_time +
                momentum.azimuth / std::sqrt(metric.g_phiphi)};
}

ConservedConstants conserved_constants(const double mass, const double spin_length,
                                       const double radius, const KerrFourMomentum& momentum) {
    const auto metric = equatorial_metric(mass, spin_length, radius);
    const double covariant_time = metric.g_tt * momentum.coordinate_time +
                                  metric.g_tphi * momentum.azimuth;
    const double covariant_azimuth = metric.g_tphi * momentum.coordinate_time +
                                      metric.g_phiphi * momentum.azimuth;
    return {-covariant_time, covariant_azimuth};
}

double coordinate_mass_shell_residual(const double mass, const double spin_length,
                                      const double radius, const KerrFourMomentum& momentum,
                                      const double rest_mass) {
    const auto metric = equatorial_metric(mass, spin_length, radius);
    const double squared_norm =
        metric.g_tt * momentum.coordinate_time * momentum.coordinate_time +
        2.0 * metric.g_tphi * momentum.coordinate_time * momentum.azimuth +
        metric.g_rr * momentum.radial * momentum.radial +
        metric.g_phiphi * momentum.azimuth * momentum.azimuth;
    return std::abs(squared_norm + rest_mass * rest_mass);
}

double coordinate_momentum_difference(const KerrFourMomentum& left,
                                      const KerrFourMomentum& right) {
    return std::max({std::abs(left.coordinate_time - right.coordinate_time),
                     std::abs(left.radial - right.radial),
                     std::abs(left.azimuth - right.azimuth)});
}

bool normalize_spacelike(LocalMomentum* vector) {
    const double norm = spacelike_norm(*vector);
    if (!std::isfinite(norm) || norm <= std::numeric_limits<double>::epsilon()) {
        return false;
    }
    *vector = scale(*vector, 1.0 / norm);
    return true;
}

LocalMomentum project_orthogonal_to_timelike(const LocalMomentum& seed,
                                             const LocalMomentum& timelike_unit) {
    return add(seed, scale(timelike_unit, minkowski_dot(timelike_unit, seed)));
}

bool build_equatorial_split_basis(const LocalMomentum& parent_unit_velocity,
                                  LocalMomentum* radial_basis,
                                  LocalMomentum* azimuth_basis) {
    *radial_basis = project_orthogonal_to_timelike({0.0, 1.0, 0.0}, parent_unit_velocity);
    if (!normalize_spacelike(radial_basis)) {
        return false;
    }

    *azimuth_basis = project_orthogonal_to_timelike({0.0, 0.0, 1.0}, parent_unit_velocity);
    *azimuth_basis = subtract(*azimuth_basis,
                              scale(*radial_basis, minkowski_dot(*radial_basis, *azimuth_basis)));
    return normalize_spacelike(azimuth_basis);
}

PenroseEventStatus failed_integration_status(const Trajectory& trajectory) {
    return trajectory.termination == TrajectoryTermination::invalid_state ||
                   trajectory.termination == TrajectoryTermination::completed
               ? PenroseEventStatus::integration_failed
               : PenroseEventStatus::captured_or_non_escaping;
}
}  // namespace

std::string_view penrose_event_status_name(const PenroseEventStatus status) {
    switch (status) {
    case PenroseEventStatus::outside_ergosphere:
        return "outside_ergosphere";
    case PenroseEventStatus::physics_invalid:
        return "physics_invalid";
    case PenroseEventStatus::captured_or_non_escaping:
        return "captured_or_non_escaping";
    case PenroseEventStatus::physically_feasible:
        return "physically_feasible";
    case PenroseEventStatus::integration_failed:
        return "integration_failed";
    }
    return "unknown";
}

PenroseEventResult evaluate_equatorial_penrose_event(
    const EquatorialPenroseScenario& scenario, const PenroseSplitParameters& split) {
    validate_scenario(scenario);

    PenroseEventResult result;
    result.split = split;
    if (!finite_split_parameters(split)) {
        result.status = PenroseEventStatus::physics_invalid;
        return result;
    }

    const double mass = scenario.black_hole_mass;
    const double spin_length = kerr_spin_length(mass, scenario.dimensionless_spin);
    const double horizon = kerr_outer_horizon(mass, spin_length);
    const double static_limit = kerr_static_limit_radius(mass, spin_length, 1.57079632679489661923);
    const double split_radius = split.split_radius_over_m * mass;
    result.horizon_radius = horizon;
    result.static_limit_radius = static_limit;
    result.split_radius = split_radius;
    if (!(split_radius > horizon && split_radius < static_limit)) {
        result.status = PenroseEventStatus::outside_ergosphere;
        return result;
    }

    const double initial_radius = scenario.initial_radius_over_m * mass;
    const double escape_radius = scenario.escape_radius_over_m * mass;
    if (initial_radius <= split_radius || escape_radius <= static_limit) {
        result.status = PenroseEventStatus::physics_invalid;
        return result;
    }

    const double input_energy = scenario.incoming_specific_energy * scenario.parent_rest_mass;
    const double input_angular_momentum =
        split.incoming_lz_over_m_m * scenario.parent_rest_mass * mass;
    result.input_energy = input_energy;

    try {
        const KerrOrbit incoming{mass, spin_length, input_energy, input_angular_momentum,
                                 scenario.parent_rest_mass, -1};
        result.incoming_trajectory = integrate_kerr_to_radius(
            incoming, initial_radius, split_radius, scenario.integration_step,
            scenario.max_integration_steps);
        if (result.incoming_trajectory.termination != TrajectoryTermination::reached_target_radius) {
            result.status = failed_integration_status(result.incoming_trajectory);
            return result;
        }

        const KerrFourMomentum incoming_coordinate_momentum =
            kerr_equatorial_four_momentum(incoming, split_radius);
        const LocalMomentum incoming_local = coordinate_to_zamo(
            mass, spin_length, split_radius, incoming_coordinate_momentum);
        const double incoming_mass_shell =
            std::abs(minkowski_dot(incoming_local, incoming_local) +
                     scenario.parent_rest_mass * scenario.parent_rest_mass);
        if (incoming_local.time <= 0.0 || incoming_mass_shell > scenario.residual_tolerance) {
            result.status = PenroseEventStatus::physics_invalid;
            result.mass_shell_residual = incoming_mass_shell;
            return result;
        }

        const LocalMomentum parent_unit_velocity =
            scale(incoming_local, 1.0 / scenario.parent_rest_mass);
        LocalMomentum radial_basis;
        LocalMomentum azimuth_basis;
        if (!build_equatorial_split_basis(parent_unit_velocity, &radial_basis, &azimuth_basis)) {
            result.status = PenroseEventStatus::physics_invalid;
            return result;
        }

        const double daughter_com_energy = scenario.parent_rest_mass / 2.0;
        const double daughter_com_momentum = std::sqrt(
            std::max(0.0, daughter_com_energy * daughter_com_energy -
                              scenario.fragment_rest_mass * scenario.fragment_rest_mass));
        const LocalMomentum direction = add(
            scale(radial_basis, std::cos(split.split_angle_rad)),
            scale(azimuth_basis, std::sin(split.split_angle_rad)));
        const LocalMomentum first_local = add(
            scale(parent_unit_velocity, daughter_com_energy),
            scale(direction, daughter_com_momentum));
        const LocalMomentum second_local = subtract(
            scale(parent_unit_velocity, daughter_com_energy),
            scale(direction, daughter_com_momentum));

        result.four_momentum_residual = max_component_abs(
            subtract(add(first_local, second_local), incoming_local));
        result.mass_shell_residual = std::max(
            std::abs(minkowski_dot(first_local, first_local) +
                     scenario.fragment_rest_mass * scenario.fragment_rest_mass),
            std::abs(minkowski_dot(second_local, second_local) +
                     scenario.fragment_rest_mass * scenario.fragment_rest_mass));
        if (first_local.time <= 0.0 || second_local.time <= 0.0 ||
            result.four_momentum_residual > scenario.residual_tolerance ||
            result.mass_shell_residual > scenario.residual_tolerance) {
            result.status = PenroseEventStatus::physics_invalid;
            return result;
        }

        const KerrFourMomentum first_coordinate =
            zamo_to_coordinate(mass, spin_length, split_radius, first_local);
        const KerrFourMomentum second_coordinate =
            zamo_to_coordinate(mass, spin_length, split_radius, second_local);
        const ConservedConstants first_constants =
            conserved_constants(mass, spin_length, split_radius, first_coordinate);
        const ConservedConstants second_constants =
            conserved_constants(mass, spin_length, split_radius, second_coordinate);
        result.energy_conservation_residual = std::abs(
            first_constants.energy + second_constants.energy - input_energy);
        result.angular_momentum_conservation_residual = std::abs(
            first_constants.angular_momentum + second_constants.angular_momentum -
            input_angular_momentum);
        const double coordinate_mass_shell = std::max(
            coordinate_mass_shell_residual(mass, spin_length, split_radius, first_coordinate,
                                           scenario.fragment_rest_mass),
            coordinate_mass_shell_residual(mass, spin_length, split_radius, second_coordinate,
                                           scenario.fragment_rest_mass));
        result.mass_shell_residual = std::max(result.mass_shell_residual, coordinate_mass_shell);
        if (result.energy_conservation_residual > scenario.residual_tolerance ||
            result.angular_momentum_conservation_residual > scenario.residual_tolerance ||
            result.mass_shell_residual > scenario.residual_tolerance) {
            result.status = PenroseEventStatus::physics_invalid;
            return result;
        }

        const KerrOrbit first_reconstructed_orbit{
            mass, spin_length, first_constants.energy, first_constants.angular_momentum,
            scenario.fragment_rest_mass, first_coordinate.radial < 0.0 ? -1 : 1};
        const KerrOrbit second_reconstructed_orbit{
            mass, spin_length, second_constants.energy, second_constants.angular_momentum,
            scenario.fragment_rest_mass, second_coordinate.radial < 0.0 ? -1 : 1};
        const KerrFourMomentum first_reconstructed =
            kerr_equatorial_four_momentum(first_reconstructed_orbit, split_radius);
        const KerrFourMomentum second_reconstructed =
            kerr_equatorial_four_momentum(second_reconstructed_orbit, split_radius);
        result.geodesic_initialization_residual = std::max(
            coordinate_momentum_difference(first_coordinate, first_reconstructed),
            coordinate_momentum_difference(second_coordinate, second_reconstructed));
        if (result.geodesic_initialization_residual > scenario.residual_tolerance) {
            result.status = PenroseEventStatus::physics_invalid;
            return result;
        }

        struct Fragment {
            KerrFourMomentum coordinate_momentum;
            ConservedConstants constants;
        };
        const std::array<Fragment, 2> fragments{{{first_coordinate, first_constants},
                                                  {second_coordinate, second_constants}}};
        for (const auto captured_index : {0U, 1U}) {
            const auto escaping_index = 1U - captured_index;
            const Fragment& captured = fragments[captured_index];
            const Fragment& escaping = fragments[escaping_index];
            if (captured.constants.energy >= -scenario.residual_tolerance ||
                escaping.constants.energy <= scenario.residual_tolerance ||
                captured.coordinate_momentum.radial >= -scenario.residual_tolerance ||
                escaping.coordinate_momentum.radial <= scenario.residual_tolerance) {
                continue;
            }

            const KerrOrbit captured_orbit{mass, spin_length, captured.constants.energy,
                                           captured.constants.angular_momentum,
                                           scenario.fragment_rest_mass, -1};
            const KerrOrbit escaping_orbit{mass, spin_length, escaping.constants.energy,
                                           escaping.constants.angular_momentum,
                                           scenario.fragment_rest_mass, 1};
            result.captured_trajectory = integrate_kerr(
                captured_orbit, split_radius, scenario.integration_step,
                scenario.max_integration_steps, escape_radius);
            result.escaping_trajectory = integrate_kerr(
                escaping_orbit, split_radius, scenario.integration_step,
                scenario.max_integration_steps, escape_radius);
            result.captured_energy = captured.constants.energy;
            result.escaping_energy = escaping.constants.energy;
            if (result.captured_trajectory.termination != TrajectoryTermination::crossed_horizon ||
                result.escaping_trajectory.termination !=
                    TrajectoryTermination::reached_escape_radius) {
                result.status =
                    result.captured_trajectory.termination == TrajectoryTermination::invalid_state ||
                            result.escaping_trajectory.termination ==
                                TrajectoryTermination::invalid_state
                        ? PenroseEventStatus::integration_failed
                        : PenroseEventStatus::captured_or_non_escaping;
                return result;
            }

            result.eta_penrose =
                (result.escaping_energy - result.input_energy) / result.input_energy;
            result.extracted_energy =
                std::max(0.0, result.escaping_energy - result.input_energy);
            result.status = PenroseEventStatus::physically_feasible;
            return result;
        }

        result.status = PenroseEventStatus::captured_or_non_escaping;
        return result;
    } catch (const std::invalid_argument&) {
        result.status = PenroseEventStatus::physics_invalid;
        return result;
    }
}
}  // namespace bh
