#include "bh/penrose_model.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

#if defined(BH_ENABLE_AVX2)
#include <immintrin.h>
#endif

namespace bh {
namespace {
#if defined(BH_ENABLE_AVX2)
__m256d load4(const DoubleBatch4& values) {
    return _mm256_loadu_pd(values.data());
}

void store4(DoubleBatch4* output, const __m256d values) {
    _mm256_storeu_pd(output->data(), values);
}

__m256d abs4(const __m256d values) {
    return _mm256_andnot_pd(_mm256_set1_pd(-0.0), values);
}

void split_component_batch4(const DoubleBatch4& parent, const DoubleBatch4& radial_basis,
                            const DoubleBatch4& azimuth_basis, const __m256d cosines,
                            const __m256d sines, const __m256d daughter_energy,
                            const __m256d daughter_momentum, DoubleBatch4* first,
                            DoubleBatch4* second) {
    const __m256d direction = _mm256_add_pd(
        _mm256_mul_pd(load4(radial_basis), cosines),
        _mm256_mul_pd(load4(azimuth_basis), sines));
    const __m256d center = _mm256_mul_pd(load4(parent), daughter_energy);
    const __m256d offset = _mm256_mul_pd(direction, daughter_momentum);
    store4(first, _mm256_add_pd(center, offset));
    store4(second, _mm256_sub_pd(center, offset));
}
#endif

bool finite_batch4(const DoubleBatch4& values) {
    return std::all_of(values.begin(), values.end(),
                       [](const double value) { return std::isfinite(value); });
}

bool finite_local_batch4(const PenroseLocalMomentumBatch4& momenta) {
    return finite_batch4(momenta.time) && finite_batch4(momenta.radial) &&
           finite_batch4(momenta.azimuth);
}

bool finite_coordinate_batch4(const KerrFourMomentumBatch4& momenta) {
    return finite_batch4(momenta.coordinate_time) && finite_batch4(momenta.radial) &&
           finite_batch4(momenta.azimuth);
}

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

#if defined(BH_ENABLE_AVX2)
struct ZamoMetricBatch4 {
    DoubleBatch4 lapse{};
    DoubleBatch4 sqrt_delta{};
    DoubleBatch4 sqrt_g_phiphi{};
    DoubleBatch4 frame_dragging_omega{};
};
#endif

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
    if (!std::isfinite(value.time) || !std::isfinite(value.radial) ||
        !std::isfinite(value.azimuth)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::max({std::abs(value.time), std::abs(value.radial), std::abs(value.azimuth)});
}

double max_component_abs(const KerrFourMomentum& value) {
    if (!std::isfinite(value.coordinate_time) || !std::isfinite(value.radial) ||
        !std::isfinite(value.azimuth)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::max(
        {std::abs(value.coordinate_time), std::abs(value.radial), std::abs(value.azimuth)});
}

double normalized_residual(const double residual, const double scale) {
    if (!std::isfinite(residual) || !std::isfinite(scale)) {
        return std::numeric_limits<double>::infinity();
    }
    return std::abs(residual) / std::max(1.0, std::abs(scale));
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
        !std::isfinite(scenario.integration_step) ||
        !std::isfinite(scenario.integration_control.absolute_tolerance) ||
        !std::isfinite(scenario.integration_control.relative_tolerance) ||
        !std::isfinite(scenario.integration_control.minimum_step) ||
        !std::isfinite(scenario.residual_tolerance) ||
        scenario.black_hole_mass <= 0.0 ||
        scenario.dimensionless_spin < 0.0 || scenario.dimensionless_spin >= 1.0 ||
        scenario.parent_rest_mass <= 0.0 || scenario.fragment_rest_mass < 0.0 ||
        scenario.fragment_rest_mass > scenario.parent_rest_mass / 2.0 ||
        scenario.incoming_specific_energy <= 0.0 || scenario.initial_radius_over_m <= 0.0 ||
        scenario.escape_radius_over_m <= scenario.initial_radius_over_m ||
        scenario.integration_step <= 0.0 || scenario.max_integration_steps == 0 ||
        scenario.integration_control.absolute_tolerance <= 0.0 ||
        scenario.integration_control.relative_tolerance <= 0.0 ||
        scenario.integration_control.minimum_step < 0.0 ||
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

#if defined(BH_ENABLE_AVX2)
ZamoMetricBatch4 zamo_metric_batch4(const PenroseZamoGeometryBatch4& geometry) {
    ZamoMetricBatch4 result;
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const EquatorialMetric metric = equatorial_metric(
            geometry.black_hole_masses[lane], geometry.spin_lengths[lane],
            geometry.radii[lane]);
        result.lapse[lane] = metric.lapse;
        result.sqrt_delta[lane] = std::sqrt(metric.delta);
        result.sqrt_g_phiphi[lane] = std::sqrt(metric.g_phiphi);
        result.frame_dragging_omega[lane] = metric.frame_dragging_omega;
    }
    return result;
}
#endif

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

PenroseEnergyBatch4Result penrose_energy_extraction_batch4(
    const PenroseEnergyBatch4Input& input) {
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        if (!std::isfinite(input.input_energies[lane]) ||
            !std::isfinite(input.escaping_energies[lane]) ||
            input.input_energies[lane] <= 0.0) {
            throw std::invalid_argument(
                "Penrose energy batch requires finite escaping energy and positive input energy");
        }
    }

    PenroseEnergyBatch4Result result;
#if defined(BH_ENABLE_AVX2)
    const __m256d input_energy = _mm256_loadu_pd(input.input_energies.data());
    const __m256d escaping_energy = _mm256_loadu_pd(input.escaping_energies.data());
    const __m256d difference = _mm256_sub_pd(escaping_energy, input_energy);
    const __m256d eta = _mm256_div_pd(difference, input_energy);
    const __m256d extracted = _mm256_max_pd(difference, _mm256_setzero_pd());
    _mm256_storeu_pd(result.eta_penrose.data(), eta);
    _mm256_storeu_pd(result.extracted_energies.data(), extracted);
#else
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const double difference =
            input.escaping_energies[lane] - input.input_energies[lane];
        result.eta_penrose[lane] = difference / input.input_energies[lane];
        result.extracted_energies[lane] = std::max(0.0, difference);
    }
#endif
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        if (!std::isfinite(result.eta_penrose[lane]) ||
            !std::isfinite(result.extracted_energies[lane])) {
            throw std::overflow_error("Penrose energy batch overflowed");
        }
    }
    return result;
}

PenroseLocalMomentumBatch4 coordinate_to_zamo_batch4(
    const PenroseZamoGeometryBatch4& geometry,
    const KerrFourMomentumBatch4& coordinate_momenta) {
    if (!finite_coordinate_batch4(coordinate_momenta)) {
        throw std::invalid_argument("ZAMO coordinate momenta must be finite");
    }

    PenroseLocalMomentumBatch4 result;
#if defined(BH_ENABLE_AVX2)
    const ZamoMetricBatch4 metric = zamo_metric_batch4(geometry);
    const __m256d coordinate_time = load4(coordinate_momenta.coordinate_time);
    store4(&result.time, _mm256_mul_pd(load4(metric.lapse), coordinate_time));
    store4(&result.radial,
           _mm256_div_pd(
               _mm256_mul_pd(load4(geometry.radii), load4(coordinate_momenta.radial)),
               load4(metric.sqrt_delta)));
    store4(&result.azimuth,
           _mm256_mul_pd(
               load4(metric.sqrt_g_phiphi),
               _mm256_sub_pd(load4(coordinate_momenta.azimuth),
                             _mm256_mul_pd(load4(metric.frame_dragging_omega),
                                           coordinate_time))));
#else
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const LocalMomentum local = coordinate_to_zamo(
            geometry.black_hole_masses[lane], geometry.spin_lengths[lane],
            geometry.radii[lane],
            {coordinate_momenta.coordinate_time[lane], coordinate_momenta.radial[lane],
             coordinate_momenta.azimuth[lane]});
        result.time[lane] = local.time;
        result.radial[lane] = local.radial;
        result.azimuth[lane] = local.azimuth;
    }
#endif
    if (!finite_local_batch4(result)) {
        throw std::overflow_error("coordinate-to-ZAMO batch overflowed");
    }
    return result;
}

KerrFourMomentumBatch4 zamo_to_coordinate_batch4(
    const PenroseZamoGeometryBatch4& geometry,
    const PenroseLocalMomentumBatch4& local_momenta) {
    if (!finite_local_batch4(local_momenta)) {
        throw std::invalid_argument("ZAMO local momenta must be finite");
    }

    KerrFourMomentumBatch4 result;
#if defined(BH_ENABLE_AVX2)
    const ZamoMetricBatch4 metric = zamo_metric_batch4(geometry);
    const __m256d coordinate_time =
        _mm256_div_pd(load4(local_momenta.time), load4(metric.lapse));
    store4(&result.coordinate_time, coordinate_time);
    store4(&result.radial,
           _mm256_div_pd(
               _mm256_mul_pd(load4(local_momenta.radial), load4(metric.sqrt_delta)),
               load4(geometry.radii)));
    store4(&result.azimuth,
           _mm256_add_pd(
               _mm256_mul_pd(load4(metric.frame_dragging_omega), coordinate_time),
               _mm256_div_pd(load4(local_momenta.azimuth),
                             load4(metric.sqrt_g_phiphi))));
#else
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const KerrFourMomentum coordinate = zamo_to_coordinate(
            geometry.black_hole_masses[lane], geometry.spin_lengths[lane],
            geometry.radii[lane],
            {local_momenta.time[lane], local_momenta.radial[lane],
             local_momenta.azimuth[lane]});
        result.coordinate_time[lane] = coordinate.coordinate_time;
        result.radial[lane] = coordinate.radial;
        result.azimuth[lane] = coordinate.azimuth;
    }
#endif
    if (!finite_coordinate_batch4(result)) {
        throw std::overflow_error("ZAMO-to-coordinate batch overflowed");
    }
    return result;
}

PenroseFragmentSplitBatch4Result split_penrose_fragments_batch4(
    const PenroseFragmentSplitBatch4Input& input) {
    if (!finite_local_batch4(input.parent_unit_velocities) ||
        !finite_local_batch4(input.radial_bases) ||
        !finite_local_batch4(input.azimuth_bases) ||
        !finite_batch4(input.split_angles_rad) ||
        !finite_batch4(input.daughter_com_energies) ||
        !finite_batch4(input.daughter_com_momenta)) {
        throw std::invalid_argument("Penrose split batch inputs must be finite");
    }
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        if (input.daughter_com_energies[lane] <= 0.0 ||
            input.daughter_com_momenta[lane] < 0.0) {
            throw std::invalid_argument("Penrose daughter energies and momenta must be physical");
        }
    }

    PenroseFragmentSplitBatch4Result result;
#if defined(BH_ENABLE_AVX2)
    DoubleBatch4 cosines{};
    DoubleBatch4 sines{};
    // AVX2 has no native trigonometric instructions; only these four lookups stay scalar.
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        cosines[lane] = std::cos(input.split_angles_rad[lane]);
        sines[lane] = std::sin(input.split_angles_rad[lane]);
    }
    const __m256d cosine = load4(cosines);
    const __m256d sine = load4(sines);
    const __m256d daughter_energy = load4(input.daughter_com_energies);
    const __m256d daughter_momentum = load4(input.daughter_com_momenta);
    split_component_batch4(input.parent_unit_velocities.time, input.radial_bases.time,
                           input.azimuth_bases.time, cosine, sine, daughter_energy,
                           daughter_momentum, &result.first.time, &result.second.time);
    split_component_batch4(input.parent_unit_velocities.radial, input.radial_bases.radial,
                           input.azimuth_bases.radial, cosine, sine, daughter_energy,
                           daughter_momentum, &result.first.radial, &result.second.radial);
    split_component_batch4(input.parent_unit_velocities.azimuth,
                           input.radial_bases.azimuth, input.azimuth_bases.azimuth,
                           cosine, sine, daughter_energy, daughter_momentum,
                           &result.first.azimuth, &result.second.azimuth);
#else
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const LocalMomentum parent{input.parent_unit_velocities.time[lane],
                                   input.parent_unit_velocities.radial[lane],
                                   input.parent_unit_velocities.azimuth[lane]};
        const LocalMomentum radial_basis{input.radial_bases.time[lane],
                                         input.radial_bases.radial[lane],
                                         input.radial_bases.azimuth[lane]};
        const LocalMomentum azimuth_basis{input.azimuth_bases.time[lane],
                                          input.azimuth_bases.radial[lane],
                                          input.azimuth_bases.azimuth[lane]};
        const LocalMomentum direction = add(
            scale(radial_basis, std::cos(input.split_angles_rad[lane])),
            scale(azimuth_basis, std::sin(input.split_angles_rad[lane])));
        const LocalMomentum first =
            add(scale(parent, input.daughter_com_energies[lane]),
                scale(direction, input.daughter_com_momenta[lane]));
        const LocalMomentum second =
            subtract(scale(parent, input.daughter_com_energies[lane]),
                     scale(direction, input.daughter_com_momenta[lane]));
        result.first.time[lane] = first.time;
        result.first.radial[lane] = first.radial;
        result.first.azimuth[lane] = first.azimuth;
        result.second.time[lane] = second.time;
        result.second.radial[lane] = second.radial;
        result.second.azimuth[lane] = second.azimuth;
    }
#endif
    if (!finite_local_batch4(result.first) || !finite_local_batch4(result.second)) {
        throw std::overflow_error("Penrose fragment split batch overflowed");
    }
    return result;
}

PenroseConservationBatch4Result penrose_conservation_residuals_batch4(
    const PenroseConservationBatch4Input& input) {
    if (!finite_local_batch4(input.parent) || !finite_local_batch4(input.first) ||
        !finite_local_batch4(input.second) ||
        !finite_batch4(input.fragment_rest_masses) ||
        !finite_batch4(input.incoming_constants.energies) ||
        !finite_batch4(input.incoming_constants.angular_momenta) ||
        !finite_batch4(input.first_constants.energies) ||
        !finite_batch4(input.first_constants.angular_momenta) ||
        !finite_batch4(input.second_constants.energies) ||
        !finite_batch4(input.second_constants.angular_momenta)) {
        throw std::invalid_argument("Penrose conservation batch inputs must be finite");
    }
    if (std::any_of(input.fragment_rest_masses.begin(), input.fragment_rest_masses.end(),
                    [](const double mass) { return mass < 0.0; })) {
        throw std::invalid_argument("Penrose fragment rest masses must be non-negative");
    }

    PenroseConservationBatch4Result result;
#if defined(BH_ENABLE_AVX2)
    const __m256d first_time = load4(input.first.time);
    const __m256d first_radial = load4(input.first.radial);
    const __m256d first_azimuth = load4(input.first.azimuth);
    const __m256d second_time = load4(input.second.time);
    const __m256d second_radial = load4(input.second.radial);
    const __m256d second_azimuth = load4(input.second.azimuth);
    const __m256d time_residual = abs4(_mm256_sub_pd(
        _mm256_add_pd(first_time, second_time), load4(input.parent.time)));
    const __m256d radial_residual = abs4(_mm256_sub_pd(
        _mm256_add_pd(first_radial, second_radial), load4(input.parent.radial)));
    const __m256d azimuth_residual = abs4(_mm256_sub_pd(
        _mm256_add_pd(first_azimuth, second_azimuth), load4(input.parent.azimuth)));
    store4(&result.four_momentum_residuals,
           _mm256_max_pd(time_residual,
                         _mm256_max_pd(radial_residual, azimuth_residual)));

    const __m256d rest_mass = load4(input.fragment_rest_masses);
    const __m256d rest_mass_squared = _mm256_mul_pd(rest_mass, rest_mass);
    const __m256d first_norm = _mm256_add_pd(
        _mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(first_radial, first_radial),
                                    _mm256_mul_pd(first_azimuth, first_azimuth)),
                      _mm256_mul_pd(first_time, first_time)),
        rest_mass_squared);
    const __m256d second_norm = _mm256_add_pd(
        _mm256_sub_pd(_mm256_add_pd(_mm256_mul_pd(second_radial, second_radial),
                                    _mm256_mul_pd(second_azimuth, second_azimuth)),
                      _mm256_mul_pd(second_time, second_time)),
        rest_mass_squared);
    store4(&result.mass_shell_residuals,
           _mm256_max_pd(abs4(first_norm), abs4(second_norm)));
    store4(&result.energy_residuals,
           abs4(_mm256_sub_pd(
               _mm256_add_pd(load4(input.first_constants.energies),
                             load4(input.second_constants.energies)),
               load4(input.incoming_constants.energies))));
    store4(&result.angular_momentum_residuals,
           abs4(_mm256_sub_pd(
               _mm256_add_pd(load4(input.first_constants.angular_momenta),
                             load4(input.second_constants.angular_momenta)),
               load4(input.incoming_constants.angular_momenta))));
#else
    for (std::size_t lane = 0; lane < avx2_double_lanes; ++lane) {
        const LocalMomentum parent{input.parent.time[lane], input.parent.radial[lane],
                                   input.parent.azimuth[lane]};
        const LocalMomentum first{input.first.time[lane], input.first.radial[lane],
                                  input.first.azimuth[lane]};
        const LocalMomentum second{input.second.time[lane], input.second.radial[lane],
                                   input.second.azimuth[lane]};
        result.four_momentum_residuals[lane] =
            max_component_abs(subtract(add(first, second), parent));
        const double mass_squared = input.fragment_rest_masses[lane] *
                                    input.fragment_rest_masses[lane];
        result.mass_shell_residuals[lane] =
            std::max(std::abs(minkowski_dot(first, first) + mass_squared),
                     std::abs(minkowski_dot(second, second) + mass_squared));
        result.energy_residuals[lane] = std::abs(
            input.first_constants.energies[lane] +
            input.second_constants.energies[lane] -
            input.incoming_constants.energies[lane]);
        result.angular_momentum_residuals[lane] = std::abs(
            input.first_constants.angular_momenta[lane] +
            input.second_constants.angular_momenta[lane] -
            input.incoming_constants.angular_momenta[lane]);
    }
#endif
    if (!finite_batch4(result.four_momentum_residuals) ||
        !finite_batch4(result.mass_shell_residuals) ||
        !finite_batch4(result.energy_residuals) ||
        !finite_batch4(result.angular_momentum_residuals)) {
        throw std::overflow_error("Penrose conservation residual batch overflowed");
    }
    return result;
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
            scenario.max_integration_steps, scenario.integration_control);
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
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(incoming_mass_shell,
                                scenario.parent_rest_mass * scenario.parent_rest_mass));
        if (incoming_local.time <= 0.0 ||
            result.maximum_normalized_residual > scenario.residual_tolerance) {
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
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(result.four_momentum_residual, max_component_abs(incoming_local)));
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(result.mass_shell_residual,
                                scenario.parent_rest_mass * scenario.parent_rest_mass));
        if (first_local.time <= 0.0 || second_local.time <= 0.0 ||
            result.maximum_normalized_residual > scenario.residual_tolerance) {
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
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(result.energy_conservation_residual, input_energy));
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(result.angular_momentum_conservation_residual,
                                input_angular_momentum));
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(coordinate_mass_shell,
                                scenario.fragment_rest_mass * scenario.fragment_rest_mass));
        if (result.maximum_normalized_residual > scenario.residual_tolerance) {
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
        result.maximum_normalized_residual = std::max(
            result.maximum_normalized_residual,
            normalized_residual(
                result.geodesic_initialization_residual,
                std::max(max_component_abs(first_coordinate), max_component_abs(second_coordinate))));
        if (result.maximum_normalized_residual > scenario.residual_tolerance) {
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
                scenario.max_integration_steps, escape_radius, scenario.integration_control);
            result.escaping_trajectory = integrate_kerr(
                escaping_orbit, split_radius, scenario.integration_step,
                scenario.max_integration_steps, escape_radius, scenario.integration_control);
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

            const double maximum_trajectory_radial_residual = std::max(
                {result.incoming_trajectory.diagnostics.maximum_normalized_radial_residual,
                 result.captured_trajectory.diagnostics.maximum_normalized_radial_residual,
                 result.escaping_trajectory.diagnostics.maximum_normalized_radial_residual});
            result.maximum_normalized_residual = std::max(
                result.maximum_normalized_residual, maximum_trajectory_radial_residual);
            if (!std::isfinite(maximum_trajectory_radial_residual)) {
                result.status = PenroseEventStatus::integration_failed;
                return result;
            }
            if (result.maximum_normalized_residual > scenario.residual_tolerance) {
                result.status = PenroseEventStatus::physics_invalid;
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
    } catch (const std::overflow_error&) {
        result.status = PenroseEventStatus::physics_invalid;
        return result;
    }
}
}  // namespace bh
