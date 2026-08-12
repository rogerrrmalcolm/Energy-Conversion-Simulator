#include "cli_json.hpp"

#include "bh/trajectory.hpp"

#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <ostream>
#include <string>
#include <string_view>

namespace bh::cli_json {
namespace {
class StreamState {
public:
    explicit StreamState(std::ostream& output)
        : output_(output), flags_(output.flags()), precision_(output.precision()) {
        output_ << std::defaultfloat
                << std::setprecision(std::numeric_limits<double>::max_digits10);
    }

    ~StreamState() {
        output_.flags(flags_);
        output_.precision(precision_);
    }

private:
    std::ostream& output_;
    std::ios::fmtflags flags_;
    std::streamsize precision_;
};

void write_string(std::ostream& output, const std::string_view value) {
    constexpr char hex[] = "0123456789abcdef";
    output << '"';
    for (const unsigned char character : value) {
        switch (character) {
        case '"':
            output << "\\\"";
            break;
        case '\\':
            output << "\\\\";
            break;
        case '\b':
            output << "\\b";
            break;
        case '\f':
            output << "\\f";
            break;
        case '\n':
            output << "\\n";
            break;
        case '\r':
            output << "\\r";
            break;
        case '\t':
            output << "\\t";
            break;
        default:
            if (character < 0x20) {
                output << "\\u00" << hex[(character >> 4) & 0x0f]
                       << hex[character & 0x0f];
            } else {
                output << static_cast<char>(character);
            }
        }
    }
    output << '"';
}

void write_number(std::ostream& output, const double value) {
    if (std::isfinite(value)) {
        output << value;
    } else {
        output << "null";
    }
}

void write_bool(std::ostream& output, const bool value) {
    output << (value ? "true" : "false");
}

std::string_view termination_name(const TrajectoryTermination termination) {
    switch (termination) {
    case TrajectoryTermination::completed:
        return "completed";
    case TrajectoryTermination::crossed_horizon:
        return "crossed_horizon";
    case TrajectoryTermination::reached_escape_radius:
        return "reached_escape_radius";
    case TrajectoryTermination::reached_target_radius:
        return "reached_target_radius";
    case TrajectoryTermination::turning_point:
        return "turning_point";
    case TrajectoryTermination::invalid_state:
        return "invalid_state";
    }
    return "unknown";
}

double elapsed_milliseconds(const PenroseDijkstraSearchDiagnostics& diagnostics) {
    return std::chrono::duration<double, std::milli>(diagnostics.elapsed).count();
}

double throughput(const PenroseDijkstraSearchDiagnostics& diagnostics) {
    const double seconds = std::chrono::duration<double>(diagnostics.elapsed).count();
    return seconds > 0.0 ? static_cast<double>(diagnostics.nodes_evaluated) / seconds : 0.0;
}

std::string_view phase_map_backend(const PenroseDijkstraSearchDiagnostics& diagnostics) {
    if (diagnostics.avx2_four_lane_batches > 0) {
        return "avx2-four-lane-single-thread";
    }
    if (diagnostics.four_lane_batches > 0) {
        return "portable-scalar-batch4-single-thread";
    }
    return "scalar-single-thread";
}

void write_split(std::ostream& output, const PenroseSplitParameters& split,
                 const std::string_view indentation) {
    output << "{\n" << indentation << "  \"split_radius_over_m\": ";
    write_number(output, split.split_radius_over_m);
    output << ",\n" << indentation << "  \"incoming_lz_over_m_m\": ";
    write_number(output, split.incoming_lz_over_m_m);
    output << ",\n" << indentation << "  \"split_angle_rad\": ";
    write_number(output, split.split_angle_rad);
    output << "\n" << indentation << "}";
}

void write_rotational_energy(std::ostream& output,
                             const RotationalEnergyResult& result,
                             const std::string_view indentation) {
    output << "{\n" << indentation << "  \"mass_energy_joules\": ";
    write_number(output, result.mass_energy_joules);
    output << ",\n" << indentation << "  \"irreducible_mass_kg\": ";
    write_number(output, result.irreducible_mass_kg);
    output << ",\n" << indentation << "  \"irreducible_mass_fraction\": ";
    write_number(output, result.irreducible_mass_fraction);
    output << ",\n" << indentation << "  \"rotational_energy_joules\": ";
    write_number(output, result.rotational_energy_joules);
    output << ",\n" << indentation << "  \"rotational_fraction\": ";
    write_number(output, result.rotational_fraction);
    output << ",\n" << indentation << "  \"spin_sensitivity_joules\": ";
    write_number(output, result.d_rotational_energy_d_spin_joules);
    output << "\n" << indentation << "}";
}

void write_integration_diagnostics(std::ostream& output, const Trajectory& trajectory,
                                   const std::string_view indentation) {
    output << "{\n" << indentation << "  \"termination\": ";
    write_string(output, termination_name(trajectory.termination));
    output << ",\n" << indentation << "  \"accepted_steps\": "
           << trajectory.diagnostics.accepted_steps
           << ",\n" << indentation << "  \"rejected_steps\": "
           << trajectory.diagnostics.rejected_steps
           << ",\n" << indentation << "  \"maximum_normalized_error\": ";
    write_number(output, trajectory.diagnostics.maximum_normalized_error);
    output << ",\n" << indentation
           << "  \"maximum_normalized_radial_residual\": ";
    write_number(output, trajectory.diagnostics.maximum_normalized_radial_residual);
    output << ",\n" << indentation << "  \"final_step\": ";
    write_number(output, trajectory.diagnostics.final_step);
    output << ",\n" << indentation << "  \"terminal_radius\": ";
    if (trajectory.points.empty()) {
        output << "null";
    } else {
        write_number(output, trajectory.points.back().radius);
    }
    output << "\n" << indentation << "}";
}

void write_event(std::ostream& output, const PenroseEventResult& result,
                 const std::string_view indentation) {
    output << "{\n" << indentation << "  \"status\": ";
    write_string(output, penrose_event_status_name(result.status));
    output << ",\n" << indentation << "  \"split\": ";
    write_split(output, result.split, std::string(indentation) + "  ");
    output << ",\n" << indentation << "  \"horizon_radius\": ";
    write_number(output, result.horizon_radius);
    output << ",\n" << indentation << "  \"static_limit_radius\": ";
    write_number(output, result.static_limit_radius);
    output << ",\n" << indentation << "  \"input_energy\": ";
    write_number(output, result.input_energy);
    output << ",\n" << indentation << "  \"captured_energy\": ";
    write_number(output, result.captured_energy);
    output << ",\n" << indentation << "  \"escaping_energy\": ";
    write_number(output, result.escaping_energy);
    output << ",\n" << indentation << "  \"net_extracted_energy\": ";
    write_number(output, result.extracted_energy);
    output << ",\n" << indentation << "  \"penrose_efficiency\": ";
    write_number(output, result.eta_penrose);
    output << ",\n" << indentation << "  \"maximum_normalized_residual\": ";
    write_number(output, result.maximum_normalized_residual);
    output << ",\n" << indentation << "  \"residuals\": {\n"
           << indentation << "    \"four_momentum\": ";
    write_number(output, result.four_momentum_residual);
    output << ",\n" << indentation << "    \"mass_shell\": ";
    write_number(output, result.mass_shell_residual);
    output << ",\n" << indentation << "    \"energy\": ";
    write_number(output, result.energy_conservation_residual);
    output << ",\n" << indentation << "    \"angular_momentum\": ";
    write_number(output, result.angular_momentum_conservation_residual);
    output << ",\n" << indentation << "    \"geodesic_initialization\": ";
    write_number(output, result.geodesic_initialization_residual);
    output << "\n" << indentation << "  },\n"
           << indentation << "  \"trajectories\": {\n"
           << indentation << "    \"incoming\": ";
    write_integration_diagnostics(output, result.incoming_trajectory,
                                  std::string(indentation) + "    ");
    output << ",\n" << indentation << "    \"captured\": ";
    write_integration_diagnostics(output, result.captured_trajectory,
                                  std::string(indentation) + "    ");
    output << ",\n" << indentation << "    \"escaping\": ";
    write_integration_diagnostics(output, result.escaping_trajectory,
                                  std::string(indentation) + "    ");
    output << "\n" << indentation << "  }\n" << indentation << "}";
}

void write_status_counts(std::ostream& output,
                         const PenroseDijkstraSearchDiagnostics& diagnostics,
                         const std::string_view indentation) {
    output << "{\n" << indentation << "  \"outside_ergosphere\": "
           << diagnostics.outside_ergosphere
           << ",\n" << indentation << "  \"physics_invalid\": "
           << diagnostics.physics_invalid
           << ",\n" << indentation << "  \"captured_or_non_escaping\": "
           << diagnostics.captured_or_non_escaping
           << ",\n" << indentation << "  \"escaping_without_target\": "
           << diagnostics.escaping_without_target
           << ",\n" << indentation << "  \"integration_failed\": "
           << diagnostics.integration_failed
           << ",\n" << indentation << "  \"goal_feasible\": "
           << diagnostics.goal_feasible << "\n" << indentation << "}";
}

void write_search_diagnostics(std::ostream& output,
                              const PenroseDijkstraSearchDiagnostics& diagnostics,
                              const std::string_view indentation) {
    output << "{\n" << indentation << "  \"candidates_in_domain\": "
           << diagnostics.candidates_in_domain
           << ",\n" << indentation << "  \"nodes_generated\": "
           << diagnostics.nodes_generated
           << ",\n" << indentation << "  \"nodes_evaluated\": "
           << diagnostics.nodes_evaluated
           << ",\n" << indentation << "  \"nodes_expanded\": "
           << diagnostics.nodes_expanded
           << ",\n" << indentation << "  \"four_lane_batches\": "
           << diagnostics.four_lane_batches
           << ",\n" << indentation << "  \"avx2_four_lane_batches\": "
           << diagnostics.avx2_four_lane_batches
           << ",\n" << indentation << "  \"four_lane_nodes\": "
           << diagnostics.four_lane_nodes
           << ",\n" << indentation << "  \"scalar_nodes\": "
           << diagnostics.scalar_nodes
           << ",\n" << indentation << "  \"duplicate_nodes_skipped\": "
           << diagnostics.duplicate_nodes_skipped
           << ",\n" << indentation << "  \"outside_search_domain_neighbors\": "
           << diagnostics.outside_search_domain_neighbors
           << ",\n" << indentation << "  \"final_verification_evaluations\": "
           << diagnostics.final_verification_evaluations
           << ",\n" << indentation << "  \"elapsed_ms\": ";
    write_number(output, elapsed_milliseconds(diagnostics));
    output << ",\n" << indentation << "  \"throughput_nodes_per_second\": ";
    write_number(output, throughput(diagnostics));
    output << ",\n" << indentation << "  \"status_counts\": ";
    write_status_counts(output, diagnostics, std::string(indentation) + "  ");
    output << "\n" << indentation << "}";
}

void write_candidate(std::ostream& output, const PenroseDijkstraNode& node,
                     const std::string_view indentation) {
    output << "{\n" << indentation << "  \"key\": [" << node.key[0] << ", "
           << node.key[1] << ", " << node.key[2] << "],\n"
           << indentation << "  \"split\": ";
    write_split(output, node.split, std::string(indentation) + "  ");
    output << ",\n" << indentation << "  \"status\": ";
    write_string(output, penrose_dijkstra_node_status_name(node.status));
    output << ",\n" << indentation << "  \"g_cost\": " << node.g_cost
           << ",\n" << indentation << "  \"h_cost\": " << node.h_cost
           << ",\n" << indentation << "  \"f_cost\": " << node.f_cost
           << ",\n" << indentation << "  \"penrose_efficiency\": ";
    write_number(output, node.eta_penrose);
    output << ",\n" << indentation << "  \"net_extracted_energy\": ";
    write_number(output, node.extracted_energy);
    output << ",\n" << indentation << "  \"maximum_normalized_residual\": ";
    write_number(output, node.maximum_normalized_residual);
    output << ",\n" << indentation << "  \"captured_termination\": ";
    write_string(output, termination_name(node.captured_termination));
    output << ",\n" << indentation << "  \"escaping_termination\": ";
    write_string(output, termination_name(node.escaping_termination));
    output << "\n" << indentation << "}";
}
}  // namespace

void write_version(std::ostream& output, const std::string_view version) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"version\",\n"
           << "  \"version\": ";
    write_string(output, version);
    output << "\n}\n";
}

void write_algebraic(std::ostream& output, const RotationalEnergyResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"algebraic\",\n"
           << "  \"result\": ";
    write_rotational_energy(output, result, "  ");
    output << "\n}\n";
}

void write_algebraic_range(std::ostream& output,
                           const RotationalEnergyRangeResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n"
           << "  \"command\": \"algebraic-range\",\n"
           << "  \"result\": {\n    \"lower\": ";
    write_rotational_energy(output, result.lower, "    ");
    output << ",\n    \"central\": ";
    write_rotational_energy(output, result.central, "    ");
    output << ",\n    \"upper\": ";
    write_rotational_energy(output, result.upper, "    ");
    output << ",\n    \"uncertainty_minus_joules\": ";
    write_number(output, result.rotational_energy_uncertainty_minus_joules);
    output << ",\n    \"uncertainty_plus_joules\": ";
    write_number(output, result.rotational_energy_uncertainty_plus_joules);
    output << "\n  }\n}\n";
}

void write_toy_plasma(std::ostream& output, const PlasmaInput& input,
                      const PlasmaResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"toy-plasma\",\n"
           << "  \"model\": \"reduced-toy-plasma\",\n  \"input\": {\n"
           << "    \"magnetic_field_tesla\": ";
    write_number(output, input.magnetic_field_tesla);
    output << ",\n    \"mass_density_kg_m3\": ";
    write_number(output, input.mass_density_kg_m3);
    output << ",\n    \"flow_area_m2\": ";
    write_number(output, input.flow_area_m2);
    output << ",\n    \"dimensionless_spin\": ";
    write_number(output, input.dimensionless_spin);
    output << ",\n    \"duration_seconds\": ";
    write_number(output, input.duration_seconds);
    output << "\n  },\n  \"result\": {\n    \"magnetization\": ";
    write_number(output, result.magnetization);
    output << ",\n    \"alfven_speed_m_s\": ";
    write_number(output, result.alfven_speed_m_s);
    output << ",\n    \"poynting_power_watts\": ";
    write_number(output, result.poynting_power_watts);
    output << ",\n    \"spin_coupling_efficiency\": ";
    write_number(output, result.spin_coupling_efficiency);
    output << ",\n    \"outward_electromagnetic_power_watts\": ";
    write_number(output, result.outward_electromagnetic_power_watts);
    output << ",\n    \"outward_electromagnetic_energy_joules\": ";
    write_number(output, result.outward_electromagnetic_energy_joules);
    output << "\n  }\n}\n";
}

void write_penrose_event(std::ostream& output,
                         const EquatorialPenroseScenario& scenario,
                         const PenroseEventResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"scenario\",\n"
           << "  \"model\": \"restricted-equatorial-penrose\",\n"
           << "  \"residual_tolerance\": ";
    write_number(output, scenario.residual_tolerance);
    output << ",\n  \"result\": ";
    write_event(output, result, "  ");
    output << "\n}\n";
}

void write_penrose_search(std::ostream& output,
                          const EquatorialPenroseDijkstraInput& input,
                          const PenroseDijkstraSearchResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"search\",\n"
           << "  \"algorithm\": ";
    write_string(output, penrose_search_algorithm_name(input.search.algorithm));
    output << ",\n  \"execution_backend\": \"scalar-single-thread\",\n"
           << "  \"target_efficiency\": ";
    write_number(output, input.search.eta_target);
    output << ",\n  \"status\": ";
    write_string(output, penrose_dijkstra_search_status_name(result.status));
    output << ",\n  \"found\": ";
    write_bool(output, result.found);
    output << ",\n  \"target_reached\": ";
    write_bool(output, result.target_reached);
    output << ",\n  \"failure_message\": ";
    write_string(output, result.failure_message);
    output << ",\n  \"diagnostics\": ";
    write_search_diagnostics(output, result.diagnostics, "  ");
    output << ",\n  \"selected_candidate\": ";
    if (result.parameter_adjustment_path.empty()) {
        output << "null";
    } else {
        write_candidate(output, result.parameter_adjustment_path.back(), "  ");
    }
    output << ",\n  \"selected_event\": ";
    if (!result.found) {
        output << "null";
    } else {
        write_event(output, result.selected_event, "  ");
    }
    output << ",\n  \"parameter_adjustment_path\": [";
    for (std::size_t index = 0; index < result.parameter_adjustment_path.size(); ++index) {
        output << (index == 0 ? "\n    " : ",\n    ");
        write_candidate(output, result.parameter_adjustment_path[index], "    ");
    }
    if (!result.parameter_adjustment_path.empty()) {
        output << '\n';
    }
    output << "  ]\n}\n";
}

void write_penrose_phase_map(std::ostream& output,
                             const EquatorialPenroseDijkstraInput& input,
                             const PenrosePhaseMapResult& result) {
    StreamState state(output);
    output << "{\n  \"schema_version\": 1,\n  \"command\": \"map\",\n"
           << "  \"model\": \"bounded-equatorial-penrose-phase-map\",\n"
           << "  \"execution_backend\": ";
    write_string(output, phase_map_backend(result.diagnostics));
    output << ",\n  \"target_efficiency\": ";
    write_number(output, input.search.eta_target);
    output << ",\n  \"status\": ";
    write_string(output, penrose_phase_map_status_name(result.status));
    output << ",\n  \"complete\": ";
    write_bool(output, result.complete);
    output << ",\n  \"failure_message\": ";
    write_string(output, result.failure_message);
    output << ",\n  \"diagnostics\": ";
    write_search_diagnostics(output, result.diagnostics, "  ");
    output << ",\n  \"best_candidate\": ";
    if (!result.best_validated_candidate) {
        output << "null";
    } else {
        write_candidate(output, *result.best_validated_candidate, "  ");
    }
    output << ",\n  \"best_event\": ";
    if (!result.best_validated_candidate || !result.complete) {
        output << "null";
    } else {
        write_event(output, result.best_event, "  ");
    }
    output << "\n}\n";
}
}  // namespace bh::cli_json
