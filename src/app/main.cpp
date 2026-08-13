#include "bh/algebraic_model.hpp"
#include "bh/dijkstra.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"
#include "bh/trajectory.hpp"
#include "cli_json.hpp"

#include <algorithm>
#include <charconv>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <optional>
#include <span>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if defined(_WIN32)
#include <io.h>
#else
#include <unistd.h>
#endif

#ifndef BH_VERSION
#define BH_VERSION "development"
#endif

namespace {
enum class OutputFormat {
    text,
    json
};

enum class ProgressMode {
    automatic,
    always,
    never
};

struct CliOptions {
    OutputFormat format{OutputFormat::text};
    ProgressMode progress{ProgressMode::automatic};
    bool verbose{};
};

struct ParsedCli {
    std::string_view command{"interactive"};
    std::vector<std::string_view> arguments;
    CliOptions options{};
};

OutputFormat parse_output_format(const std::string_view value) {
    if (value == "text") {
        return OutputFormat::text;
    }
    if (value == "json") {
        return OutputFormat::json;
    }
    throw std::invalid_argument("--format must be text or json");
}

std::string_view normalize_command(const std::string_view command) {
    if (command == "--algebraic") {
        return "algebraic";
    }
    if (command == "--algebraic-range") {
        return "algebraic-range";
    }
    if (command == "--toy-plasma") {
        return "toy-plasma";
    }
    if (command == "--scenario") {
        return "scenario";
    }
    if (command == "--search-penrose") {
        return "search";
    }
    if (command == "--map-penrose") {
        return "map";
    }
    if (command == "--interactive") {
        return "interactive";
    }
    if (command == "--version") {
        return "version";
    }
    if (command == "--help" || command == "-h") {
        return "help";
    }
    return command;
}

ParsedCli parse_cli(const int argc, char* argv[]) {
    ParsedCli parsed;
    bool has_command = false;
    for (int index = 1; index < argc; ++index) {
        const std::string_view token = argv[index];
        if (token == "--verbose") {
            parsed.options.verbose = true;
            continue;
        }
        if (token == "--summary") {
            parsed.options.verbose = false;
            continue;
        }
        if (token == "--progress") {
            parsed.options.progress = ProgressMode::always;
            continue;
        }
        if (token == "--no-progress") {
            parsed.options.progress = ProgressMode::never;
            continue;
        }
        if (token == "--format") {
            if (++index >= argc) {
                throw std::invalid_argument("--format requires text or json");
            }
            parsed.options.format = parse_output_format(argv[index]);
            continue;
        }
        if (token.starts_with("--format=")) {
            parsed.options.format =
                parse_output_format(token.substr(std::string_view("--format=").size()));
            continue;
        }

        if (!has_command) {
            parsed.command = normalize_command(token);
            has_command = true;
        } else {
            parsed.arguments.push_back(token);
        }
    }
    return parsed;
}

double evaluated_nodes_per_second(const bh::PenroseDijkstraSearchDiagnostics& diagnostics) {
    const double seconds = std::chrono::duration<double>(diagnostics.elapsed).count();
    return seconds > 0.0 ? static_cast<double>(diagnostics.nodes_evaluated) / seconds : 0.0;
}

bool standard_error_is_terminal() {
#if defined(_WIN32)
    return _isatty(_fileno(stderr)) != 0;
#else
    return isatty(fileno(stderr)) != 0;
#endif
}

bool running_in_ci() {
    const char* const value = std::getenv("CI");
    return value != nullptr && std::string_view(value) != "" &&
           std::string_view(value) != "false";
}

bool should_report_progress(const CliOptions& options) {
    if (options.progress == ProgressMode::always) {
        return true;
    }
    if (options.progress == ProgressMode::never) {
        return false;
    }
    return options.format == OutputFormat::text && !running_in_ci() &&
           standard_error_is_terminal();
}

class TerminalProgress {
public:
    void update(const bh::PenroseSearchProgress& progress) {
        if (last_nodes_ && *last_nodes_ == progress.nodes_evaluated) {
            return;
        }
        last_nodes_ = progress.nodes_evaluated;
        const double percent = progress.candidates_in_domain == 0
                                   ? 0.0
                                   : 100.0 * static_cast<double>(progress.nodes_evaluated) /
                                         static_cast<double>(progress.candidates_in_domain);
        const double seconds = std::chrono::duration<double>(progress.elapsed).count();
        const double nodes_per_second =
            seconds > 0.0 ? static_cast<double>(progress.nodes_evaluated) / seconds : 0.0;
        constexpr std::size_t progress_bar_width = 24;
        const std::size_t completed_bar_cells = progress.candidates_in_domain == 0
                                                    ? 0
                                                    : std::min(
                                                          progress_bar_width,
                                                          progress.nodes_evaluated *
                                                              progress_bar_width /
                                                              progress.candidates_in_domain);

        std::ostringstream line;
        line << "Progress: " << progress.nodes_evaluated << "/"
             << progress.candidates_in_domain << " (" << std::fixed << std::setprecision(1)
             << percent << "%) [" << std::string(completed_bar_cells, '#')
             << std::string(progress_bar_width - completed_bar_cells, '-') << "] | "
             << std::setprecision(2) << nodes_per_second << " nodes/s"
             << " | best ";
        if (progress.best_eta_penrose) {
            line << std::setprecision(4) << 100.0 * *progress.best_eta_penrose << "%";
        } else {
            line << "n/a";
        }

        const std::string text = line.str();
        std::cerr << '\r' << text;
        if (text.size() < previous_width_) {
            std::cerr << std::string(previous_width_ - text.size(), ' ');
        }
        std::cerr << std::flush;
        previous_width_ = text.size();
        active_ = true;
    }

    void finish() {
        if (active_) {
            std::cerr << '\n';
            active_ = false;
            previous_width_ = 0;
            last_nodes_.reset();
        }
    }

    ~TerminalProgress() {
        finish();
    }

private:
    std::size_t previous_width_{};
    bool active_{};
    std::optional<std::size_t> last_nodes_{};
};

bh::PenroseSearchProgressObserver make_progress_observer(
    const std::size_t candidates, const CliOptions& options, TerminalProgress& display) {
    if (!should_report_progress(options)) {
        return {};
    }
    bh::PenroseSearchProgressObserver observer;
    observer.report_every_nodes = std::max<std::size_t>(10, candidates / 100);
    observer.callback = [&display](const bh::PenroseSearchProgress& progress) {
        display.update(progress);
    };
    return observer;
}

double parse_double(const std::string_view text, const std::string_view argument_name) {
    double value{};
    const auto [end, error] = std::from_chars(
        text.data(), text.data() + text.size(), value, std::chars_format::general);
    if (error != std::errc{} || end != text.data() + text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("expected a finite number for " + std::string(argument_name));
    }
    return value;
}

std::string trim_interactive_input(std::string text) {
    const std::size_t first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) {
        return {};
    }
    const std::size_t last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

std::string read_interactive_line() {
    std::string text;
    if (!std::getline(std::cin, text)) {
        throw std::runtime_error("interactive input stream closed");
    }
    return trim_interactive_input(text);
}

double prompt_double(const std::string_view label, const double default_value) {
    while (true) {
        std::cout << "  " << label << " [" << default_value << "]: " << std::flush;
        const std::string text = read_interactive_line();
        if (text.empty()) {
            return default_value;
        }
        try {
            return parse_double(text, label);
        } catch (const std::invalid_argument& error) {
            std::cout << "  " << error.what() << ". Try again.\n";
        }
    }
}

std::size_t prompt_step_count(const std::string_view label, const std::size_t default_value) {
    while (true) {
        const double value = prompt_double(label, static_cast<double>(default_value));
        if (std::floor(value) == value && value > 0.0 &&
            value <= static_cast<double>(std::numeric_limits<std::size_t>::max())) {
            return static_cast<std::size_t>(value);
        }
        std::cout << "  Enter a positive whole number.\n";
    }
}

std::string_view termination_name(const bh::TrajectoryTermination termination) {
    switch (termination) {
    case bh::TrajectoryTermination::completed:
        return "completed";
    case bh::TrajectoryTermination::crossed_horizon:
        return "crossed_horizon";
    case bh::TrajectoryTermination::reached_escape_radius:
        return "reached_escape_radius";
    case bh::TrajectoryTermination::reached_target_radius:
        return "reached_target_radius";
    case bh::TrajectoryTermination::turning_point:
        return "turning_point";
    case bh::TrajectoryTermination::invalid_state:
        return "invalid_state";
    }
    return "unknown";
}

void print_integration_diagnostics(const std::string_view label,
                                   const bh::Trajectory& trajectory) {
    std::cout << "  " << label << " accepted/rejected: "
              << trajectory.diagnostics.accepted_steps << "/"
              << trajectory.diagnostics.rejected_steps << "\n"
              << "  " << label << " max error: "
              << trajectory.diagnostics.maximum_normalized_error << "\n"
              << "  " << label << " max radial residual: "
              << trajectory.diagnostics.maximum_normalized_radial_residual << "\n"
              << "  " << label << " final step: "
              << trajectory.diagnostics.final_step << "\n";
    if (!trajectory.points.empty()) {
        std::cout << "  " << label << " terminal radius: "
                  << trajectory.points.back().radius << " M\n";
    }
}

void print_usage() {
    std::cout << "Usage:\n"
              << "  black-hole-sim\n"
              << "  black-hole-sim algebraic <mass_kg> <a_star>\n"
              << "  black-hole-sim algebraic-range <mass_low_kg> <mass_mid_kg> <mass_high_kg>"
                 " <spin_low> <spin_mid> <spin_high>\n"
              << "  black-hole-sim toy-plasma <magnetic_field_tesla> <mass_density_kg_m3>"
                 " <flow_area_m2>"
                 " <a_star> <duration_seconds>\n"
              << "  black-hole-sim scenario <event.cfg>\n"
              << "  black-hole-sim search <search.cfg>\n"
              << "  black-hole-sim map <search.cfg>\n"
              << "  black-hole-sim interactive\n"
              << "  black-hole-sim version\n\n"
        << "Output options:\n"
        << "  --summary          Print the concise result summary (default)\n"
        << "  --verbose          Include full physics and integration diagnostics\n"
        << "  --format text|json Select human-readable or structured output\n"
        << "  --progress         Show search progress even when output is redirected\n"
        << "  --no-progress      Disable search progress\n\n"
              << "Running without arguments opens the sequential guided pipeline.\n"
              << "Legacy --algebraic, --scenario, --search-penrose, and --map-penrose flags"
                 " remain supported.\n"
              << "Search runs Dijkstra over at most 25000 parameter nodes. Map evaluates every"
                 " node and reports the greatest validated extraction in that bounded grid.\n"
              << "The toy-plasma command is a reduced, ideal-MHD-inspired transport scaling,"
                 " not an MHD or GRMHD simulation.\n"
              << "The interactive command collects inputs in physics order, then runs the"
                  " bounded exhaustive phase map.\n";
}

void print_version() {
    std::cout << "black-hole-sim " << BH_VERSION << '\n';
}

void print_algebraic_result(const bh::RotationalEnergyResult& result,
                            const bool verbose) {
    std::cout << "Algebraic Kerr rotational-energy reservoir\n"
              << "  total mass-energy:       " << result.mass_energy_joules << " J\n"
              << "  irreducible mass:        " << result.irreducible_mass_kg << " kg\n"
              << "  rotational reservoir:    " << result.rotational_energy_joules << " J\n"
              << "  rotational fraction:     " << 100.0 * result.rotational_fraction << " %\n";
    if (verbose) {
        std::cout << "  irreducible fraction:    "
                  << 100.0 * result.irreducible_mass_fraction << " %\n"
                  << "  dE_rot / da_star:        "
                  << result.d_rotational_energy_d_spin_joules << " J\n";
    }
}

void print_algebraic_range_result(const bh::RotationalEnergyRangeResult& result,
                                  const bool verbose) {
    std::cout << "Algebraic Kerr rotational-energy range\n";
    std::cout << "  lower reservoir:         " << result.lower.rotational_energy_joules << " J\n";
    std::cout << "  central reservoir:       " << result.central.rotational_energy_joules
              << " J\n";
    std::cout << "  upper reservoir:         " << result.upper.rotational_energy_joules << " J\n";
    std::cout << "  lower uncertainty:       -" << result.rotational_energy_uncertainty_minus_joules
              << " J\n";
    std::cout << "  upper uncertainty:       +" << result.rotational_energy_uncertainty_plus_joules
              << " J\n";
    if (verbose) {
        std::cout << "  lower rotational fraction:   "
                  << 100.0 * result.lower.rotational_fraction << " %\n"
                  << "  central rotational fraction: "
                  << 100.0 * result.central.rotational_fraction << " %\n"
                  << "  upper rotational fraction:   "
                  << 100.0 * result.upper.rotational_fraction << " %\n";
    }
}

void print_toy_plasma_result(const bh::PlasmaResult& result,
                             const bool verbose) {
    if (!verbose) {
        std::cout << "Reduced toy-plasma transport estimate\n"
                  << "  magnetization sigma:           " << result.magnetization << "\n"
                  << "  outward electromagnetic power: "
                  << result.outward_electromagnetic_power_watts << " W\n"
                  << "  outward electromagnetic energy: "
                  << result.outward_electromagnetic_energy_joules << " J\n"
                  << "  model scope: educational 0-D scaling, not MHD or GRMHD\n";
        return;
    }
    std::cout << "Reduced toy-plasma transport estimate\n"
              << "  model: 0-D ideal-MHD-inspired scaling, not an MHD solution or GRMHD\n"
              << "  magnetization sigma:                 " << result.magnetization << "\n"
              << "  relativistic Alfven speed:           " << result.alfven_speed_m_s
              << " m/s\n"
              << "  raw Poynting-flux scaling:           " << result.poynting_power_watts
              << " W\n"
              << "  heuristic spin coupling factor:      "
              << result.spin_coupling_efficiency << "\n"
              << "  outward electromagnetic power:       "
              << result.outward_electromagnetic_power_watts << " W\n"
              << "  outward electromagnetic energy:      "
              << result.outward_electromagnetic_energy_joules << " J\n\n"
              << "Assumes a field perpendicular to representative flow. It does not calculate matter"
                 " flux, collector capture, conversion, storage, transmission, or usable energy.\n";
}

void print_penrose_event_result(const bh::EquatorialPenroseScenario& scenario,
                                const bh::PenroseEventResult& result,
                                const bool verbose) {
    if (!verbose) {
        std::cout << "Idealized equatorial Penrose event\n"
                  << "  status: "
                  << bh::penrose_event_status_name(result.status) << "\n"
                  << "  split (r/M, Lz, angle):    (" << result.split.split_radius_over_m << ", "
                  << result.split.incoming_lz_over_m_m << ", "
                  << result.split.split_angle_rad << ")\n"
                  << "  extraction / efficiency:  " << result.extracted_energy << " / "
                  << 100.0 * result.eta_penrose << " %\n"
                  << "  captured / escaping:      "
                  << termination_name(result.captured_trajectory.termination) << " / "
                  << termination_name(result.escaping_trajectory.termination) << "\n"
                  << "  maximum residual:         "
                  << result.maximum_normalized_residual << "\n";
        return;
    }
    std::cout << "Restricted equatorial Penrose event\n"
              << "  model: neutral test particle, local two-body split, G = c = 1\n"
              << "  status: " << bh::penrose_event_status_name(result.status) << "\n\n"
              << "Geometry\n"
              << "  outer horizon r+:        " << result.horizon_radius << " M\n"
              << "  equatorial static limit: " << result.static_limit_radius << " M\n"
              << "  split radius:            " << result.split_radius << " M\n\n"
              << "Integration control\n"
              << "  maximum step:            " << scenario.integration_step << "\n"
              << "  absolute tolerance:      "
              << scenario.integration_control.absolute_tolerance << "\n"
              << "  relative tolerance:      "
              << scenario.integration_control.relative_tolerance << "\n"
              << "  minimum step:            "
              << scenario.integration_control.minimum_step << "\n\n"
              << "Energy ledger (normalized geometric units)\n"
              << "  input energy:            " << result.input_energy << "\n"
              << "  captured-fragment energy: " << result.captured_energy << "\n"
              << "  escaping-fragment energy: " << result.escaping_energy << "\n"
              << "  net extracted energy:    " << result.extracted_energy << "\n"
              << "  Penrose efficiency:      " << 100.0 * result.eta_penrose << " %\n\n"
              << "Validation\n"
              << "  four-momentum residual:  " << result.four_momentum_residual << "\n"
              << "  mass-shell residual:     " << result.mass_shell_residual << "\n"
              << "  energy residual:         " << result.energy_conservation_residual << "\n"
              << "  angular-momentum residual: "
              << result.angular_momentum_conservation_residual << "\n"
              << "  geodesic initialization: " << result.geodesic_initialization_residual
              << "\n"
              << "  maximum normalized residual: " << result.maximum_normalized_residual
              << "\n"
              << "  incoming trajectory:     "
              << termination_name(result.incoming_trajectory.termination) << "\n"
              << "  captured trajectory:     "
              << termination_name(result.captured_trajectory.termination) << "\n"
              << "  escaping trajectory:     "
              << termination_name(result.escaping_trajectory.termination) << "\n\n"
              << "Integration diagnostics\n";
    print_integration_diagnostics("incoming", result.incoming_trajectory);
    print_integration_diagnostics("captured", result.captured_trajectory);
    print_integration_diagnostics("escaping", result.escaping_trajectory);
    std::cout << "\nThis is an idealized event calculation, not a measurement of extracted"
                 " astrophysical energy or a delivery model.\n";
}

void print_new_parameter_window_guidance(
    const bh::PenroseDijkstraSearchConfig& search) {
    std::cout << "The complete discrete window contains no candidate reaching "
              << 100.0 * search.eta_target << "%. Run a new window with at least one"
                 " split-radius, Lz, or angle interval below or above the current interval;"
                 " split-radius bounds must remain inside the ergosphere. Keep the new window"
                 " at or below "
              << bh::max_penrose_search_nodes
              << " nodes. Values between grid points require smaller steps and a"
                 " correspondingly smaller window.\n";
}

void print_penrose_dijkstra_result(const bh::EquatorialPenroseDijkstraInput& input,
                                   const bh::PenroseDijkstraSearchResult& result,
                                   const bool verbose,
                                   const bool include_physical_event = true) {
    const bh::PenroseDijkstraSearchConfig& search = input.search;
    const bh::PenroseDijkstraSearchDiagnostics& diagnostics = result.diagnostics;
    if (!verbose) {
        std::cout << "Dijkstra Penrose search\n"
                  << "  status:                    "
                  << bh::penrose_dijkstra_search_status_name(result.status) << "\n"
                  << "  work:                      " << diagnostics.nodes_evaluated << "/"
                  << diagnostics.candidates_in_domain << " nodes, "
                  << evaluated_nodes_per_second(diagnostics) << " nodes/s\n";
        if (!result.failure_message.empty()) {
            std::cout << "  detail:                    " << result.failure_message << "\n";
        }
        if (result.found && !result.parameter_adjustment_path.empty()) {
            const bh::PenroseDijkstraNode& selected = result.parameter_adjustment_path.back();
            std::cout << "  target / result:           " << 100.0 * search.eta_target
                      << " % / " << 100.0 * selected.eta_penrose << " %"
                      << (result.target_reached ? "\n" : " (bounded fallback)\n")
                      << "  adjustment cost:           " << selected.g_cost << "\n";
            if (include_physical_event) {
                std::cout << "  split (r/M, Lz, angle):     ("
                          << selected.split.split_radius_over_m << ", "
                          << selected.split.incoming_lz_over_m_m << ", "
                          << selected.split.split_angle_rad << ")\n"
                          << "  captured / escaping:       "
                          << termination_name(selected.captured_termination) << " / "
                          << termination_name(selected.escaping_termination) << "\n"
                          << "  maximum residual:          "
                          << selected.maximum_normalized_residual << "\n";
            }
        }
        if (result.status == bh::PenroseDijkstraSearchStatus::best_feasible_below_target) {
            std::cout << "  note: target not reached; showing the best validated candidate"
                         " in this grid\n";
        } else if (result.status ==
                   bh::PenroseDijkstraSearchStatus::no_solution_within_bounds) {
            std::cout << "  note: no positive, physically validated extraction in this grid\n";
        }
        return;
    }
    std::cout << "Dijkstra Penrose parameter search\n"
              << "  model: bounded parameter graph; each node calls the restricted"
                 " equatorial Penrose evaluator\n"
              << "  algorithm: " << bh::penrose_search_algorithm_name(search.algorithm)
              << " (h = 0, so f = g)\n"
              << "  execution backend: scalar-single-thread\n"
              << "  black-hole mass scale:   " << input.scenario.black_hole_mass << "\n"
              << "  dimensionless spin:      " << input.scenario.dimensionless_spin << "\n"
              << "  target Penrose efficiency: " << 100.0 * search.eta_target << " %\n"
              << "  start split:              ("
              << search.start.split_radius_over_m << ", "
              << search.start.incoming_lz_over_m_m << ", "
              << search.start.split_angle_rad << ")\n"
              << "  grid steps:               ("
              << search.step.split_radius_over_m << ", "
              << search.step.incoming_lz_over_m_m << ", "
              << search.step.split_angle_rad << ")\n"
              << "  lower bounds:             ("
              << search.lower_bound.split_radius_over_m << ", "
              << search.lower_bound.incoming_lz_over_m_m << ", "
              << search.lower_bound.split_angle_rad << ")\n"
              << "  upper bounds:             ("
              << search.upper_bound.split_radius_over_m << ", "
              << search.upper_bound.incoming_lz_over_m_m << ", "
              << search.upper_bound.split_angle_rad << ")\n"
              << "  edge costs:               (" << search.edge_costs[0] << ", "
              << search.edge_costs[1] << ", " << search.edge_costs[2] << ")\n"
              << "  maximum evaluations:      " << search.max_evaluations << "\n"
              << "  window node limit:        " << bh::max_penrose_search_nodes << "\n"
              << "  time budget:              " << search.time_budget.count()
              << " ms (0 = disabled)\n\n"
              << "Search diagnostics\n"
              << "  candidates in domain:     " << diagnostics.candidates_in_domain << "\n"
              << "  nodes generated:          " << diagnostics.nodes_generated << "\n"
              << "  nodes evaluated:          " << diagnostics.nodes_evaluated << "\n"
              << "  nodes expanded:           " << diagnostics.nodes_expanded << "\n"
              << "  duplicate nodes skipped:  " << diagnostics.duplicate_nodes_skipped << "\n"
              << "  out-of-domain neighbors:  "
              << diagnostics.outside_search_domain_neighbors << "\n"
              << "  fresh final checks:       "
              << diagnostics.final_verification_evaluations << "\n"
              << "  outside ergosphere:       " << diagnostics.outside_ergosphere << "\n"
              << "  physics invalid:          " << diagnostics.physics_invalid << "\n"
              << "  captured/non-escaping:    " << diagnostics.captured_or_non_escaping << "\n"
              << "  escaped below target:     " << diagnostics.escaping_without_target << "\n"
              << "  integration failed:       " << diagnostics.integration_failed << "\n"
              << "  goal feasible:            " << diagnostics.goal_feasible << "\n"
              << "  elapsed:                  "
              << std::chrono::duration<double, std::milli>(diagnostics.elapsed).count()
              << " ms\n"
              << "  throughput:               " << evaluated_nodes_per_second(diagnostics)
              << " nodes/s\n\n"
              << "Search outcome\n"
              << "  status: " << bh::penrose_dijkstra_search_status_name(result.status) << "\n";

    if (!result.failure_message.empty()) {
        std::cout << "  detail: " << result.failure_message << "\n";
    }

    if (!result.found) {
        std::cout << "\nNo validated candidate is available for this terminated search.\n";
        if (result.status == bh::PenroseDijkstraSearchStatus::no_solution_within_bounds) {
            print_new_parameter_window_guidance(search);
        }
        return;
    }

    const bh::PenroseDijkstraNode& selected = result.parameter_adjustment_path.back();
    std::cout << "\nSelected candidate\n"
              << "  target reached:           "
              << (result.target_reached ? "yes" : "no; returning bounded fallback") << "\n"
              << "  status: " << bh::penrose_dijkstra_node_status_name(selected.status) << "\n"
              << "  parameter-adjustment cost: " << selected.g_cost << "\n"
              << "  h / f:                    " << selected.h_cost << " / " << selected.f_cost
              << "\n"
              << "  split parameters:         ("
              << selected.split.split_radius_over_m << ", "
              << selected.split.incoming_lz_over_m_m << ", "
              << selected.split.split_angle_rad << ")\n"
              << "  Penrose efficiency:       " << 100.0 * selected.eta_penrose << " %\n"
              << "  net extracted energy:     " << selected.extracted_energy << "\n"
              << "  maximum residual:         " << selected.maximum_normalized_residual << "\n\n"
              << "Parameter-adjustment trace (not a particle trajectory)\n";
    if (result.target_reached) {
        std::cout << "  Goal contract: this parameter set passed the Penrose evaluator,"
                     " capture/escape checks, conservation tolerance, and the configured"
                     " efficiency threshold.\n";
    }
    for (const bh::PenroseDijkstraNode& node : result.parameter_adjustment_path) {
        std::cout << "  g/h/f=" << node.g_cost << "/" << node.h_cost << "/" << node.f_cost
                  << " key=(" << node.key[0] << ", " << node.key[1] << ", " << node.key[2]
                  << ") delta=(" << node.local_change[0] << ", " << node.local_change[1]
                  << ", " << node.local_change[2] << ") split=("
                  << node.split.split_radius_over_m << ", "
                  << node.split.incoming_lz_over_m_m << ", " << node.split.split_angle_rad
                  << ") status=" << bh::penrose_dijkstra_node_status_name(node.status)
                  << " eta=" << 100.0 * node.eta_penrose << "% extracted="
                  << node.extracted_energy << " residual=" << node.maximum_normalized_residual
                  << " capture=" << termination_name(node.captured_termination)
                  << " escape=" << termination_name(node.escaping_termination) << "\n";
    }
    if (result.status == bh::PenroseDijkstraSearchStatus::best_feasible_below_target) {
        std::cout << "\n";
        print_new_parameter_window_guidance(search);
    }
    if (include_physical_event) {
        std::cout << "\nPhysical selected-event diagnostics\n";
        print_penrose_event_result(input.scenario, result.selected_event, true);
    }
}

void print_penrose_phase_map_result(const bh::EquatorialPenroseDijkstraInput& input,
                                     const bh::PenrosePhaseMapResult& result,
                                     const bool verbose,
                                     const bool include_physical_event = true) {
    const bh::PenroseDijkstraSearchConfig& search = input.search;
    const bh::PenroseDijkstraSearchDiagnostics& diagnostics = result.diagnostics;
    const std::string_view execution_backend =
        diagnostics.avx2_four_lane_batches > 0
            ? "avx2-four-lane-single-thread"
            : diagnostics.four_lane_batches > 0
                  ? "portable-scalar-batch4-single-thread"
                  : "scalar-single-thread";
    if (!verbose) {
        std::cout << "Bounded Penrose phase-space map\n"
                  << "  status / backend:         "
                  << bh::penrose_phase_map_status_name(result.status) << " / "
                  << execution_backend << "\n"
                  << "  work:                     " << diagnostics.nodes_evaluated << "/"
                  << diagnostics.candidates_in_domain << " nodes, "
                  << evaluated_nodes_per_second(diagnostics) << " nodes/s\n"
                  << "  AVX2 batches:             "
                  << diagnostics.avx2_four_lane_batches << "\n";
        if (!result.failure_message.empty()) {
            std::cout << "  detail:                   " << result.failure_message << "\n";
        }
        if (result.best_validated_candidate) {
            const bh::PenroseDijkstraNode& best = *result.best_validated_candidate;
            std::cout << "  target / best:            " << 100.0 * search.eta_target
                      << " % / " << 100.0 * best.eta_penrose << " %\n"
                      << "  split (r/M, Lz, angle):   (" << best.split.split_radius_over_m << ", "
                      << best.split.incoming_lz_over_m_m << ", "
                      << best.split.split_angle_rad << ")\n"
                      << "  maximum residual:         "
                      << best.maximum_normalized_residual << "\n";
        } else {
            std::cout << "  best validated candidate: none\n";
        }
        if (result.complete && diagnostics.goal_feasible == 0) {
            if (result.best_validated_candidate) {
                std::cout << "  note: target not reached; showing the best validated candidate"
                             " in this grid\n";
            } else {
                std::cout << "  note: no positive, physically validated extraction in this grid\n";
            }
        }
        return;
    }
    std::cout << "Bounded Penrose phase-space map\n"
              << "  model: four-lane exhaustive evaluation of the declared parameter grid\n"
              << "  execution backend: " << execution_backend << "\n"
              << "  lower bounds:             ("
              << search.lower_bound.split_radius_over_m << ", "
              << search.lower_bound.incoming_lz_over_m_m << ", "
              << search.lower_bound.split_angle_rad << ")\n"
              << "  upper bounds:             ("
              << search.upper_bound.split_radius_over_m << ", "
              << search.upper_bound.incoming_lz_over_m_m << ", "
              << search.upper_bound.split_angle_rad << ")\n"
              << "  grid steps:               ("
              << search.step.split_radius_over_m << ", "
              << search.step.incoming_lz_over_m_m << ", "
              << search.step.split_angle_rad << ")\n"
              << "  candidates in domain:     " << diagnostics.candidates_in_domain << "\n"
              << "  nodes evaluated:          " << diagnostics.nodes_evaluated << "\n"
              << "  four-lane batches:        " << diagnostics.four_lane_batches << "\n"
              << "  AVX2 four-lane batches:   " << diagnostics.avx2_four_lane_batches << "\n"
              << "  four-lane nodes:          " << diagnostics.four_lane_nodes << "\n"
              << "  scalar tail nodes:        " << diagnostics.scalar_nodes << "\n"
              << "  retained map entries:     " << result.candidates.size() << "\n"
              << "  maximum evaluations:      " << search.max_evaluations << "\n"
              << "  window node limit:        " << bh::max_penrose_search_nodes << "\n"
              << "  elapsed:                  "
              << std::chrono::duration<double, std::milli>(diagnostics.elapsed).count()
              << " ms\n"
              << "  throughput:               " << evaluated_nodes_per_second(diagnostics)
              << " nodes/s\n"
              << "  status: " << bh::penrose_phase_map_status_name(result.status) << "\n";
    if (!result.failure_message.empty()) {
        std::cout << "  detail: " << result.failure_message << "\n";
    }
    std::cout << "\nStatus counts\n"
              << "  physics invalid:          " << diagnostics.physics_invalid << "\n"
              << "  captured/non-escaping:    " << diagnostics.captured_or_non_escaping << "\n"
              << "  escaped below target:     " << diagnostics.escaping_without_target << "\n"
              << "  integration failed:       " << diagnostics.integration_failed << "\n"
              << "  goal feasible:            " << diagnostics.goal_feasible << "\n";

    if (!result.best_validated_candidate) {
        std::cout << "\nNo positive, physically validated extraction was found in the"
                     " evaluated portion of this grid.\n";
        if (result.complete) {
            print_new_parameter_window_guidance(search);
        }
        return;
    }

    const bh::PenroseDijkstraNode& best = *result.best_validated_candidate;
    std::cout << "\nBest validated candidate "
              << (result.complete ? "within the fully evaluated bounded grid" :
                                    "encountered before termination")
              << "\n"
              << "  split parameters:         (" << best.split.split_radius_over_m << ", "
              << best.split.incoming_lz_over_m_m << ", " << best.split.split_angle_rad
              << ")\n"
              << "  net extracted energy:     " << best.extracted_energy << "\n"
              << "  Penrose efficiency:       " << 100.0 * best.eta_penrose << " %\n"
              << "  maximum residual:         " << best.maximum_normalized_residual << "\n";
    if (!result.complete) {
        std::cout << "\nThis is not a maximum over the full declared grid because the map"
                     " terminated early.\n";
        return;
    }

    if (diagnostics.goal_feasible == 0) {
        std::cout << "\n";
        print_new_parameter_window_guidance(search);
    }

    if (include_physical_event) {
        std::cout << "\nPhysical best-event diagnostics\n";
        print_penrose_event_result(input.scenario, result.best_event, true);
    }
}

void run_algebraic(const std::span<const std::string_view> arguments,
                   const CliOptions& options) {
    if (arguments.size() != 2) {
        throw std::invalid_argument("algebraic requires <mass_kg> <a_star>");
    }
    const bh::RotationalEnergyResult result = bh::rotational_energy(
        parse_double(arguments[0], "mass_kg"), parse_double(arguments[1], "a_star"));
    if (options.format == OutputFormat::json) {
        bh::cli_json::write_algebraic(std::cout, result);
    } else {
        print_algebraic_result(result, options.verbose);
    }
}

void run_algebraic_range(const std::span<const std::string_view> arguments,
                         const CliOptions& options) {
    if (arguments.size() != 6) {
        throw std::invalid_argument(
            "algebraic-range requires three mass values followed by three spin values");
    }

    const bh::RotationalEnergyRangeResult result = bh::rotational_energy_range(
        {{parse_double(arguments[0], "mass_low_kg"),
          parse_double(arguments[1], "mass_mid_kg"),
          parse_double(arguments[2], "mass_high_kg")},
         {parse_double(arguments[3], "spin_low"),
          parse_double(arguments[4], "spin_mid"),
          parse_double(arguments[5], "spin_high")}});
    if (options.format == OutputFormat::json) {
        bh::cli_json::write_algebraic_range(std::cout, result);
    } else {
        print_algebraic_range_result(result, options.verbose);
    }
}

void run_toy_plasma(const std::span<const std::string_view> arguments,
                    const CliOptions& options) {
    if (arguments.size() != 5) {
        throw std::invalid_argument(
            "toy-plasma requires <magnetic_field_tesla> <mass_density_kg_m3> <flow_area_m2>"
            " <a_star> <duration_seconds>");
    }

    const bh::PlasmaInput input{
        parse_double(arguments[0], "magnetic_field_tesla"),
        parse_double(arguments[1], "mass_density_kg_m3"),
        parse_double(arguments[2], "flow_area_m2"),
        parse_double(arguments[3], "a_star"),
        parse_double(arguments[4], "duration_seconds")};
    const bh::PlasmaResult result = bh::estimate_toy_plasma_transport(input);

    if (options.format == OutputFormat::json) {
        bh::cli_json::write_toy_plasma(std::cout, input, result);
    } else {
        print_toy_plasma_result(result, options.verbose);
    }
}

void run_penrose_event(const std::span<const std::string_view> arguments,
                       const CliOptions& options) {
    if (arguments.size() != 1) {
        throw std::invalid_argument("scenario requires a path to a .cfg scenario file");
    }

    const bh::EquatorialPenroseEventInput input =
        bh::load_equatorial_penrose_event_input(std::string(arguments[0]));
    const bh::PenroseEventResult result =
        bh::evaluate_equatorial_penrose_event(input.scenario, input.split);
    if (options.format == OutputFormat::json) {
        bh::cli_json::write_penrose_event(std::cout, input.scenario, result);
    } else {
        print_penrose_event_result(input.scenario, result, options.verbose);
    }
}

void run_penrose_dijkstra(const std::span<const std::string_view> arguments,
                          const CliOptions& options) {
    if (arguments.size() != 1) {
        throw std::invalid_argument("search requires a path to a .cfg search scenario");
    }

    const bh::EquatorialPenroseDijkstraInput input =
        bh::load_equatorial_penrose_dijkstra_input(std::string(arguments[0]));
    bh::require_complete_penrose_search_window(input.scenario, input.search);
    const bh::PenroseSearchWindowSummary window =
        bh::describe_penrose_search_window(input.scenario, input.search);
    TerminalProgress display;
    const bh::PenroseDijkstraSearchResult result = bh::find_penrose_dijkstra_path(
        input.scenario, input.search, {},
        make_progress_observer(window.candidates, options, display));
    display.finish();
    if (options.format == OutputFormat::json) {
        bh::cli_json::write_penrose_search(std::cout, input, result);
    } else {
        print_penrose_dijkstra_result(input, result, options.verbose);
    }
}

void run_penrose_phase_map(const std::span<const std::string_view> arguments,
                           const CliOptions& options) {
    if (arguments.size() != 1) {
        throw std::invalid_argument("map requires a path to a .cfg search scenario");
    }

    const bh::EquatorialPenroseDijkstraInput input =
        bh::load_equatorial_penrose_dijkstra_input(std::string(arguments[0]));
    bh::require_complete_penrose_search_window(input.scenario, input.search);
    const bh::PenroseSearchWindowSummary window =
        bh::describe_penrose_search_window(input.scenario, input.search);
    TerminalProgress display;
    const bh::PenrosePhaseMapResult result = bh::evaluate_penrose_phase_map(
        input.scenario, input.search, {},
        make_progress_observer(window.candidates, options, display));
    display.finish();
    if (options.format == OutputFormat::json) {
        bh::cli_json::write_penrose_phase_map(std::cout, input, result);
    } else {
        print_penrose_phase_map_result(input, result, options.verbose);
    }
}

constexpr double normalized_kerr_mass = 1.0;

struct SharedBlackHoleState {
    double mass_kg{1.98847e31};
    double dimensionless_spin{0.999};
};

struct SharedTrajectoryControl {
    double maximum_step{0.002};
    std::size_t maximum_steps{50'000};
    bh::KerrIntegrationControl integration_control{};
};

struct PenroseEventState {
    double parent_rest_mass{1.0};
    double fragment_rest_mass{0.0};
    double incoming_specific_energy{1.0};
    double initial_radius_over_m{10.0};
    double escape_radius_over_m{20.0};
    double residual_tolerance{1.0e-7};
};

struct PenroseSearchState {
    bh::PenroseSplitParameters start{1.095, 2.07, -2.0};
    bh::PenroseSplitParameters lower_bound{1.05, 2.03, -2.12};
    bh::PenroseSplitParameters upper_bound{1.941, 2.12, -1.88};
    bh::PenroseSplitParameters step{0.009, 0.01, 0.01};
};

struct InteractiveSession {
    SharedBlackHoleState black_hole{};
    SharedTrajectoryControl trajectory_control{};
    PenroseEventState penrose{};
    PenroseSearchState penrose_search{};
};

bh::EquatorialPenroseScenario make_interactive_penrose_scenario(
    const InteractiveSession& session) {
    const PenroseEventState& input = session.penrose;
    return {normalized_kerr_mass,
            session.black_hole.dimensionless_spin,
            input.parent_rest_mass,
            input.fragment_rest_mass,
            input.incoming_specific_energy,
            input.initial_radius_over_m,
            input.escape_radius_over_m,
            session.trajectory_control.maximum_step,
            session.trajectory_control.maximum_steps,
            session.trajectory_control.integration_control,
            input.residual_tolerance};
}

void prompt_penrose_scenario_inputs(InteractiveSession& session) {
    PenroseEventState& input = session.penrose;
    input.parent_rest_mass = prompt_double("Parent rest mass", input.parent_rest_mass);
    input.fragment_rest_mass = prompt_double("Fragment rest mass", input.fragment_rest_mass);
    input.incoming_specific_energy = prompt_double("Incoming specific energy",
                                                    input.incoming_specific_energy);
    input.initial_radius_over_m = prompt_double("Initial radius r / M",
                                                input.initial_radius_over_m);
    input.escape_radius_over_m = prompt_double("Escape radius r / M",
                                               input.escape_radius_over_m);
    input.residual_tolerance = prompt_double("Conservation residual tolerance",
                                             input.residual_tolerance);
}

void configure_shared_black_hole(InteractiveSession& session) {
    std::cout << "\nShared black-hole configuration\n"
              << "  These values feed the algebraic result and fixed Kerr/Penrose scenario.\n";
    session.black_hole.mass_kg = prompt_double("Black-hole mass [kg]", session.black_hole.mass_kg);
    session.black_hole.dimensionless_spin = prompt_double(
        "Dimensionless spin a_star", session.black_hole.dimensionless_spin);
}

void configure_shared_trajectory_control(InteractiveSession& session) {
    std::cout << "\nShared Kerr integration controls\n"
              << "  Kerr and Penrose trajectory validation use these same settings.\n";
    session.trajectory_control.maximum_step = prompt_double(
        "Maximum RK4 step", session.trajectory_control.maximum_step);
    session.trajectory_control.maximum_steps = prompt_step_count(
        "Maximum integration steps", session.trajectory_control.maximum_steps);
    session.trajectory_control.integration_control.absolute_tolerance = prompt_double(
        "Absolute integration tolerance",
        session.trajectory_control.integration_control.absolute_tolerance);
    session.trajectory_control.integration_control.relative_tolerance = prompt_double(
        "Relative integration tolerance",
        session.trajectory_control.integration_control.relative_tolerance);
    session.trajectory_control.integration_control.minimum_step = prompt_double(
        "Minimum integration step (0 = automatic)",
        session.trajectory_control.integration_control.minimum_step);
}

void print_shared_session_header(const InteractiveSession& session) {
    std::cout << "\nFixed scenario\n"
              << "  black-hole mass:         " << session.black_hole.mass_kg << " kg\n"
              << "  dimensionless spin:       " << session.black_hole.dimensionless_spin << "\n"
              << "  Kerr/Penrose mass scale:  M = 1 (radii and angular momentum normalized by M)\n";
}

void run_interactive_algebraic(const InteractiveSession& session, const bool verbose) {
    print_algebraic_result(
        bh::rotational_energy(session.black_hole.mass_kg, session.black_hole.dimensionless_spin),
        verbose);
}

bh::EquatorialPenroseDijkstraInput prompt_interactive_phase_map_input(
    InteractiveSession& session) {
    const double spin_length = bh::kerr_spin_length(
        normalized_kerr_mass, session.black_hole.dimensionless_spin);
    const double horizon = bh::kerr_outer_horizon(normalized_kerr_mass, spin_length);
    const double static_limit = bh::kerr_static_limit_radius(
        normalized_kerr_mass, spin_length, 1.57079632679489661923);
    std::cout << "  Configure the bounded parameter grid used by the exhaustive phase map.\n"
              << "  Every candidate is passed to Penrose, which calls Kerr for the"
                  " incoming, captured, and escaping paths.\n"
              << "  Four adjacent angle candidates are batched when the backend supports it.\n"
              << "  Maximum map size: " << bh::max_penrose_search_nodes << " nodes.\n"
              << "  Valid split-radius interval: " << horizon << " < r_split / M < "
              << static_limit << "\n";

    PenroseSearchState& input = session.penrose_search;
    input.start.split_radius_over_m = prompt_double(
        "Starting split radius r_split / M", input.start.split_radius_over_m);
    input.start.incoming_lz_over_m_m = prompt_double(
        "Starting angular momentum Lz / (m M)", input.start.incoming_lz_over_m_m);
    input.start.split_angle_rad = prompt_double(
        "Starting split angle [radians]", input.start.split_angle_rad);

    input.lower_bound.split_radius_over_m = prompt_double(
        "Split-radius lower bound", input.lower_bound.split_radius_over_m);
    input.upper_bound.split_radius_over_m = prompt_double(
        "Split-radius upper bound", input.upper_bound.split_radius_over_m);
    input.step.split_radius_over_m = prompt_double(
        "Split-radius step", input.step.split_radius_over_m);
    input.lower_bound.incoming_lz_over_m_m = prompt_double(
        "Lz lower bound", input.lower_bound.incoming_lz_over_m_m);
    input.upper_bound.incoming_lz_over_m_m = prompt_double(
        "Lz upper bound", input.upper_bound.incoming_lz_over_m_m);
    input.step.incoming_lz_over_m_m = prompt_double(
        "Lz step", input.step.incoming_lz_over_m_m);
    input.lower_bound.split_angle_rad = prompt_double(
        "Split-angle lower bound [radians]", input.lower_bound.split_angle_rad);
    input.upper_bound.split_angle_rad = prompt_double(
        "Split-angle upper bound [radians]", input.upper_bound.split_angle_rad);
    input.step.split_angle_rad = prompt_double(
        "Split-angle step [radians]", input.step.split_angle_rad);

    bh::PenroseDijkstraSearchConfig search;
    search.start = input.start;
    search.lower_bound = input.lower_bound;
    search.upper_bound = input.upper_bound;
    search.step = input.step;
    search.eta_target = 0.15;
    search.max_evaluations = bh::max_penrose_search_nodes;
    search.time_budget = std::chrono::milliseconds(0);

    return {make_interactive_penrose_scenario(session), search};
}

void run_interactive(const CliOptions& options) {
    InteractiveSession session;
    std::cout << "Guided black-hole simulation pipeline\n"
              << "  Inputs are collected once in execution order.\n"
              << "  Candidate parameters are configured before evaluation; the phase map"
                 " exists only after Penrose and Kerr evaluate those candidates.\n";

    std::cout << "\nStage 1/7 - Black-hole inputs\n";
    configure_shared_black_hole(session);

    std::cout << "\nStage 2/7 - Kerr integration controls\n";
    configure_shared_trajectory_control(session);
    print_shared_session_header(session);

    std::cout << "\nStage 3/7 - Algebraic Kerr reservoir\n";
    run_interactive_algebraic(session, options.verbose);

    std::cout << "\nStage 4/7 - Particle and event inputs\n"
              << "  Energies and particle masses use normalized geometrized units, not joules.\n";
    prompt_penrose_scenario_inputs(session);

    std::cout << "\nStage 5/7 - Candidate parameter grid\n";
    const bh::EquatorialPenroseDijkstraInput request =
        prompt_interactive_phase_map_input(session);
    const bh::PenroseSearchWindowSummary window =
        bh::describe_penrose_search_window(request.scenario, request.search);
    bh::require_complete_penrose_search_window(request.scenario, request.search);
    std::cout << "\nAccepted phase-map window: "
              << window.dimension_sizes[0] << " x "
              << window.dimension_sizes[1] << " x "
              << window.dimension_sizes[2] << " = "
              << window.candidates << " nodes.\n";

    std::cout << "\nStage 6/7 - Exhaustive phase-map evaluation\n" << std::flush;
    TerminalProgress display;
    const bh::PenrosePhaseMapResult result = bh::evaluate_penrose_phase_map(
        request.scenario, request.search, {},
        make_progress_observer(window.candidates, options, display));
    display.finish();
    print_penrose_phase_map_result(request, result, options.verbose, false);

    std::cout << "\nStage 7/7 - Selected Penrose event\n";
    if (!result.complete) {
        std::cout << "The phase map did not complete, so no bounded-grid result is shown.\n";
        return;
    }
    if (!result.best_validated_candidate) {
        std::cout << "No validated candidate was selected, so there is no final event to show.\n";
        return;
    }
    print_penrose_event_result(request.scenario, result.best_event, options.verbose);
}
}  // namespace

int main(const int argc, char* argv[]) {
    std::cout << std::scientific << std::setprecision(6);
    try {
        const ParsedCli cli = parse_cli(argc, argv);
        const std::span<const std::string_view> arguments(cli.arguments);

        if (cli.command == "help") {
            if (!arguments.empty()) {
                throw std::invalid_argument("help does not accept positional arguments");
            }
            print_usage();
        } else if (cli.command == "version") {
            if (!arguments.empty()) {
                throw std::invalid_argument("version does not accept positional arguments");
            }
            if (cli.options.format == OutputFormat::json) {
                bh::cli_json::write_version(std::cout, BH_VERSION);
            } else {
                print_version();
            }
        } else if (cli.command == "interactive") {
            if (!arguments.empty()) {
                throw std::invalid_argument("interactive does not accept positional arguments");
            }
            if (cli.options.format == OutputFormat::json) {
                throw std::invalid_argument(
                    "interactive mode requires text output; use a non-interactive command for JSON");
            }
            run_interactive(cli.options);
        } else if (cli.command == "algebraic") {
            run_algebraic(arguments, cli.options);
        } else if (cli.command == "algebraic-range") {
            run_algebraic_range(arguments, cli.options);
        } else if (cli.command == "toy-plasma") {
            run_toy_plasma(arguments, cli.options);
        } else if (cli.command == "scenario") {
            run_penrose_event(arguments, cli.options);
        } else if (cli.command == "search") {
            run_penrose_dijkstra(arguments, cli.options);
        } else if (cli.command == "map") {
            run_penrose_phase_map(arguments, cli.options);
        } else {
            std::cerr << "Error: unknown command '" << cli.command << "'\n\n";
            print_usage();
            return 2;
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
