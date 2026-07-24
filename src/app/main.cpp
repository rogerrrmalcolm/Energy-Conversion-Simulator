#include "bh/algebraic_model.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/trajectory.hpp"

#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

namespace {
double parse_double(const std::string_view text, const std::string_view argument_name) {
    double value{};
    const auto [end, error] = std::from_chars(
        text.data(), text.data() + text.size(), value, std::chars_format::general);
    if (error != std::errc{} || end != text.data() + text.size() || !std::isfinite(value)) {
        throw std::invalid_argument("expected a finite number for " + std::string(argument_name));
    }
    return value;
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
              << "  " << label << " final step: "
              << trajectory.diagnostics.final_step << "\n";
}

void print_usage(const std::string_view program_name) {
    std::cout << "Usage:\n"
              << "  " << program_name << " --algebraic <mass_kg> <a_star>\n"
              << "  " << program_name
              << " --algebraic-range <mass_low_kg> <mass_mid_kg> <mass_high_kg>"
                 " <spin_low> <spin_mid> <spin_high>\n"
              << "  " << program_name << " --scenario <path-to-event.cfg>\n\n"
              << "The scenario command evaluates one declared equatorial Penrose event."
                 " It does not search or optimize parameters.\n";
}

void print_algebraic_result(const bh::RotationalEnergyResult& result) {
    std::cout << "Algebraic Kerr rotational-energy reservoir\n"
              << "  total mass-energy:       " << result.mass_energy_joules << " J\n"
              << "  irreducible mass:        " << result.irreducible_mass_kg << " kg\n"
              << "  rotational reservoir:    " << result.rotational_energy_joules << " J\n"
              << "  rotational fraction:     " << 100.0 * result.rotational_fraction << " %\n";
}

void run_algebraic(const int argc, char* argv[]) {
    if (argc != 4) {
        throw std::invalid_argument("--algebraic requires <mass_kg> <a_star>");
    }
    print_algebraic_result(
        bh::rotational_energy(parse_double(argv[2], "mass_kg"), parse_double(argv[3], "a_star")));
}

void run_algebraic_range(const int argc, char* argv[]) {
    if (argc != 8) {
        throw std::invalid_argument(
            "--algebraic-range requires three mass values followed by three spin values");
    }

    const bh::RotationalEnergyRangeResult result = bh::rotational_energy_range(
        {{parse_double(argv[2], "mass_low_kg"), parse_double(argv[3], "mass_mid_kg"),
          parse_double(argv[4], "mass_high_kg")},
         {parse_double(argv[5], "spin_low"), parse_double(argv[6], "spin_mid"),
          parse_double(argv[7], "spin_high")}});
    std::cout << "Algebraic Kerr rotational-energy range\n";
    std::cout << "  lower reservoir:         " << result.lower.rotational_energy_joules << " J\n";
    std::cout << "  central reservoir:       " << result.central.rotational_energy_joules
              << " J\n";
    std::cout << "  upper reservoir:         " << result.upper.rotational_energy_joules << " J\n";
    std::cout << "  lower uncertainty:       -" << result.rotational_energy_uncertainty_minus_joules
              << " J\n";
    std::cout << "  upper uncertainty:       +" << result.rotational_energy_uncertainty_plus_joules
              << " J\n";
}

void run_penrose_event(const int argc, char* argv[]) {
    if (argc != 3) {
        throw std::invalid_argument("--scenario requires a path to a .cfg scenario file");
    }

    const bh::EquatorialPenroseEventInput input =
        bh::load_equatorial_penrose_event_input(argv[2]);
    const bh::PenroseEventResult result =
        bh::evaluate_equatorial_penrose_event(input.scenario, input.split);

    std::cout << "Restricted equatorial Penrose event\n"
              << "  model: neutral test particle, local two-body split, G = c = 1\n"
              << "  status: " << bh::penrose_event_status_name(result.status) << "\n\n"
              << "Geometry\n"
              << "  outer horizon r+:        " << result.horizon_radius << " M\n"
              << "  equatorial static limit: " << result.static_limit_radius << " M\n"
              << "  split radius:            " << result.split_radius << " M\n\n"
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
}  // namespace

int main(const int argc, char* argv[]) {
    std::cout << std::scientific << std::setprecision(6);
    try {
        if (argc == 2 && (std::string_view(argv[1]) == "--help" ||
                          std::string_view(argv[1]) == "-h")) {
            print_usage(argv[0]);
            return 0;
        }
        if (argc >= 2 && std::string_view(argv[1]) == "--algebraic") {
            run_algebraic(argc, argv);
            return 0;
        }
        if (argc >= 2 && std::string_view(argv[1]) == "--algebraic-range") {
            run_algebraic_range(argc, argv);
            return 0;
        }
        if (argc >= 2 && std::string_view(argv[1]) == "--scenario") {
            run_penrose_event(argc, argv);
            return 0;
        }

        print_usage(argv[0]);
        return 2;
    } catch (const std::exception& error) {
        std::cerr << "Error: " << error.what() << '\n';
        return 1;
    }
}
