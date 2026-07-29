#include "bh/algebraic_model.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"
#include "bh/trajectory.hpp"

#include <charconv>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
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

int prompt_choice(const std::string_view label, const int default_value, const int minimum,
                  const int maximum) {
    while (true) {
        const double value = prompt_double(label, static_cast<double>(default_value));
        if (std::floor(value) == value && value >= minimum && value <= maximum) {
            return static_cast<int>(value);
        }
        std::cout << "  Enter an integer from " << minimum << " to " << maximum << ".\n";
    }
}

int prompt_radial_direction() {
    while (true) {
        const int direction = prompt_choice("Radial direction (-1 = inward, 1 = outward)", -1,
                                            -1, 1);
        if (direction != 0) {
            return direction;
        }
        std::cout << "  Radial direction must be -1 or 1.\n";
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
              << "  " << label << " final step: "
              << trajectory.diagnostics.final_step << "\n";
    if (!trajectory.points.empty()) {
        std::cout << "  " << label << " terminal radius: "
                  << trajectory.points.back().radius << " M\n";
    }
}

void print_usage(const std::string_view program_name) {
    std::cout << "Usage:\n"
              << "  " << program_name << " --algebraic <mass_kg> <a_star>\n"
              << "  " << program_name
              << " --algebraic-range <mass_low_kg> <mass_mid_kg> <mass_high_kg>"
                 " <spin_low> <spin_mid> <spin_high>\n"
              << "  " << program_name
              << " --toy-plasma <magnetic_field_tesla> <mass_density_kg_m3> <flow_area_m2>"
                 " <a_star> <duration_seconds>\n"
              << "  " << program_name << " --scenario <path-to-event.cfg>\n\n"
              << "  " << program_name << " --interactive\n\n"
              << "The scenario command evaluates one declared equatorial Penrose event."
                 " It does not search or optimize parameters.\n"
              << "The toy-plasma command is a reduced 0-D electromagnetic transport estimate,"
                 " not an MHD or GRMHD simulation.\n"
              << "The interactive command opens a shared session and retains black-hole state"
                 " across engines.\n";
}

void print_algebraic_result(const bh::RotationalEnergyResult& result) {
    std::cout << "Algebraic Kerr rotational-energy reservoir\n"
              << "  total mass-energy:       " << result.mass_energy_joules << " J\n"
              << "  irreducible mass:        " << result.irreducible_mass_kg << " kg\n"
              << "  rotational reservoir:    " << result.rotational_energy_joules << " J\n"
              << "  rotational fraction:     " << 100.0 * result.rotational_fraction << " %\n";
}

void print_algebraic_range_result(const bh::RotationalEnergyRangeResult& result) {
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

void print_toy_plasma_result(const bh::PlasmaResult& result) {
    std::cout << "Reduced toy-plasma transport estimate\n"
              << "  model: 0-D ideal-MHD scaling, not MHD or GRMHD\n"
              << "  magnetization sigma:                 " << result.magnetization << "\n"
              << "  relativistic Alfven speed:           " << result.alfven_speed_m_s
              << " m/s\n"
              << "  raw electromagnetic Poynting power:  " << result.poynting_power_watts
              << " W\n"
              << "  heuristic spin coupling:             "
              << result.spin_coupling_efficiency << "\n"
              << "  outward electromagnetic power:       "
              << result.outward_electromagnetic_power_watts << " W\n"
              << "  outward electromagnetic energy:      "
              << result.outward_electromagnetic_energy_joules << " J\n\n"
              << "This model does not calculate matter flux, collector capture, conversion,"
                 " storage, transmission, or usable delivered energy.\n";
}

void print_penrose_event_result(const bh::EquatorialPenroseScenario& scenario,
                                const bh::PenroseEventResult& result) {
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
    print_algebraic_range_result(result);
}

void run_toy_plasma(const int argc, char* argv[]) {
    if (argc != 7) {
        throw std::invalid_argument(
            "--toy-plasma requires <magnetic_field_tesla> <mass_density_kg_m3> <flow_area_m2>"
            " <a_star> <duration_seconds>");
    }

    const bh::PlasmaInput input{
        parse_double(argv[2], "magnetic_field_tesla"),
        parse_double(argv[3], "mass_density_kg_m3"),
        parse_double(argv[4], "flow_area_m2"),
        parse_double(argv[5], "a_star"),
        parse_double(argv[6], "duration_seconds")};
    const bh::PlasmaResult result = bh::estimate_plasma_extraction(input);

    print_toy_plasma_result(result);
}

void run_penrose_event(const int argc, char* argv[]) {
    if (argc != 3) {
        throw std::invalid_argument("--scenario requires a path to a .cfg scenario file");
    }

    const bh::EquatorialPenroseEventInput input =
        bh::load_equatorial_penrose_event_input(argv[2]);
    const bh::PenroseEventResult result =
        bh::evaluate_equatorial_penrose_event(input.scenario, input.split);
    print_penrose_event_result(input.scenario, result);
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

struct KerrPathState {
    double energy{1.0};
    double angular_momentum{0.0};
    double rest_mass{1.0};
    int radial_direction{-1};
    double initial_radius_over_m{10.0};
    double escape_radius_over_m{20.0};
};

struct PenroseEventState {
    double parent_rest_mass{1.0};
    double fragment_rest_mass{0.0};
    double incoming_specific_energy{1.0};
    double initial_radius_over_m{10.0};
    double escape_radius_over_m{20.0};
    double residual_tolerance{1.0e-7};
    double split_radius_over_m{1.10};
    double incoming_lz_over_m_m{2.07};
    double split_angle_rad{-2.0};
};

struct ToyPlasmaState {
    double magnetic_field_tesla{1.0};
    double mass_density_kg_m3{1.0e-8};
    double flow_area_m2{10.0};
    double duration_seconds{2.0};
};

struct InteractiveSession {
    SharedBlackHoleState black_hole{};
    SharedTrajectoryControl trajectory_control{};
    KerrPathState kerr{};
    PenroseEventState penrose{};
    ToyPlasmaState toy_plasma{};
};

void configure_shared_black_hole(InteractiveSession& session) {
    std::cout << "\nShared black-hole configuration\n"
              << "  These values are retained for Algebraic, Kerr, Penrose, and plasma runs.\n";
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
    std::cout << "\nShared simulation session\n"
              << "  black-hole mass:         " << session.black_hole.mass_kg << " kg\n"
              << "  dimensionless spin:       " << session.black_hole.dimensionless_spin << "\n"
              << "  Kerr/Penrose mass scale:  M = 1 (radii and angular momentum normalized by M)\n";
}

void run_interactive_algebraic(const InteractiveSession& session) {
    std::cout << "\nAlgebraic Kerr reservoir using the shared black-hole state\n";
    print_algebraic_result(
        bh::rotational_energy(session.black_hole.mass_kg, session.black_hole.dimensionless_spin));
}

void run_interactive_kerr(InteractiveSession& session) {
    std::cout << "\nEquatorial Kerr trajectory input\n"
              << "  The shared black-hole mass and spin are already loaded.\n"
              << "  Enter only this orbit's conserved quantities and boundaries.\n";
    KerrPathState& input = session.kerr;
    input.energy = prompt_double("Conserved energy E", input.energy);
    input.angular_momentum = prompt_double("Conserved axial angular momentum Lz",
                                            input.angular_momentum);
    input.rest_mass = prompt_double("Rest mass mu (1 = timelike, 0 = null)", input.rest_mass);
    input.radial_direction = prompt_radial_direction();
    input.initial_radius_over_m = prompt_double("Initial radius r / M", input.initial_radius_over_m);
    input.escape_radius_over_m = prompt_double(
        "Configured escape radius r / M (must exceed initial radius)",
        input.escape_radius_over_m);

    const double spin_length = bh::kerr_spin_length(
        normalized_kerr_mass, session.black_hole.dimensionless_spin);
    const bh::KerrOrbit orbit{normalized_kerr_mass, spin_length, input.energy,
                              input.angular_momentum, input.rest_mass, input.radial_direction};
    const bh::Trajectory trajectory = bh::integrate_kerr(
        orbit, input.initial_radius_over_m * normalized_kerr_mass,
        session.trajectory_control.maximum_step, session.trajectory_control.maximum_steps,
        input.escape_radius_over_m * normalized_kerr_mass,
        session.trajectory_control.integration_control);

    std::cout << "Equatorial Kerr trajectory\n"
              << "  outer horizon r+:        "
              << bh::kerr_outer_horizon(normalized_kerr_mass, spin_length) << " M\n"
              << "  equatorial static limit: "
              << bh::kerr_static_limit_radius(
                     normalized_kerr_mass, spin_length, 1.57079632679489661923)
              << " M\n"
              << "  termination:             " << termination_name(trajectory.termination) << "\n";
    print_integration_diagnostics("trajectory", trajectory);
}

void run_interactive_penrose(InteractiveSession& session) {
    std::cout << "\nRestricted equatorial Penrose event input\n"
              << "  The shared black-hole mass and spin are already loaded.\n"
              << "  Enter only particle, split, and event-boundary parameters.\n";
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
    input.split_radius_over_m = prompt_double("Split radius r_split / M",
                                              input.split_radius_over_m);
    input.incoming_lz_over_m_m = prompt_double("Incoming angular momentum Lz / (m M)",
                                               input.incoming_lz_over_m_m);
    input.split_angle_rad = prompt_double("Split angle [radians]", input.split_angle_rad);

    const bh::EquatorialPenroseScenario scenario{
        normalized_kerr_mass,
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
    const bh::PenroseSplitParameters split{
        input.split_radius_over_m, input.incoming_lz_over_m_m, input.split_angle_rad};
    const bh::PenroseEventResult result = bh::evaluate_equatorial_penrose_event(scenario, split);
    print_penrose_event_result(scenario, result);
}

void run_interactive_toy_plasma(InteractiveSession& session) {
    std::cout << "\nReduced toy-plasma transport input\n"
              << "  This model uses the shared black-hole spin. Its reduced formula does not use mass.\n";
    ToyPlasmaState& input = session.toy_plasma;
    input.magnetic_field_tesla = prompt_double("Magnetic field [tesla]",
                                                input.magnetic_field_tesla);
    input.mass_density_kg_m3 = prompt_double("Mass density [kg / m^3]",
                                              input.mass_density_kg_m3);
    input.flow_area_m2 = prompt_double("Flow area [m^2]", input.flow_area_m2);
    input.duration_seconds = prompt_double("Duration [seconds]", input.duration_seconds);
    print_toy_plasma_result(bh::estimate_plasma_extraction(
        {input.magnetic_field_tesla, input.mass_density_kg_m3, input.flow_area_m2,
         session.black_hole.dimensionless_spin, input.duration_seconds}));
}

void run_interactive() {
    InteractiveSession session;
    std::cout << "Unified interactive simulation session\n"
              << "  Configure the black hole once. Every engine reuses that in-memory state.\n"
              << "  Kerr validates paths; Penrose evaluates split energy and delegates path checks to Kerr.\n";
    configure_shared_black_hole(session);

    while (true) {
        print_shared_session_header(session);
        std::cout << "  1. Update shared black-hole inputs\n"
                  << "  2. Update shared Kerr integration controls\n"
                  << "  3. Run Algebraic rotational-energy reservoir\n"
                  << "  4. Run Kerr trajectory validation\n"
                  << "  5. Run Penrose energy-extraction event\n"
                  << "  6. Run reduced toy-plasma transport\n"
                  << "  7. End session\n";
        switch (prompt_choice("Choose an action", 5, 1, 7)) {
        case 1:
            configure_shared_black_hole(session);
            break;
        case 2:
            configure_shared_trajectory_control(session);
            break;
        case 3:
            run_interactive_algebraic(session);
            break;
        case 4:
            run_interactive_kerr(session);
            break;
        case 5:
            run_interactive_penrose(session);
            break;
        case 6:
            run_interactive_toy_plasma(session);
            break;
        case 7:
            return;
        }
    }
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
        if (argc >= 2 && std::string_view(argv[1]) == "--toy-plasma") {
            run_toy_plasma(argc, argv);
            return 0;
        }
        if (argc == 2 && std::string_view(argv[1]) == "--interactive") {
            run_interactive();
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
