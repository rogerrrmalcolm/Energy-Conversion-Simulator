#include "bh/penrose_scenario_io.hpp"

#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <unordered_set>

namespace bh {
namespace {
bool is_space(const char value) {
    return value == ' ' || value == '\t' || value == '\r';
}

std::string_view trim(std::string_view value) {
    while (!value.empty() && is_space(value.front())) {
        value.remove_prefix(1);
    }
    while (!value.empty() && is_space(value.back())) {
        value.remove_suffix(1);
    }
    return value;
}

std::string line_error(const std::size_t line_number, const std::string_view message) {
    return "scenario line " + std::to_string(line_number) + ": " + std::string(message);
}

double parse_double(const std::string_view value, const std::size_t line_number,
                    const std::string_view key) {
    double parsed{};
    const auto [end, error] = std::from_chars(
        value.data(), value.data() + value.size(), parsed, std::chars_format::general);
    if (error != std::errc{} || end != value.data() + value.size() || !std::isfinite(parsed)) {
        throw std::invalid_argument(line_error(
            line_number, "expected a finite numeric value for '" + std::string(key) + "'"));
    }
    return parsed;
}

std::size_t parse_size(const std::string_view value, const std::size_t line_number,
                       const std::string_view key) {
    std::size_t parsed{};
    const auto [end, error] = std::from_chars(value.data(), value.data() + value.size(), parsed);
    if (error != std::errc{} || end != value.data() + value.size()) {
        throw std::invalid_argument(
            line_error(line_number, "expected a whole-number value for '" + std::string(key) + "'"));
    }
    return parsed;
}
}  // namespace

EquatorialPenroseEventInput load_equatorial_penrose_event_input(
    const std::filesystem::path& path) {
    std::ifstream file(path);
    if (!file) {
        throw std::invalid_argument("could not open scenario file: " + path.string());
    }

    EquatorialPenroseEventInput input;
    std::unordered_set<std::string> seen_keys;
    bool has_version = false;
    std::size_t scenario_version{};
    std::string line;
    std::size_t line_number = 0;

    while (std::getline(file, line)) {
        ++line_number;
        const std::size_t comment = line.find('#');
        std::string_view content = line;
        if (comment != std::string::npos) {
            content = content.substr(0, comment);
        }
        content = trim(content);
        if (content.empty()) {
            continue;
        }

        const std::size_t separator = content.find('=');
        if (separator == std::string_view::npos) {
            throw std::invalid_argument(line_error(line_number, "expected key = value"));
        }

        const std::string_view key = trim(content.substr(0, separator));
        const std::string_view value = trim(content.substr(separator + 1));
        if (key.empty() || value.empty()) {
            throw std::invalid_argument(line_error(line_number, "key and value must both be present"));
        }
        if (!seen_keys.emplace(key).second) {
            throw std::invalid_argument(
                line_error(line_number, "duplicate key '" + std::string(key) + "'"));
        }

        if (key == "scenario_version") {
            scenario_version = parse_size(value, line_number, key);
            has_version = true;
        } else if (key == "black_hole_mass") {
            input.scenario.black_hole_mass = parse_double(value, line_number, key);
        } else if (key == "dimensionless_spin") {
            input.scenario.dimensionless_spin = parse_double(value, line_number, key);
        } else if (key == "parent_rest_mass") {
            input.scenario.parent_rest_mass = parse_double(value, line_number, key);
        } else if (key == "fragment_rest_mass") {
            input.scenario.fragment_rest_mass = parse_double(value, line_number, key);
        } else if (key == "incoming_specific_energy") {
            input.scenario.incoming_specific_energy = parse_double(value, line_number, key);
        } else if (key == "initial_radius_over_m") {
            input.scenario.initial_radius_over_m = parse_double(value, line_number, key);
        } else if (key == "escape_radius_over_m") {
            input.scenario.escape_radius_over_m = parse_double(value, line_number, key);
        } else if (key == "integration_step") {
            input.scenario.integration_step = parse_double(value, line_number, key);
        } else if (key == "max_integration_steps") {
            input.scenario.max_integration_steps = parse_size(value, line_number, key);
        } else if (key == "integration_absolute_tolerance") {
            input.scenario.integration_control.absolute_tolerance =
                parse_double(value, line_number, key);
        } else if (key == "integration_relative_tolerance") {
            input.scenario.integration_control.relative_tolerance =
                parse_double(value, line_number, key);
        } else if (key == "integration_minimum_step") {
            input.scenario.integration_control.minimum_step =
                parse_double(value, line_number, key);
        } else if (key == "residual_tolerance") {
            input.scenario.residual_tolerance = parse_double(value, line_number, key);
        } else if (key == "split_radius_over_m") {
            input.split.split_radius_over_m = parse_double(value, line_number, key);
        } else if (key == "incoming_lz_over_m_m") {
            input.split.incoming_lz_over_m_m = parse_double(value, line_number, key);
        } else if (key == "split_angle_rad") {
            input.split.split_angle_rad = parse_double(value, line_number, key);
        } else {
            throw std::invalid_argument(
                line_error(line_number, "unknown scenario key '" + std::string(key) + "'"));
        }
    }

    constexpr std::array<std::string_view, 17> required_keys{
        "scenario_version",
        "black_hole_mass",
        "dimensionless_spin",
        "parent_rest_mass",
        "fragment_rest_mass",
        "incoming_specific_energy",
        "initial_radius_over_m",
        "escape_radius_over_m",
        "integration_step",
        "max_integration_steps",
        "integration_absolute_tolerance",
        "integration_relative_tolerance",
        "integration_minimum_step",
        "residual_tolerance",
        "split_radius_over_m",
        "incoming_lz_over_m_m",
        "split_angle_rad"};
    for (const std::string_view key : required_keys) {
        if (!seen_keys.contains(std::string(key))) {
            throw std::invalid_argument("scenario is missing required key '" + std::string(key) + "'");
        }
    }
    if (!has_version || scenario_version != 1) {
        throw std::invalid_argument("unsupported scenario_version; expected 1");
    }

    return input;
}
}  // namespace bh
