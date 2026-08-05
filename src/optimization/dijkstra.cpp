#include "bh/dijkstra.hpp"

#include "bh/kerr_geodesic.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <limits>
#include <map>
#include <queue>
#include <stdexcept>
#include <string>
#include <utility>

namespace bh {
namespace {
struct GridQueueEntry {
    int cost{};
    DijkstraGridKey key{};
};

struct GridLater {
    bool operator()(const GridQueueEntry& left, const GridQueueEntry& right) const {
        if (left.cost != right.cost) {
            return left.cost > right.cost;
        }
        return left.key > right.key;
    }
};

struct SearchQueueEntry {
    std::size_t g_cost{};
    DijkstraGridKey key{};
    std::size_t discovery_order{};
};

struct SearchLater {
    bool operator()(const SearchQueueEntry& left, const SearchQueueEntry& right) const {
        // h = 0 in this Dijkstra mode, so f and g are identical. The key and
        // discovery order make every equal-cost decision reproducible.
        if (left.g_cost != right.g_cost) {
            return left.g_cost > right.g_cost;
        }
        if (left.key != right.key) {
            return left.key > right.key;
        }
        return left.discovery_order > right.discovery_order;
    }
};

constexpr std::array<DijkstraGridKey, 6> neighbor_changes{{{{-1, 0, 0}}, {{1, 0, 0}},
                                                             {{0, -1, 0}}, {{0, 1, 0}},
                                                             {{0, 0, -1}}, {{0, 0, 1}}}};

bool within(const DijkstraGridKey& key, const DijkstraGridKey& lower,
            const DijkstraGridKey& upper) {
    for (std::size_t index = 0; index < key.size(); ++index) {
        if (key[index] < lower[index] || key[index] > upper[index]) {
            return false;
        }
    }
    return true;
}

using ParameterValues = std::array<double, 3>;

struct CandidateGrid {
    ParameterValues lower{};
    ParameterValues step{};
    DijkstraGridKey upper{};
    DijkstraGridKey start{};
    std::size_t candidates_in_domain{};
};

struct CompactCandidateEvaluation {
    PenroseDijkstraNodeStatus status{PenroseDijkstraNodeStatus::physics_invalid};
    PenroseEventStatus event_status{PenroseEventStatus::physics_invalid};
    double input_energy{};
    double captured_energy{};
    double escaping_energy{};
    double eta_penrose{};
    double extracted_energy{};
    double maximum_normalized_residual{};
    TrajectoryTermination captured_termination{TrajectoryTermination::completed};
    TrajectoryTermination escaping_termination{TrajectoryTermination::completed};
};

ParameterValues parameter_values(const PenroseSplitParameters& split) {
    return {split.split_radius_over_m, split.incoming_lz_over_m_m, split.split_angle_rad};
}

PenroseSplitParameters split_from_values(const ParameterValues& values) {
    return {values[0], values[1], values[2]};
}

bool close_to_integer(const double value) {
    const double nearest = std::round(value);
    return std::abs(value - nearest) <= 1.0e-9 * std::max(1.0, std::abs(value));
}

std::size_t checked_multiply(const std::size_t left, const std::size_t right,
                             const char* const message) {
    if (right != 0 && left > std::numeric_limits<std::size_t>::max() / right) {
        throw std::invalid_argument(message);
    }
    return left * right;
}

CandidateGrid make_candidate_grid(const PenroseDijkstraSearchConfig& config) {
    if (!std::isfinite(config.eta_target) || config.eta_target < 0.0) {
        throw std::invalid_argument("eta_target must be finite and non-negative");
    }
    if (config.algorithm != PenroseSearchAlgorithm::dijkstra_h_zero) {
        throw std::invalid_argument("only dijkstra_h_zero is supported by the baseline search");
    }
    if (config.max_evaluations == 0) {
        throw std::invalid_argument("max_evaluations must be positive");
    }
    if (config.time_budget.count() < 0) {
        throw std::invalid_argument("time_budget must not be negative");
    }
    for (const std::size_t edge_cost : config.edge_costs) {
        if (edge_cost != 1) {
            throw std::invalid_argument(
                "Dijkstra baseline requires every one-step parameter edge cost to equal one");
        }
    }

    const ParameterValues start = parameter_values(config.start);
    const ParameterValues lower = parameter_values(config.lower_bound);
    const ParameterValues upper = parameter_values(config.upper_bound);
    const ParameterValues step = parameter_values(config.step);
    CandidateGrid grid{lower, step, {}, {}, 1};

    for (std::size_t index = 0; index < lower.size(); ++index) {
        if (!std::isfinite(start[index]) || !std::isfinite(lower[index]) ||
            !std::isfinite(upper[index]) || !std::isfinite(step[index])) {
            throw std::invalid_argument("Dijkstra search parameters must be finite");
        }
        if (lower[index] > upper[index] || step[index] <= 0.0) {
            throw std::invalid_argument("Dijkstra search bounds must be ordered with positive steps");
        }

        const double upper_index = (upper[index] - lower[index]) / step[index];
        const double start_index = (start[index] - lower[index]) / step[index];
        if (!close_to_integer(upper_index) || !close_to_integer(start_index)) {
            throw std::invalid_argument(
                "Dijkstra start and upper bounds must align with the declared grid step");
        }

        const double rounded_upper = std::round(upper_index);
        const double rounded_start = std::round(start_index);
        if (rounded_upper > static_cast<double>(std::numeric_limits<int>::max()) ||
            rounded_start < 0.0 || rounded_start > rounded_upper) {
            throw std::invalid_argument("Dijkstra start must be inside the declared grid bounds");
        }
        grid.upper[index] = static_cast<int>(rounded_upper);
        grid.start[index] = static_cast<int>(rounded_start);
        const std::size_t dimension_size = static_cast<std::size_t>(grid.upper[index]) + 1;
        grid.candidates_in_domain = checked_multiply(
            grid.candidates_in_domain, dimension_size,
            "Dijkstra candidate grid is too large to count safely");
    }
    if (grid.candidates_in_domain > max_penrose_search_nodes) {
        throw std::invalid_argument(
            "Penrose search window contains " +
            std::to_string(grid.candidates_in_domain) +
            " nodes; the scalar limit is " +
            std::to_string(max_penrose_search_nodes) +
            ". Increase a step size or submit a separate smaller parameter window");
    }
    return grid;
}

void validate_radius_bounds(const EquatorialPenroseScenario& scenario,
                            const CandidateGrid& grid) {
    if (!std::isfinite(scenario.black_hole_mass) || !std::isfinite(scenario.dimensionless_spin) ||
        scenario.black_hole_mass <= 0.0 || scenario.dimensionless_spin < 0.0 ||
        scenario.dimensionless_spin >= 1.0) {
        throw std::invalid_argument("search scenario has an invalid Kerr mass or spin");
    }

    const double mass = scenario.black_hole_mass;
    const double spin_length = kerr_spin_length(mass, scenario.dimensionless_spin);
    const double horizon = kerr_outer_horizon(mass, spin_length);
    const double static_limit = kerr_static_limit_radius(
        mass, spin_length, 1.57079632679489661923);
    const double lower_radius = grid.lower[0] * mass;
    const double upper_radius =
        (grid.lower[0] + static_cast<double>(grid.upper[0]) * grid.step[0]) * mass;
    if (!std::isfinite(lower_radius) || !std::isfinite(upper_radius) ||
        !(lower_radius > horizon && upper_radius < static_limit)) {
        throw std::invalid_argument(
            "search split-radius bounds must remain strictly inside the equatorial ergosphere");
    }
}

CandidateGrid validated_candidate_grid(const EquatorialPenroseScenario& scenario,
                                       const PenroseDijkstraSearchConfig& config) {
    CandidateGrid grid = make_candidate_grid(config);
    validate_radius_bounds(scenario, grid);
    return grid;
}

PenroseSplitParameters split_at(const CandidateGrid& grid, const DijkstraGridKey& key) {
    ParameterValues values{};
    for (std::size_t index = 0; index < values.size(); ++index) {
        values[index] = grid.lower[index] + static_cast<double>(key[index]) * grid.step[index];
    }
    return split_from_values(values);
}

PenroseDijkstraNodeStatus classify_node(const PenroseEventResult& event,
                                        const EquatorialPenroseScenario& scenario,
                                        const double eta_target) {
    switch (event.status) {
    case PenroseEventStatus::outside_ergosphere:
        return PenroseDijkstraNodeStatus::outside_ergosphere;
    case PenroseEventStatus::physics_invalid:
        return PenroseDijkstraNodeStatus::physics_invalid;
    case PenroseEventStatus::captured_or_non_escaping:
        return PenroseDijkstraNodeStatus::captured_or_non_escaping;
    case PenroseEventStatus::integration_failed:
        return PenroseDijkstraNodeStatus::integration_failed;
    case PenroseEventStatus::physically_feasible:
        return event.extracted_energy > 0.0 && event.eta_penrose >= eta_target &&
                       event.maximum_normalized_residual <= scenario.residual_tolerance
                   ? PenroseDijkstraNodeStatus::goal_feasible
                   : PenroseDijkstraNodeStatus::escaping_without_target;
    }
    return PenroseDijkstraNodeStatus::physics_invalid;
}

bool has_validated_positive_extraction(const PenroseEventResult& event,
                                       const EquatorialPenroseScenario& scenario) {
    return event.status == PenroseEventStatus::physically_feasible &&
           event.extracted_energy > 0.0 &&
           event.maximum_normalized_residual <= scenario.residual_tolerance;
}

CompactCandidateEvaluation compact_evaluation(const PenroseEventResult& event,
                                              const EquatorialPenroseScenario& scenario,
                                              const double eta_target) {
    return {classify_node(event, scenario, eta_target),
            event.status,
            event.input_energy,
            event.captured_energy,
            event.escaping_energy,
            event.eta_penrose,
            event.extracted_energy,
            event.maximum_normalized_residual,
            event.captured_trajectory.termination,
            event.escaping_trajectory.termination};
}

PenroseDijkstraNode make_node(const CandidateGrid& grid, const DijkstraGridKey& key,
                              const CompactCandidateEvaluation& evaluation,
                              const std::size_t g_cost,
                              const std::size_t discovery_order) {
    PenroseDijkstraNode node;
    node.key = key;
    node.split = split_at(grid, key);
    node.status = evaluation.status;
    node.event_status = evaluation.event_status;
    node.g_cost = g_cost;
    node.h_cost = 0;
    node.f_cost = g_cost;
    node.discovery_order = discovery_order;
    node.input_energy = evaluation.input_energy;
    node.captured_energy = evaluation.captured_energy;
    node.escaping_energy = evaluation.escaping_energy;
    node.eta_penrose = evaluation.eta_penrose;
    node.extracted_energy = evaluation.extracted_energy;
    node.maximum_normalized_residual = evaluation.maximum_normalized_residual;
    node.captured_termination = evaluation.captured_termination;
    node.escaping_termination = evaluation.escaping_termination;
    return node;
}

void record_status(PenroseDijkstraSearchDiagnostics& diagnostics,
                   const PenroseDijkstraNodeStatus status) {
    switch (status) {
    case PenroseDijkstraNodeStatus::outside_ergosphere:
        ++diagnostics.outside_ergosphere;
        return;
    case PenroseDijkstraNodeStatus::physics_invalid:
        ++diagnostics.physics_invalid;
        return;
    case PenroseDijkstraNodeStatus::captured_or_non_escaping:
        ++diagnostics.captured_or_non_escaping;
        return;
    case PenroseDijkstraNodeStatus::escaping_without_target:
        ++diagnostics.escaping_without_target;
        return;
    case PenroseDijkstraNodeStatus::integration_failed:
        ++diagnostics.integration_failed;
        return;
    case PenroseDijkstraNodeStatus::goal_feasible:
        ++diagnostics.goal_feasible;
        return;
    }
}

bool approximately_equal(const double left, const double right) {
    if (!std::isfinite(left) || !std::isfinite(right)) {
        return left == right;
    }
    const double scale = std::max({1.0, std::abs(left), std::abs(right)});
    return std::abs(left - right) <= 1.0e-12 * scale;
}

bool matches_fresh_selection(const CompactCandidateEvaluation& cached,
                             const PenroseEventResult& fresh,
                             const EquatorialPenroseScenario& scenario,
                             const double eta_target) {
    const CompactCandidateEvaluation rechecked = compact_evaluation(fresh, scenario, eta_target);
    return cached.status == rechecked.status &&
           cached.event_status == rechecked.event_status &&
           cached.captured_termination == rechecked.captured_termination &&
           cached.escaping_termination == rechecked.escaping_termination &&
           approximately_equal(cached.input_energy, rechecked.input_energy) &&
           approximately_equal(cached.captured_energy, rechecked.captured_energy) &&
           approximately_equal(cached.escaping_energy, rechecked.escaping_energy) &&
           approximately_equal(cached.eta_penrose, rechecked.eta_penrose) &&
           approximately_equal(cached.extracted_energy, rechecked.extracted_energy) &&
           approximately_equal(cached.maximum_normalized_residual,
                               rechecked.maximum_normalized_residual);
}

bool better_fallback(const CompactCandidateEvaluation& candidate,
                     const std::size_t candidate_cost,
                     const DijkstraGridKey& candidate_key,
                     const CompactCandidateEvaluation& current_best,
                     const std::size_t current_best_cost,
                     const DijkstraGridKey& current_best_key) {
    if (!approximately_equal(candidate.eta_penrose, current_best.eta_penrose)) {
        return candidate.eta_penrose > current_best.eta_penrose;
    }
    if (!approximately_equal(candidate.extracted_energy, current_best.extracted_energy)) {
        return candidate.extracted_energy > current_best.extracted_energy;
    }
    if (candidate_cost != current_best_cost) {
        return candidate_cost < current_best_cost;
    }
    return candidate_key < current_best_key;
}

bool better_extraction(const PenroseDijkstraNode& candidate,
                       const PenroseDijkstraNode& current_best) {
    if (candidate.extracted_energy != current_best.extracted_energy) {
        return candidate.extracted_energy > current_best.extracted_energy;
    }
    if (candidate.eta_penrose != current_best.eta_penrose) {
        return candidate.eta_penrose > current_best.eta_penrose;
    }
    return candidate.key < current_best.key;
}

bool time_budget_exhausted(const std::chrono::steady_clock::time_point start,
                           const PenroseDijkstraSearchConfig& config) {
    return config.time_budget.count() > 0 &&
           std::chrono::steady_clock::now() - start >= config.time_budget;
}

template <typename Result>
Result with_elapsed(Result result, const std::chrono::steady_clock::time_point start) {
    result.diagnostics.elapsed = std::chrono::steady_clock::now() - start;
    return result;
}

PenroseDijkstraNode reconstruct_node(
    const CandidateGrid& grid, const DijkstraGridKey& key,
    const std::map<DijkstraGridKey, CompactCandidateEvaluation>& evaluations,
    const std::map<DijkstraGridKey, std::size_t>& best_cost,
    const std::map<DijkstraGridKey, std::size_t>& discovery_order,
    const std::map<DijkstraGridKey, DijkstraGridKey>& parent) {
    PenroseDijkstraNode node = make_node(
        grid, key, evaluations.at(key), best_cost.at(key), discovery_order.at(key));
    const auto parent_entry = parent.find(key);
    if (parent_entry != parent.end()) {
        node.parent_key = parent_entry->second;
        for (std::size_t index = 0; index < node.local_change.size(); ++index) {
            node.local_change[index] = key[index] - parent_entry->second[index];
        }
    }
    return node;
}

std::vector<PenroseDijkstraNode> reconstruct_path(
    const CandidateGrid& grid, const DijkstraGridKey& selected,
    const std::map<DijkstraGridKey, CompactCandidateEvaluation>& evaluations,
    const std::map<DijkstraGridKey, std::size_t>& best_cost,
    const std::map<DijkstraGridKey, std::size_t>& discovery_order,
    const std::map<DijkstraGridKey, DijkstraGridKey>& parent) {
    std::vector<PenroseDijkstraNode> path;
    for (DijkstraGridKey key = selected;; key = parent.at(key)) {
        path.push_back(reconstruct_node(
            grid, key, evaluations, best_cost, discovery_order, parent));
        if (key == grid.start) {
            break;
        }
    }
    std::reverse(path.begin(), path.end());
    return path;
}
}  // namespace

DijkstraGridResult find_dijkstra_grid_path(const DijkstraGridKey& start,
                                            const DijkstraGridKey& goal,
                                            const DijkstraGridKey& lower,
                                            const DijkstraGridKey& upper) {
    if (!within(start, lower, upper) || !within(goal, lower, upper)) {
        return {};
    }

    std::priority_queue<GridQueueEntry, std::vector<GridQueueEntry>, GridLater> open;
    std::map<DijkstraGridKey, int> best_cost;
    std::map<DijkstraGridKey, DijkstraGridKey> parent;
    open.push({0, start});
    best_cost.emplace(start, 0);
    while (!open.empty()) {
        const GridQueueEntry current = open.top();
        open.pop();
        if (best_cost.at(current.key) != current.cost) {
            continue;
        }
        if (current.key == goal) {
            DijkstraGridResult result{true, {}};
            for (DijkstraGridKey key = goal;; key = parent.at(key)) {
                result.parameter_adjustment_path.push_back(key);
                if (key == start) {
                    break;
                }
            }
            std::reverse(result.parameter_adjustment_path.begin(),
                         result.parameter_adjustment_path.end());
            return result;
        }
        for (const DijkstraGridKey& change : neighbor_changes) {
            DijkstraGridKey next = current.key;
            for (std::size_t index = 0; index < next.size(); ++index) {
                next[index] += change[index];
            }
            if (!within(next, lower, upper)) {
                continue;
            }
            const int next_cost = current.cost + 1;
            const auto [entry, inserted] = best_cost.try_emplace(next, next_cost);
            if (!inserted && next_cost >= entry->second) {
                continue;
            }
            entry->second = next_cost;
            parent.insert_or_assign(next, current.key);
            open.push({next_cost, next});
        }
    }
    return {};
}

PenroseSearchWindowSummary describe_penrose_search_window(
    const EquatorialPenroseScenario& scenario,
    const PenroseDijkstraSearchConfig& config) {
    const CandidateGrid grid = validated_candidate_grid(scenario, config);
    return {{{static_cast<std::size_t>(grid.upper[0]) + 1,
              static_cast<std::size_t>(grid.upper[1]) + 1,
              static_cast<std::size_t>(grid.upper[2]) + 1}},
            grid.candidates_in_domain};
}

void require_complete_penrose_search_window(
    const EquatorialPenroseScenario& scenario,
    const PenroseDijkstraSearchConfig& config) {
    const PenroseSearchWindowSummary window =
        describe_penrose_search_window(scenario, config);
    if (config.max_evaluations < window.candidates) {
        throw std::invalid_argument(
            "complete Penrose search requires max_evaluations >= " +
            std::to_string(window.candidates) +
            " so no candidate node in the declared window can be skipped");
    }
    if (config.time_budget.count() != 0) {
        throw std::invalid_argument(
            "complete Penrose search requires time_budget = 0; use cancellation to stop an"
            " intentionally abandoned run");
    }
}

std::string_view penrose_dijkstra_node_status_name(const PenroseDijkstraNodeStatus status) {
    switch (status) {
    case PenroseDijkstraNodeStatus::outside_ergosphere:
        return "outside_ergosphere";
    case PenroseDijkstraNodeStatus::physics_invalid:
        return "physics_invalid";
    case PenroseDijkstraNodeStatus::captured_or_non_escaping:
        return "captured_or_non_escaping";
    case PenroseDijkstraNodeStatus::escaping_without_target:
        return "escaping_without_target";
    case PenroseDijkstraNodeStatus::integration_failed:
        return "integration_failed";
    case PenroseDijkstraNodeStatus::goal_feasible:
        return "goal_feasible";
    }
    return "unknown";
}

std::string_view penrose_dijkstra_search_status_name(const PenroseDijkstraSearchStatus status) {
    switch (status) {
    case PenroseDijkstraSearchStatus::found_goal:
        return "found_goal";
    case PenroseDijkstraSearchStatus::best_feasible_below_target:
        return "best_feasible_below_target";
    case PenroseDijkstraSearchStatus::no_solution_within_bounds:
        return "no_solution_within_bounds";
    case PenroseDijkstraSearchStatus::target_unattainable_under_model:
        return "target_unattainable_under_model";
    case PenroseDijkstraSearchStatus::node_budget_exhausted:
        return "node_budget_exhausted";
    case PenroseDijkstraSearchStatus::time_budget_exhausted:
        return "time_budget_exhausted";
    case PenroseDijkstraSearchStatus::cancelled:
        return "cancelled";
    case PenroseDijkstraSearchStatus::evaluation_failure:
        return "evaluation_failure";
    }
    return "unknown";
}

std::string_view penrose_phase_map_status_name(const PenrosePhaseMapStatus status) {
    switch (status) {
    case PenrosePhaseMapStatus::completed:
        return "completed";
    case PenrosePhaseMapStatus::node_budget_exhausted:
        return "node_budget_exhausted";
    case PenrosePhaseMapStatus::time_budget_exhausted:
        return "time_budget_exhausted";
    case PenrosePhaseMapStatus::cancelled:
        return "cancelled";
    case PenrosePhaseMapStatus::evaluation_failure:
        return "evaluation_failure";
    }
    return "unknown";
}

std::string_view penrose_search_algorithm_name(const PenroseSearchAlgorithm algorithm) {
    switch (algorithm) {
    case PenroseSearchAlgorithm::dijkstra_h_zero:
        return "dijkstra_h_zero";
    }
    return "unknown";
}

PenroseDijkstraSearchResult find_penrose_dijkstra_path(
    const EquatorialPenroseScenario& scenario, const PenroseDijkstraSearchConfig& config,
    const std::stop_token stop_token) {
    const auto start_time = std::chrono::steady_clock::now();
    PenroseDijkstraSearchResult result;
    const CandidateGrid grid = validated_candidate_grid(scenario, config);
    result.diagnostics.candidates_in_domain = grid.candidates_in_domain;

    if (config.eta_target >= classical_penrose_efficiency_limit()) {
        result.status = PenroseDijkstraSearchStatus::target_unattainable_under_model;
        result.failure_message =
            "the requested efficiency reaches the unattainable extremal single-split limit";
        return with_elapsed(std::move(result), start_time);
    }
    if (stop_token.stop_requested()) {
        result.status = PenroseDijkstraSearchStatus::cancelled;
        return with_elapsed(std::move(result), start_time);
    }

    constexpr DijkstraGridKey lower{0, 0, 0};
    std::priority_queue<SearchQueueEntry, std::vector<SearchQueueEntry>, SearchLater> open;
    std::map<DijkstraGridKey, std::size_t> best_cost;
    std::map<DijkstraGridKey, DijkstraGridKey> parent;
    std::map<DijkstraGridKey, std::size_t> discovery_order;
    std::map<DijkstraGridKey, CompactCandidateEvaluation> evaluations;
    std::optional<DijkstraGridKey> best_fallback_key;
    std::size_t next_discovery_order = 1;
    open.push({0, grid.start, 0});
    best_cost.emplace(grid.start, 0);
    discovery_order.emplace(grid.start, 0);
    result.diagnostics.nodes_generated = 1;

    while (!open.empty()) {
        if (stop_token.stop_requested()) {
            result.status = PenroseDijkstraSearchStatus::cancelled;
            return with_elapsed(std::move(result), start_time);
        }
        if (time_budget_exhausted(start_time, config)) {
            result.status = PenroseDijkstraSearchStatus::time_budget_exhausted;
            return with_elapsed(std::move(result), start_time);
        }

        const SearchQueueEntry current = open.top();
        open.pop();
        if (best_cost.at(current.key) != current.g_cost) {
            continue;
        }

        auto evaluation = evaluations.find(current.key);
        if (evaluation == evaluations.end()) {
            if (result.diagnostics.nodes_evaluated >= config.max_evaluations) {
                result.status = PenroseDijkstraSearchStatus::node_budget_exhausted;
                return with_elapsed(std::move(result), start_time);
            }
            try {
                const PenroseEventResult event =
                    evaluate_equatorial_penrose_event(scenario, split_at(grid, current.key));
                const auto [entry, inserted] = evaluations.emplace(
                    current.key, compact_evaluation(event, scenario, config.eta_target));
                (void)inserted;
                evaluation = entry;
                ++result.diagnostics.nodes_evaluated;
                record_status(result.diagnostics, evaluation->second.status);
                if (has_validated_positive_extraction(event, scenario) &&
                    evaluation->second.status ==
                        PenroseDijkstraNodeStatus::escaping_without_target &&
                    (!best_fallback_key ||
                     better_fallback(evaluation->second, current.g_cost, current.key,
                                     evaluations.at(*best_fallback_key),
                                     best_cost.at(*best_fallback_key),
                                     *best_fallback_key))) {
                    best_fallback_key = current.key;
                }
            } catch (const std::exception& error) {
                result.status = PenroseDijkstraSearchStatus::evaluation_failure;
                result.failure_message = "candidate evaluator threw: " + std::string(error.what());
                return with_elapsed(std::move(result), start_time);
            }
        }

        if (evaluation->second.status == PenroseDijkstraNodeStatus::goal_feasible) {
            PenroseEventResult fresh_goal;
            try {
                fresh_goal = evaluate_equatorial_penrose_event(
                    scenario, split_at(grid, current.key));
                ++result.diagnostics.final_verification_evaluations;
            } catch (const std::exception& error) {
                result.status = PenroseDijkstraSearchStatus::evaluation_failure;
                result.failure_message =
                    "fresh goal verification threw: " + std::string(error.what());
                return with_elapsed(std::move(result), start_time);
            }
            if (evaluation->second.status != PenroseDijkstraNodeStatus::goal_feasible ||
                !matches_fresh_selection(
                    evaluation->second, fresh_goal, scenario, config.eta_target)) {
                result.status = PenroseDijkstraSearchStatus::evaluation_failure;
                result.failure_message =
                    "fresh goal verification did not match the cached candidate evaluation";
                return with_elapsed(std::move(result), start_time);
            }

            result.status = PenroseDijkstraSearchStatus::found_goal;
            result.found = true;
            result.target_reached = true;
            result.selected_event = std::move(fresh_goal);
            result.parameter_adjustment_path = reconstruct_path(
                grid, current.key, evaluations, best_cost, discovery_order, parent);
            return with_elapsed(std::move(result), start_time);
        }

        ++result.diagnostics.nodes_expanded;
        for (const DijkstraGridKey& change : neighbor_changes) {
            DijkstraGridKey next = current.key;
            for (std::size_t index = 0; index < next.size(); ++index) {
                next[index] += change[index];
            }
            if (!within(next, lower, grid.upper)) {
                ++result.diagnostics.outside_search_domain_neighbors;
                continue;
            }

            const std::size_t next_cost = current.g_cost + 1;
            const auto [entry, inserted] = best_cost.try_emplace(next, next_cost);
            if (!inserted && next_cost >= entry->second) {
                ++result.diagnostics.duplicate_nodes_skipped;
                continue;
            }
            entry->second = next_cost;
            parent.insert_or_assign(next, current.key);
            if (inserted) {
                discovery_order.emplace(next, next_discovery_order++);
                ++result.diagnostics.nodes_generated;
            }
            open.push({next_cost, next, discovery_order.at(next)});
        }
    }

    if (!best_fallback_key) {
        result.status = PenroseDijkstraSearchStatus::no_solution_within_bounds;
        return with_elapsed(std::move(result), start_time);
    }

    PenroseEventResult fresh_fallback;
    try {
        fresh_fallback = evaluate_equatorial_penrose_event(
            scenario, split_at(grid, *best_fallback_key));
        ++result.diagnostics.final_verification_evaluations;
    } catch (const std::exception& error) {
        result.status = PenroseDijkstraSearchStatus::evaluation_failure;
        result.failure_message =
            "fresh fallback verification threw: " + std::string(error.what());
        return with_elapsed(std::move(result), start_time);
    }
    const CompactCandidateEvaluation& cached_fallback = evaluations.at(*best_fallback_key);
    if (cached_fallback.status != PenroseDijkstraNodeStatus::escaping_without_target ||
        !has_validated_positive_extraction(fresh_fallback, scenario) ||
        !matches_fresh_selection(
            cached_fallback, fresh_fallback, scenario, config.eta_target)) {
        result.status = PenroseDijkstraSearchStatus::evaluation_failure;
        result.failure_message =
            "fresh fallback verification did not match the cached candidate evaluation";
        return with_elapsed(std::move(result), start_time);
    }

    result.status = PenroseDijkstraSearchStatus::best_feasible_below_target;
    result.found = true;
    result.target_reached = false;
    result.selected_event = std::move(fresh_fallback);
    result.parameter_adjustment_path = reconstruct_path(
        grid, *best_fallback_key, evaluations, best_cost, discovery_order, parent);
    return with_elapsed(std::move(result), start_time);
}

PenrosePhaseMapResult evaluate_penrose_phase_map(
    const EquatorialPenroseScenario& scenario, const PenroseDijkstraSearchConfig& config,
    const std::stop_token stop_token) {
    const auto start_time = std::chrono::steady_clock::now();
    PenrosePhaseMapResult result;
    const CandidateGrid grid = validated_candidate_grid(scenario, config);
    result.diagnostics.candidates_in_domain = grid.candidates_in_domain;
    const std::size_t reserve_size = std::min(
        {grid.candidates_in_domain, config.max_evaluations, static_cast<std::size_t>(1'000'000)});
    result.candidates.reserve(reserve_size);

    std::size_t discovery_order = 0;
    for (std::size_t radius_index = 0;
         radius_index <= static_cast<std::size_t>(grid.upper[0]); ++radius_index) {
        for (std::size_t angular_momentum_index = 0;
             angular_momentum_index <= static_cast<std::size_t>(grid.upper[1]);
             ++angular_momentum_index) {
            for (std::size_t angle_index = 0;
                 angle_index <= static_cast<std::size_t>(grid.upper[2]); ++angle_index) {
                if (stop_token.stop_requested()) {
                    result.status = PenrosePhaseMapStatus::cancelled;
                    return with_elapsed(std::move(result), start_time);
                }
                if (time_budget_exhausted(start_time, config)) {
                    result.status = PenrosePhaseMapStatus::time_budget_exhausted;
                    return with_elapsed(std::move(result), start_time);
                }
                if (result.diagnostics.nodes_evaluated >= config.max_evaluations) {
                    result.status = PenrosePhaseMapStatus::node_budget_exhausted;
                    return with_elapsed(std::move(result), start_time);
                }

                const DijkstraGridKey key{
                    static_cast<int>(radius_index),
                    static_cast<int>(angular_momentum_index),
                    static_cast<int>(angle_index)};
                PenroseEventResult event;
                try {
                    event = evaluate_equatorial_penrose_event(scenario, split_at(grid, key));
                } catch (const std::exception& error) {
                    result.status = PenrosePhaseMapStatus::evaluation_failure;
                    result.failure_message = "candidate evaluator threw: " + std::string(error.what());
                    return with_elapsed(std::move(result), start_time);
                }

                const CompactCandidateEvaluation evaluation =
                    compact_evaluation(event, scenario, config.eta_target);
                PenroseDijkstraNode node = make_node(grid, key, evaluation, 0, discovery_order++);
                ++result.diagnostics.nodes_generated;
                ++result.diagnostics.nodes_evaluated;
                record_status(result.diagnostics, node.status);
                if (has_validated_positive_extraction(event, scenario) &&
                    (!result.best_validated_candidate ||
                     better_extraction(node, *result.best_validated_candidate))) {
                    result.best_validated_candidate = node;
                }
                result.candidates.push_back(std::move(node));
            }
        }
    }

    result.complete = true;
    result.status = PenrosePhaseMapStatus::completed;
    if (!result.best_validated_candidate) {
        return with_elapsed(std::move(result), start_time);
    }

    try {
        result.best_event = evaluate_equatorial_penrose_event(
            scenario, result.best_validated_candidate->split);
        ++result.diagnostics.final_verification_evaluations;
    } catch (const std::exception& error) {
        result.complete = false;
        result.status = PenrosePhaseMapStatus::evaluation_failure;
        result.failure_message = "fresh best-candidate verification threw: " + std::string(error.what());
        return with_elapsed(std::move(result), start_time);
    }
    if (!has_validated_positive_extraction(result.best_event, scenario) ||
        !approximately_equal(result.best_event.extracted_energy,
                             result.best_validated_candidate->extracted_energy) ||
        !approximately_equal(result.best_event.eta_penrose,
                             result.best_validated_candidate->eta_penrose)) {
        result.complete = false;
        result.status = PenrosePhaseMapStatus::evaluation_failure;
        result.failure_message =
            "fresh best-candidate verification did not match the phase-map evaluation";
    }
    return with_elapsed(std::move(result), start_time);
}
}  // namespace bh
