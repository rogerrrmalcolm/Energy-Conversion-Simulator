#pragma once

#include "bh/penrose_model.hpp"

#include <array>
#include <chrono>
#include <cstddef>
#include <limits>
#include <optional>
#include <stop_token>
#include <string>
#include <string_view>
#include <vector>

namespace bh {
using DijkstraGridKey = std::array<int, 3>;

struct DijkstraGridResult {
    bool found{};
    std::vector<DijkstraGridKey> parameter_adjustment_path;
};

// Basic deterministic Dijkstra search on a bounded three-coordinate integer grid.
[[nodiscard]] DijkstraGridResult find_dijkstra_grid_path(
    const DijkstraGridKey& start, const DijkstraGridKey& goal,
    const DijkstraGridKey& lower_bound, const DijkstraGridKey& upper_bound);

enum class PenroseSearchAlgorithm {
    dijkstra_h_zero
};

// A fixed Penrose scenario with a bounded, three-coordinate candidate grid.
// All values use the normalized units of PenroseSplitParameters. The baseline
// requires each edge cost to remain one and uses h = 0, so f = g.
struct PenroseDijkstraSearchConfig {
    PenroseSplitParameters start{};
    PenroseSplitParameters lower_bound{};
    PenroseSplitParameters upper_bound{};
    PenroseSplitParameters step{};
    double eta_target{};
    std::array<std::size_t, 3> edge_costs{1, 1, 1};
    std::size_t max_evaluations{std::numeric_limits<std::size_t>::max()};
    std::chrono::milliseconds time_budget{};
    PenroseSearchAlgorithm algorithm{PenroseSearchAlgorithm::dijkstra_h_zero};
};

enum class PenroseDijkstraSearchStatus {
    found_goal,
    best_feasible_below_target,
    no_solution_within_bounds,
    target_unattainable_under_model,
    node_budget_exhausted,
    time_budget_exhausted,
    cancelled,
    evaluation_failure
};

enum class PenroseDijkstraNodeStatus {
    outside_ergosphere,
    physics_invalid,
    captured_or_non_escaping,
    escaping_without_target,
    integration_failed,
    goal_feasible
};

// A parameter-space candidate, not a physical trajectory point.
struct PenroseDijkstraNode {
    DijkstraGridKey key{};
    PenroseSplitParameters split{};
    PenroseDijkstraNodeStatus status{PenroseDijkstraNodeStatus::physics_invalid};
    PenroseEventStatus event_status{PenroseEventStatus::physics_invalid};
    std::optional<DijkstraGridKey> parent_key{};
    DijkstraGridKey local_change{};
    std::size_t g_cost{};
    std::size_t h_cost{};
    std::size_t f_cost{};
    std::size_t discovery_order{};
    double input_energy{};
    double captured_energy{};
    double escaping_energy{};
    double eta_penrose{};
    double extracted_energy{};
    double maximum_normalized_residual{};
    TrajectoryTermination captured_termination{TrajectoryTermination::completed};
    TrajectoryTermination escaping_termination{TrajectoryTermination::completed};
};

struct PenroseDijkstraSearchDiagnostics {
    std::size_t candidates_in_domain{};
    std::size_t nodes_generated{};
    std::size_t nodes_evaluated{};
    std::size_t nodes_expanded{};
    std::size_t duplicate_nodes_skipped{};
    std::size_t outside_search_domain_neighbors{};
    std::size_t final_verification_evaluations{};
    std::size_t outside_ergosphere{};
    std::size_t physics_invalid{};
    std::size_t captured_or_non_escaping{};
    std::size_t escaping_without_target{};
    std::size_t integration_failed{};
    std::size_t goal_feasible{};
    std::chrono::nanoseconds elapsed{};
};

struct PenroseDijkstraSearchResult {
    PenroseDijkstraSearchStatus status{PenroseDijkstraSearchStatus::no_solution_within_bounds};
    // True when either a target-reaching goal or a completed-search fallback
    // candidate was selected.
    bool found{};
    bool target_reached{};
    // Ordered minimum-cost parameter-adjustment trace from the requested start
    // to the selected candidate.
    std::vector<PenroseDijkstraNode> parameter_adjustment_path;
    PenroseDijkstraSearchDiagnostics diagnostics{};
    // The only physical trajectories returned by the search: those of the
    // freshly verified selected event.
    PenroseEventResult selected_event{};
    std::string failure_message{};
};

enum class PenrosePhaseMapStatus {
    completed,
    node_budget_exhausted,
    time_budget_exhausted,
    cancelled,
    evaluation_failure
};

// A bounded, exhaustive deterministic scan. A completed result can report the
// greatest validated net extraction found within its declared grid, never a
// global physical maximum.
struct PenrosePhaseMapResult {
    PenrosePhaseMapStatus status{PenrosePhaseMapStatus::completed};
    bool complete{};
    std::vector<PenroseDijkstraNode> candidates;
    std::optional<PenroseDijkstraNode> best_validated_candidate{};
    PenroseEventResult best_event{};
    PenroseDijkstraSearchDiagnostics diagnostics{};
    std::string failure_message{};
};

[[nodiscard]] std::string_view penrose_dijkstra_node_status_name(
    PenroseDijkstraNodeStatus status);
[[nodiscard]] std::string_view penrose_dijkstra_search_status_name(
    PenroseDijkstraSearchStatus status);
[[nodiscard]] std::string_view penrose_phase_map_status_name(PenrosePhaseMapStatus status);
[[nodiscard]] std::string_view penrose_search_algorithm_name(PenroseSearchAlgorithm algorithm);

// Finds the minimum number of one-grid-step parameter adjustments needed to reach
// a physically feasible event whose Penrose efficiency reaches eta_target. If a
// completed search has no such goal, returns the greatest validated extraction
// in the grid, breaking extraction ties by minimum adjustment cost.
[[nodiscard]] PenroseDijkstraSearchResult find_penrose_dijkstra_path(
    const EquatorialPenroseScenario& scenario, const PenroseDijkstraSearchConfig& config,
    std::stop_token stop_token = {});

// Evaluates every candidate in deterministic key order. It is a phase-space
// diagnostic and a bounded maximum-within-grid check, not a graph path search.
[[nodiscard]] PenrosePhaseMapResult evaluate_penrose_phase_map(
    const EquatorialPenroseScenario& scenario, const PenroseDijkstraSearchConfig& config,
    std::stop_token stop_token = {});
}
