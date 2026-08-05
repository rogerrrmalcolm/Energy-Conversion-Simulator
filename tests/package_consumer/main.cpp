#include "bh/algebraic_model.hpp"
#include "bh/dijkstra.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"

int main() {
    const auto reservoir = bh::rotational_energy(1.0, 0.5);
    const double horizon = bh::kerr_outer_horizon(1.0, 0.5);
    const auto grid_path = bh::find_dijkstra_grid_path({0, 0, 0}, {1, 0, 0},
                                                        {0, 0, 0}, {1, 0, 0});
    return reservoir.rotational_energy_joules > 0.0 && horizon > 1.0 &&
                   bh::classical_penrose_efficiency_limit() > 0.2 && grid_path.found &&
                   bh::penrose_dijkstra_node_status_name(
                       bh::PenroseDijkstraNodeStatus::goal_feasible) == "goal_feasible" &&
                   bh::penrose_dijkstra_search_status_name(
                       bh::PenroseDijkstraSearchStatus::found_goal) == "found_goal" &&
                   bh::penrose_phase_map_status_name(
                       bh::PenrosePhaseMapStatus::completed) == "completed"
               ? 0
               : 1;
}
