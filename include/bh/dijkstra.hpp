#pragma once

#include <array>
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
}
