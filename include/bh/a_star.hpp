#pragma once

#include <array>
#include <vector>

namespace bh {
using AStarGridKey = std::array<int, 3>;

struct AStarGridResult {
    bool found{};
    std::vector<AStarGridKey> parameter_adjustment_path;
};

// Basic deterministic A* baseline on a bounded three-coordinate integer grid.
// It deliberately uses h = 0 until a physics-specific admissible heuristic exists.
[[nodiscard]] AStarGridResult find_a_star_grid_path(
    const AStarGridKey& start, const AStarGridKey& goal,
    const AStarGridKey& lower_bound, const AStarGridKey& upper_bound);
}
