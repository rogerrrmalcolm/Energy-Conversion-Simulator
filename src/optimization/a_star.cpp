#include "bh/a_star.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <queue>

namespace bh {
namespace {
struct QueueEntry {
    int f{};
    int g{};
    AStarGridKey key{};
};

struct Later {
    bool operator()(const QueueEntry& left, const QueueEntry& right) const {
        if (left.f != right.f) return left.f > right.f;
        if (left.g != right.g) return left.g > right.g;
        return left.key > right.key;
    }
};

bool within(const AStarGridKey& key, const AStarGridKey& lower,
            const AStarGridKey& upper) {
    for (std::size_t index = 0; index < key.size(); ++index) {
        if (key[index] < lower[index] || key[index] > upper[index]) return false;
    }
    return true;
}
}  // namespace

AStarGridResult find_a_star_grid_path(const AStarGridKey& start, const AStarGridKey& goal,
                                       const AStarGridKey& lower, const AStarGridKey& upper) {
    if (!within(start, lower, upper) || !within(goal, lower, upper)) return {};

    std::priority_queue<QueueEntry, std::vector<QueueEntry>, Later> open;
    std::map<AStarGridKey, int> best_cost;
    std::map<AStarGridKey, AStarGridKey> parent;
    open.push({0, 0, start});
    best_cost.emplace(start, 0);
    constexpr std::array<AStarGridKey, 6> changes{{{{-1, 0, 0}}, {{1, 0, 0}},
                                                     {{0, -1, 0}}, {{0, 1, 0}},
                                                     {{0, 0, -1}}, {{0, 0, 1}}}};

    while (!open.empty()) {
        const QueueEntry current = open.top();
        open.pop();
        if (best_cost.at(current.key) != current.g) continue;
        if (current.key == goal) {
            AStarGridResult result{true, {}};
            for (AStarGridKey key = goal;; key = parent.at(key)) {
                result.parameter_adjustment_path.push_back(key);
                if (key == start) break;
            }
            std::reverse(result.parameter_adjustment_path.begin(),
                         result.parameter_adjustment_path.end());
            return result;
        }
        for (const AStarGridKey& change : changes) {
            AStarGridKey next = current.key;
            for (std::size_t index = 0; index < next.size(); ++index) next[index] += change[index];
            if (!within(next, lower, upper)) continue;
            const int next_g = current.g + 1;
            const auto [entry, inserted] = best_cost.try_emplace(next, next_g);
            if (!inserted && next_g >= entry->second) continue;
            entry->second = next_g;
            parent.insert_or_assign(next, current.key);
            open.push({next_g, next_g, next});  // h = 0, therefore f = g.
        }
    }
    return {};
}
}  // namespace bh
