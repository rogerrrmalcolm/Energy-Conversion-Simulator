#include "bh/dijkstra.hpp"

#include <algorithm>
#include <array>
#include <cstddef>
#include <map>
#include <queue>

namespace bh {
namespace {
struct QueueEntry {
    int cost{};
    DijkstraGridKey key{};
};

struct Later {
    bool operator()(const QueueEntry& left, const QueueEntry& right) const {
        if (left.cost != right.cost) return left.cost > right.cost;
        return left.key > right.key;
    }
};

bool within(const DijkstraGridKey& key, const DijkstraGridKey& lower,
            const DijkstraGridKey& upper) {
    for (std::size_t index = 0; index < key.size(); ++index) {
        if (key[index] < lower[index] || key[index] > upper[index]) return false;
    }
    return true;
}
}  // namespace

DijkstraGridResult find_dijkstra_grid_path(const DijkstraGridKey& start,
                                            const DijkstraGridKey& goal,
                                            const DijkstraGridKey& lower,
                                            const DijkstraGridKey& upper) {
    if (!within(start, lower, upper) || !within(goal, lower, upper)) return {};

    std::priority_queue<QueueEntry, std::vector<QueueEntry>, Later> open;
    std::map<DijkstraGridKey, int> best_cost;
    std::map<DijkstraGridKey, DijkstraGridKey> parent;
    open.push({0, start});
    best_cost.emplace(start, 0);
    constexpr std::array<DijkstraGridKey, 6> changes{{{{-1, 0, 0}}, {{1, 0, 0}},
                                                        {{0, -1, 0}}, {{0, 1, 0}},
                                                        {{0, 0, -1}}, {{0, 0, 1}}}};

    while (!open.empty()) {
        const QueueEntry current = open.top();
        open.pop();
        if (best_cost.at(current.key) != current.cost) continue;
        if (current.key == goal) {
            DijkstraGridResult result{true, {}};
            for (DijkstraGridKey key = goal;; key = parent.at(key)) {
                result.parameter_adjustment_path.push_back(key);
                if (key == start) break;
            }
            std::reverse(result.parameter_adjustment_path.begin(),
                         result.parameter_adjustment_path.end());
            return result;
        }
        for (const DijkstraGridKey& change : changes) {
            DijkstraGridKey next = current.key;
            for (std::size_t index = 0; index < next.size(); ++index) next[index] += change[index];
            if (!within(next, lower, upper)) continue;
            const int next_cost = current.cost + 1;
            const auto [entry, inserted] = best_cost.try_emplace(next, next_cost);
            if (!inserted && next_cost >= entry->second) continue;
            entry->second = next_cost;
            parent.insert_or_assign(next, current.key);
            open.push({next_cost, next});
        }
    }
    return {};
}
}  // namespace bh
