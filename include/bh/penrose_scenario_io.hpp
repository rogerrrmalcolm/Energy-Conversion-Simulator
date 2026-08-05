#pragma once

#include "bh/dijkstra.hpp"
#include "bh/penrose_model.hpp"

#include <filesystem>

namespace bh {
struct EquatorialPenroseEventInput {
    EquatorialPenroseScenario scenario{};
    PenroseSplitParameters split{};
};

struct EquatorialPenroseDijkstraInput {
    EquatorialPenroseScenario scenario{};
    PenroseDijkstraSearchConfig search{};
};

[[nodiscard]] EquatorialPenroseEventInput load_equatorial_penrose_event_input(
    const std::filesystem::path& path);
[[nodiscard]] EquatorialPenroseDijkstraInput load_equatorial_penrose_dijkstra_input(
    const std::filesystem::path& path);
}
