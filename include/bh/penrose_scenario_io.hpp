#pragma once

#include "bh/penrose_model.hpp"

#include <filesystem>

namespace bh {
struct EquatorialPenroseEventInput {
    EquatorialPenroseScenario scenario{};
    PenroseSplitParameters split{};
};

[[nodiscard]] EquatorialPenroseEventInput load_equatorial_penrose_event_input(
    const std::filesystem::path& path);
}
