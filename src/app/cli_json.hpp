#pragma once

#include "bh/algebraic_model.hpp"
#include "bh/dijkstra.hpp"
#include "bh/penrose_scenario_io.hpp"
#include "bh/plasma_model.hpp"

#include <iosfwd>
#include <string_view>

namespace bh::cli_json {
void write_version(std::ostream& output, std::string_view version);
void write_algebraic(std::ostream& output, const RotationalEnergyResult& result);
void write_algebraic_range(std::ostream& output,
                           const RotationalEnergyRangeResult& result);
void write_toy_plasma(std::ostream& output, const PlasmaInput& input,
                      const PlasmaResult& result);
void write_penrose_event(std::ostream& output,
                         const EquatorialPenroseScenario& scenario,
                         const PenroseEventResult& result);
void write_penrose_search(std::ostream& output,
                          const EquatorialPenroseDijkstraInput& input,
                          const PenroseDijkstraSearchResult& result);
void write_penrose_phase_map(std::ostream& output,
                             const EquatorialPenroseDijkstraInput& input,
                             const PenrosePhaseMapResult& result);
}
