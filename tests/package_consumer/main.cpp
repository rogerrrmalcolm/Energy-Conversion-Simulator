#include "bh/algebraic_model.hpp"
#include "bh/kerr_geodesic.hpp"
#include "bh/penrose_model.hpp"
#include "bh/penrose_scenario_io.hpp"

int main() {
    const auto reservoir = bh::rotational_energy(1.0, 0.5);
    const double horizon = bh::kerr_outer_horizon(1.0, 0.5);
    return reservoir.rotational_energy_joules > 0.0 && horizon > 1.0 &&
                   bh::classical_penrose_efficiency_limit() > 0.2
               ? 0
               : 1;
}
