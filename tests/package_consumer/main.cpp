#include "bh/algebraic_model.hpp"

int main() {
    const auto result = bh::rotational_energy(1.0, 0.5);
    return result.rotational_energy_joules > 0.0 ? 0 : 1;
}
