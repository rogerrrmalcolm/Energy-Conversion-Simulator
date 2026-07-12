# Algebraic Model Upgrade Guide

This guide keeps the algebraic model at milestone 1: a scalar, deterministic Kerr
rotational-energy bound. It does not model a Penrose event, extraction hardware,
or deliverable energy.

## 1. Fix library/header resolution in the editor

Configure the project with CMake before trusting IntelliSense:

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

The public headers live in `include/bh`. The CMake target already publishes that
directory through `target_include_directories`, so code should include headers as:

```cpp
#include "bh/algebraic_model.hpp"
```

If VS Code still underlines those includes, point the C/C++ extension at:

```text
build/compile_commands.json
```

or use the CMake Tools extension and select the configured `build` directory.

## 2. Keep the physics formula small and reusable

Implement these pure dimensionless helpers first:

```cpp
double irreducible_mass_fraction(double a_star);
double rotational_energy_fraction(double a_star);
```

They should validate `0 <= a_star < 1` and then calculate:

```text
M_irr / M = sqrt((1 + sqrt(1 - a_star^2)) / 2)
E_rot / E_mass = 1 - M_irr / M
```

This makes the SI boundary easier to test because the fractions do not depend on
kilograms, joules, or constants.

## 3. Use a named input instead of loose doubles

Prefer:

```cpp
struct RotationalEnergyInput {
    double mass_kg;
    double dimensionless_spin;
    SpinRange spin_uncertainty;
};
```

This makes it harder to mix up mass, spin, and uncertainty when the catalogue
loader is added. Require the uncertainty central value to match the selected
spin so observed estimates and scenario values stay auditable.

## 4. Return an auditable result

A useful milestone-1 result should include:

- total mass-energy in joules;
- irreducible mass in kilograms;
- irreducible mass fraction;
- theoretical rotational-energy reservoir in joules;
- rotational-energy fraction;
- lower and upper rotational-energy bounds from spin uncertainty;
- local sensitivity `dE_rot / da_star`.

Do not label any of these as extracted, captured, transmitted, or delivered
energy. Those belong to later model stages.

## 5. Test limits before adding features

Add tests for:

- `a_star = 0` gives zero rotational reservoir;
- near-extremal spin approaches `1 - 1 / sqrt(2)`;
- invalid mass and spin values throw;
- uncertainty bounds are ordered;
- higher spin gives a higher rotational-energy bound;
- spin sensitivity is positive for positive spin.

Run this after each change:

```powershell
cmake --build build
ctest --test-dir build --output-on-failure
```

## 6. Next improvement

The next clean step is catalogue validation: preserve source mass and spin values,
normalize accepted values, and create labelled spin scenarios only when the
catalogue lacks an observed spin. Do not invent observed spin values.
