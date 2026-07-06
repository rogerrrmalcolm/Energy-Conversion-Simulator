# Black Hole Energy Simulation — First Pass

This C++20 foundation separates four levels of physical modeling so that a cheap reference calculation is never confused with a trajectory proof or a plasma simulation. Trajectories use geometrized units (`G = c = 1`), metric signature `(-,+,+,+)`, Boyer–Lindquist coordinates for Kerr, and distances normalized by black-hole mass `M`; SI values appear only in the algebraic and toy plasma boundaries. No component claims that energy has been extracted from an observed black hole or that astrophysical energy delivery is technologically feasible.

## Model reports (three sentences each)

### 1. Algebraic rotational-energy model

This model converts an observed mass and dimensionless spin into total mass-energy, irreducible mass, and the theoretical Kerr rotational-energy reservoir. It is fast, deterministic, and useful for validating the zero-spin and near-extremal limits, but it does not model a Penrose event or imply that the reservoir can all be recovered. Its central tradeoff is simplicity versus physical detail: one calculation gives a rigorous theoretical bound, while saying nothing about trajectory feasibility, collection, or delivery losses.

### 2. Schwarzschild geodesic model

This reference integrator evolves an equatorial timelike test particle around a non-rotating black hole with a fourth-order Runge–Kutta stepper. It establishes understandable orbital behavior and a baseline against which rotating Kerr trajectories can be checked, including the `a -> 0` limit. Its tradeoff is that fixed steps and coordinate time near the horizon are easy to inspect but less robust than adaptive integration and horizon-regular coordinates.

### 3. Kerr geodesic model

This model evolves an equatorial timelike or null trajectory using the separated Kerr radial potential and conserved energy and axial angular momentum. It exposes frame dragging, the spin-dependent horizon, capture, and escape, which are essential foundations for a later deterministic Penrose split optimizer. This first pass intentionally does not switch direction at radial turning points, integrate off-equatorial motion, or adapt its step size, so it is an educational kernel rather than evidence of a 20.7% feasible event.

### 4. Plasma model

This transparent zero-dimensional model estimates magnetization, relativistic Alfvén speed, Poynting power, and an explicitly heuristic spin-coupling energy. It demonstrates how magnetic field strength, density, area, spin, and duration could enter a continuous electromagnetic extraction ledger without pretending to solve plasma dynamics. Its tradeoff is interpretability versus realism: the outputs are useful scaling scenarios, but a defensible extraction prediction requires general-relativistic magnetohydrodynamics (GRMHD), boundary conditions, field geometry, and numerical convergence studies.

## Build and run

```powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
.\build\black_hole_demo.exe
```

Those commands target single-configuration generators such as Ninja and place the executable directly in `build/`; use `./build/black_hole_demo` in that case. With a multi-configuration generator such as Visual Studio, omit `-DCMAKE_BUILD_TYPE=Release` and pass `--config Release` to the build, test, and install commands.

## Use as a CMake library

Install the library into a local prefix:

```powershell
cmake --install build --prefix install
```

Then consume it from another CMake project:

```cmake
find_package(BlackHoleModels 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE bh::models)
```

Configure the consuming project with `-DCMAKE_PREFIX_PATH="<repository>/install"`. The exported target supplies the public include directory and C++20 requirement automatically.

## Current numerical boundaries

- The Schwarzschild implementation is equatorial, timelike, fixed-step RK4 and stops just outside the coordinate-singular horizon.
- The Kerr implementation is equatorial (`Q=0`), uses a fixed radial direction and first-order stepping, and rejects a forbidden radial potential rather than locating a turning point.
- Both trajectory implementations use test particles and ignore radiation reaction and backreaction on the black hole.
- The plasma model is a scenario calculator, not particle-in-cell, MHD, or GRMHD software.

The next correctness milestone is adaptive error control plus invariant monitoring, followed by a hand-constructed momentum-conserving split with one captured negative-energy fragment and one verified escaping fragment.
