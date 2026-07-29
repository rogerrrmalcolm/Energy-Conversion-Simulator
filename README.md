# Black Hole Energy Simulation

A C++20 scientific-computing foundation for three separate questions:

1. Given a black-hole mass and spin, what is its theoretical rotational-energy reservoir?
2. Given a declared equatorial Kerr test-particle event, does a local Penrose split conserve four-momentum, capture a negative-energy fragment, and let a positive-energy fragment escape?
3. Given a clearly labelled toy plasma scenario, how do a few magnetic and density inputs scale a reduced-model energy estimate?

The program does not claim to measure energy extracted from an observed black hole, and it does not claim that delivery of energy to Earth is feasible.

## Current physics

### Algebraic Kerr reservoir

At the SI boundary, the algebraic model accepts mass in kilograms and a dimensionless Kerr spin named a_star, where 0 <= a_star < 1.

~~~text
M_irr = M * sqrt((1 + sqrt(1 - a_star^2)) / 2)
E_mass = M * c^2
E_rot  = (M - M_irr) * c^2
~~~

E_rot is the theoretical rotational-energy reservoir, not energy guaranteed to be extracted by one event. The model supports ordered lower, central, and upper mass and spin values. It evaluates the same formula at each bound and reports asymmetric uncertainty relative to the central result.

### Equatorial Kerr geodesics

The trajectory kernel uses geometrized units, G = c = 1, Boyer-Lindquist coordinates, metric signature (-,+,+,+), and the equatorial restriction Q = 0. The user supplies a black-hole mass scale M, a Kerr spin length a = a_star * M, conserved energy E, axial angular momentum L_z, rest mass mu, and an initial radial direction.

~~~text
Delta(r) = r^2 - 2*M*r + a^2
r_plus   = M + sqrt(M^2 - a^2)

P(r) = E*(r^2 + a^2) - a*L_z
R(r) = P(r)^2 - Delta(r) * (mu^2*r^2 + (L_z - a*E)^2)
~~~

R(r) >= 0 is required for the radial state to be allowed. The integrator uses adaptive RK4 step doubling: it compares one full step with two half steps, accepts the higher-accuracy state only when the normalized radial and azimuthal error is within tolerance, and records accepted/rejected-step diagnostics. Boyer-Lindquist coordinate time is reported but excluded from step control because it becomes coordinate-singular at the horizon. The result returns explicit termination states for a horizon crossing, a configured escape radius, a requested target radius, a turning point, or an invalid state.

A turning point is localized at the first radial-potential root, R(r) = 0, in the selected radial direction. This restricted baseline stops at that event rather than reversing the particle's radial direction: an inward branch that turns before the horizon is not captured, and an outward branch that turns before the configured escape radius is not escaping.

### Restricted Penrose event

The event model evaluates one declared split; it does not search or optimize parameters. The split radius must be inside the equatorial ergosphere:

~~~text
r_plus < r_split < r_static_limit
r_static_limit = 2*M at the equator
~~~

At the split, the engine converts the incoming four-momentum into a local ZAMO frame, constructs an idealized neutral two-body split, and checks local four-momentum conservation and both mass-shell constraints. It converts each fragment back to Kerr conserved quantities, then integrates the captured and escaping geodesics separately.

An event is physically feasible in this restricted model only when:

- the captured fragment has negative conserved energy and crosses the horizon;
- the escaping fragment has positive conserved energy and reaches the configured large escape_radius;
- four-momentum, mass-shell, energy, angular-momentum, and geodesic-initialization residuals are within the declared tolerance.

The escape radius is separate from the ergosphere boundary. The ergosphere determines where a Penrose split can produce a negative-energy fragment. Escape is confirmed only when the outward fragment reaches the independently configured large radius. The coordinate integrator stops immediately outside the Boyer-Lindquist horizon and records that event as crossed_horizon.

~~~text
eta_penrose = (E_escape - E_input) / E_input
E_extracted = max(0, E_escape - E_input)
~~~

Those event energies are normalized geometrized quantities. They are not joules and are not a claim about an observed astrophysical extraction.

The familiar approximately 20.71 percent classical Penrose efficiency is a restrictive ideal-limit check, not a result this engine claims for an arbitrary event. The reference scenario is deliberately below that limit.

## Run the CLI

Configure and build:

~~~powershell
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
~~~

Those commands target a single-configuration generator such as Ninja. With a multi-configuration generator such as Visual Studio, omit CMAKE_BUILD_TYPE and add --config Release to the build, test, and install commands.

Evaluate an algebraic reservoir:

~~~powershell
.\build\black_hole_demo.exe --algebraic 1.98847e31 0.9
~~~

Evaluate an algebraic mass-and-spin range:

~~~powershell
.\build\black_hole_demo.exe --algebraic-range 1.590776e31 1.98847e31 2.386164e31 0.8 0.9 0.99
~~~

Evaluate one explicit Penrose event:

~~~powershell
.\build\black_hole_demo.exe --scenario .\scenarios\equatorial_penrose_reference.cfg
~~~

The reference file is a reproducible normalized test-particle case. It is intentionally labelled idealized; it is not an astrophysical observation.

## Scenario file

The scenarios/equatorial_penrose_reference.cfg file is a versioned key = value input. The loader rejects unknown keys, duplicate keys, malformed values, missing keys, and unsupported versions.

~~~text
scenario_version
black_hole_mass
dimensionless_spin
parent_rest_mass
fragment_rest_mass
incoming_specific_energy
initial_radius_over_m
escape_radius_over_m
integration_step
max_integration_steps
integration_absolute_tolerance
integration_relative_tolerance
integration_minimum_step
residual_tolerance
split_radius_over_m
incoming_lz_over_m_m
split_angle_rad
~~~

All values in this file are normalized geometrized inputs except the dimensionless spin and angle, which is measured in radians. The integration tolerances govern adaptive RK4 error control; residual_tolerance is the maximum normalized conservation or geodesic-initialization residual allowed for a physically feasible Penrose event.

## Validation

model_tests checks:

- zero-spin and near-extremal algebraic limits;
- mass and spin range ordering and uncertainty reporting;
- analytic spin sensitivity, floating-point uncertainty roundoff, and non-finite-result rejection;
- Schwarzschild timelike radial mass-shell validation and analytic radial free-fall time;
- Kerr horizon and static-limit formulas;
- allowed radial potential, horizon capture, and configured-radius escape;
- a local Penrose split with conservation residual checks;
- an analytic radial Schwarzschild-limit Kerr infall and adaptive-step rejection/refinement;
- scale invariance of the normalized Penrose efficiency;
- rejection of a split outside the ergosphere;
- loading, rejecting malformed scenarios, and executing the versioned reference scenario;
- CLI help, algebraic, uncertainty-range, valid-Penrose, and rejected-scenario workflows;
- installation, downstream CMake configuration, compilation, and execution of a package consumer;
- toy-plasma input validation, overflow rejection, and causal Alfven speed.

## Deliberate limits and next work

- Kerr motion is equatorial only, has Q = 0, uses adaptive RK4, and does not yet continue through radial turning points.
- The Penrose model is a neutral, idealized two-body test-particle split. It has no charged-particle fields, radiation reaction, backreaction, collision model, nuclear fragmentation model, or GRMHD.
- The plasma component is a reduced toy model, not MHD or GRMHD.
- A* parameter search is deliberately not implemented yet. Any future search must call this evaluator for each parameter set and must never be presented as a particle trajectory.
- Multithreading and SIMD are deliberately deferred until scalar correctness, integration error control, and invariant monitoring are expanded and profiled.

## CMake library

Install the library:

~~~powershell
cmake --install build --prefix install
~~~

Consume it from another CMake project:

~~~cmake
find_package(BlackHoleModels 0.1 CONFIG REQUIRED)
target_link_libraries(your_target PRIVATE bh::models)
~~~

Set CMAKE_PREFIX_PATH to the local install prefix when configuring the consuming project.
