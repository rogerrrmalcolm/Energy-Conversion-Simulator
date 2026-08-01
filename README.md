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

E_rot is the theoretical rotational-energy reservoir, not energy guaranteed to be extracted by one event. The model supports ordered lower, central, and upper mass and spin values. Because the formula is monotonic in both positive mass and sub-extremal spin, it evaluates the two bounding corners `(M_low, a_low)` and `(M_high, a_high)`. The reported range is a bound calculation, not a statistical confidence interval unless the supplied input bounds already have that interpretation. Near extremality, the implementation evaluates `1 - a_star^2` in a cancellation-resistant form; exact `a_star = 1` remains outside this sub-extremal model.

### Schwarzschild reference trajectory

The code also retains a separate, fixed-step RK4 solver for equatorial timelike Schwarzschild motion. It validates the initial radial mass-shell relation, reports the radial-first-integral residual, and is tested against circular motion and analytic radial free-fall. It is a zero-spin validation reference; the Penrose evaluator uses the adaptive Kerr solver, not this component. It does not provide the full turning-point handling or adaptive control of the Kerr path.

### Equatorial Kerr geodesics

The trajectory kernel uses geometrized units, G = c = 1, Boyer-Lindquist coordinates, metric signature (-,+,+,+), and the equatorial restriction Q = 0. The user supplies a black-hole mass scale M, a Kerr spin length a = a_star * M, conserved energy E, axial angular momentum L_z, rest mass mu, and an initial radial direction.

~~~text
Delta(r) = r^2 - 2*M*r + a^2
r_plus   = M + sqrt(M^2 - a^2)

P(r) = E*(r^2 + a^2) - a*L_z
R(r) = P(r)^2 - Delta(r) * (mu^2*r^2 + (L_z - a*E)^2)
~~~

R(r) >= 0 is required for the radial state to be allowed. The solver rejects a state whose Boyer-Lindquist time component is not future-directed. The integrator uses adaptive RK4 step doubling: it compares one full step with two half steps, accepts the higher-accuracy state only when the normalized radial and azimuthal error is within tolerance, and records accepted/rejected-step diagnostics. It also reports a normalized residual of the radial first integral, `(Sigma * dr/dlambda)^2 = R(r)`, at every reported trajectory point. Boyer-Lindquist coordinate time is reported but excluded from step control because it becomes coordinate-singular at the horizon. The result returns explicit termination states for a horizon crossing, a configured escape radius, a requested target radius, a turning point, or an invalid state.

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
- four-momentum, mass-shell, energy, angular-momentum, geodesic-initialization, and radial-first-integral residuals are within the declared tolerance.

The escape radius is separate from the ergosphere boundary. The ergosphere determines where a Penrose split can produce a negative-energy fragment. Escape is confirmed only when the outward fragment reaches the independently configured large radius. The coordinate integrator stops immediately outside the Boyer-Lindquist horizon and records that event as crossed_horizon. This first baseline accepts only a direct split: the selected captured fragment must initially point inward and the selected escaping fragment must initially point outward. It does not yet model branches that reverse direction at a later turning point.

~~~text
eta_penrose = (E_escape - E_input) / E_input
E_extracted = max(0, E_escape - E_input)
~~~

Those event energies are normalized geometrized quantities. They are not joules and are not a claim about an observed astrophysical extraction.

The familiar approximately 20.71 percent classical Penrose efficiency is a restrictive ideal-limit check, not a result this engine claims for an arbitrary event. Because this model requires `a_star < 1`, the exact extremal limit is not an attainable target in this baseline. The reference scenario is deliberately below that limit.

### Reduced toy-plasma transport

The plasma component is a transparent 0-D, ideal-MHD-inspired transport scaling, not an MHD or GRMHD simulation. Given a magnetic field, mass density, flow area, black-hole spin, and duration, it calculates magnetization, a causal relativistic Alfven speed, a raw electromagnetic flux scaling through the declared surface, a visible heuristic spin-coupling factor, and the resulting outward electromagnetic power and energy.

~~~text
sigma = B^2 / (mu_0 * rho * c^2)
v_A   = c * sqrt(sigma / (1 + sigma))
P_EM  = (B^2 / mu_0) * v_A * area
P_out = P_EM * eta_spin_coupling
E_out = P_out * duration
~~~

The raw flux expression assumes a magnetic field perpendicular to a representative flow. `E_out` is only the toy model's outward electromagnetic energy through the configured surface; it is not energy proven to have been extracted from the black hole. The model has no field geometry, black-hole mass scaling, accretion solution, matter flux, collector capture, conversion, storage, transmission, or usable delivered energy.

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

Evaluate one reduced toy-plasma transport estimate:

~~~powershell
.\build\black_hole_demo.exe --toy-plasma 1.0 1e-8 10.0 0.9 2.0
~~~

The five values are magnetic field in tesla, mass density in kg/m^3, flow area in m^2, dimensionless spin, and duration in seconds.

### Interactive input mode

Run the interactive CLI when you want one in-memory simulation session without preparing a command line or Penrose scenario file:

~~~powershell
.\build\black_hole_demo.exe --interactive
~~~

The session asks for black-hole mass in kilograms and dimensionless spin once at startup, then keeps them in memory for every engine. The menu can explicitly update that shared black-hole state, configure shared Kerr integration controls, run the Algebraic reservoir, validate a Kerr trajectory, evaluate a Penrose event, run toy-plasma transport, or end the session.

Kerr and Penrose use the same selected black hole in normalized geometrized coordinates with `M = 1`; radii are entered as `r / M` and Penrose angular momentum as `Lz / (m M)`. This is a coordinate normalization of the shared physical black-hole mass, not a second black-hole input.

The Kerr menu asks only for an orbit's `E`, `Lz`, rest mass, direction, and boundaries. The Penrose menu asks only for the parent/fragment properties and split parameters; it reuses the shared spin and Kerr integration controls, then delegates incoming, capture, and escape path validation to the Kerr integrator. The toy-plasma baseline reuses shared spin; its reduced formula does not use mass.

A basic deterministic Dijkstra grid baseline now exists in `src/optimization/dijkstra.cpp`. It searches three bounded integer coordinates with six one-step neighbors and unit adjustment costs. It is deliberately only a tested graph-algorithm baseline: it does not yet call the Penrose evaluator or appear in the CLI. The future physics connection will keep the black-hole session fixed and vary only `split_radius_over_m`, `incoming_lz_over_m_m`, and `split_angle_rad`.

The original argument-based commands and versioned Penrose scenario files remain the reproducible interface for benchmarks and tests.

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
- allowed radial potential, future-directed momentum, radial-invariant diagnostics, horizon capture, and configured-radius escape;
- a local Penrose split with conservation residual checks;
- an analytic radial Schwarzschild-limit Kerr infall and adaptive-step rejection/refinement;
- scale invariance of the normalized Penrose efficiency;
- a deterministic bounded three-coordinate Dijkstra baseline, including shortest-path, tie-break, and bounds tests;
- rejection of a split outside the ergosphere;
- loading, rejecting malformed scenarios, and executing the versioned reference scenario;
- CLI help, algebraic, uncertainty-range, valid-Penrose, and rejected-scenario workflows;
- installation, downstream CMake configuration, compilation, and execution of a package consumer;
- toy-plasma input validation, zero-spin coupling, overflow rejection, and causal Alfven speed.

## Deliberate limits and next work

- Kerr motion is equatorial only, has Q = 0, uses adaptive RK4, and does not yet continue through radial turning points. E and L_z are fixed input constants rather than separately evolved variables, so the reported radial-first-integral residual is the invariant check available to this restricted solver.
- The Penrose model is a neutral, idealized two-body test-particle split with identical daughter rest masses and an immediate inward/outward branch requirement. It has no charged-particle fields, radiation reaction, backreaction, collision model, nuclear fragmentation model, unequal daughter masses, or GRMHD.
- The plasma component is a reduced, ideal-MHD-inspired toy transport model with a heuristic spin factor, not an MHD or GRMHD solution.
- The Dijkstra layer currently searches only a toy integer grid. Connecting it to Penrose candidates, physical goal states, scenario files, and CLI output is the next optimization milestone; its returned path must never be presented as a particle trajectory.
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
