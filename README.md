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

Search a bounded Penrose parameter grid with Dijkstra:

~~~powershell
.\build\black_hole_demo.exe --search-penrose .\scenarios\equatorial_penrose_dijkstra_reference.cfg
~~~

The reference search starts at `r_split/M = 1.09`, changes one radius-grid step to `1.10`, and reaches a validated event above its declared 4 percent efficiency target. Its trace is a record of parameter adjustments, not a particle trajectory.

Search for a physically valid event reaching 15 percent net efficiency across
the full `0.01 M` split-radius grid inside the equatorial ergosphere:

~~~powershell
.\build\black_hole_demo.exe --search-penrose .\scenarios\equatorial_penrose_dijkstra_15_percent.cfg
~~~

For the scenario's `a_star = 0.999`, the outer horizon is approximately
`1.04471018 M` and the equatorial static limit is `2 M`; therefore, the declared
radius grid is `1.05 M` through `1.99 M`. These are open physical boundaries,
so neither boundary itself is a candidate. The angular-momentum and split-angle
bounds remain independent exploratory choices. A failed search means no valid
15 percent event exists in this declared 2,375-candidate grid, not that such an
event is impossible throughout the continuous physical domain.

With the current restricted evaluator, the completed map of this grid finds no
15 percent goal. Its best validated candidate is `(1.14, 2.09, -1.98)` at
`9.065063%` net extraction. Because both non-radius coordinates lie on their
upper bounds, the next domain study should expand `incoming_lz_over_m_m` and
`split_angle_rad` while retaining explicit finite bounds.

Exhaustively map the same bounded candidate grid:

~~~powershell
.\build\black_hole_demo.exe --map-penrose .\scenarios\equatorial_penrose_dijkstra_reference.cfg
~~~

`--map-penrose` is intentionally different from Dijkstra. It evaluates every candidate within the declared bounds, then reports the greatest *validated net extraction found within that fully evaluated grid*. It does not claim a global physical maximum. The reference map evaluates 125 candidates; use its printed elapsed time as the local performance measurement for your machine.

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

### Deterministic Penrose parameter search

The Dijkstra layer keeps one declared `EquatorialPenroseScenario` fixed: black-hole mass and spin, particle masses, incoming specific energy, integration controls, residual tolerance, and escape radius do not change during a run. It varies only these normalized candidate parameters:

- `split_radius_over_m`;
- `incoming_lz_over_m_m`;
- `split_angle_rad`.

The search stores a node as three canonical integer grid indices, then derives the three floating-point split values from the configured lower bounds and steps. This prevents equivalent floating-point values from producing duplicate graph nodes. A node evaluation calls `evaluate_equatorial_penrose_event`; Dijkstra itself does not calculate a geodesic or Penrose split.

The only supported algorithm is `dijkstra_h_zero`: `h = 0`, `f = g`, and every one-coordinate, one-grid-step edge has a declared cost of exactly `1`. Therefore, the returned goal has the fewest parameter adjustments from the declared start among the searched candidates. The goal predicate requires a physically feasible event, positive net extracted energy, residuals within tolerance, and `eta_penrose >= search_eta_target`. A returned `parameter_adjustment_path` is not the path of a particle through spacetime. The physical incoming, captured, and escaping trajectories shown in the CLI come only from the freshly re-evaluated selected goal event.

The search has explicit terminal statuses: `found_goal`, `no_solution_within_bounds`, `target_unattainable_under_model`, `node_budget_exhausted`, `time_budget_exhausted`, `cancelled`, and `evaluation_failure`. A target at or above the ideal 20.710678 percent extremal single-split limit returns `target_unattainable_under_model` before evaluating candidates. This is separate from the roughly 29 percent theoretical rotational-reservoir limit.

The search is a threshold solver, not proof that its first answer has globally maximum extraction. A 15 percent target means `E_escape >= 1.15 * E_input`; if no bounded candidate reaches that target, `no_solution_within_bounds` is the correct result. Use `--map-penrose` only to make a bounded-grid maximum claim, or use a later constrained optimizer before making a broader maximum-energy claim.

### Resolution and performance guide

Use this sequence when improving a search scenario:

1. Run a coarse Dijkstra search and record its target, selected parameters, `g` cost, residual, node counts, and elapsed time.
2. Halve one or more grid steps while keeping the start and bounds exactly aligned to the new grid. Increase `search_max_evaluations` to cover the larger domain.
3. Repeat the search. A stable conclusion means both grids find a feasible target with similar validated parameters and energy; it is not proof of a global optimum.
4. Run `--map-penrose` only when the full grid fits the declared budget. Its maximum is valid only when the result status is `completed`.

The implementation is intentionally scalar and single-threaded. Geodesic integration dominates the workload, so the search first reduces memory pressure by caching compact status and ledger values keyed by canonical integer grid key. It does not retain full trajectory vectors for every evaluated node; it retains them only for the fresh final verification. Do not parallelize the priority queue or add SIMD until profiling a representative larger workload identifies a real bottleneck.

The original argument-based commands and versioned Penrose scenario files remain the reproducible interface for benchmarks and tests.

## Scenario files

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

`scenarios/equatorial_penrose_dijkstra_reference.cfg` contains the same fixed event fields plus this bounded-grid contract:

~~~text
search_eta_target
search_split_radius_lower
search_split_radius_upper
search_split_radius_step
search_incoming_lz_lower
search_incoming_lz_upper
search_incoming_lz_step
search_split_angle_lower
search_split_angle_upper
search_split_angle_step
search_radius_step_cost
search_incoming_lz_step_cost
search_split_angle_step_cost
search_max_evaluations
search_time_budget_ms
search_algorithm
~~~

The declared `split_*` values are the Dijkstra start candidate. Each lower-to-upper interval must be an exact multiple of its declared positive step, and the start must lie on that grid. The split-radius bounds must remain strictly inside the equatorial ergosphere. All three `search_*_step_cost` values must be `1`; this project intentionally rejects weighted costs until a defensible preference is defined. `search_max_evaluations` is required and positive. `search_time_budget_ms = 0` disables a wall-clock limit for a reproducible reference run. `search_algorithm` must be `dijkstra_h_zero`.

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
- deterministic bounded Dijkstra graph behavior, parent reconstruction, unit costs, duplicate suppression, fixed tie breaking, cancellation, explicit budgets, unreachable targets, and unattainable-limit rejection;
- physics-backed Penrose search with a fresh final-event verification, plus a bounded exhaustive phase map and coarse-versus-fine grid-resolution check;
- rejection of a split outside the ergosphere;
- loading, rejecting malformed scenarios, and executing the versioned reference scenario;
- CLI help, algebraic, uncertainty-range, valid-Penrose, and rejected-scenario workflows;
- installation, downstream CMake configuration, compilation, and execution of a package consumer;
- toy-plasma input validation, zero-spin coupling, overflow rejection, and causal Alfven speed.

## Deliberate limits and next work

- Kerr motion is equatorial only, has Q = 0, uses adaptive RK4, and does not yet continue through radial turning points. E and L_z are fixed input constants rather than separately evolved variables, so the reported radial-first-integral residual is the invariant check available to this restricted solver.
- The Penrose model is a neutral, idealized two-body test-particle split with identical daughter rest masses and an immediate inward/outward branch requirement. It has no charged-particle fields, radiation reaction, backreaction, collision model, nuclear fragmentation model, unequal daughter masses, or GRMHD.
- The plasma component is a reduced, ideal-MHD-inspired toy transport model with a heuristic spin factor, not an MHD or GRMHD solution.
- Dijkstra assigns unit cost to every one-step parameter adjustment. That is a transparent baseline for minimum-change search, not a physical impulse, elapsed-time, or energy cost. A heuristic is deliberately disabled because no admissible lower bound on efficiency improvement per step has been proven.
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
