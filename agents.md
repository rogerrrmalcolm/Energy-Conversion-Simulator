# Black Hole Energy Simulation - Agent Guide

## 1. Project mission

Build a high-performance C++20 and Python application that ingests a catalogue of observed black holes, validates the available measurements, estimates each object's theoretical rotational-energy reservoir, solves a deterministic Penrose-process optimization problem for configurable atomic-nucleus projectiles, and visualizes how hypothetical deliverable energy compares with regional and worldwide energy demand.

The project is also a systems-programming and scientific-computing platform. Multithreading, concurrency, SIMD, data-oriented design, L1/L2 cache behavior, and measured latency are first-class requirements, not optional polish. Every optimized implementation must be compared with a correct scalar reference implementation, and Python orchestration must call the same validated C++ engine rather than reimplementing the physics.

The intended audience includes quantitative developers, performance-oriented C++ engineers, scientific-computing teams, and big-technology infrastructure engineers. Optimize the project for mathematical transparency, reproducibility, data quality, numerical correctness, concurrency safety, and defensible performance measurements.

## 2. Required scientific claim

Use this description consistently:

> Given observed mass and spin estimates, calculate a black hole's theoretical rotational-energy reservoir, optimize an idealized Penrose extraction event for a labelled particle species, estimate a hypothetical usable-energy chain, and compare it with cited regional and worldwide energy demand.

Do not claim that the program measures energy that has actually been extracted from an observed black hole. A catalogue normally contains measured or inferred properties such as mass, spin, distance, and uncertainty. It does not contain a measurement of energy removed by a Penrose process.

Do not claim that delivering energy from an astrophysical black hole to Earth is currently feasible. The collection and delivery stages are transparent scenario models whose assumptions can be changed by the user.

Do not say that a particle is shot out of the event horizon. The incident particle splits or interacts inside the ergosphere; one negative-energy fragment may cross the horizon, while the positive-energy fragment must remain outside the horizon and escape to a configured large radius.

Atomic nuclei are configurable incident-particle scenarios, not claims about naturally observed Penrose events. Support known elements and isotopes through a versioned particle-species dataset containing atomic number, mass number, charge, rest mass, provenance, and stability metadata; do not invent nuclear properties or treat a bare nucleus as neutral matter.

Keep these values separate throughout the engine, result schema, and user interface:

1. total black-hole mass-energy;
2. theoretical rotational-energy reservoir;
3. energy carried to infinity by the optimized escaping fragment;
4. energy intercepted by a hypothetical collector;
5. energy remaining after conversion, storage, and transmission;
6. energy hypothetically available at the selected comparison boundary;
7. cited regional or worldwide demand that the delivered energy is compared against.

Never label one of these quantities as another.

## 3. Sequential end-to-end workflow

The application must implement and present the following stages in order.

### Stage 1: ingest a black-hole catalogue

Read a versioned CSV or Parquet dataset through a schema-aware loader. Support these fields when available:

- stable object identifier and display name;
- black-hole mass and mass unit;
- dimensionless spin estimate `a_star`;
- distance and distance unit;
- lower and upper uncertainty bounds for mass and spin;
- object class, such as stellar-mass or supermassive;
- measurement method, source citation, and catalogue version.

Mass is essential for an energy estimate. Spin is essential for a specific rotational-energy estimate. Distance is useful for context and delivery scenarios but does not change the intrinsic rotational-energy reservoir.

Do not invent a spin value when the catalogue does not provide one. Represent missing values explicitly. The program may calculate clearly labelled spin scenarios, such as `a_star = 0.5`, `0.9`, and `0.99`, but must not present those scenarios as observations.

The loader must report malformed rows, unsupported units, missing required fields, duplicate identifiers, and provenance. Preserve the source values alongside normalized values so conversions can be audited.

### Stage 2: validate and normalize the records

Convert each accepted record into a strongly typed internal representation. Validate at least:

- `M > 0`;
- `0 <= a_star < 1` for the classical sub-extremal Kerr model;
- finite numeric values;
- ordered and physically valid uncertainty bounds;
- supported and dimensionally correct units;
- non-negative distance when supplied.

Use geometrized units (`G = c = 1`) inside the relativistic physics engine. Convert to and from SI only at explicit boundaries. Document the metric signature, coordinate system, constants, and conversion rules.

Each input row must finish this stage in exactly one state:

- accepted as an observation;
- accepted as a labelled scenario because an observed value is missing;
- rejected with structured validation errors.

### Stage 3: calculate the rotational-energy reservoir

For black-hole mass `M` and dimensionless spin `a_star`, calculate the irreducible mass:

```text
M_irr = M * sqrt((1 + sqrt(1 - a_star^2)) / 2)
```

Then calculate the theoretical rotational-energy reservoir:

```text
E_rot = (M - M_irr) * c^2
```

Also calculate total mass-energy:

```text
E_mass = M * c^2
```

`E_rot` is a global theoretical reservoir, not energy automatically obtainable in one event. For a nearly extremal Kerr black hole, the ideal rotational reservoir approaches about 29% of its mass-energy. Treat that as a limiting validation case, not an engineering prediction.

Propagate input uncertainty into a range or local sensitivity for `E_rot`. Keep an observed estimate distinct from a hypothetical spin scenario.

### Stage 4: solve a deterministic Penrose-process problem

Model an incident particle entering the Kerr ergosphere and splitting into two fragments. Enforce local four-momentum conservation. A physically successful event requires:

- one negative-energy fragment, relative to infinity, crossing the horizon;
- one positive-energy fragment escaping to a configured large radius;
- timelike or null four-momenta as required by the chosen model;
- valid Kerr geodesics;
- conserved energy `E`, axial angular momentum `L_z`, and Carter constant `Q` within documented tolerances;
- a small four-momentum and numerical integration residual.

Define event efficiency as:

```text
eta_penrose = (E_escape - E_input) / E_input
E_extracted = max(0, E_escape - E_input)
```

`E_escape` includes the incident energy supplied to the event. `E_extracted` is the net gain attributed to loss of black-hole rotational energy. Use `E_extracted`, not total `E_escape`, as the input to the usable-energy calculation.

The main research problem is deterministic constrained trajectory design:

> Maximize escaping energy subject to four-momentum conservation, causal momenta, capture of the negative-energy fragment, escape of the positive-energy fragment, and the selected Kerr geodesic equations.

Begin with equatorial geodesics in Boyer-Lindquist coordinates. Add full `r`, `theta`, `phi`, and `t` evolution only after the restricted solver is validated. Candidate solution methods include multiple shooting, direct collocation, sequential quadratic programming, continuation, adaptive subdivision, and a bounded A* search over a discretized candidate-parameter graph. An A* search is an optimization-layer tool: it searches combinations of input and split parameters, then calls the physics evaluator for each combination. It must never replace geodesic integration or be described as the particle's physical path through spacetime. Cross-check local solutions using a second method or bounded search over the low-dimensional feasible domain.

Report feasibility, primal and dual residuals, Karush-Kuhn-Tucker residuals where applicable, integration error, termination reason, and whether the result is known only to be a local optimum. A visually convincing trajectory is not proof of correctness.

The projectile selector may instantiate a proton, alpha particle, or a nucleus from the versioned isotope table. Nuclear species determine rest mass and charge, but the first validated geodesic baseline must remain a test-particle model; electromagnetic forces, ionization state, fragmentation, radiation reaction, and plasma coupling must be enabled only by models that explicitly implement and report them.

Implement the physics ladder in this order and never present a lower rung as a higher one:

1. algebraic Kerr reservoir and local energy-ledger model;
2. equatorial geodesics and event detection;
3. full Kerr geodesics with conserved invariants;
4. local four-momentum-conserving Penrose split optimization;
5. charged-particle motion in documented electromagnetic fields;
6. relativistic plasma and, ultimately, validated GRMHD modeling.

### Stage 5: convert extracted energy into a usable-energy scenario

Apply an explicit efficiency chain rather than one unexplained percentage:

```text
E_delivered = E_extracted
         * eta_capture
         * eta_conversion
         * eta_storage
         * eta_transmission
```

Where:

- `E_extracted` is the net rotational-energy gain from the modeled event or sequence of events, excluding supplied incident energy;
- `eta_capture` is the fraction intercepted by a hypothetical collector;
- `eta_conversion` is the fraction converted into a usable carrier;
- `eta_storage` accounts for storage losses;
- `eta_transmission` accounts for losses to the selected hypothetical delivery boundary.

All efficiencies must be in `[0, 1]`, individually configurable, and visible in the output. Do not hide them in source code. Provide best-case, nominal, and conservative scenario files, but label all three as hypothetical.

Distance must not be represented by a fabricated efficiency formula. If a delivery model is implemented, state its transport mechanism, geometry, travel time, collector area, beam divergence or inverse-square assumptions, and losses. Otherwise, use distance only as context and explicitly say that the delivery physics is not modeled.

#### Plasma fidelity contract

A coefficient-only plasma estimate is permitted as an educational baseline but must be named `toy-plasma` or `reduced-plasma`, not realistic plasma. A model may be labelled realistic only when it solves documented relativistic MHD or GRMHD equations on a stated grid, includes magnetic-field geometry and boundary/initial conditions, demonstrates conservation and resolution convergence, and is validated against analytic tests or a trusted independent code.

The plasma layer estimates energy transported outward through a defined surface; it does not by itself prove energy was intercepted or made usable. Keep these measurements distinct:

```text
E_outward_flux     = integral of outward matter and electromagnetic energy flux
E_collector_input  = E_outward_flux * eta_capture
E_usable           = E_collector_input * eta_conversion * eta_storage * eta_transmission
```

Report matter-energy flux and electromagnetic Poynting flux separately, along with numerical dissipation, boundary flux balance, grid resolution, timestep policy, and divergence-control error. Python may configure, analyze, and visualize plasma runs, but the performance-critical finite-volume or finite-difference update belongs in the tested C++ engine.

### Stage 6: compare the result with regional and worldwide demand

Maintain a separate, versioned energy-demand dataset containing:

- metric name, such as electricity demand or total final energy use;
- value and unit;
- geographic scope;
- reporting year;
- authoritative source URL or citation;
- retrieval date.

Do not mix electricity consumption, primary energy, and total final energy consumption. Let the user choose a metric and geographic scope, including a country, region, or the world.

Calculate perspectives such as:

```text
years_of_demand = E_delivered / annual_selected_demand
average_power    = E_delivered / delivery_duration
```

Also show the assumed delivery duration, equivalent continuous power, fraction of annual demand, number of households only when supported by a cited definition, and the difference between theoretical reservoir and usable energy. Time-sensitive demand values belong in versioned data files, not hard-coded UI strings.

### Stage 7: present auditable results

For every black hole, display:

- original and normalized catalogue values;
- observed values versus assumed scenarios;
- mass-energy and rotational-energy reservoir;
- incident energy, total escaping energy, net extracted energy, and Penrose efficiency;
- solver status and conservation residuals;
- every collection and delivery efficiency;
- hypothetical usable delivered energy;
- selected geographic scope, comparison metric, year, source, and result;
- uncertainty or sensitivity ranges;
- runtime and execution backend.

The visualization should show the horizon and ergosphere, incoming/captured/escaping trajectories, selected element or isotope, a live energy ledger, constraint residuals, and an efficiency waterfall from `E_rot` to `E_delivered`. It must visually distinguish the ergosphere from the event horizon and never animate a fragment emerging from inside the horizon. Phase-space maps and Pareto frontiers are more useful than distributions of randomly generated events.

### Stage 8: process the catalogue efficiently

After the scalar implementation is correct, execute the workflow through three distinct performance layers:

1. **Multithreading:** distribute independent black holes, scenarios, or parameter tiles across a fixed worker pool.
2. **Concurrency:** overlap catalogue decoding, validation, computation, and result serialization through bounded queues and explicit backpressure.
3. **SIMD:** evaluate batches of identical arithmetic operations within each worker using a structure-of-arrays layout.

The application must allow each layer to be enabled independently so its behavior and speedup can be measured.

## 4. Multithreading requirements

Use a fixed-size thread pool rather than creating one thread per black hole. A reasonable starting point is C++20 `std::jthread` with cooperative cancellation through `std::stop_token`.

Parallel work units may include:

- separate catalogue records;
- separate observed/scenario variants of one record;
- independent regions of a deterministic parameter sweep;
- phase-space tiles;
- finite-difference or sensitivity evaluations when dependencies permit.

Follow these rules:

- solver instances and scratch buffers are thread-local;
- immutable catalogue and configuration data may be shared;
- results are written to preassigned indices or merged deterministically;
- exceptions are captured and reported with their work item;
- cancellation does not leave partially published results;
- thread count is configurable, including `1` for reference testing;
- avoid nested parallelism and oversubscription;
- separate frequently written per-thread data to reduce false sharing;
- do not use a mutex around the whole numerical kernel.

Test the same workload with `1`, `2`, `4`, and up to the available hardware thread count. Report throughput, latency, CPU utilization, parallel efficiency, and scaling. Explain where Amdahl's law limits speedup.

## 5. Concurrency requirements

Implement a bounded producer/consumer pipeline:

```text
catalogue reader -> validator/batcher -> compute workers -> ordered result writer
```

The stages may overlap, but the output must remain deterministic. Assign every input a monotonically increasing sequence number and reorder completed work before final serialization when necessary.

Use concurrency constructs deliberately:

- `std::mutex` and `std::condition_variable` for an initial bounded queue;
- `std::jthread` and `std::stop_token` for lifetime and cancellation;
- `std::latch` or `std::barrier` only where phase coordination is genuinely needed;
- atomics for small counters or state flags, not as a default replacement for ownership;
- RAII for threads, locks, queue shutdown, and output resources.

The bounded queue must provide backpressure so a fast reader cannot exhaust memory while compute workers are busy. Define shutdown behavior for normal completion, cancellation, parser failure, solver failure, and writer failure. Avoid detached threads.

Concurrency tests must cover empty input, one record, queue capacity one, more workers than records, cancellation, worker exceptions, slow output, and repeated execution under ThreadSanitizer on a supported toolchain.

## 6. SIMD requirements

Keep a scalar kernel as the correctness oracle. Only vectorize after profiling identifies a suitable hot loop.

Good SIMD candidates include:

- Kerr metric and radial-potential evaluation for several states;
- irreducible-mass and rotational-energy calculations across records;
- integration-stage arithmetic for trajectories sharing an execution shape;
- constraint and residual evaluation across a parameter tile;
- unit conversions and efficiency-chain calculations.

Use a structure-of-arrays representation for hot batches:

```text
masses[]
spins[]
radii[]
energies[]
angular_momenta[]
```

This is preferable to an array of large heterogeneous objects when adjacent SIMD lanes need the same field.

Implementation progression:

1. write and test a scalar kernel;
2. inspect compiler vectorization reports;
3. help auto-vectorization with contiguous data, simple loops, alignment, and non-aliasing where valid;
4. add an explicit portable SIMD abstraction if the toolchain supports one reliably;
5. use platform intrinsics only when measurement justifies the portability cost;
6. handle remainder lanes and non-multiple batch sizes;
7. compare scalar and SIMD results using documented floating-point tolerances.

Adaptive geodesics may cause lane divergence because different trajectories need different step counts. Do not force poorly matched trajectories into one vector. Group work by integration regime, use fixed-size kernels where scientifically acceptable, compact active lanes, or restrict SIMD to uniform subkernels such as metric and residual evaluation.

Record ISA, vector width, compiler, flags, alignment, batch size, and numerical mode in benchmark output. Do not enable unsafe floating-point transformations without testing their effect on invariants and optimizer decisions.

## 7. Combined execution model

The preferred CPU design is hierarchical:

```text
bounded input pipeline
        |
        v
fixed thread pool
        |
        v
one structure-of-arrays batch per worker
        |
        v
SIMD numerical kernels
        |
        v
ordered result aggregation
```

Thread-level parallelism operates across independent batches. SIMD operates within a batch. The pipeline overlaps I/O and compute. Avoid assigning the same parallelism to multiple layers because that causes oversubscription and noisy benchmarks.

Support these execution modes from the CLI:

- `scalar-single-thread`;
- `scalar-multithread`;
- `simd-single-thread`;
- `simd-multithread`;
- `concurrent-pipeline`, using the selected compute backend.

Every mode must consume the same validated inputs and produce equivalent ordered results within documented tolerances.

## 8. C++ and Python architecture

Target C++20 or newer with CMake for the numerical engine, and a supported modern Python version for orchestration, analysis, notebooks, and visualization. Keep physics, numerical algorithms, scheduling, I/O, demand comparisons, and rendering independent; expose the stable C++ API to Python through a thin binding such as pybind11 only after the C++ result schema and ownership model are tested.

```text
src/
  domain/        Strong types, catalogue records, scenarios, result schemas
  io/            Catalogue, isotope, demand, configuration, and result I/O
  physics/       Algebraic model, Kerr geometry, particles, plasma, units
  integrators/   Scalar geodesics, adaptive stepping, event detection
  optimization/  Candidate graphs, A* search, constraints, shooting/collocation, continuation, sensitivity
  simd/          Structure-of-arrays batches and vectorized kernels
  concurrency/   Bounded queues, thread pool, cancellation, pipeline
  analytics/     Phase maps, Pareto frontiers, uncertainty and sensitivity
  rendering/     Trajectories, energy waterfall, comparison views
  app/           C++ CLI and service entry points
python/
  bindings/      Thin access to versioned C++ APIs
  analysis/      Reproducible analysis and sensitivity workflows
  visualization/ Interactive trajectories, ledgers, and demand views
tests/
benchmarks/
data/catalogues/
data/particles/
data/demand/
scenarios/
docs/
```

Prefer value types, RAII, explicit ownership, `std::span` or ranges for non-owning views, and strong quantity types at SI boundaries. Avoid global mutable state. Keep solver state, thread count, batch width, tolerances, and scenario assumptions explicit.

Potential dependencies must be justified. Suitable C++ candidates include Eigen, Boost.Odeint, Catch2 or GoogleTest, Google Benchmark, `nlohmann/json`, pybind11, and Arrow/Parquet or HDF5 when dataset scale warrants them; suitable Python tools include NumPy, SciPy, pandas or Polars, Plotly, and pytest. Python must not sit inside a per-step geodesic or optimization hot loop, and crossing the language boundary must operate on batches or coarse-grained jobs.

### Cache and latency engineering

Treat sub-second latency as a workload-specific service-level objective, never as an unsupported claim about full-fidelity optimization or GRMHD. Define named latency classes, for example: algebraic catalogue query, cached trajectory lookup, single validated trajectory, bounded parameter tile, and full plasma simulation; report median, p95, and p99 latency separately for each class.

Use data-oriented layouts for hot kernels:

- structure-of-arrays batches for mass, spin, radius, energy, angular momentum, charge, and species mass;
- compact hot records separated from citations, names, strings, and visualization metadata;
- contiguous, aligned buffers sized through measurement rather than assumed cache sizes;
- per-worker scratch arenas reused between jobs to avoid allocator traffic;
- tiles chosen to keep the active metric, state, and residual arrays within measured L1 or L2 working sets where possible;
- padding or ownership partitioning for frequently written counters to prevent false sharing;
- precomputed immutable constants and lookup tables only when versioned and numerically verified.

Measure L1/L2 miss rates, branch misses, instructions per cycle, allocation counts, and memory bandwidth with platform-appropriate profilers. A cache optimization is accepted only when it improves a representative benchmark without changing scalar-reference results outside documented tolerances.

## 9. Deterministic optimization and analysis

Do not use Monte Carlo as the central analysis method. Prefer:

- constrained nonlinear optimization;
- continuation across spin, split radius, and incident orbit;
- phase-space maps of escape, capture, and forbidden regions;
- bifurcation detection near separatrices and unstable photon orbits;
- analytic, automatic-differentiation, or adjoint sensitivities;
- Pareto frontiers for escaping energy, escape time, required impulse, and delivery assumptions;
- interval bounds or adaptive deterministic subdivision for low-dimensional global checks.

An advanced version may formulate the event as an optimal-control problem. The decision variables are the incident state, split event, and fragment momenta; the dynamics are Kerr geodesics; the objective is energy delivered to infinity. This provides genuine overlap with optimization and calibration systems used by quantitative firms.

### A* parameter-space search: primary algorithmic design

The project may include an A* search layer as a deterministic way to explore a bounded, discretized space of Penrose-event candidate parameters. This is an educational and engineering feature that adds graph algorithms, state-space design, priority queues, heuristic design, reproducibility, and search diagnostics to the project.

The A* layer has one strict boundary:

> A* searches a graph of candidate parameter sets. The Kerr physics engine evaluates each set. The resulting A* path is a record of parameter adjustments, not the physical trajectory of a particle through the ergosphere.

The physical trajectory is obtained only by integrating the selected candidate's geodesic(s). Do not render the A* path as a line moving through the black-hole scene. A visualization may show the physical geodesic separately and may show the A* path as a parameter-space table, plot, or graph.

#### A* search objective

The baseline A* objective is:

> Find a physically feasible candidate that reaches a configured net-energy-extraction target with the lowest parameter-adjustment cost from a declared starting candidate.

For every candidate, calculate:

```text
eta_penrose = (E_escape - E_input) / E_input
E_extracted = max(0, E_escape - E_input)
```

The goal predicate must be explicit:

```text
goal(candidate) =
    candidate is in the configured parameter domain
    AND the split conserves four-momentum within tolerance
    AND the captured fragment has valid negative energy at infinity and crosses the horizon
    AND the escaping fragment remains outside the horizon and reaches the configured escape radius
    AND geodesic invariants and integration residuals are within tolerance
    AND eta_penrose >= eta_target
```

Use the term **net extracted energy**, not "contracted energy." `E_escape` includes energy supplied by the incident particle; `E_extracted` is the net gain attributed to the idealized extraction event.

This objective deliberately gives A* a well-defined minimum-cost problem. A* does not automatically find the globally highest-energy candidate. If the user asks for "the maximum energy," the program must use one of these accurate descriptions:

1. Run a bounded exhaustive search and return the highest validated `E_extracted` found.
2. Run a deterministic optimizer or branch-and-bound method and report its termination and optimality status.
3. Run A* toward one or more increasing `eta_target` values, then report the highest target reached and the best validated candidate encountered. This is a threshold search, not proof of the global maximum.

Never claim that a first A* goal is the global physical optimum unless a separate exhaustive or mathematically justified global-optimality argument supports that claim.

#### Fixed scenario inputs versus A* search variables

One A* run represents one fully declared hypothetical scenario. The black-hole configuration and the particle model are fixed during that run; A* changes only the bounded candidate variables. This prevents the algorithm from mixing physically different black holes or particle species into a single search path.

| Fixed for one search run | May vary from node to node |
| --- | --- |
| black-hole mass and source/scenario label | split radius |
| dimensionless spin `a_star` | incoming angular momentum or impact parameter |
| coordinate convention and physical constants | incident direction or inward angle |
| particle species, rest mass, and charge model | split direction angle |
| model tier: neutral test particle, equatorial Kerr baseline | permitted momentum-share parameter |
| escape radius, integration tolerances, and residual tolerances | later, a tightly bounded fragment-mass or energy-share variable |
| `eta_target`, parameter bounds, and step schedule | no other value unless documented in the scenario |

Mass and spin are scenario inputs, not ordinary A* coordinates. For an observed catalogue record, preserve the original mass and spin values, their units, uncertainty, provenance, and whether the value is observed or a labelled scenario. Do not let a search silently adjust an observation merely to produce a better result.

#### Start with a small, named parameter vector

The first A* implementation must vary no more than three parameters. A recommended initial vector is:

```text
x = (r_split_over_M, incoming_Lz_over_mM, split_angle_rad)
```

where:

- `r_split_over_M` is the equatorial split radius normalized by the black-hole mass scale in geometrized units;
- `incoming_Lz_over_mM` is the incident axial angular momentum normalized by the particle rest mass and black-hole mass scale;
- `split_angle_rad` is the escaping fragment's local split direction in radians.

For example, a candidate must be written with field names and units, not as an unexplained tuple:

```text
(r_split_over_M = 3.0,
 incoming_Lz_over_mM = 1.9,
 split_angle_rad = 0.28)
```

Only introduce a fourth variable, such as a normalized momentum-share parameter, after the three-variable graph, evaluator, and tests are correct. A grid with 21 possible values in each of six independent dimensions already contains more than 85 million candidates. This state-space growth is the reason to start small, bound every variable, and record every discretization choice.

Use geometrized and dimensionless quantities inside the search. Convert user-facing SI values only at the scenario-input boundary and report both original and normalized values in the result. Angles are radians internally. Do not mix kilograms, metres, joules, and geometrized coordinates inside a graph key.

Not every momentum component is an independent search variable. The candidate parameterization must leave enough quantities for the physics evaluator to derive from the mass-shell condition and local four-momentum conservation. Do not let users type arbitrary independent momentum and angle values and then call an algebraically inconsistent state a valid node.

#### Candidate-node contract

Each graph node represents one candidate parameter set and its deterministic evaluation. It should conceptually contain:

```text
node id and canonical quantized parameter key
candidate parameter values in normalized units
parent node id and the one local change used to reach this node
g cost, h heuristic, and f = g + h priority
physics-evaluation status
E_input, E_escape, E_extracted, and eta_penrose
capture, escape, invariant, and residual diagnostics
search depth and deterministic discovery order
```

The canonical key must be based on quantized integer grid indices, not raw floating-point bit patterns. For example, store a split-radius index, angular-momentum index, and angle index; derive the displayed floating-point values from the fixed grid origin and step. This prevents numerically equivalent candidates from being evaluated repeatedly because of floating-point roundoff.

Every node must have one of these explicit statuses:

- `outside_search_domain`: violates configured parameter bounds; do not enqueue it.
- `physics_invalid`: fails a required conservation, mass-shell, causal, coordinate, or finite-value check.
- `captured_or_non_escaping`: evaluates cleanly but the intended escaping fragment does not reach the escape radius.
- `escaping_without_target`: escapes but does not reach `eta_target`.
- `goal_feasible`: satisfies every physical condition and the energy target.
- `integration_failed`: integrator or event detector did not produce a trustworthy result; preserve its failure reason.

A node is a point in parameter space, not automatically a feasible physical event. The returned A* path may therefore include non-goal parameter states. Label it `parameter_adjustment_path` or `search_trace`, never `particle_path`. Each entry must show its status so the user can see why the search continued.

#### Graph construction and neighbor generation

Generate a deterministic local graph. From a three-parameter node, the baseline graph has up to six direct neighbors:

```text
r_split index:            -1, +1
incoming L_z index:       -1, +1
split-angle index:        -1, +1
```

Each edge changes exactly one grid index by one configured step. Define a fixed neighbor order, for example radius minus, radius plus, angular momentum minus, angular momentum plus, angle minus, angle plus. Stable ordering is required when priorities tie so repeated runs produce the same result and the same returned path.

Bounds are mandatory:

- the split radius must remain outside the outer horizon and inside the selected ergosphere region;
- angular momentum, energy, and angles must remain within explicitly documented scenario limits;
- all values must be finite;
- the step size must be positive and recorded;
- a candidate already visited at an equal or lower `g` cost must not be re-expanded;
- parameter changes must never silently alter fixed black-hole, particle, tolerance, or model-tier inputs.

The initial version uses a fixed coarse grid. A later two-stage version may first search a coarse grid, then construct a smaller, finer grid around the best validated candidate. Record both stages, their bounds, step sizes, node counts, and the parent result. Do not describe coarse-to-fine refinement as a proof of global optimality.

#### A* cost and heuristic contract

For A*, every edge cost must be non-negative. The simple and auditable baseline is a weighted normalized adjustment cost:

```text
edge_cost =
    w_r * abs(delta_r_split_over_M) / step_r
  + w_L * abs(delta_incoming_Lz_over_mM) / step_L
  + w_a * abs(delta_split_angle_rad) / step_angle
```

With one-coordinate, one-step neighbors, each edge normally has one positive component. Use `g(node)` for the sum of all edge costs from the declared start node. The weights express a search preference, not a physical force; store them in the scenario and show them in the output.

The baseline heuristic is valid only when it is an admissible lower bound on remaining adjustment cost. If a proven upper bound on achievable efficiency improvement per one-step change is available, a conservative heuristic may be:

```text
h(node) = minimum_possible_adjustment_cost
          * ceil(max(0, eta_target - eta_penrose(node)) / proven_eta_gain_per_step)
```

Do not invent `proven_eta_gain_per_step`. The energy surface can be non-monotonic near capture and escape boundaries. Until such a bound is justified and tested, set:

```text
h(node) = 0
```

This makes the baseline search Dijkstra's algorithm, which is slower but retains the correct minimum-cost guarantee. A non-admissible guessed heuristic is permitted only when explicitly labelled `weighted_a_star` or `heuristic_best_first`; it may find a useful candidate faster but cannot claim an optimal adjustment path.

The priority rule is always:

```text
f(node) = g(node) + h(node)
```

When two nodes have the same `f`, break ties deterministically using, in order: lower `h`, lower `g`, lexicographic grid key, then discovery order. Do not let pointer addresses, hash-table iteration order, or thread timing choose the result.

#### Physics evaluation sequence for every candidate

The A* module must call a pure, deterministic candidate evaluator owned by the physics/integrator layer. The evaluator must not read a global priority queue, mutate search state, print UI output, or depend on the order in which A* discovers nodes.

For an equatorial neutral test-particle baseline, each evaluation must proceed in this order:

1. Validate finite normalized candidate parameters and convert them to the documented Kerr coordinate/state convention.
2. Verify that the proposed split event is outside the outer horizon and within the allowed ergosphere domain.
3. Construct the incident state and integrate or otherwise validate the incident trajectory to the split event.
4. Construct the two fragment states by enforcing local four-momentum conservation and the applicable mass-shell/causality constraints.
5. Verify that the candidate negative-energy fragment has the required energy relative to infinity and is captured across the horizon.
6. Integrate the positive-energy fragment outward; it must remain outside the horizon and reach the configured large escape radius.
7. Check conservation of `E`, `L_z`, and `Q` where the selected model defines them, plus all integration and split residuals.
8. Calculate `E_input`, `E_escape`, `E_extracted`, and `eta_penrose`, then return a structured status rather than a bare boolean.

If any required calculation fails, return a structured failure status and numerical diagnostic. Do not turn an integrator failure into zero energy or treat it as a captured particle without recording why.

#### Meaning and reporting of the 20.7 percent limit

The approximately 20.7 percent value belongs to a restricted idealized classical Penrose-process baseline: a neutral single split around an extremal Kerr black hole under its stated assumptions. In the conventional efficiency definition it is:

```text
eta_penrose_limit = (sqrt(2) - 1) / 2 approximately 0.20710678
```

It is not the approximately 29 percent total rotational-energy-reservoir limit, it is not a universal result for all Penrose variants, and it is not a delivery efficiency. Do not combine these quantities.

The catalogue model requires sub-extremal `0 <= a_star < 1`, so an exact 20.7 percent result is generally an unattainable limiting target in the baseline model. The UI and CLI must say `near classical single-split limit` rather than `maximum obtainable energy` when a result approaches it. A target at or above the applicable theoretical bound must produce a clear `target_unattainable_under_model` status, not an endless search.

For early demonstrations, choose a lower labelled target such as `eta_target = 0.05` or `0.10` only when the underlying scenario and model can support it. Never choose an energy target merely because it produces a visually pleasing A* path.

#### Required A* result and user-visible output

When A* finds a goal, return a structured result containing at least:

```text
search status: found_goal, no_solution_within_bounds, target_unattainable_under_model,
               node_budget_exhausted, time_budget_exhausted, cancelled, or evaluation_failure

fixed scenario: black-hole data/provenance, particle model, model tier, units, tolerances,
                parameter bounds, grid steps, weights, eta_target, and escape radius

final candidate: named normalized values and user-facing converted values
energy ledger: E_input, E_escape, E_extracted, eta_penrose, and the separate E_rot reservoir
physical result: capture status, escape status, event radii, invariant residuals, split residual,
                 integration error, and termination reason
search result: total g cost, h value, f value, nodes generated, nodes evaluated, nodes expanded,
               duplicate nodes skipped, invalid nodes, elapsed time, and deterministic seed/version
parameter_adjustment_path: ordered nodes from the declared start candidate to the returned goal,
                           including each named parameter value, local change, status, g/h/f,
                           eta_penrose, E_extracted, and evaluator diagnostic
```

For example, the output may truthfully say:

```text
final candidate:
  (r_split_over_M = 3.1,
   incoming_Lz_over_mM = 2.0,
   split_angle_rad = 0.30)

parameter_adjustment_path:
  (3.0, 1.9, 0.28) -> (3.1, 1.9, 0.28)
  -> (3.1, 2.0, 0.28) -> (3.1, 2.0, 0.30)
```

The example values are a graph-search illustration, not validated astrophysical inputs or an assertion that this candidate escapes. Results must never imply otherwise.

The result must show both the A* parameter-adjustment path and the selected candidate's physically integrated incident/captured/escaping trajectories. They are different artifacts with different meanings.

#### Input, configuration, and reproducibility rules

Start with a versioned scenario file and headless CLI, not a free-form UI. An interactive UI may later create the same scenario file, but it must not introduce hidden defaults. The scenario must explicitly state:

- observed black-hole identifier or a clearly labelled black-hole scenario;
- original mass/spin values, units, uncertainty, provenance, and normalized values;
- selected particle species and what force model is enabled;
- coordinate system, metric signature, physical constants, and unit convention;
- search-variable names, lower and upper bounds, grid origins, step sizes, and canonical quantization rule;
- declared start candidate;
- `eta_target`, edge-cost weights, node/time budgets, and heuristic mode;
- escape radius, integrator settings, residual tolerances, and model tier;
- result-schema version, code revision, and deterministic ordering policy.

Reject a scenario with inconsistent units, an invalid starting candidate, a spin outside the model domain, a non-positive step size, reversed bounds, a target outside the declared model's permitted range, or an unlabelled assumed value.

#### Required tests and validation evidence

Test the graph algorithm independently from Kerr physics before connecting them. Use a tiny hand-constructed graph with known edge costs and a known shortest goal path to prove that A*, Dijkstra mode (`h = 0`), parent reconstruction, duplicate handling, and deterministic tie breaking work correctly.

Then test the candidate evaluator independently using analytic limits and fixed hand-constructed scenarios. Only after both layers pass independently may an integration test call A* with the real evaluator.

The A* test suite must cover at least:

- one known toy graph where A* returns the correct minimum-cost path and exact parent sequence;
- A* with `h = 0` returning the same cost as Dijkstra on the same graph;
- deterministic repeated runs returning identical final candidate, path, counts, and statuses;
- invalid bounds, zero or negative step sizes, non-finite inputs, and invalid start candidates;
- a start candidate that already satisfies the goal;
- no feasible goal within the bounded graph;
- a theoretical target rejected as unattainable before unnecessary expansion;
- duplicate candidate keys reached by different paths, retaining only the lower-cost parent;
- a tie between candidates, verifying the documented deterministic tie break;
- physics-invalid, captured, non-escaping, escaping-but-below-target, and goal-feasible evaluation outcomes;
- error propagation from a failed integrator without hiding the failure;
- path reconstruction where every returned node contains named values and diagnostics;
- comparison of the returned goal by a fresh independent evaluator call;
- evidence that no returned escaping fragment crosses the horizon.

For every search result, re-evaluate the final candidate after the search finishes. A cache hit is not sufficient verification. The returned energy ledger and feasibility status must match within documented floating-point tolerances.

#### Performance and scaling rules

The first A* version is scalar, single-threaded, bounded, and deterministic. Profile the evaluator before optimizing the priority queue. The expensive work is expected to be geodesic and constraint evaluation, not the heap operation.

Cache evaluations by canonical integer grid key only when the cache value includes the full evaluator status and diagnostics. The cache key must include every fixed scenario field that can affect physics; never reuse a candidate result across different black holes, particle models, tolerances, or escape radii.

Do not parallelize one priority queue until the scalar semantics and deterministic result are proven. Later safe parallel work may include independent scenario runs, independent target thresholds, independent bounded tiles, or independent coarse-grid searches. The combined result must still state whether it is an A* result, a tile comparison, or a non-global parallel search.

Kubernetes, when introduced, must run independent deterministic experiment jobs after the local CLI is validated. A Kubernetes Job may vary one explicitly recorded scenario, grid resolution, target, or bounded tile and write a separate result file. Kubernetes must not be presented as part of the A* algorithm or used to hide nondeterministic search behavior.

#### Implementation order for this layer

Implement the algorithmic layer in this order:

Define strongly typed scenario, candidate, node-status, evaluation-result, and search-result schemas in domain/.
Implement a deterministic scalar candidate evaluator in physics/ and integrators/; verify that it works independently of A*.
Define the A* search contract: state representation, goal condition, transition cost, heuristic, bounds, and canonical state keys.
Implement the generic A* engine with a priority queue, best-known costs, closed-state handling, parent links, and path reconstruction.
Test the generic A* engine on a toy integer-grid graph with known optimal answers.
Add deterministic neighbor generation for the real parameter space, initially using no more than three search variables.
Connect A* to one fixed black-hole scenario and the verified scalar evaluator.
Report the final named parameters, net extracted energy, efficiency, physical diagnostics, total search cost, and complete parameter-adjustment path.
Establish a reproducible baseline, then add threshold sweeps, coarse-to-fine refinement, and comparisons against Dijkstra and greedy best-first search.
Add profiling, evaluator caching, independent experiment parallelism, and Kubernetes jobs only after correctness is demonstrated.

Keep the source boundary explicit:

```text
domain/        search scenario, candidate, evaluation, node, and result value types
physics/       Kerr formulas and local split constraints; no priority queue
integrators/   trajectory propagation and event detection; no graph ownership
optimization/  A* frontier, neighbor generation, costs, heuristics, and path reconstruction
app/           CLI/UI parsing, scenario persistence, and result presentation
```

This division keeps the A* design suitable for algorithmic and systems-oriented portfolio work without overstating its scientific meaning.

## 10. Correctness and validation

Every physics feature requires an analytic limit, invariant, or independent reference test. Include:

- Schwarzschild behavior as spin approaches zero;
- horizon and static-limit locations;
- the near-extremal rotational-energy limit;
- normalization of timelike and null momenta;
- conservation of `E`, `L_z`, and `Q` along geodesics;
- four-momentum conservation at the split;
- convergence under decreasing integration tolerances;
- robust horizon-crossing and escape event detection;
- dimensional-analysis tests for all SI conversions;
- derivative checks against finite differences;
- fixed-input regression tests for optimized scenarios;
- solver feasibility and optimality residual checks.

Performance implementations additionally require:

- scalar versus SIMD numerical comparisons;
- one-thread versus many-thread result comparisons;
- sequential versus concurrent-pipeline output comparisons;
- deterministic output-order tests;
- race, deadlock, cancellation, and exception-path tests;
- odd-size and remainder-lane SIMD tests.

Report absolute and relative residuals. Never hide rejected records, infeasible states, numerical instability, conservation drift, failed optimizations, or incomplete cancellation.

## 11. Benchmarking policy

Correctness precedes optimization:

1. establish the scalar single-thread reference;
2. create representative datasets and solver workloads;
3. profile to find measured bottlenecks;
4. add SIMD and measure kernel speedup;
5. add multithreading and measure scaling;
6. add the concurrent pipeline and measure overlap, backpressure, and memory use;
7. test the combined implementation for oversubscription and numerical equivalence.

Every benchmark report must include:

- code revision;
- CPU model, core count, and SIMD ISA;
- compiler, version, flags, and build type;
- operating system;
- record count and solver workload;
- thread count, queue capacity, and batch size;
- accuracy target and tolerances;
- wall time, throughput, latency, CPU utilization, and peak memory;
- median, p95, and p99 latency for the named workload class;
- L1/L2 cache misses and branch misses when supported by the measurement platform;
- reference result and measured speedup.

Do not quote a speedup without its baseline. Do not benchmark only trivial energy formulas while implying that the result represents full geodesic optimization.

## 12. Run metadata and reproducibility

Each run must capture:

- input catalogue name, version, checksum, and citations;
- regional/world energy-demand dataset version and citations;
- particle/isotope dataset version and citations;
- code revision and result-schema version;
- compiler and build configuration;
- physical constants and unit convention;
- all observed inputs and hypothetical assumptions;
- integrator, tolerances, solver, initial guess, and stopping criteria;
- execution mode, thread count, queue capacity, SIMD backend, and batch width;
- counts of accepted, scenario-derived, rejected, feasible, and infeasible records;
- conservation and optimization residuals;
- runtime and hardware metadata.

The CLI must support headless execution for CI and compute environments. The UI must call the same engine and consume the same result schema as the CLI.

## 13. Sequential implementation milestones

### Milestone 1: data and analytic baseline

- define strong types and schemas;
- ingest and validate a small catalogue;
- ingest versioned particle-species and regional/world demand datasets;
- calculate `E_mass`, `M_irr`, and `E_rot` with a scalar single-thread kernel;
- add unit, limit, and dimensional-analysis tests.

### Milestone 2: scalar trajectory engine

- implement and validate equatorial geodesics before full Kerr geodesics;
- add adaptive error control and event detection;
- monitor invariants;
- validate one hand-constructed Penrose event;
- establish end-to-end scalar reference outputs.

### Milestone 3: deterministic optimization

- formulate the constrained split problem;
- define a bounded, normalized two- or three-variable candidate-parameter graph with named units, fixed black-hole inputs, a declared start candidate, and deterministic neighbor order;
- implement and test scalar A* mechanics separately on a hand-constructed graph before connecting it to Kerr physics;
- connect A* to the pure candidate evaluator so it finds a minimum parameter-adjustment-cost path to a feasible configured energy target, never a substitute for a physical geodesic;
- return the final named candidate, net extracted energy, `eta_penrose`, full physical diagnostics, and the ordered parameter-adjustment path used to reach it;
- compare the A* result with Dijkstra mode (`h = 0`) and report whether a heuristic is admissible, weighted, or disabled;
- solve and verify feasible capture/escape trajectories;
- add continuation, sensitivities, and phase-space mapping;
- record feasibility, residuals, and termination reasons.

### Milestone 4: charged particles and plasma foundation

- add isotope-backed mass and charge to incident-particle scenarios;
- implement charged-particle motion in a documented electromagnetic field;
- build a scalar relativistic plasma reference with explicit assumptions;
- add conservation, convergence, and analytic-wave tests before any GRMHD claim;
- treat realistic GRMHD as a separate validated solver tier, not a coefficient attached to extracted energy.

### Milestone 5: multithreading and cache-aware kernels

- add a fixed worker pool;
- parallelize independent records and parameter tiles;
- make solver scratch state thread-local;
- prove deterministic result equivalence;
- benchmark scaling and investigate contention and false sharing.
- profile working sets and tile uniform kernels for measured L1/L2 behavior.

### Milestone 6: concurrent pipeline

- add bounded queues and backpressure;
- overlap ingestion, validation, computation, and ordered output;
- implement cancellation and failure propagation;
- stress-test lifecycle, queue, and exception behavior.

### Milestone 7: SIMD

- convert profiled hot data to structure-of-arrays batches;
- validate compiler auto-vectorization;
- implement explicit SIMD only where beneficial;
- handle divergent trajectories and tail lanes;
- benchmark scalar, SIMD, threaded, and combined modes.

### Milestone 8: Python integration, global perspective, and visualization

- implement the explicit efficiency waterfall;
- expose coarse-grained, versioned C++ APIs to Python without duplicating physics;
- compare with a selected regional or worldwide metric and year;
- show assumptions and citations in every result;
- render the horizon, ergosphere, trajectories, residuals, phase maps, and energy perspective;
- create reproducible demonstration scenarios.

## 14. Agent working rules

When modifying this repository:

1. Follow the milestone order unless a dependency requires a documented exception.
2. State every physical, numerical, delivery, and data assumption introduced.
3. Keep observations distinct from scenarios and theoretical bounds.
4. Keep physics independent from I/O, concurrency, SIMD, comparisons, and rendering.
5. Retain a scalar single-thread reference path permanently.
6. Add equivalence tests before adding each optimized execution mode.
7. Use multithreading across independent work and SIMD within uniform batches.
8. Preserve deterministic output ordering regardless of task completion order.
9. Do not invent catalogue values, isotope properties, demand values, efficiencies, citations, or benchmark results.
10. Run the relevant tests and benchmarks and report exactly what was run.
11. Flag idealizations and current technological infeasibility in documentation and UI.
12. Treat races, deadlocks, instability, infeasibility, conservation drift, and solver failures as defects or explicit results, never values to suppress.

## 15. Definition of success

The completed project should demonstrate one auditable chain:

> Read measured black-hole properties and cited particle-species data -> validate and normalize them -> calculate the rotational-energy bound -> evaluate physically admissible Penrose-event candidates -> use a deterministic A* parameter-space search to find and report a bounded minimum-adjustment path to a feasible energy target, without confusing that search trace with a particle trajectory -> optimize and cross-check the selected event whose escaping fragment never crosses the horizon -> model explicitly scoped plasma and collection losses -> compare hypothetical delivered energy with a cited regional or worldwide demand metric -> reproduce equivalent ordered results through scalar, multithreaded, concurrent, SIMD, and Python-facing execution modes -> report correctness, latency, cache, and performance evidence.

That chain, rather than visual effects alone, is the core portfolio result.
