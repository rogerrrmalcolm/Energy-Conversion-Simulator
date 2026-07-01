# Black Hole Energy Simulation - Agent Guide

## 1. Project mission

Build a high-performance C++20 application that ingests a catalogue of observed black holes, validates the available measurements, estimates each object's theoretical rotational-energy reservoir, solves a deterministic Penrose-process optimization problem, and compares a hypothetical amount of deliverable energy with Canadian energy demand.

The project is also a systems-programming learning platform. Multithreading, concurrency, and SIMD are first-class requirements, not optional polish. Every optimized implementation must be compared with a correct scalar reference implementation.

The intended audience includes quantitative developers, performance-oriented C++ engineers, scientific-computing teams, and big-technology infrastructure engineers. Optimize the project for mathematical transparency, reproducibility, data quality, numerical correctness, concurrency safety, and defensible performance measurements.

## 2. Required scientific claim

Use this description consistently:

> Given observed mass and spin estimates, calculate a black hole's theoretical rotational-energy reservoir, optimize an idealized Penrose extraction event, and compare hypothetical deliverable energy with Canadian energy demand.

Do not claim that the program measures energy that has actually been extracted from an observed black hole. A catalogue normally contains measured or inferred properties such as mass, spin, distance, and uncertainty. It does not contain a measurement of energy removed by a Penrose process.

Do not claim that delivering energy from an astrophysical black hole to Canada is currently feasible. The delivery stage is a transparent scenario model whose assumptions can be changed by the user.

Keep these values separate throughout the engine, result schema, and user interface:

1. total black-hole mass-energy;
2. theoretical rotational-energy reservoir;
3. energy carried to infinity by the optimized escaping fragment;
4. energy intercepted by a hypothetical collector;
5. energy remaining after conversion, storage, and transmission;
6. energy hypothetically available to Canada.

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

Begin with equatorial geodesics in Boyer-Lindquist coordinates. Add full `r`, `theta`, `phi`, and `t` evolution only after the restricted solver is validated. Candidate solution methods include multiple shooting, direct collocation, sequential quadratic programming, continuation, and adaptive subdivision. Cross-check local solutions using a second method or bounded search over the low-dimensional feasible domain.

Report feasibility, primal and dual residuals, Karush-Kuhn-Tucker residuals where applicable, integration error, termination reason, and whether the result is known only to be a local optimum. A visually convincing trajectory is not proof of correctness.

### Stage 5: convert extracted energy into a usable-energy scenario

Apply an explicit efficiency chain rather than one unexplained percentage:

```text
E_canada = E_extracted
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
- `eta_transmission` accounts for delivery losses to Canada.

All efficiencies must be in `[0, 1]`, individually configurable, and visible in the output. Do not hide them in source code. Provide best-case, nominal, and conservative scenario files, but label all three as hypothetical.

Distance must not be represented by a fabricated efficiency formula. If a delivery model is implemented, state its transport mechanism, geometry, travel time, collector area, beam divergence or inverse-square assumptions, and losses. Otherwise, use distance only as context and explicitly say that the delivery physics is not modeled.

### Stage 6: compare the result with Canadian demand

Maintain a separate, versioned Canadian energy dataset containing:

- metric name, such as electricity demand or total final energy use;
- value and unit;
- geographic scope;
- reporting year;
- authoritative source URL or citation;
- retrieval date.

Do not mix electricity consumption with total energy consumption. Let the user choose the comparison metric.

Calculate perspectives such as:

```text
years_of_demand = E_canada / annual_canadian_demand
average_power    = E_canada / delivery_duration
```

Also show the assumed delivery duration, equivalent continuous power, fraction of annual demand, and the difference between theoretical reservoir and usable energy. Time-sensitive Canadian values belong in data files, not hard-coded UI strings.

### Stage 7: present auditable results

For every black hole, display:

- original and normalized catalogue values;
- observed values versus assumed scenarios;
- mass-energy and rotational-energy reservoir;
- incident energy, total escaping energy, net extracted energy, and Penrose efficiency;
- solver status and conservation residuals;
- every collection and delivery efficiency;
- hypothetical usable Canadian energy;
- Canadian comparison metric, year, source, and result;
- uncertainty or sensitivity ranges;
- runtime and execution backend.

The visualization should show the horizon and ergosphere, incoming/captured/escaping trajectories, a live energy ledger, constraint residuals, and an efficiency waterfall from `E_rot` to `E_canada`. Phase-space maps and Pareto frontiers are more useful than distributions of randomly generated events.

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

## 8. C++ architecture

Target C++20 or newer and use CMake. Keep physics, numerical algorithms, scheduling, I/O, Canadian comparisons, and rendering independent.

```text
src/
  domain/        Strong types, catalogue records, scenarios, result schemas
  io/            Catalogue readers, Canadian datasets, configuration, output
  physics/       Kerr geometry, invariants, units, rotational energy
  integrators/   Scalar geodesics, adaptive stepping, event detection
  optimization/  Constraints, shooting/collocation, continuation, sensitivity
  simd/          Structure-of-arrays batches and vectorized kernels
  concurrency/   Bounded queues, thread pool, cancellation, pipeline
  analytics/     Phase maps, Pareto frontiers, uncertainty and sensitivity
  rendering/     Trajectories, energy waterfall, comparison views
  app/           CLI and interactive entry points
tests/
benchmarks/
data/catalogues/
data/canada/
scenarios/
docs/
```

Prefer value types, RAII, explicit ownership, `std::span` or ranges for non-owning views, and strong quantity types at SI boundaries. Avoid global mutable state. Keep solver state, thread count, batch width, tolerances, and scenario assumptions explicit.

Potential dependencies must be justified. Suitable candidates include Eigen, Boost.Odeint, Catch2 or GoogleTest, Google Benchmark, `nlohmann/json`, and Arrow/Parquet or HDF5 when dataset scale warrants them.

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
- reference result and measured speedup.

Do not quote a speedup without its baseline. Do not benchmark only trivial energy formulas while implying that the result represents full geodesic optimization.

## 12. Run metadata and reproducibility

Each run must capture:

- input catalogue name, version, checksum, and citations;
- Canadian comparison dataset version and citation;
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
- ingest a versioned Canadian comparison dataset;
- calculate `E_mass`, `M_irr`, and `E_rot` with a scalar single-thread kernel;
- add unit, limit, and dimensional-analysis tests.

### Milestone 2: scalar trajectory engine

- implement equatorial Kerr geodesics;
- add adaptive error control and event detection;
- monitor invariants;
- validate one hand-constructed Penrose event;
- establish end-to-end scalar reference outputs.

### Milestone 3: deterministic optimization

- formulate the constrained split problem;
- solve and verify feasible capture/escape trajectories;
- add continuation, sensitivities, and phase-space mapping;
- record feasibility, residuals, and termination reasons.

### Milestone 4: multithreading

- add a fixed worker pool;
- parallelize independent records and parameter tiles;
- make solver scratch state thread-local;
- prove deterministic result equivalence;
- benchmark scaling and investigate contention and false sharing.

### Milestone 5: concurrent pipeline

- add bounded queues and backpressure;
- overlap ingestion, validation, computation, and ordered output;
- implement cancellation and failure propagation;
- stress-test lifecycle, queue, and exception behavior.

### Milestone 6: SIMD

- convert profiled hot data to structure-of-arrays batches;
- validate compiler auto-vectorization;
- implement explicit SIMD only where beneficial;
- handle divergent trajectories and tail lanes;
- benchmark scalar, SIMD, threaded, and combined modes.

### Milestone 7: Canadian perspective and visualization

- implement the explicit efficiency waterfall;
- compare with a selected Canadian metric and year;
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
9. Do not invent catalogue values, Canadian demand values, efficiencies, citations, or benchmark results.
10. Run the relevant tests and benchmarks and report exactly what was run.
11. Flag idealizations and current technological infeasibility in documentation and UI.
12. Treat races, deadlocks, instability, infeasibility, conservation drift, and solver failures as defects or explicit results, never values to suppress.

## 15. Definition of success

The completed project should demonstrate one auditable chain:

> Read measured black-hole properties -> validate and normalize them -> calculate the rotational-energy bound -> optimize a physically admissible Penrose event -> apply explicit hypothetical collection and delivery losses -> compare with a cited Canadian demand metric -> reproduce the same ordered results through scalar, multithreaded, concurrent, and SIMD execution modes -> report correctness and performance evidence.

That chain, rather than visual effects alone, is the core portfolio result.
