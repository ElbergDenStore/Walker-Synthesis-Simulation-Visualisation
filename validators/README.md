# validators/

Correctness checks and performance benchmarks. Not on the default path —
add with `addpath('validators')` before running.

| File | Purpose |
|---|---|
| `smoke_test.m` | 3-layer repo health check: checkcode sweep, tiny runtime sim, dependency resolution. |
| `validate_simulators.m` | Cache integrity + fast-math vs toolbox equivalence (bit-identical). |
| `validate_propagator_accuracy.m` | Pure-math `fast_walker_ecef` vs toolbox propagator over 3 h / 24 h. |
| `benchmark_aer_methods.m` | `aer()` workflow vs vectorised ENU geometry benchmark. |
| `benchmark_simulator_speed.m` | Parallelism / worker-scaling benchmark. |
| `diagnose_coverage_gap.m` | Reproduce a specific coverage failure for debugging. |
