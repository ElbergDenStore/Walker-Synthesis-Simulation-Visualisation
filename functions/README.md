# functions/

Engine library. Every file here is a **function** (called by the top-level
scripts, never run directly). `path_setup` adds this folder recursively, so
nested subfolders resolve automatically by filename.

## Core
| File | Purpose |
|---|---|
| `path_setup.m` | Recursive path bootstrap (`addpath(genpath(...))`). |
| `gridsearch.m` | Parallel gridsearch engine over (planes, sats/plane, inclination). |
| `fast_walker_ecef.m` | Pure-math ECEF propagator. Walker **Star** + **Delta** via the `WalkerStar` flag. |
| `generate_walker_states.m` | Frozen-orbit ECI initial states (Star + Delta), Newton-solved. |
| `calculate_walker_star.m` | Analytical minimum-sats (Street-of-Coverage); polar + inclined (4th arg). |
| `generate_walker_star_scenario.m` | Build a Walker Star in a `satelliteScenario` (toolbox propagators). |
| `synthesis_preset.m` | Named `Master_config` presets for `numerical_walker_synthesis`. |
| `run_altitude_sweep.m` | Checkpoint/resume altitude-loop + save engine behind the numerical synthesis. |
| `get_adjusted_tx_gain.m` | Tx gain scaled for a constant ground footprint vs altitude. |
| `generate_p618_lookup_table.m` | Precompute the ITU-R P.618 atmospheric-attenuation lookup table. |

## Subfolders
- `link_budget/` — vectorised link budget, interference, throughput, PFD.
- `beams/` — beam-grid geometry models (hexagonal, OneWeb).
- `ue_generation/` — user-equipment grid generators (equal-area, population, random).
- `utils/` — infrastructure helpers (cancellation, progress).
