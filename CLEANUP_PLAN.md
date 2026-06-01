# Cleanup Plan — inventory + rename/merge map

Transient planning artifact. Edit the **Suggested name** and **Action** columns,
then I execute (renames are smoke-tested, merges are validator-tested). Once done,
the `KEEP` rows get distilled into a short README per folder and this file is deleted.

**Action codes:** `KEEP` · `RENAME` · `MERGE` (into target) · `MOVE` (to subfolder) · `DELETE` · `REVIEW` (needs your decision)

---

## ROOT — main entry points

> **Design rule — "press Run":** every top-level main is a **SCRIPT** (no `function`
> line) so it runs from the MATLAB/VS Code Run button with NO typed arguments. Config
> is edited in a `USER CONFIG` block at the top (a single `preset` selector variable +
> a few overridable fields). Everything in `functions/` is a **FUNCTION** (the engine,
> called by scripts, never run directly). `constellation_simulator.m` stays a function.

| Current | Purpose | Suggested name | Action |
|---|---|---|---|
| constellation_simulator.m | Core simulator: propagate → visibility → coverage → optional link budget | constellation_simulator.m | KEEP (engine fn) |
| get_cfg.m | Build a Cfg struct (delta/star + link-budget params) for examples | get_cfg.m | KEEP (engine fn) |
| show_constellation.m | 3D scenario viewer of a constellation + UE overlay | show_constellation.m | KEEP |
| plot_simulation.m | Publication plots from a metrics struct | plot_simulation.m | KEEP |
| run_single_coverage_test.m | Example: run ONE simulation end-to-end (script) | run_constellation_simulation_example.m | RENAMED ✓ |
| Find_valid_constellations.m | Multi-altitude **delta** gridsearch (regional) | numerical_walker_synthesis.m | MERGE-TARGET (script) |
| Find_valid_walker_delta_global.m | Same driver, **global** delta preset | → numerical_walker_synthesis.m | MERGE (preset `"global_delta"`) |
| Find_valid_walker_star_gridsearch.m | Same driver, **star** preset (i=87°) | → numerical_walker_synthesis.m | MERGE (preset `"star"`) |
| run_single_gridsearch.m | Same driver, **single altitude** | → numerical_walker_synthesis.m | MERGE (scalar `heights_km`) |
| find_valid_walker_star_inclined.m | Analytical star sizing sweep over altitude (SOC) | analytical_walker_star_synthesis.m | MERGE-TARGET (script); resolves prior REVIEW |

> **Two synthesis mains (both SCRIPTS, press Run):**
> - **`numerical_walker_synthesis.m`** — gridsearch/simulation-based sizing. `USER CONFIG`
>   block: `preset ∈ {"regional_delta","global_delta","star"}` + `heights_km`
>   (vector = sweep, scalar = single altitude). Calls helper `synthesis_preset(preset)`
>   (new fn in `functions/`) for the `Master_config` block, then `gridsearch`.
>   Collapses 4 scripts → 1.
> - **`analytical_walker_star_synthesis.m`** — closed-form Street-of-Coverage sizing.
>   `USER CONFIG` block: altitude range, inclination, min elevation/latitude. Calls the
>   merged `calculate_walker_star` (polar + inclined). Absorbs `find_valid_walker_star_inclined`.
> - **New helper:** `functions/synthesis_preset.m` — holds the named `Master_config`
>   blocks once, so the numerical script stays short and duplication-free.

---

## functions/ — core library (walker pairs to MERGE per your choice)

| Current | Purpose | Suggested name | Action |
|---|---|---|---|
| path_setup.m | Recursive path bootstrap | path_setup.m | KEEP |
| gridsearch.m | Parallel gridsearch engine over (P, S, inc) | gridsearch.m | KEEP |
| constellation_simulator's backends ↓ | | | |
| fast_walker_ecef.m | Pure-math ECEF propagator, **delta** | fast_walker_ecef.m | MERGE-TARGET |
| fast_walker_star_ecef.m | Pure-math ECEF propagator, **star** | → fast_walker_ecef.m | MERGE (add `WalkerStar` dispatch) |
| generate_walker_delta_states.m | Frozen-orbit ECI initial states, **delta** | generate_walker_states.m | MERGE-TARGET (rename) |
| generate_walker_star_states.m | Frozen-orbit ECI initial states, **star** | → generate_walker_states.m | MERGE |
| calculate_walker_star.m | Analytical min-sats, **polar** (SOC) | calculate_walker_star.m | MERGE-TARGET |
| calculate_walker_star_inclined.m | Analytical min-sats, **inclined** (general) | → calculate_walker_star.m | MERGE (inclined generalises polar; your README TODO) |
| asymmetrical_walker_star_generation.m | Build star in a satelliteScenario (toolbox) | generate_walker_star_scenario.m | RENAMED ✓ |
| calculate_hexagonal_beams.m | Hexagonal beam-grid model | calculate_hexagonal_beams.m | MOVE → functions/beams/ |
| calculate_OneWeb_beams.m | OneWeb phased-array beam geometry | calculate_oneweb_beams.m | MOVE → functions/beams/ + RENAME (case) |
| generate_equal_ish_area_UEs.m | Equal-area UE grid | generate_equal_area_ues.m | MOVE → functions/ue_generation/ + RENAME |
| generate_population_based_UEs.m | Population-weighted UE grid (WorldPop) | generate_population_ues.m | MOVE → functions/ue_generation/ + RENAME |
| generate_random_UEs.m | Uniform-random UE grid | generate_random_ues.m | MOVE → functions/ue_generation/ + RENAME |
| get_adjusted_tx_gain.m | Tx gain scaled for constant footprint | get_adjusted_tx_gain.m | KEEP (or MOVE → link_budget/) |
| generate_p618_lookup_table.m | Precompute ITU P.618 attenuation LUT | generate_p618_lookup_table.m | KEEP (or MOVE → link_budget/) |

> **Walker merge note:** you chose flag-dispatched merges. Each merged function gets a `WalkerStar` (logical) argument and validated, so callers like `constellation_simulator`'s `propagate_fast_math` call one function instead of branching. `validate_simulators` + `compare_toolbox_vs_math` guard correctness.

### functions/link_budget/ — KEEP all (good cohesive folder)
| Current | Purpose | Action |
|---|---|---|
| link_calc_matrix.m | Vectorised link budget (FSPL, P.618, gain, SINR, PFD) | KEEP |
| beam_gain_and_interference.m | Per-UE carrier/interference PSDs via KD-tree | KEEP |
| calculate_throughput_matrix.m | SINR → throughput (modified Shannon, bw sharing) | KEEP |
| PFD_calc.m | Reverse-solve Tx power for ITU PFD target | KEEP |

### functions/non-human_functions/ — RESOLVED ✓ rename folder → functions/utils/
| Current | Purpose | Action |
|---|---|---|
| CancelToken.m | Cancellation handle for the simulator | MOVE → functions/utils/ |
| updateLiveScriptProgress.m | Live-script progress bar | MOVE → functions/utils/ |

> Suggest renaming `non-human_functions/` → `utils/` (the hyphen + label is unclear; also a hyphen in a path is fine but the name is cryptic).

---

## plotting_scripts/ — propose subfolders (with genpath, zero resolution risk)

Suggested grouping:

**beams/**
| Current | Suggested name | Action |
|---|---|---|
| plot_beamgrid.m | plot_beamgrid.m | MOVE |
| plot_oneweb_beamgrid_example.m | plot_oneweb_beamgrid_example.m | MOVE |
| beams_on_earth_matlab_viewer.m | plot_beams_on_earth.m | MOVE + RENAME |
| IllustrateBeamUtilization.m | plot_beam_utilization.m | MOVE + RENAME (case) |

**orbits/** (propagator/altitude validation plots)
| Current | Suggested name | Action |
|---|---|---|
| plot_altitude_vs_latitude.m | plot_altitude_vs_latitude.m | MOVE |
| plot_circular_orbit_altitude.m | plot_circular_orbit_altitude.m | MOVE |
| constellation_altitude_spread.m | plot_star_altitude_spread.m | MOVE + RENAME |
| walker_delta_altitude_spread.m | plot_delta_altitude_spread.m | MOVE + RENAME |
| validate_satellite_initialization.m | plot_orbit_asymmetry.m | MOVE + RENAME (RESOLVED ✓ plot) |
| survey_frozen_eccentricity.m | plot_frozen_eccentricity.m | MOVE + RENAME |
| plot_seam_geometry.m | plot_seam_geometry.m | MOVE |
| plot_constellation_comparison.m | plot_constellation_comparison.m | MOVE |

**sweeps/** (gridsearch result plots)
| Current | Suggested name | Action |
|---|---|---|
| plot_sweep_delta.m | plot_sweep_delta.m | MOVE |
| plot_sweep_star.m | plot_sweep_star.m | MOVE |
| plot_walker_star_inclined_sweep.m | plot_sweep_star_inclined.m | MOVE + RENAME |
| summarize_sweep.m | summarize_sweep.m | MOVE |
| plot_gridsearch.m | plot_gridsearch.m | MOVE |
| plot_coverage_area_vs_altitude.m | plot_coverage_area_vs_altitude.m | MOVE |

**link_budget/**
| Current | Suggested name | Action |
|---|---|---|
| steering_loss.m | plot_steering_loss.m | MOVE + RENAME |
| plot_tx_gain_vs_altitude.m | plot_tx_gain_vs_altitude.m | MOVE |
| plot_payload_power_model.m | plot_payload_power_model.m | MOVE |
| interference.m | plot_interference.m | MOVE + RENAME |

**p618/**
| Current | Suggested name | Action |
|---|---|---|
| plot_p618_lookup_max_absorption_3d.m | plot_p618_lookup_map.m | MOVE + RENAME (file/function name mismatch!) |
| p618_aalborg_freq_elevation_study.m | plot_p618_aalborg_study.m | MOVE + RENAME |

**ue_distribution/**
| Current | Suggested name | Action |
|---|---|---|
| compare_ue_distribution.m | compare_ue_distribution.m | MOVE |
| userdensitymodelling.m | plot_population_density.m | MOVE + RENAME |
| plot_roi.m | plot_roi.m | MOVE |
| elevation_angle_histogram.m | plot_elevation_histogram.m | MOVE + RENAME |
| coverage_population_comparison/ (existing) | coverage_population_comparison/ | KEEP |

**misc / REVIEW**
| Current | Note | Action |
|---|---|---|
| compare_walker_star_delta.m | 48h delta-vs-star link comparison | RESOLVED ✓ plot → orbits/plot_star_vs_delta_comparison.m |
| plot_oneweb_tle_altitude.m | OneWeb TLE altitude trace | MOVE → orbits/ |
| plot_oneweb_tle_altitude.py | Python sibling | RESOLVED ✓ KEEP → orbits/ (alongside .m) |

> ⚠️ **Name/function mismatch:** `plot_p618_lookup_max_absorption_3d.m` defines `function plot_p618_lookup_map_20deg(...)`. File name ≠ function name — MATLAB calls the **file** name to dispatch but errors if you call the function name. Fix during rename.

---

## validators/ — mostly KEEP, two REVIEW

| Current | Purpose | Suggested name | Action |
|---|---|---|---|
| smoke_test.m | 3-layer repo health check | smoke_test.m | KEEP |
| validate_simulators.m | Cache integrity + fast-vs-toolbox equivalence | validate_simulators.m | KEEP |
| compare_toolbox_vs_math.m | Pure-math vs toolbox over 3/24h | validate_propagator_accuracy.m | RENAMED ✓ |
| AER_comparison.m | aer() vs vectorised ENU benchmark | benchmark_aer_methods.m | RENAMED ✓ |
| compare_simulator_speed.m | Parallelism / worker-scaling benchmark | benchmark_simulator_speed.m | RENAMED ✓ |
| diagnose_coverage_gap.m | Reproduce a specific coverage failure | diagnose_coverage_gap.m | KEEP |

---

## Suggested execution order
1. **DELETE/REVIEW pass** — resolve REVIEW rows, delete confirmed dead files.
2. **Synthesis mains merge** — collapse the 4 `Find_valid_*` / `run_single_gridsearch`
   scripts → `numerical_walker_synthesis.m` (script + `synthesis_preset.m` helper);
   create `analytical_walker_star_synthesis.m`. Both stay SCRIPTS (press Run). Smoke-test.
3. **Walker function merges** — 3 pairs → 3 functions. Run `validate_simulators` + `compare_toolbox_vs_math` after each.
4. **Renames** — apply RENAME rows (file + `function` line + every caller + `@handle`). Smoke-test per batch.
5. **Moves** — apply MOVE rows into subfolders. Smoke-test (genpath covers them).
6. **Folder READMEs** — distil KEEP rows into a short README.md per folder; delete this plan.

