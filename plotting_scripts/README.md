# plotting_scripts/

Standalone visualisation scripts. Each is run directly (most via the Run
button; some take optional arguments). Output figures are written under
`figures/`. Scripts bootstrap the engine via `path_setup`, so they can be run
from anywhere.

## beams/
| File | Purpose |
|---|---|
| `plot_beamgrid.m` | Plot a beam grid. |
| `plot_oneweb_beamgrid_example.m` | Example OneWeb beam-grid plot. |
| `plot_beams_on_earth.m` | Beams projected onto the Earth (3D viewer). |
| `plot_beam_utilization.m` | Beam-utilisation illustration. |

## orbits/
Propagator / altitude validation and geometry plots.
| File | Purpose |
|---|---|
| `plot_altitude_vs_latitude.m` | Geodetic altitude vs latitude over a pass. |
| `plot_circular_orbit_altitude.m` | Circular-orbit altitude behaviour. |
| `plot_star_altitude_spread.m` | Walker-Star altitude spread. |
| `plot_delta_altitude_spread.m` | Walker-Delta altitude spread. |
| `plot_orbit_asymmetry.m` | Orbit non-circularity (asc/desc + N/S asymmetry) across propagators. |
| `plot_frozen_eccentricity.m` | Frozen-orbit eccentricity survey. |
| `plot_seam_geometry.m` | Walker-Star seam geometry. |
| `plot_constellation_comparison.m` | Compare optimal constellations. |
| `plot_star_vs_delta_comparison.m` | 48 h Walker Star-vs-Delta link comparison. |
| `plot_oneweb_tle_altitude.m` | OneWeb TLE altitude trace (`.py` sibling cross-validates via SGP4). |

## sweeps/
Gridsearch result plots.
| File | Purpose |
|---|---|
| `plot_sweep_delta.m` | Walker-Delta sweep results. |
| `plot_sweep_star.m` | Walker-Star sweep results. |
| `plot_sweep_star_inclined.m` | Inclined Walker-Star sweep results. |
| `summarize_sweep.m` | Summarise a sweep run. |
| `plot_gridsearch.m` | Plot a single gridsearch result. |
| `plot_coverage_area_vs_altitude.m` | Coverage area vs altitude. |

## link_budget/
| File | Purpose |
|---|---|
| `plot_steering_loss.m` | Beam-steering loss vs altitude. |
| `plot_tx_gain_vs_altitude.m` | Adjusted Tx gain vs altitude. |
| `plot_payload_power_model.m` | Payload power model. |
| `plot_interference.m` | Interference visualisation. |

## p618/
| File | Purpose |
|---|---|
| `plot_p618_lookup_map.m` | Plot P.618 attenuation (20° elevation) on a geographic map. |
| `plot_p618_aalborg_study.m` | P.618 attenuation vs frequency and elevation for a fixed location. |

## ue_distribution/
| File | Purpose |
|---|---|
| `compare_ue_distribution.m` | Compare UE distributions (population vs uniform). |
| `plot_population_density.m` | Population-density modelling plot. |
| `plot_roi.m` | Plot the region of interest. |
| `plot_elevation_histogram.m` | Elevation-angle histogram. |
| `coverage_population_comparison/` | Population-weighted coverage comparison (data + tables). |
