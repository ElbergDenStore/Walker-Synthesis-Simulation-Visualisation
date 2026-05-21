#!/usr/bin/env python3
"""
plot_oneweb_tle_altitude.py

Python equivalent of plot_oneweb_tle_altitude.m for cross-validation.

Propagates 4 OneWeb satellites in the same orbital plane from their
TLEs using SGP4 and plots geodetic altitude vs latitude for one
orbital pass per satellite.

Key implementation notes vs MATLAB:
  - MATLAB's satelliteScenario returns positions in ECEF directly.
  - Python's sgp4 library returns positions in TEME.
  - astropy is used here to convert TEME -> ITRS (ECEF) to match MATLAB.
  - WGS-84 geodetic conversion is implemented to match MATLAB's ecef2lla.

Dependencies:
    pip install sgp4 astropy matplotlib numpy
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timezone, timedelta

from sgp4.api import Satrec, jday
from astropy.coordinates import TEME, ITRS, CartesianRepresentation
from astropy.time import Time
from astropy import units as u


def plot_oneweb_tle_altitude():
    # TLE definitions (identical to MATLAB script)
    # sats_def = [
    #     ('OW-61609',
    #      '1 61609U 24188R   26140.19895427 -.00000120  00000-0 -37973-3 0  9990',
    #      '2 61609  87.8870  41.8227 0000915 106.0390 254.0838 13.10369549 79625'),
    #     ('OW-51629',
    #      '1 51629U 22012H   26140.14084356 -.00000497  00000-0 -14613-2 0  9993',
    #      '2 51629  87.8856  41.8567 0001890  77.6988 282.4351 13.10373994205610'),
    #     ('OW-45149',
    #      '1 45149U 20008U   26140.18071808  .00000634  00000-0  17791-2 0  9999',
    #      '2 45149  87.8854  41.8403 0002488  80.3957 279.7451 13.10380156304515'),
    #     ('OW-51651',
    #      '1 51651U 22012AF  26140.16244155  .00000302  00000-0  82907-3 0  9990',
    #      '2 51651  87.8855  41.8592 0002020  71.7990 288.3357 13.10375373205654'),
    # ]
    sats_def = [
    ('starlink-1',
        '1 68155U 26049P   26140.01867369  .00020070  00000-0  62221-3 0  9994',
        '2 68155  97.2844 161.6747 0001303  94.3542 265.7850 15.33452098 10573'),
    ('starlink-2',
        '1 68155U 26049P   26140.01867369  .00020070  00000-0  62221-3 0  9994',
        '2 68155  97.2844 161.6747 0001303  94.3542 265.7850 15.33452098 10573')
    ]


    # Orbital period from mean motion of first satellite (matches MATLAB)
    n_revday = 13.10369549
    T_s = 86400.0 / n_revday
    print(f'\nOrbital period: {T_s/60:.2f} min')

    # Parse TLE epochs and use the latest as scenario start
    epoch_strs = ['26140.19895427', '26140.14084356', '26140.18071808', '26140.16244155']
    epoch_dts = [parse_tle_epoch(s) for s in epoch_strs]
    start_dt = max(epoch_dts)
    print(f'Scenario start: {start_dt.strftime("%Y-%m-%d %H:%M:%S")} UTC\n')

    # Time array: 5-second steps for one full orbit + 30 s (matches MATLAB SampleTime=5)
    dt_step = 5.0
    n_steps = int((T_s + 30.0) / dt_step) + 1
    time_offsets_s = np.arange(n_steps) * dt_step

    times_list = [start_dt + timedelta(seconds=float(s)) for s in time_offsets_s]

    # Pre-build astropy Time array for vectorised TEME->ITRS conversion
    iso_strings = [t.strftime('%Y-%m-%dT%H:%M:%S.%f') for t in times_list]
    astropy_times = Time(iso_strings, format='isot', scale='utc')

    # Pre-compute Julian dates for sgp4_array
    jds = np.zeros(n_steps)
    frs = np.zeros(n_steps)
    for i, t in enumerate(times_list):
        jds[i], frs[i] = jday(t.year, t.month, t.day,
                              t.hour, t.minute,
                              t.second + t.microsecond / 1e6)

    # Colours matching MATLAB (colourblind-friendly)
    cols = [
        (0.00, 0.45, 0.74),   # blue
        (0.85, 0.33, 0.10),   # orange
        (0.47, 0.67, 0.19),   # green
        (0.49, 0.18, 0.56),   # purple
    ]

    # -----------------------------------------------------------------------
    # Propagate all satellites; collect ECEF-based and TEME-based LLA data
    # -----------------------------------------------------------------------
    sat_data = []  # list of dicts, one per satellite

    for k, (name, tle1, tle2) in enumerate(sats_def):
        sat = Satrec.twoline2rv(tle1, tle2)

        # Propagate all timesteps at once (r in km, TEME frame)
        e_arr, r_arr, _v_arr = sat.sgp4_array(jds, frs)

        valid = e_arr == 0
        r_teme_m = r_arr[valid] * 1e3  # km -> m

        # --- TEME-based geodetic (apply WGS-84 directly to TEME x,y,z) ---
        # TEME differs from ECEF only by a rotation around z (GMST angle).
        # Rotation around z preserves z and |x,y|, so geodetic latitude and
        # altitude are virtually identical; only longitude is offset.
        teme_lats, _teme_lons, teme_alts = ecef2lla(
            r_teme_m[:, 0], r_teme_m[:, 1], r_teme_m[:, 2]
        )
        teme_alts_km = teme_alts / 1e3

        # --- ECEF-based geodetic (proper TEME -> ITRS via astropy) ---
        at_valid = astropy_times[valid]
        teme_frame = TEME(
            CartesianRepresentation(
                r_teme_m[:, 0] * u.m,
                r_teme_m[:, 1] * u.m,
                r_teme_m[:, 2] * u.m,
            ),
            obstime=at_valid,
        )
        itrs_frame = teme_frame.transform_to(ITRS(obstime=at_valid))
        ecef = itrs_frame.cartesian

        ecef_lats, _ecef_lons, ecef_alts = ecef2lla(
            ecef.x.value, ecef.y.value, ecef.z.value
        )
        ecef_alts_km = ecef_alts / 1e3

        print(f'{name}')
        print(f'  ECEF  min={np.min(ecef_alts_km):.2f} km  max={np.max(ecef_alts_km):.2f} km  '
              f'p-p={(np.max(ecef_alts_km) - np.min(ecef_alts_km)) * 1e3:.1f} m')
        print(f'  TEME  min={np.min(teme_alts_km):.2f} km  max={np.max(teme_alts_km):.2f} km  '
              f'p-p={(np.max(teme_alts_km) - np.min(teme_alts_km)) * 1e3:.1f} m')
        print(f'  Max |alt diff| = {np.max(np.abs(ecef_alts_km - teme_alts_km)) * 1e3:.4f} m')
        print(f'  Max |lat diff| = {np.max(np.abs(ecef_lats   - teme_lats  )) * 3600:.4f} arcsec\n')

        sat_data.append(dict(
            name=name,
            ecef_lats=ecef_lats, ecef_alts_km=ecef_alts_km,
            teme_lats=teme_lats, teme_alts_km=teme_alts_km,
        ))

    script_dir = os.path.dirname(os.path.abspath(__file__))
    out_dir = os.path.normpath(os.path.join(script_dir, '..', 'figures'))
    os.makedirs(out_dir, exist_ok=True)

    # -----------------------------------------------------------------------
    # Figure 1 – ECEF-based (matches MATLAB)
    # -----------------------------------------------------------------------
    fig1, ax1 = plt.subplots(figsize=(10.24, 5.75), facecolor='white')
    for k, d in enumerate(sat_data):
        ax1.plot(*_ascending_pass(d['ecef_lats'], d['ecef_alts_km']),
                 '-', color=cols[k], linewidth=1.8, label=d['name'])
    ax1.set_xlabel('Geodetic latitude (deg)')
    ax1.set_ylabel('Geodetic altitude (km)')
    ax1.set_title(
        'OneWeb \u2013 altitude vs latitude  [ECEF/WGS-84, matches MATLAB]',
        fontweight='bold',
    )
    ax1.legend(loc='best')
    ax1.grid(True)
    ax1.set_xlim([-90, 90])
    fig1.tight_layout()
    fname1 = os.path.join(out_dir, 'oneweb_same_plane_altitude_vs_latitude_python_ecef.png')
    fig1.savefig(fname1, dpi=300, bbox_inches='tight')
    print(f'Figure 1 saved -> {fname1}')

    # -----------------------------------------------------------------------
    # Figure 2 – Raw TEME positions fed directly into WGS-84 converter
    # -----------------------------------------------------------------------
    fig2, ax2 = plt.subplots(figsize=(10.24, 5.75), facecolor='white')
    for k, d in enumerate(sat_data):
        ax2.plot(*_ascending_pass(d['teme_lats'], d['teme_alts_km']),
                 '-', color=cols[k], linewidth=1.8, label=d['name'])
    ax2.set_xlabel('TEME pseudo-latitude (deg)')
    ax2.set_ylabel('Geodetic altitude (km)')
    ax2.set_title(
        'OneWeb \u2013 altitude vs latitude  [raw TEME coords into WGS-84, no GMST rotation]',
        fontweight='bold',
    )
    ax2.legend(loc='best')
    ax2.grid(True)
    ax2.set_xlim([-90, 90])
    fig2.tight_layout()
    fname2 = os.path.join(out_dir, 'oneweb_same_plane_altitude_vs_latitude_python_teme.png')
    fig2.savefig(fname2, dpi=300, bbox_inches='tight')
    print(f'Figure 2 saved -> {fname2}')

    # -----------------------------------------------------------------------
    # Figure 3 – Difference between ECEF and TEME altitudes (per satellite)
    # -----------------------------------------------------------------------
    fig3, axes = plt.subplots(2, 2, figsize=(12, 7), facecolor='white', sharex=True)
    axes = axes.flatten()
    for k, d in enumerate(sat_data):
        lat_e, alt_e = _ascending_pass(d['ecef_lats'], d['ecef_alts_km'])
        lat_t, alt_t = _ascending_pass(d['teme_lats'], d['teme_alts_km'])
        # Interpolate TEME onto ECEF latitude grid for a clean diff
        alt_t_interp = np.interp(lat_e, lat_t, alt_t)
        diff_m = (alt_e - alt_t_interp) * 1e3  # km -> m
        axes[k].plot(lat_e, diff_m, '-', color=cols[k], linewidth=1.5)
        axes[k].axhline(0, color='k', linewidth=0.6, linestyle='--')
        axes[k].set_title(d['name'])
        axes[k].set_ylabel('\u0394alt (m)')
        axes[k].set_xlabel('Geodetic latitude (deg)')
        axes[k].grid(True)
        axes[k].set_xlim([-90, 90])
    fig3.suptitle(
        'Altitude difference: ECEF \u2212 TEME  (effect of GMST rotation)',
        fontweight='bold',
    )
    fig3.tight_layout()
    fname3 = os.path.join(out_dir, 'oneweb_altitude_ecef_minus_teme_diff.png')
    fig3.savefig(fname3, dpi=300, bbox_inches='tight')
    print(f'Figure 3 saved -> {fname3}')

    plt.show()


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def _ascending_pass(lats: np.ndarray, alts_km: np.ndarray):
    """Return (lat_sorted, alt_sorted) for the ascending half of the orbit."""
    asc_mask = np.concatenate([[False], np.diff(lats) >= 0])
    lat_asc  = lats[asc_mask]
    alt_asc  = alts_km[asc_mask]
    sort_ix  = np.argsort(lat_asc)
    return lat_asc[sort_ix], alt_asc[sort_ix]


def parse_tle_epoch(epoch_str: str) -> datetime:
    """
    Parse a TLE epoch string 'YYDDD.FFFFFFFF' into a UTC datetime.
    Matches MATLAB's parse_tle_epoch local function.
    """
    yy = int(epoch_str[:2])
    year_full = 2000 + yy
    day_frac = float(epoch_str[2:])
    dt = datetime(year_full, 1, 1, tzinfo=timezone.utc) + timedelta(days=day_frac - 1)
    return dt


def ecef2lla(x: np.ndarray, y: np.ndarray, z: np.ndarray):
    """
    Convert ECEF (metres) to geodetic latitude (deg), longitude (deg),
    and altitude above the WGS-84 ellipsoid (metres).

    Uses the iterative Bowring method, which matches MATLAB's ecef2lla
    default behaviour (WGS-84, altitude in metres).

    Parameters
    ----------
    x, y, z : array_like  [metres]

    Returns
    -------
    lat : ndarray  [degrees]
    lon : ndarray  [degrees]
    alt : ndarray  [metres]
    """
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)

    # WGS-84 constants
    a  = 6378137.0           # semi-major axis (m)
    f  = 1.0 / 298.257223563
    b  = a * (1.0 - f)       # semi-minor axis
    e2 = 2.0 * f - f ** 2    # first eccentricity squared

    lon = np.degrees(np.arctan2(y, x))
    p   = np.sqrt(x ** 2 + y ** 2)

    # Iterative latitude solution
    lat = np.degrees(np.arctan2(z, p * (1.0 - e2)))  # initial approximation
    for _ in range(10):
        sin_lat = np.sin(np.radians(lat))
        N       = a / np.sqrt(1.0 - e2 * sin_lat ** 2)
        lat_new = np.degrees(np.arctan2(z + e2 * N * sin_lat, p))
        if np.max(np.abs(lat_new - lat)) < 1e-12:
            break
        lat = lat_new

    lat_r   = np.radians(lat)
    sin_lat = np.sin(lat_r)
    cos_lat = np.cos(lat_r)
    N       = a / np.sqrt(1.0 - e2 * sin_lat ** 2)

    # Altitude: equatorial formula; switch to polar formula near poles
    pole = np.abs(cos_lat) < 1e-10
    alt  = np.where(pole,
                    np.abs(z) / np.abs(sin_lat) - N * (1.0 - e2),
                    p / cos_lat - N)

    return lat, lon, alt


if __name__ == '__main__':
    plot_oneweb_tle_altitude()
