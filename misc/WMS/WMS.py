"""
Synthetic WMS Diagnostic - User Data Version (v4 - corrected)
=============================================================

Reads two data files (temperature and H2O mole fraction) in the format:
    time  space  value

For each unique timestep, extracts the spatial profiles along the beam,
computes the synthetic absorbance, and retrieves the effective T and c
that the WMS sensor would report.

Usage:
    python synthetic_wms_diagnostic.py --temp <T_file> --h2o <H2O_file> \
        --L <path_length_cm> [--partfun <Q_file>] [options]

Based on:
  - Qu et al., Optics Express 23(12), 16492 (2015)
  - Qu & Schmidt, Appl. Phys. B 119, 45-53 (2015)

v4 changes:
  - FIXED: Partition function updated from user-provided Partfun_H16OH.txt
    (old TIPS-2011 had 10-21% errors at combustion temperatures)
  - Added --partfun option to load Q(T) from external file
  - Embedded fallback: 201-point table (every 25 K, 25-5000 K)
v3 changes:
  - FIXED: Doppler HWHM used molar mass instead of molecular mass
v2 changes:
  - FIXED: Partition function polynomial replaced by spline interpolation
  - Added numerical guards throughout to prevent segfaults
"""
import matplotlib
matplotlib.use('Agg')

import numpy as np
from scipy.optimize import least_squares
from scipy.special import wofz
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import argparse
import os
import sys

# =============================================================================
# CONSTANTS
# =============================================================================
c2 = 1.4388          # second radiation constant hc/k [cm*K]
T_ref = 296.0         # HITRAN reference temperature [K]
kB = 1.380649e-23       # Boltzmann constant [J/K per molecule]
c_light = 2.99792458e8  # speed of light [m/s]
N_A = 6.02214076e23     # Avogadro's number [mol^-1]
m_H2O = 18.015e-3 / N_A # molecular mass of H2O [kg/molecule]

# =============================================================================
# PARTITION FUNCTION for H2-16-O
# =============================================================================
# Default: embedded subsample (every 25 K, 25-5000 K, plus T_ref=296 K)
# from user-provided Partfun_H16OH.txt. Override with --partfun <file>.
# Spline on log(Q) for smooth interpolation.
# =============================================================================
_TIPS_T_DEFAULT = np.array([
    25, 50, 75, 100, 125, 150, 175, 200, 225, 250,
    275, 296, 300, 325, 350, 375, 400, 425, 450, 475,
    500, 525, 550, 575, 600, 625, 650, 675, 700, 725,
    750, 775, 800, 825, 850, 875, 900, 925, 950, 975,
    1000, 1025, 1050, 1075, 1100, 1125, 1150, 1175, 1200, 1225,
    1250, 1275, 1300, 1325, 1350, 1375, 1400, 1425, 1450, 1475,
    1500, 1525, 1550, 1575, 1600, 1625, 1650, 1675, 1700, 1725,
    1750, 1775, 1800, 1825, 1850, 1875, 1900, 1925, 1950, 1975,
    2000, 2025, 2050, 2075, 2100, 2125, 2150, 2175, 2200, 2225,
    2250, 2275, 2300, 2325, 2350, 2375, 2400, 2425, 2450, 2475,
    2500, 2525, 2550, 2575, 2600, 2625, 2650, 2675, 2700, 2725,
    2750, 2775, 2800, 2825, 2850, 2875, 2900, 2925, 2950, 2975,
    3000, 3025, 3050, 3075, 3100, 3125, 3150, 3175, 3200, 3225,
    3250, 3275, 3300, 3325, 3350, 3375, 3400, 3425, 3450, 3475,
    3500, 3525, 3550, 3575, 3600, 3625, 3650, 3675, 3700, 3725,
    3750, 3775, 3800, 3825, 3850, 3875, 3900, 3925, 3950, 3975,
    4000, 4025, 4050, 4075, 4100, 4125, 4150, 4175, 4200, 4225,
    4250, 4275, 4300, 4325, 4350, 4375, 4400, 4425, 4450, 4475,
    4500, 4525, 4550, 4575, 4600, 4625, 4650, 4675, 4700, 4725,
    4750, 4775, 4800, 4825, 4850, 4875, 4900, 4925, 4950, 4975,
    5000], dtype=float)

_TIPS_Q_DEFAULT = np.array([
    4.7114, 12.9615, 23.1702, 35.1531, 48.7080, 63.6774,
    79.9435, 97.4150, 116.0196, 135.7001, 156.4116, 174.5810,
    178.1203, 200.8022, 224.4419, 249.0313, 274.5686, 301.0574,
    328.5055, 356.9241, 386.3275, 416.7322, 448.1571, 480.6225,
    514.1505, 548.7646, 584.4895, 621.3514, 659.3775, 698.5966,
    739.0384, 780.7338, 823.7148, 868.0148, 913.6681, 960.7101,
    1009.1775, 1059.1077, 1110.5395, 1163.5124, 1218.0672, 1274.2454,
    1332.0897, 1391.6436, 1452.9515, 1516.0588, 1581.0117, 1647.8574,
    1716.6439, 1787.4200, 1860.2355, 1935.1409, 2012.1877, 2091.4280,
    2172.9151, 2256.7028, 2342.8458, 2431.3999, 2522.4213, 2615.9675,
    2712.0964, 2810.8672, 2912.3395, 3016.5741, 3123.6324, 3233.5770,
    3346.4710, 3462.3786, 3581.3647, 3703.4955, 3828.8376, 3957.4587,
    4089.4276, 4224.8138, 4363.6877, 4506.1207, 4652.1853, 4801.9548,
    4955.5033, 5112.9062, 5274.2396, 5439.5808, 5609.0079, 5782.6002,
    5960.4377, 6142.6017, 6329.1744, 6520.2391, 6715.8799, 6916.1824,
    7121.2326, 7331.1182, 7545.9275, 7765.7501, 7990.6765, 8220.7984,
    8456.2086, 8697.0008, 8943.2698, 9195.1118, 9452.6236, 9715.9035,
    9985.0506, 10260.1652, 10541.3486, 10828.7034, 11122.3331, 11422.3421,
    11728.8363, 12041.9222, 12361.7077, 12688.3016, 13021.8137, 13362.3550,
    13710.0374, 14064.9737, 14427.2780, 14797.0651, 15174.4509, 15559.5523,
    15952.4871, 16353.3740, 16762.3327, 17179.4836, 17604.9483, 18038.8489,
    18481.3085, 18932.4511, 19392.4015, 19861.2850, 20339.2278, 20826.3571,
    21322.8002, 21828.6857, 22344.1424, 22869.2999, 23404.2883, 23949.2384,
    24504.2814, 25069.5490, 25645.1734, 26231.2873, 26828.0236, 27435.5159,
    28053.8978, 28683.3035, 29323.8672, 29975.7235, 30639.0072, 31313.8532,
    32000.3967, 32698.7727, 33409.1167, 34131.5637, 34866.2493, 35613.3085,
    36372.8765, 37145.0884, 37930.0792, 38727.9835, 39538.9360, 40363.0708,
    41200.5221, 42051.4233, 42915.9080, 43794.1089, 44686.1587, 45592.1895,
    46512.3326, 47446.7194, 48395.4801, 49358.7448, 50336.6428, 51329.3026,
    52336.8523, 53359.4191, 54397.1295, 55450.1092, 56518.4833, 57602.3757,
    58701.9097, 59817.2078, 60948.3912, 62095.5806, 63258.8954, 64438.4542,
    65634.3745, 66846.7728, 68075.7645, 69321.4638, 70583.9840, 71863.4372,
    73159.9342, 74473.5848, 75804.4974, 77152.7795, 78518.5370, 79901.8748,
    81302.8963, 82721.7039, 84158.3984], dtype=float)

# Module-level spline (initialized with default, can be overridden)
_Q_spline = CubicSpline(_TIPS_T_DEFAULT, np.log(_TIPS_Q_DEFAULT))
_Q_Tref = np.exp(_Q_spline(T_ref))
_Q_T_min = _TIPS_T_DEFAULT[0]
_Q_T_max = _TIPS_T_DEFAULT[-1]


def load_partition_function(filepath):
    """
    Load partition function from file (two columns: T[K]  Q(T)).
    Rebuilds the global spline interpolator.
    """
    global _Q_spline, _Q_Tref, _Q_T_min, _Q_T_max
    data = np.loadtxt(filepath, comments='#')
    T_data, Q_data = data[:, 0], data[:, 1]
    _Q_spline = CubicSpline(T_data, np.log(Q_data))
    _Q_Tref = np.exp(_Q_spline(T_ref))
    _Q_T_min = T_data[0]
    _Q_T_max = T_data[-1]
    print(f"  Loaded partition function: {len(T_data)} points, "
          f"T = [{T_data[0]:.0f}, {T_data[-1]:.0f}] K, "
          f"Q(296) = {_Q_Tref:.4f}")


def partition_function(T):
    """Q(T) for H2-16-O via cubic spline on log(Q)."""
    T_clamped = np.clip(T, _Q_T_min, _Q_T_max)
    return np.exp(_Q_spline(T_clamped))


def partition_function_ratio(T):
    """Q(T_ref) / Q(T) for H2-16-O."""
    return _Q_Tref / partition_function(T)


# =============================================================================
# H2O LINE PARAMETERS (HITRAN 2012, Table 1 from the papers)
# =============================================================================
lines = np.array([
    # [nu0 (cm-1), S(T_ref) (cm-2/atm), E_lower (cm-1),
    #  gamma_air (cm-1/atm), gamma_self (cm-1/atm), n_air, n_self]
    [7153.7207, 1.9e-6,  2552.88, 0.0553, 0.371, 0.62, 0.82],  # Peak 1
    [7153.7484, 5.5e-6,  2552.86, 0.0537, 0.371, 0.62, 0.82],  # Peak 1
    [7154.3533, 9.3e-5,  1789.04, 0.0210, 0.215, 0.57, 0.65],  # Peak 2
    [7154.3534, 2.8e-4,  1789.04, 0.0325, 0.333, 0.57, 0.65],  # Peak 2
    [7154.5950, 8.0e-6,  2142.60, 0.0414, 0.281, 0.41, 0.50],  # weak
])

nu0_lines  = lines[:, 0]
S_ref      = lines[:, 1]
E_lower    = lines[:, 2]
gamma_air  = lines[:, 3]
gamma_self = lines[:, 4]
n_air      = lines[:, 5]
n_self     = lines[:, 6]
N_lines    = len(nu0_lines)

# =============================================================================
# SPECTROSCOPIC FUNCTIONS
# =============================================================================

def line_strength(T, j):
    """
    Line strength S_j(T) [cm-2/atm] scaled from HITRAN reference value.

    S(T) = S(T_ref) * [Q(T_ref)/Q(T)]
           * exp[-c2*E''*(1/T - 1/T_ref)]
           * [1 - exp(-c2*nu0/T)] / [1 - exp(-c2*nu0/T_ref)]
    """
    if T < 200:
        return 0.0
    try:
        QR = partition_function_ratio(T)
        boltz = np.exp(-c2 * E_lower[j] * (1.0/T - 1.0/T_ref))
        stim = (1.0 - np.exp(-c2 * nu0_lines[j] / T)) / \
               (1.0 - np.exp(-c2 * nu0_lines[j] / T_ref))
        result = S_ref[j] * QR * (T_ref / T) * boltz * stim
        return result if np.isfinite(result) else 0.0
    except (OverflowError, FloatingPointError):
        return 0.0


def voigt_profile(nu, nu0, gL, gD):
    """
    Normalized Voigt profile [cm].
    Convolution of Lorentzian (HWHM=gL) and Gaussian (HWHM=gD).
    Falls back to pure Lorentzian if gD is negligible.
    """
    if gL < 0:
        gL = 0.0
    # If Doppler width is negligible, use pure Lorentzian
    if gD < 1e-6 * max(gL, 1e-20):
        if gL < 1e-20:
            return np.zeros_like(nu)
        return (gL / np.pi) / ((nu - nu0)**2 + gL**2)
    sigma = gD / np.sqrt(2.0 * np.log(2.0))
    z = ((nu - nu0) + 1j * gL) / (sigma * np.sqrt(2.0))
    result = np.real(wofz(z)) / (sigma * np.sqrt(2.0 * np.pi))
    return np.where(np.isfinite(result), result, 0.0)


def doppler_hwhm(nu0, T):
    """Doppler HWHM [cm-1] for H2O."""
    return nu0 / c_light * np.sqrt(2.0 * kB * T * np.log(2.0) / m_H2O)


def lorentz_hwhm(T, P, c_h2o, j):
    """
    Pressure-broadened Lorentzian HWHM [cm-1].

    gamma_L = P * [(1 - x_H2O)*gamma_air*(T_ref/T)^n_air
                    + x_H2O*gamma_self*(T_ref/T)^n_self]
    """
    ga = gamma_air[j]  * (T_ref / T) ** n_air[j]
    gs = gamma_self[j] * (T_ref / T) ** n_self[j]
    return P * ((1.0 - c_h2o) * ga + c_h2o * gs)


# =============================================================================
# CORE COMPUTATION
# =============================================================================

def compute_absorbance_spectrum(nu, x_path, T_prof, c_prof, P=1.0):
    """
    Integrate absorbance along the beam path using CFD profiles.

    alpha(nu) = integral_0^L  x_H2O(s) * P * SUM_j[S_j(T(s)) * phi_j(nu,s)] ds

    Parameters
    ----------
    nu      : wavenumber grid [cm-1], shape (N_nu,)
    x_path  : positions along beam [cm], shape (N_x,)
    T_prof  : temperature at each position [K], shape (N_x,)
    c_prof  : H2O mole fraction at each position [-], shape (N_x,)
    P       : total pressure [atm]

    Returns
    -------
    alpha   : absorbance spectrum [dimensionless], shape (N_nu,)
    """
    N_nu = len(nu)
    alpha = np.zeros(N_nu)

    if len(x_path) < 2:
        return alpha  # need at least 2 points for integration

    for i in range(len(x_path)):
        T_loc = T_prof[i]
        c_loc = c_prof[i]
        if T_loc < 200 or c_loc < 1e-10:
            continue

        local = np.zeros(N_nu)
        for j in range(N_lines):
            Sj = line_strength(T_loc, j)
            if Sj < 1e-30:
                continue
            gL = lorentz_hwhm(T_loc, P, c_loc, j)
            gD = doppler_hwhm(nu0_lines[j], T_loc)
            local += Sj * voigt_profile(nu, nu0_lines[j], gL, gD)
        local *= c_loc * P

        # Trapezoidal weight
        if i == 0:
            dx = 0.5 * (x_path[1] - x_path[0])
        elif i == len(x_path) - 1:
            dx = 0.5 * (x_path[-1] - x_path[-2])
        else:
            dx = 0.5 * (x_path[i+1] - x_path[i-1])

        alpha += local * dx

    return alpha


def simulate_uniform(nu, T, c, L, P=1.0):
    """Absorbance spectrum for uniform T, c over path L."""
    if T < 200 or c < 1e-10 or L < 1e-10:
        return np.zeros_like(nu)
    alpha = np.zeros_like(nu)
    for j in range(N_lines):
        Sj = line_strength(T, j)
        if Sj < 1e-30:
            continue
        gL = lorentz_hwhm(T, P, c, j)
        gD = doppler_hwhm(nu0_lines[j], T)
        alpha += Sj * voigt_profile(nu, nu0_lines[j], gL, gD)
    return alpha * c * P * L


def extract_effective_T_c(nu, alpha_synth, L_assumed, P=1.0,
                          T_guess=1500.0, c_guess=0.15):
    """
    Fit synthetic absorbance with uniform-medium model
    to recover T_eff and c_eff (what the sensor would report).
    """
    def residuals(params):
        T, c = params
        try:
            model = simulate_uniform(nu, T, c, L_assumed, P)
            res = model - alpha_synth
            res = np.where(np.isfinite(res), res, 1e6)
            return res
        except Exception:
            return np.ones_like(nu) * 1e6

    try:
        # Tight tolerances + generous eval budget. The fit normally converges in
        # ~7 evals; this only guards weak-signal timesteps from stopping early.
        # NOTE: on well-resolved spectra it shifts T_eff by <0.05 K vs the scipy
        # defaults (1e-8) and leaves the RMS unchanged -- the residual floor is
        # uniform-vs-nonuniform model mismatch, not optimizer tolerance.
        result = least_squares(residuals, [T_guess, c_guess],
                               bounds=([300, 1e-6], [3500, 1.0]),
                               method='trf',
                               ftol=1e-12, xtol=1e-12, gtol=1e-12,
                               max_nfev=2000)
        return result.x[0], result.x[1], result
    except Exception as e:
        print(f"    WARNING: fit failed ({e})")
        class Dummy:
            fun = np.zeros_like(nu)
        return np.nan, np.nan, Dummy()


def compute_column_density(x_path, c_prof):
    """Compute the integrated column density of H2O along the path."""
    c_density = 0.0
    for i in range(len(x_path) - 1):
        if x_path[i] <= 0 <= x_path[i+1]:
            c_interp = np.interp(0, [x_path[i], x_path[i+1]], [c_prof[i], c_prof[i+1]])
            c_density += c_interp * (x_path[i+1] - x_path[i])

    return c_density

# =============================================================================
# DATA I/O
# =============================================================================

def load_data(filepath):
    """
    Load file with columns: time  space  value
    Lines starting with # are comments.
    """
    print(f"  Reading {filepath} ...")
    raw = np.loadtxt(filepath, comments='#')
    if raw.ndim == 1:
        raw = raw.reshape(1, -1)
    assert raw.shape[1] >= 3, \
        f"Expected >= 3 columns (time space value), got {raw.shape[1]}"

    t_col, x_col, v_col = raw[:, 0], raw[:, 1], raw[:, 2]
    unique_t = np.unique(np.round(t_col, decimals=10))
    print(f"    Found {len(unique_t)} unique timesteps, "
          f"t = [{unique_t[0]:.6g} ... {unique_t[-1]:.6g}]")

    spaces, values = {}, {}
    for t in unique_t:
        mask = np.abs(t_col - t) < 1e-10
        order = np.argsort(x_col[mask])
        spaces[t] = x_col[mask][order]
        values[t] = v_col[mask][order]

    print(f"    {len(spaces[unique_t[0]])} spatial points per timestep")
    print(f"    Space range: [{spaces[unique_t[0]][0]:.6g}, "
          f"{spaces[unique_t[0]][-1]:.6g}]")
    return unique_t, spaces, values


# =============================================================================
# MAIN
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="Synthetic WMS diagnostic from CFD profiles (v4)",
        formatter_class=argparse.RawTextHelpFormatter,
        epilog="""
EXAMPLE:
  python synthetic_wms_diagnostic.py \\
      --temp T_profile_2mm.dat \\
      --h2o xH2O_profile_2mm.dat \\
      --L 1.5 \\
      --space-unit m \\
      --partfun Partfun_H16OH.txt \\
      --output results_2mm

DATA FILE FORMAT (whitespace-separated, # for comments):
  # time   space   value
  0   0.000   1123.0
  0   0.001   1125.3
  ...
        """)

    parser.add_argument('--temp', required=True,
                        help='Temperature file (time space T[K])')
    parser.add_argument('--h2o', required=True,
                        help='H2O mole fraction file (time space x_H2O[-])')
    parser.add_argument('--L', type=float, required=True,
                        help='Assumed absorption path length [cm]')
    parser.add_argument('--partfun', default=None,
                        help='Partition function file (T[K] Q(T)), '
                             'e.g. Partfun_H16OH.txt. If not given, uses '
                             'built-in 201-point table.')
    parser.add_argument('--P', type=float, default=1.0,
                        help='Total pressure [atm] (default: 1.0)')
    parser.add_argument('--space-unit', choices=['cm', 'mm', 'm'],
                        default='cm',
                        help='Unit of the space column (default: cm)')
    parser.add_argument('--output', default='results',
                        help='Output directory (default: results)')
    parser.add_argument('--T-guess', type=float, default=1500.0,
                        help='Initial T guess for fit [K]')
    parser.add_argument('--c-guess', type=float, default=0.15,
                        help='Initial c guess for fit')
    parser.add_argument('--nu-min', type=float, default=7153.0,
                        help='Wavenumber grid start [cm-1]')
    parser.add_argument('--nu-max', type=float, default=7155.2,
                        help='Wavenumber grid end [cm-1]')
    parser.add_argument('--nu-points', type=int, default=2000,
                        help='Number of wavenumber grid points')
    parser.add_argument('--plot-every', type=int, default=1,
                        help='Plot every N timesteps (0=no plots)')
    parser.add_argument('--T-threshold', type=float, default=500.0,
                        help='Hot-region threshold [K] for arithmetic avg')
    parser.add_argument('--xmin', type=float, default=0.01,
                        help='H2O mole-fraction threshold defining the flame '
                             'body for the absorbance integration. The beam is '
                             'integrated from the first sample outward only up '
                             'to the first point where x_H2O < xmin, so the '
                             'contiguous flame (incl. its shoulder) is kept but '
                             'detached far-field H2O pockets are excluded. '
                             'Set 0 to integrate the full provided path '
                             '(legacy behaviour).')

    args = parser.parse_args()

    scale = {'cm': 1.0, 'mm': 0.1, 'm': 100.0}[args.space_unit]

    if args.partfun:
        print("\nLoading partition function:")
        load_partition_function(args.partfun)

    print("\nLoading data files:")
    times_T, spaces_T, vals_T = load_data(args.temp)
    times_c, spaces_c, vals_c = load_data(args.h2o)

    common_times = np.intersect1d(
        np.round(times_T, 10), np.round(times_c, 10))
    if len(common_times) == 0:
        print("\nERROR: No common timesteps found.")
        sys.exit(1)
    print(f"\n  {len(common_times)} common timesteps found.")

    os.makedirs(args.output, exist_ok=True)
    nu = np.linspace(args.nu_min, args.nu_max, args.nu_points)

    results = []
    T_g, c_g = args.T_guess, args.c_guess

    print(f"\nProcessing {len(common_times)} timesteps ...\n")
    hdr = (f"{'time':>12s}  {'T_eff[K]':>10s}  {'c_eff[%]':>10s}  "
           f"{'T_arith[K]':>10s}  {'c_arith[%]':>10s}  "
           f"{'T_peak[K]':>10s}  {'c_peak[%]':>10s}  "
           f"{'L_hot[cm]':>10s}  {'RMS_res':>10s}")
    print(hdr)
    print("-" * len(hdr))

    for i, t in enumerate(common_times):
      try:
        x_T = spaces_T[t] * scale
        T_prof = vals_T[t].copy()
        x_c = spaces_c[t] * scale
        c_prof = vals_c[t].copy()

        if len(x_T) != len(x_c) or not np.allclose(x_T, x_c, atol=1e-8):
            x_common = np.union1d(x_T, x_c)
            T_prof = np.interp(x_common, x_T, T_prof)
            c_prof = np.interp(x_common, x_c, c_prof)
            x_path = x_common
        else:
            x_path = x_T

        # Restrict the synthesized beam to the contiguous flame body: integrate
        # from the first sample outward, stopping at the first point where
        # x_H2O drops below args.xmin. This keeps the flame and its real shoulder
        # but excludes detached far-field H2O pockets (small x_H2O, but long and
        # hot path) that otherwise contaminate the spectrum and make the
        # retrieved T_eff/c_eff oscillate. Only the absorbance path is masked;
        # the arithmetic/peak/column-density diagnostics below are unchanged.
        x_abs, T_abs, c_abs = x_path, T_prof, c_prof
        if args.xmin > 0:
            below = np.where(c_prof < args.xmin)[0]
            if len(below):
                cut = below[0]
                x_abs, T_abs, c_abs = x_path[:cut], T_prof[:cut], c_prof[:cut]

        alpha_synth = compute_absorbance_spectrum(
            nu, x_abs, T_abs, c_abs, P=args.P)

        # Only the sensor-equivalent fit depends on the (masked) absorbance. If
        # the beam carries no significant flame signal (before ignition / near
        # burnout, or the whole window is below xmin), leave T_eff/c_eff
        # undefined but still report the arith/peak diagnostics computed below,
        # so columns 4-8 do not depend on xmin.
        if np.max(alpha_synth) < 1e-10:
            T_eff, c_eff, rms, fit = np.nan, np.nan, 0.0, None
        else:
            T_eff, c_eff, fit = extract_effective_T_c(
                nu, alpha_synth, args.L, P=args.P,
                T_guess=T_g, c_guess=c_g)
            rms = np.sqrt(np.mean(fit.fun**2))
            if np.isfinite(T_eff) and np.isfinite(c_eff):
                T_g, c_g = T_eff, c_eff

        # hot = T_prof > args.T_threshold
        # if np.any(hot):
        #     T_arith = np.mean(T_prof[hot])
        #     c_arith = np.mean(c_prof[hot])
        #     x_hot = x_path[hot]
        #     L_hot = x_hot[-1] - x_hot[0] if len(x_hot) > 1 else 0.0
        # else:
        #     T_arith = np.mean(T_prof)
        #     c_arith = np.mean(c_prof)
        #     L_hot = 0.0

        # C_mean: mean over the L region
        # T_mean: h2o-weighted mean over the L region
        inside = (x_path >= 0) & (x_path <= args.L)
        if np.any(inside):
            c_arith = np.mean(c_prof[inside])
            T_arith = np.sum(T_prof[inside] * c_prof[inside]) / np.sum(c_prof[inside]) \
                      if np.sum(c_prof[inside]) > 1e-10 else np.mean(T_prof[inside])
            L_hot = args.L
        else:
            c_arith = np.mean(c_prof)
            T_arith = np.mean(T_prof)
            L_hot = 0.0

        results.append([t, T_eff, c_eff*100, T_arith, c_arith*100,
                        np.max(T_prof), np.max(c_prof)*100, L_hot, rms])

        print(f"{t:12.6g}  {T_eff:10.1f}  {c_eff*100:10.2f}  "
              f"{T_arith:10.1f}  {c_arith*100:10.2f}  "
              f"{np.max(T_prof):10.1f}  {np.max(c_prof)*100:10.2f}  "
              f"{L_hot:10.2f}  {rms:10.2e}")

        if args.plot_every > 0 and i % args.plot_every == 0 and np.isfinite(T_eff):
            alpha_fit = simulate_uniform(nu, T_eff, c_eff, args.L, args.P)
            fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(16, 4))

            ax1b = ax1.twinx()
            ax1.plot(x_path / scale, T_prof, 'r-', lw=2, label='T(x)')
            ax1b.plot(x_path / scale, c_prof * 100, 'b--', lw=2,
                      label='c(x)')
            ax1.axhline(T_eff, color='r', ls=':', alpha=0.5,
                        label=f'T_eff={T_eff:.0f}K')
            ax1b.axhline(c_eff*100, color='b', ls=':', alpha=0.5,
                        label=f'c_eff={c_eff*100:.1f}%')
            ax1.set_xlabel(f'Position [{args.space_unit}]')
            ax1.set_ylabel('T [K]', color='r')
            ax1b.set_ylabel('H$_2$O [%]', color='b')
            ax1.set_title(f't = {t:.4g}')
            ax1.legend(loc='lower left', fontsize=7)
            ax1b.legend(loc='upper right', fontsize=7)

            ax2.plot(nu, alpha_synth, 'k-', lw=2, label='Synthetic')
            ax2.plot(nu, alpha_fit, 'r--', lw=1.5,
                     label=f'Fit T={T_eff:.0f}K c={c_eff*100:.1f}%')
            ax2.set_xlabel('Wavenumber [cm$^{-1}$]')
            ax2.set_ylabel('Absorbance')
            ax2.legend(fontsize=7)
            ax2.set_title('Absorbance spectrum')

            ax3.plot(nu, alpha_synth - alpha_fit, 'k-', lw=1)
            ax3.axhline(0, color='gray', ls='--', lw=0.5)
            ax3.set_xlabel('Wavenumber [cm$^{-1}$]')
            ax3.set_ylabel('Residual')
            ax3.set_title(f'Fit residual (RMS={rms:.2e})')

            plt.tight_layout()
            plt.savefig(os.path.join(args.output, f'spectrum_t{i:05d}.png'),
                        dpi=120, bbox_inches='tight')
            plt.close()

      except Exception as e:
        print(f"{t:12.6g}  ERROR: {e}")
        results.append([t]+[np.nan]*6+[0.0, 0.0])
        continue

    # --- Save results ---
    results = np.array(results)
    file_header = (
        f"Synthetic WMS diagnostic results (v4 - corrected partition function)\n"
        f"T file: {args.temp}\n"
        f"H2O file: {args.h2o}\n"
        f"L_assumed = {args.L} cm, P = {args.P} atm\n"
        f"T_threshold = {args.T_threshold} K\n"
        f"space_unit = {args.space_unit}\n"
        f"\n"
        f"Columns:\n"
        f"  1: time\n"
        f"  2: T_eff [K]       - effective temperature (sensor equivalent)\n"
        f"  3: c_eff [%]       - effective H2O mole fraction (sensor equiv)\n"
        f"  4: T_arith [K]     - arithmetic mean T (hot region)\n"
        f"  5: c_arith [%]     - arithmetic mean c (hot region)\n"
        f"  6: T_peak [K]      - peak T in CFD profile\n"
        f"  7: c_peak [%]      - peak c in CFD profile\n"
        f"  8: L_hot [cm]      - extent of hot region (T > threshold)\n"
        f"  9: RMS_residual    - fit quality"
    )
    outfile = os.path.join(args.output, 'effective_values.dat')
    np.savetxt(outfile, results, fmt='%14.6g', header=file_header)
    print(f"\nResults saved to: {outfile}")

    # --- Time history plot ---
    if len(results) > 1 and not np.all(np.isnan(results[:, 1])):
        valid = ~np.isnan(results[:, 1])
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)
        ax1.plot(results[valid, 0], results[valid, 1], 'ro-', ms=3,
                 label='T$_{eff}$ (sensor equivalent)')
        ax1.plot(results[valid, 0], results[valid, 3], 'b^--', ms=3,
                 label=f'T$_{{arith}}$ (T>{args.T_threshold:.0f}K)')
        ax1.plot(results[valid, 0], results[valid, 5], 'g+:', ms=4,
                 label='T$_{peak}$ (CFD max)')
        ax1.set_ylabel('Temperature [K]')
        ax1.legend()
        ax1.set_title('Comparison: sensor-equivalent vs CFD values')
        ax1.grid(True, alpha=0.3)

        ax2.plot(results[valid, 0], results[valid, 2], 'ro-', ms=3,
                 label='c$_{eff}$ (sensor equivalent)')
        ax2.plot(results[valid, 0], results[valid, 4], 'b^--', ms=3,
                 label=f'c$_{{arith}}$ (T>{args.T_threshold:.0f}K)')
        ax2.plot(results[valid, 0], results[valid, 6], 'g+:', ms=4,
                 label='c$_{peak}$ (CFD max)')
        ax2.set_ylabel('H$_2$O [%]')
        ax2.set_xlabel('Time')
        ax2.legend()
        ax2.grid(True, alpha=0.3)

        plt.tight_layout()
        figfile = os.path.join(args.output, 'time_history.png')
        plt.savefig(figfile, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Time history plot: {figfile}")

        # Integrated column density plot
        fig, ax = plt.subplots(figsize=(10, 4))
        column_densities = []
        for t in common_times:
            x_c = spaces_c[t] * scale
            c_prof = vals_c[t].copy()
            if len(x_c) != len(x_T) or not np.allclose(x_c, x_T, atol=1e-8):
                x_common = np.union1d(x_T, x_c)
                c_prof = np.interp(x_common, x_c, c_prof)
                x_path = x_common
            else:
                x_path = x_c
            column_density = compute_column_density(x_path, c_prof)
            column_densities.append(column_density)
        
        # Load the experimental data
        exp_times, exp_c = np.loadtxt("~/basilisk/basilisk-sandbox-rcaraccio/data/fatehi/yH2O-2mm", unpack=True)

        ax.plot(common_times, column_densities, 'ms-', ms=4,
                label='Integrated column density of H$_2$O')
        ax.plot(exp_times, exp_c*args.L, 'ko', ms=4, label='Experimental data')
        ax.set_xlabel('Time')
        ax.set_ylabel('Column density [cm]')
        ax.legend()
        ax.grid(True, alpha=0.3)
        figfile = os.path.join(args.output, 'column_density.png')
        plt.savefig(figfile, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Column density plot: {figfile}")

    print("\nDone.")


if __name__ == "__main__":
    main()
