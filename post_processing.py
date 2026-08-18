"""
This script is responsible for the lift, moment, and pressure coefficient plotting and data exporting. 
"""

import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd

# Plotting methods used across different scripts to display coefficient data for different computational results

def plot_coeffs(a, a_xfVisc=None, cL_TAT=None, cL_KJ=None, cL_P=None, cL_xf=None, cL_xfVisc=None, cmLE_TAT=None, cmLE_VPM=None, cmLE_xf=None, cmLE_xfVisc=None, cmqc_TAT=None, cmqc_VPM=None, cmqc_xf=None, cmqc_xfVisc=None, title=None, moment=False):
    """
    Plots any combination of lift and moment coefficients against angle of attack on a single set of axes.

    Inputs:
    Angle of Attack sweep (a) in DEGREES (size N), used as the x-axis for every series except the viscous ones.
    Angle of Attack sweep from XFOIL viscous runs (a_xfVisc) in DEGREES, passed separately because XFOIL
    omits angles where the boundary layer solution fails to converge, giving a shorter array than the input sweep.
    Any combination of coefficient arrays passed by keyword. Series left as None are skipped.

    Naming convention: TAT for thin airfoil theory, KJ and P for the two VPM lift routes,
    xf for XFOIL inviscid, and xfVisc for XFOIL viscous.

    Optional plot title (title) and a legend placement flag (moment) which moves the legend
    to the lower left for moment plots to avoid overlapping the data.

    Outputs:
    None. Displays a matplotlib figure.
    Raises ValueError if no coefficient series are provided.
    """
    # (y-series, x-series, label, colour) order sets the legend order
    series = [
        (cL_TAT, a, r'$C_L$ (Thin Airfoil)', 'seagreen'),
        (cL_KJ, a, r'$C_L$ (Kutta-Joukowski)', 'steelblue'),
        (cL_P,  a, r'$C_L$ (pressure integration)', 'firebrick'),
        (cL_xf, a, r'$C_L$ (XFOIL, Inviscid)', 'orange'),
        (cL_xfVisc, a_xfVisc, r'$C_L$ (XFOIL, Viscous, Re = $10^6$)', 'black'),
        (cmLE_TAT, a, r'$C_{M,LE}$ (Thin Airfoil)',  'seagreen'),
        (cmLE_VPM, a, r'$C_{M,LE}$ (VPM)',  'indigo'),
        (cmLE_xf, a, r'$C_{M,LE}$ (XFOIL, Inviscid)', 'darkorange'),
        (cmLE_xfVisc, a_xfVisc, r'$C_{M,LE}$ (XFOIL, Viscous, Re = $10^6$)', 'red'),
        (cmqc_TAT, a, r'$C_{M,c/4}$ (Thin Airfoil)', 'navy'),
        (cmqc_VPM, a, r'$C_{M,c/4}$ (VPM)',  'goldenrod'),
        (cmqc_xf, a, r'$C_{M,c/4}$ (XFOIL, Inviscid)', 'chartreuse'),
        (cmqc_xfVisc, a_xfVisc, r'$C_{M,c/4}$ (XFOIL, Viscous, Re = $10^6$)', 'indigo'),
    ]
 
    provided = [(y, x, lbl, col) for (y, x, lbl, col) in series if y is not None]

    if not provided:
        raise ValueError("plot_coeffs: no coefficient series provided.")

    fig, ax = plt.subplots(figsize=(9, 6))

    for y, x, lbl, col in provided:
        ax.plot(x, y, label=lbl, color=col, linewidth=2)  

    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)')
    ax.set_ylabel(r'Coefficient' if len(provided) > 1 else provided[0][2])
    ax.set_title(title if title else 'Constant Strength VPM — Aerodynamic Coefficients')

    if moment:
        ax.legend(loc='lower left')
    else:
        ax.legend(loc='lower right')
    ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.show()

def plot_pressure(coeff_P_matrix, midpoints, alphas, method='VPM'):
    """
    Interactively plots the coefficient of pressure distribution for a user-selected angle of attack.

    Inputs:
    Coefficient of Pressure distribution (coeff_P_matrix) across each panel across each angle of attack (size M x N).
    Set of panel midpoint coordinates (midpoints) in reverse-selig panel order (size M x 2).
    Angle of Attack sweep (alphas) in DEGREES (size N).
    Solver name (method) used only for the plot title.

    Outputs:
    None. Loops on user input, displaying one figure per requested angle until the user quits.
    Requested angles that do not appear in the sweep snap to the nearest available angle.
    Surfaces are split at the leading edge (minimum x/c) and coloured separately.
    Axis limits are set from finite values only, so masked trailing edge panels cannot rescale the plot.
    """
    x_coords = midpoints[:, 0]

    # Reverse Selig order is TE -> lower -> LE -> upper -> TE.
    # Split at the leading edge (minimum x/c) to colour the surfaces separately for easier visuals.
    le_idx = np.argmin(x_coords)

    # Show AOA from user defined range and move to nearest valid AOA if bad user input.
    print(f"Available angles of attack: {alphas}")

    while True:
        user_input = input("Enter an angle of attack to plot Cp for (or 'q' to quit): ").strip()
        if user_input.lower() == 'q':
            break

        # If not a float value, force user to enter float or quit.
        try:
            requested_aoa = float(user_input)
        except ValueError:
            print("Invalid input -- please enter a numeric angle or 'q'.")
            continue

        # For floats not in angle range, find the nearest AOA in user defined range.
        k = np.argmin(np.abs(alphas - requested_aoa))
        matched_aoa = alphas[k]

        if not np.isclose(matched_aoa, requested_aoa, atol=1e-6):
            print(f"No exact match for {requested_aoa}°; showing closest available angle: {matched_aoa:.2f}°.")

        # Isolate the kth set of coefficient of pressure values corresponding to the kth AOA in range
        cp = coeff_P_matrix[:, k]

        fig, ax = plt.subplots(figsize=(9, 6))

        ax.plot(x_coords[:le_idx+1], cp[:le_idx+1],
                color='blue', linewidth=2, marker='o', markersize=3, label='Lower surface')
        ax.plot(x_coords[le_idx:], cp[le_idx:],
                color='red', linewidth=2, marker='o', markersize=3, label='Upper surface')

        # Set limits from the finite data only, so a stray value can never rescale the axes.
        lo, hi = np.nanmin(cp), np.nanmax(cp)
        pad = 0.1 * (hi - lo)
        # already inverted, so no invert_yaxis()
        ax.set_ylim(hi + pad, lo - pad)

        ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
        ax.set_xlabel(r'$x/c$')
        ax.set_ylabel(r'$C_p$')
        ax.set_title(rf'{method} — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)')
        ax.legend(loc='lower right')
        ax.grid(True, linestyle=':', alpha=0.6)
        plt.tight_layout()
        plt.show()

def plot_cl_with_residuals(alphas, curves, ref_label, airfoil_name, panel_count):
    """
    Plots lift coefficient curves alongside their residuals against a chosen reference method.

    Inputs:
    Angle of Attack sweep (alphas) in DEGREES (size N). All curves must share this x-axis.
    Dictionary of lift coefficient arrays (curves) keyed by display label, each of size N.
    Key of the reference method (ref_label) that residuals are measured against; must exist in curves.
    Airfoil designation (airfoil_name) and panel count (panel_count), used for the plot title.

    Outputs:
    None. Displays a two-panel figure: absolute curves on top, residuals below on a tighter scale.
    The residual panel exists because the curves overlap almost completely at full scale.
    The reference series is omitted from the residual panel, since it would be identically zero.
    Raises KeyError if ref_label is not present in curves.
    """

    if ref_label not in curves:
        raise KeyError(f"Reference '{ref_label}' not found in curves: {list(curves.keys())}")

    ref = np.asarray(curves[ref_label])

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, figsize=(9, 8),
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.08}
    )

    # --- Top panel: absolute Cl curves ---
    for label, cl in curves.items():
        ax_top.plot(alphas, cl, linewidth=2, label=label)

    ax_top.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax_top.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax_top.set_ylabel(r'$C_L$')
    ax_top.set_title(fr'$C_L$ VPM to XFOIL Residual (NACA {airfoil_name}, {panel_count} panels)')
    ax_top.legend()
    ax_top.grid(True, linestyle=':', alpha=0.6)

    # --- Bottom panel: residuals relative to the reference method ---
    for label, cl in curves.items():
        if label == ref_label:
            continue   # reference minus itself is identically zero
        ax_bot.plot(alphas, np.asarray(cl) - ref, linewidth=2, label=label)

    ax_bot.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax_bot.set_xlabel(r'$\alpha$ (°)')
    ax_bot.set_ylabel(rf'$\Delta C_L$ vs. {ref_label}')
    ax_bot.grid(True, linestyle=':', alpha=0.6)

    plt.tight_layout()
    plt.show()

def plot_cp_comparison(coeff_P_matrix, midpoints, alphas, alpha_target, xfoil_cp_df, title=None):
    """
    Overlays the VPM and XFOIL inviscid pressure distributions for a single angle of attack.

    Inputs:
    Coefficient of Pressure distribution (coeff_P_matrix) across each panel across each angle of attack (size M x N),
    with np.nan at masked trailing edge panels.
    Set of panel midpoint coordinates (midpoints) in reverse-selig panel order (size M x 2).
    Angle of Attack sweep (alphas) in DEGREES (size N).
    Target angle to plot (alpha_target) in DEGREES; the nearest available sweep angle is used.
    XFOIL pressure data (xfoil_cp_df) from run_xfoil_cp() with columns x, y, and Cp, taken at the SAME angle.
    Optional plot title (title).

    Outputs:
    None. Displays a matplotlib figure.

    Note on ordering: this project's solver arrays run in reverse-selig order
    (TE -> lower -> LE -> upper -> TE), while XFOIL's CPWR output runs standard selig
    (TE -> upper -> LE -> lower -> TE). Since the two traversals are opposite, surfaces are
    identified by the SIGN of y rather than by array position. This assumes the lower surface
    stays below y = 0 along the full chord, which holds for the NACA sections used here but
    would misclassify points on a highly cambered airfoil near the trailing edge.

    Axis limits are scaled to the XFOIL data so that VPM oscillation at the finest panels
    cannot flatten the plot.
    """

    k = np.argmin(np.abs(alphas - alpha_target))
    matched = alphas[k]
    if not np.isclose(matched, alpha_target, atol=1e-6):
        print(f"No exact VPM match for {alpha_target}°; using {matched:.2f}°.")

    # --- VPM data, split by surface via y-sign ---
    x_v  = midpoints[:, 0]
    y_v  = midpoints[:, 1]
    cp_v = coeff_P_matrix[:, k]
    up_v = y_v >= 0

    # --- XFOIL data, split the same way ---
    x_x  = xfoil_cp_df["x"].to_numpy()
    y_x  = xfoil_cp_df["y"].to_numpy()
    cp_x = xfoil_cp_df["Cp"].to_numpy()
    up_x = y_x >= 0

    fig, ax = plt.subplots(figsize=(9, 6))

    # VPM on top, with markers so the panel-to-panel behaviour is visible.
    ax.plot(x_v[up_v],  cp_v[up_v],  color='firebrick', linewidth=1.2,
            marker='o', markersize=2.5, label='VPM, upper', alpha=0.4)
    ax.plot(x_v[~up_v], cp_v[~up_v], color='steelblue', linewidth=1.2,
            marker='o', markersize=2.5, label='VPM, lower', alpha=0.4)

    # XFOIL first, as a solid reference line underneath.
    ax.plot(x_x[up_x],  cp_x[up_x],  color='black', linewidth=2,
            label='XFOIL (inviscid), upper')
    ax.plot(x_x[~up_x], cp_x[~up_x], color='black', linewidth=2,
            linestyle='--', label='XFOIL (inviscid), lower')

    # Scale the axes to the XFOIL data, so a VPM outlier cannot flatten the plot.
    lo = min(np.nanmin(cp_x), -1.0)
    hi = max(np.nanmax(cp_x),  1.0)
    pad = 0.15 * (hi - lo)
    ax.set_ylim(hi + pad, lo - pad)       # inverted; no invert_yaxis() needed

    ax.axhline(0, color='grey', linewidth=0.8, linestyle=':')
    ax.set_xlabel(r'$x/c$')
    ax.set_ylabel(r'$C_p$')
    ax.set_title(title if title else
                 rf'Constant VPM vs XFOIL inviscid ($\alpha$ = {matched:.1f}°)')
    ax.legend(loc='upper right', fontsize=9)
    ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout()
    plt.show()

# --------------------------------
# Exporting methods used across different scripts to save coefficient data from different computation methods

def export_VPM_pressure(coeff_P_matrix, midpoints, alphas, input_file_name, method="VPM"):
    """
    Writes the vortex panel method pressure distribution to a .csv file in the computation_results folder.

    Inputs:
    Coefficient of Pressure distribution (coeff_P_matrix) across each panel across each angle of attack (size M x N).
    Set of panel midpoint coordinates (midpoints) in reverse-selig panel order (size M x 2).
    Angle of Attack sweep (alphas) in DEGREES (size N), used to label the output columns.
    Source .dat filename (input_file_name), stripped of its extension to name the output file.
    Solver name (method) used as a prefix in the output filename.

    Outputs:
    None. Writes a .csv with x/c in the first column and one Cp column per angle of attack.
    Masked trailing edge panels are dropped so the file contains no empty cells.
    """

    output_dir = os.path.join(os.path.dirname(__file__), "computation_results")
    os.makedirs(output_dir, exist_ok=True)

    identifier  = input_file_name.replace(".dat", "")
    output_path = os.path.join(output_dir, f"{method}_Cp_{identifier}_results.csv")

    x_coords = midpoints[:, 0]
    x_min = np.min(x_coords)
    chord = np.max(x_coords) - x_min
    norm_x = (x_coords - x_min) / chord

    column_names = [f"Cp_alpha_{a:.2f}" for a in alphas]
    df = pd.DataFrame(coeff_P_matrix, columns=column_names)
    df.insert(0, "x_over_c", norm_x)

    # Masked panels export as empty cells; drop them so the CSV has no gaps.
    df = df.dropna()

    df.to_csv(output_path, index=False, float_format="%.6f")
    print(f"{method} pressure distribution results exported: {output_path}")

def export_tat_results(alphas, coeffs, input_file_name, zl_ang):
    """
    Writes thin airfoil theory coefficients to a .csv file in the computation_results folder.

    Inputs:
    Angle of Attack sweep (alphas) in DEGREES (size N).
    Coefficient matrix (coeffs) with columns ordered as [c_l, c_mLE, c_mqc] (size N x 3).
    Source .dat filename (input_file_name), stripped of its extension to name the output file.
    Zero-lift angle of attack (zl_ang) in DEGREES.

    Outputs:
    None. Writes a .csv with the airfoil name and zero-lift angle as commented header lines,
    followed by one row per angle of attack.
    """
    
    output_dir = os.path.join(os.path.dirname(__file__), "computation_results")
    os.makedirs(output_dir, exist_ok=True)
    identifier = input_file_name.replace(".dat", "")
    output_path = os.path.join(output_dir, f"TAT_{identifier}_results.csv")

    # Write airfoil name & zero-lift angle as metadata comment header
    with open(output_path, "w") as f:
        f.write(f"# Thin Airfoil Theory Results: {identifier}\n")
        f.write(f"# Zero-lift angle of attack: {zl_ang:.4f} deg\n")

    # Use Pandas to format information in coefficient matrix via dataframe
    df = pd.DataFrame({
        "alpha_deg"  : alphas,
        "c_L"        : coeffs[:, 0],
        "c_M_LE"     : coeffs[:, 1],
        "c_M_QC"     : coeffs[:, 2]
    })

    # Convert coefficient matrix into .csv file and inform user of the output
    df.to_csv(output_path, index=False, float_format="%.6f", mode="a")
    print(f"TAT results exported: {output_path}")

def export_VPM_coeffs(alphas, c_l_KJ, c_l_P, c_mLE, c_mqc, input_file_name, method="CVPM"):
    """
    Writes vortex panel method lift and moment coefficients to a .csv file in the computation_results folder.

    Inputs:
    Angle of Attack sweep (alphas) in DEGREES (size N).
    Coefficient of Lift from Kutta-Joukowski circulation (c_l_KJ) and from pressure integration (c_l_P), each size N.
    Coefficient of Moment about the Leading Edge (c_mLE) and about the Quarter Chord (c_mqc), each size N.
    Source .dat filename (input_file_name), stripped of its extension to name the output file.
    Solver name (method) used as a prefix in the output filename.

    Outputs:
    None. Writes a .csv with one row per angle of attack.
    Note: c_l_KJ and c_l_P are computed over slightly different panel sets, since the pressure
    integration excludes the masked trailing edge panels while the circulation sum does not.
    A small residual difference between the two columns is therefore expected.
    """

    output_dir = os.path.join(os.path.dirname(__file__), "computation_results")
    os.makedirs(output_dir, exist_ok=True)

    identifier = input_file_name.replace(".dat", "")
    output_path = os.path.join(output_dir, f"{method}_coeffs_{identifier}_results.csv")

    df = pd.DataFrame({
        "alpha_deg" : alphas,
        "c_L_KJ"    : c_l_KJ,
        "c_L_P"     : c_l_P,
        "c_M_LE"    : c_mLE,
        "c_M_QC"    : c_mqc
    })

    df.to_csv(output_path, index=False, float_format="%.6f")
    print(f"{method} coefficient results exported: {output_path}")