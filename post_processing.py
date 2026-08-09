# This script is responsible for the lift, moment, and pressure coefficient
# plotting and data exporting.

import matplotlib.pyplot as plt;
import numpy as np;
import os;
import pandas as pd;

# Plotting methods used across different scripts to display coefficient data for different computational results

# Plots any combination of lift and moment coefficients against angle of attack.
# Pass only the series you have; omitted ones are simply not drawn.
def plot_coeffs(a, cL_TAT=None, cL_KJ=None, cL_P=None, cL_xf=None, cL_xfV=None, cmLE_TAT=None, cmLE_VPM=None, cmLE_xf=None, cmLE_xfV=None, cmqc_TAT=None, cmqc_VPM=None, cmqc_xf=None, cmqc_xfV=None, title=None):
    # (series, label, colour) -- order sets the legend order
    series = [
        (cL_TAT, r'$C_L$ (Thin Airfoil)', 'indigo'),
        (cL_KJ, r'$C_L$ (Kutta-Joukowski)', 'steelblue'),
        (cL_P,  r'$C_L$ (pressure integration)', 'firebrick'),
        (cL_xf,  r'$C_L$ (XFOIL, Inviscid)', 'fuchsia'),
        (cL_xfV,  r'$C_L$ (XFOIL, Viscous)', 'blue'),
        (cmLE_TAT,  r'$C_{M,LE} (Thin Airfoil)$',  'seagreen'),
        (cmLE_VPM,  r'$C_{M,LE}$ (VPM)',  'slategrey'),
        (cmLE_xf,  r'$C_{M,LE}$ (XFOIL, Invsicid)', 'darkorange'),
        (cmLE_xfV,  r'$C_{M,LE}$ (XFOIL, Viscous)', 'red'),
        (cmqc_TAT,  r'$C_{M,c/4}$ (Thin Airfoil)', 'navy'),
        (cmqc_VPM,  r'$C_{M,c/4}$ (VPM)',  'goldenrod'),
        (cmqc_xf,  r'$C_{M,c/4}$ (XFOIL, Inviscid)', 'chartreuse'),
        (cmqc_xfV,  r'$C_{M,c/4}$ (XFOIL, Viscous)', 'indigo'),
    ];

    provided = [(y, lbl, col) for (y, lbl, col) in series if y is not None];

    if not provided:
        raise ValueError("plot_coeffs: no coefficient series provided.");

    fig, ax = plt.subplots(figsize=(9, 6));

    for y, lbl, col in provided:
        ax.plot(a, y, label=lbl, color=col, linewidth=2);

    ax.axhline(0, color='black', linewidth=0.8, linestyle='--');
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--');
    ax.set_xlabel(r'$\alpha$ (°)');
    ax.set_ylabel('Coefficient' if len(provided) > 1 else provided[0][1]);
    ax.set_title(title if title else 'Constant Strength VPM — Aerodynamic Coefficients');
    ax.legend();
    ax.grid(True, linestyle=':', alpha=0.6);
    plt.tight_layout();
    plt.show();

# Plots pressure against a specified angle of attack until user decides to quit
def plot_pressure(coeff_P_matrix, midpoints, alphas, method='VPM'):
    x_coords = midpoints[:, 0];

    # Reverse Selig order is TE -> lower -> LE -> upper -> TE.
    # Split at the leading edge (minimum x/c) to colour the surfaces separately for easier visuals.
    le_idx = np.argmin(x_coords);

    # Show AOA from user defined range and move to nearest valid AOA if bad user input.
    print(f"Available angles of attack: {alphas}");

    while True:
        user_input = input("Enter an angle of attack to plot Cp for (or 'q' to quit): ").strip();
        if user_input.lower() == 'q':
            break;

        # If not a float value, force user to enter float or quit.
        try:
            requested_aoa = float(user_input);
        except ValueError:
            print("Invalid input -- please enter a numeric angle or 'q'.");
            continue;

        # For floats not in angle range, find the nearest AOA in user defined range.
        k = np.argmin(np.abs(alphas - requested_aoa));
        matched_aoa = alphas[k];

        if not np.isclose(matched_aoa, requested_aoa, atol=1e-6):
            print(f"No exact match for {requested_aoa}°; showing closest available angle: {matched_aoa:.2f}°.");

        # Isolate the kth set of coefficient of pressure values corresponding to the kth AOA in range
        cp = coeff_P_matrix[:, k];

        fig, ax = plt.subplots(figsize=(9, 6));

        ax.plot(x_coords[:le_idx+1], cp[:le_idx+1],
                color='blue', linewidth=2, marker='o', markersize=3, label='Lower surface');
        ax.plot(x_coords[le_idx:], cp[le_idx:],
                color='red', linewidth=2, marker='o', markersize=3, label='Upper surface');

        # Set limits from the finite data only, so a stray value can never rescale the axes.
        lo, hi = np.nanmin(cp), np.nanmax(cp);
        pad = 0.1 * (hi - lo);
        ax.set_ylim(hi + pad, lo - pad);          # already inverted, so no invert_yaxis()

        ax.axhline(0, color='black', linewidth=0.8, linestyle='--');
        ax.set_xlabel(r'$x/c$');
        ax.set_ylabel(r'$C_p$');
        ax.set_title(rf'{method} — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)');
        ax.legend();
        ax.grid(True, linestyle=':', alpha=0.6);
        plt.tight_layout();
        plt.show();

# --------------------------------
# Exporting methods used across different scripts to save coefficient data from different computation methods


def export_VPM_pressure(coeff_P_matrix, midpoints, alphas, input_file_name, method="VPM"):
    output_dir = os.path.join(os.path.dirname(__file__), "computation_results");
    os.makedirs(output_dir, exist_ok=True);

    identifier  = input_file_name.replace(".dat", "");
    output_path = os.path.join(output_dir, f"{method}_Cp_{identifier}_results.csv");

    x_coords = midpoints[:, 0];
    chord    = np.max(x_coords) - np.min(x_coords);
    norm_x   = x_coords / chord;

    column_names = [f"Cp_alpha_{a:.2f}" for a in alphas];
    df = pd.DataFrame(coeff_P_matrix, columns=column_names);
    df.insert(0, "x_over_c", norm_x);

    # Masked panels export as empty cells; drop them so the CSV has no gaps.
    df = df.dropna();

    df.to_csv(output_path, index=False, float_format="%.6f");
    print(f"{method} pressure distribution results exported: {output_path}");

# Function that saves aerodynamics coefficients to a .csv file for master script to plot
def export_tat_results(alphas, coeffs, input_file_name, zl_ang):
    
    # Define the output folder relative to the directory of this script
    output_dir = os.path.join(os.path.dirname(__file__), "computation_results")
    os.makedirs(output_dir, exist_ok=True)

    # Strip .dat extension for clean output naming
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

def plot_cl_with_residuals(alphas, curves, ref_label, airfoil_name):
    if ref_label not in curves:
        raise KeyError(f"Reference '{ref_label}' not found in curves: {list(curves.keys())}");

    ref = np.asarray(curves[ref_label]);

    fig, (ax_top, ax_bot) = plt.subplots(
        2, 1, sharex=True, figsize=(9, 8),
        gridspec_kw={'height_ratios': [3, 1], 'hspace': 0.08}
    );

    # --- Top panel: absolute Cl curves ---
    for label, cl in curves.items():
        ax_top.plot(alphas, cl, linewidth=2, label=label);

    ax_top.axhline(0, color='black', linewidth=0.8, linestyle='--');
    ax_top.axvline(0, color='black', linewidth=0.8, linestyle='--');
    ax_top.set_ylabel(r'$C_L$');
    ax_top.set_title(f'Lift Coefficient Results — NACA {airfoil_name}');
    ax_top.legend();
    ax_top.grid(True, linestyle=':', alpha=0.6);

    # --- Bottom panel: residuals relative to the reference method ---
    for label, cl in curves.items():
        if label == ref_label:
            continue;   # reference minus itself is identically zero
        ax_bot.plot(alphas, np.asarray(cl) - ref, linewidth=2, label=label);

    ax_bot.axhline(0, color='black', linewidth=0.8, linestyle='--');
    ax_bot.set_xlabel(r'$\alpha$ (°)');
    ax_bot.set_ylabel(rf'$\Delta C_L$ vs. {ref_label}');
    ax_bot.grid(True, linestyle=':', alpha=0.6);

    plt.tight_layout();
    plt.show();

# Overlays VPM and XFOIL inviscid Cp distributions for a single angle of attack.
#
# NOTE ON ORDERING: this project's solver arrays run in reverse Selig order
# (TE -> lower -> LE -> upper -> TE), while XFOIL's CPWR output runs standard
# Selig (TE -> upper -> LE -> lower -> TE). The two traversals are opposite, so
# each source is split at its own leading edge and the surfaces are identified
# by the SIGN of y rather than by array position.
def plot_cp_comparison(coeff_P_matrix, midpoints, alphas, alpha_target,
                       xfoil_cp_df, title=None):
    """
    Parameters
    ----------
    coeff_P_matrix : (M, N) ndarray
        VPM Cp, full length M with np.nan at invalid panels.
    midpoints : (M, 2) ndarray
        VPM collocation points, full length M.
    alphas : (N,) ndarray
        Angles of attack in DEGREES (the VPM sweep).
    alpha_target : float
        Angle to plot. The closest available VPM angle is used.
    xfoil_cp_df : DataFrame
        From run_xfoil_cp(), with columns x, y, Cp -- taken at the SAME angle.
    """
    k = np.argmin(np.abs(alphas - alpha_target));
    matched = alphas[k];
    if not np.isclose(matched, alpha_target, atol=1e-6):
        print(f"No exact VPM match for {alpha_target}°; using {matched:.2f}°.");

    # --- VPM data, split by surface via y-sign ---
    x_v  = midpoints[:, 0];
    y_v  = midpoints[:, 1];
    cp_v = coeff_P_matrix[:, k];
    up_v = y_v >= 0;

    # --- XFOIL data, split the same way ---
    x_x  = xfoil_cp_df["x"].to_numpy();
    y_x  = xfoil_cp_df["y"].to_numpy();
    cp_x = xfoil_cp_df["Cp"].to_numpy();
    up_x = y_x >= 0;

    fig, ax = plt.subplots(figsize=(9, 6));

    # XFOIL first, as a solid reference line underneath.
    ax.plot(x_x[up_x],  cp_x[up_x],  color='black', linewidth=2,
            label='XFOIL (inviscid), upper');
    ax.plot(x_x[~up_x], cp_x[~up_x], color='black', linewidth=2,
            linestyle='--', label='XFOIL (inviscid), lower');

    # VPM on top, with markers so the panel-to-panel behaviour is visible.
    ax.plot(x_v[up_v],  cp_v[up_v],  color='firebrick', linewidth=1.2,
            marker='o', markersize=2.5, label='VPM, upper');
    ax.plot(x_v[~up_v], cp_v[~up_v], color='steelblue', linewidth=1.2,
            marker='o', markersize=2.5, label='VPM, lower');

    # Scale the axes to the XFOIL data, so a VPM outlier cannot flatten the plot.
    lo = min(np.nanmin(cp_x), -1.0);
    hi = max(np.nanmax(cp_x),  1.0);
    pad = 0.15 * (hi - lo);
    ax.set_ylim(hi + pad, lo - pad);       # inverted; no invert_yaxis() needed

    ax.axhline(0, color='grey', linewidth=0.8, linestyle=':');
    ax.set_xlabel(r'$x/c$');
    ax.set_ylabel(r'$C_p$');
    ax.set_title(title if title else
                 rf'Constant VPM vs XFOIL inviscid ($\alpha$ = {matched:.1f}°)');
    ax.legend(loc='lower right', fontsize=9);
    ax.grid(True, linestyle=':', alpha=0.6);
    plt.tight_layout();
    plt.show();

def plot_cond():
    print();