# This script is responsible for the lift, moment, and pressure coefficient
# plotting and data exporting.

import matplotlib.pyplot as plt;
import numpy as np;
import os;
import pandas as pd;

# Plots any combination of lift and moment coefficients against angle of attack.
# Pass only the series you have; omitted ones are simply not drawn.
def plot_coeffs(a, c_l_TAT=None, c_l_KJ=None, c_l_P=None, c_mLE=None, c_mqc=None, title=None):
    # (series, label, colour) -- order sets the legend order
    series = [
        (c_l_TAT, r'$C_L$ (Thin Airfoil)', 'indigo'),
        (c_l_KJ, r'$C_L$ (Kutta-Joukowski)', 'steelblue'),
        (c_l_P,  r'$C_L$ (pressure integration)', 'firebrick'),
        (c_mLE,  r'$C_{M,LE}$',  'seagreen'),
        (c_mqc,  r'$C_{M,c/4}$', 'darkorange'),
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

def plot_pressure(coeff_P_matrix, midpoints, alphas):
    x_coords = midpoints[:, 0];

    # Reverse Selig order is TE -> lower -> LE -> upper -> TE.
    # Split at the leading edge (minimum x/c) to colour the surfaces separately.
    le_idx = np.argmin(x_coords);

    print(f"Available angles of attack: {alphas}");

    while True:
        user_input = input("Enter an angle of attack to plot Cp for (or 'q' to quit): ").strip();
        if user_input.lower() == 'q':
            break;

        try:
            requested_aoa = float(user_input);
        except ValueError:
            print("Invalid input -- please enter a numeric angle or 'q'.");
            continue;

        k = np.argmin(np.abs(alphas - requested_aoa));
        matched_aoa = alphas[k];

        if not np.isclose(matched_aoa, requested_aoa, atol=1e-6):
            print(f"No exact match for {requested_aoa}°; showing closest available angle: {matched_aoa:.2f}°.");

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
        ax.set_title(rf'VPM — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)');
        ax.legend();
        ax.grid(True, linestyle=':', alpha=0.6);
        plt.tight_layout();
        plt.show();

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