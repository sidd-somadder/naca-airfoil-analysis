# This script will run Thin Airfoil Theory derivations and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# The chosen .csv file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.dat" for this script to execute successfully

# Currently, this thin airfoil theory solver's scope is to only handle NACA 4-digit series 

import sys;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.integrate as scpi;
import pandas as pd;
import os;


def run_tat_solver(input_file_name, alphas):
    alphas_rad = alphas * (np.pi / 180);

    # Extract the numeric code from the required "NACA_XXXX_N#.csv" naming convention.
    code_name = input_file_name.split("_")[1];

    try:
        code = int(code_name);
    except ValueError:
        print(f"TAT solver skipped: '{code_name}' is not a valid NACA numeric code.");
        return;

    if len(code_name) != 4:
        print(f"TAT solver skipped: '{code_name}' is not a 4-digit NACA code.");
        return;

    symmetric = code // 100 == 0;
    if symmetric:
        coeffs, angle_zero_lift = sym_4digit_solver(alphas_rad);
    else:
        coeffs, angle_zero_lift = asym_4digit_solver(code, alphas_rad);

    c_l = coeffs[:,0];
    c_mLE = coeffs[:,1];
    c_mqc = coeffs[:,2];

    printvals(alphas, c_l, c_mLE, c_mqc, angle_zero_lift);
    coeff_visualizer(alphas, c_l, c_mLE, c_mqc);
    export_tat_results(alphas, coeffs, input_file_name, angle_zero_lift);

# temp print statements to verify values w/ calculator
def printvals(a, l, mLE, mqc, zla):
    print(f"alphas : {a}");
    print("---");
    print(f"Zero-Lift Angle = {zla} deg.");
    print("---");
    print(f"lift coeff : {l}");
    print("---");
    print(f"LE moment coeff : {mLE}");
    print("---");
    print(f"QC moment coeff : {mqc}");

# Temporary plotting to confirm visually
def coeff_visualizer(a, l, mLE, mqc):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(a, l,   label=r'$C_L$',       color='steelblue', linewidth=2)
    ax.plot(a, mLE, label=r'$C_{M,LE}$',  color='firebrick', linewidth=2)
    ax.plot(a, mqc, label=r'$C_{M,c/4}$', color='seagreen',  linewidth=2)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)'); ax.set_ylabel('Coefficient')
    ax.set_title('TAT — Quick Check'); ax.legend(); ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout(); plt.show()

def sym_4digit_solver(alphas_rad):
    # Use known thin airfoil theory results for symmetric airfoils.
    c_l = 2 * np.pi * alphas_rad;
    c_mLE = (-1)*(c_l)/4;
    c_mqc = np.zeros_like(alphas_rad);
    zero_lift_angle = 0;

    return np.column_stack((c_l, c_mLE, c_mqc)), zero_lift_angle;

def asym_4digit_solver(code, alphas_rad):
    # Extract digits from 4-digit series code
    m_dig = code // 1000;
    p_dig = ((code % 1000) // 100); 

    # Define mean camberline slope equations forward and aft of max camber location
    # Apply Glauert transformation (x/c = 0.5(1-cos(theta_0))) 
    def dz1(t0,m,p):
        return (2*m)/(p**2) * (p - 0.5*(1-np.cos(t0)));
            
    def dz2(t0,m,p):
        return (2*m)/((1-p)**2) * (p - 0.5*(1-np.cos(t0)));
                
    # Convert NACA digits into chordwise % values
    m = m_dig / 100;
    p = p_dig / 10;

    # Initialize Fourier coefficients 
    A_0 = np.zeros_like(alphas_rad);
    A_1 = np.zeros(1);
    A_2 = np.zeros(1);

    # Use Glauert transformation to find angle of max camber position
    t_c = np.arccos(1-2*p);

    # Integrate for Fourier coefficients
    A_0 = alphas_rad - (1/np.pi) * (scpi.quad(dz1, 0, t_c, args=(m,p))[0] + scpi.quad(dz2, t_c, np.pi, args=(m,p))[0]); 

    A_1 = (2/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * np.cos(t0), 0, t_c)[0] 
                        + scpi.quad(lambda t0: dz2(t0, m, p) * np.cos(t0), t_c, np.pi)[0]);
    
    A_2 = (2/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * np.cos(2*t0), 0, t_c)[0] 
                        + scpi.quad(lambda t0: dz2(t0, m, p) * np.cos(2*t0), t_c, np.pi)[0]);

    # Use Fourier coefficients to derive lift and moments coefficients of asymmetric airfoils
    c_l = np.pi * (2*A_0 + A_1);
    c_mLE = (-1)*(c_l/4) - (np.pi/4)*(A_1 - A_2);
    # Note, quarter-chord moment coefficient is theoretically constant
    c_mqc = np.full_like(alphas_rad, (np.pi/4) * (A_2 - A_1), dtype=float);

    # Compute zero lift angle of attack by expanding c_l and setting to 0
    zero_lift_angle = (1/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * (1-np.cos(t0)), 0, t_c)[0] 
                                + scpi.quad(lambda t0: dz2(t0, m, p) * (1-np.cos(t0)), t_c, np.pi)[0]);

    # express zero lift angle in degrees
    zero_lift_angle = (180/np.pi)*zero_lift_angle;

    # Return coefficients in matrix form and zero lift angle of attack as a tuple
    return np.column_stack((c_l, c_mLE, c_mqc)), zero_lift_angle;
    
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
