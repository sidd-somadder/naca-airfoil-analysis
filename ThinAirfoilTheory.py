# This script will run Thin Airfoil Theory derivations and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# The chosen .csv file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.csv" for this script to execute successfully

# Currently, this thin airfoil theory generator can only handle NACA 4-digit series 

# For future implementation: if airfoil coordinates are provided and are not NACA series, 
# warn user that mean camber line will be approximated numerically which introduces uncertainty and noise to results

import sys;
import numpy as np;
import matplotlib.pyplot as plt;
from scipy.integrate import quad;


def run_tat_solver(input_file_name):
    target_code = "NACA";

    # Currently using fixed values for AOA range & resolution; in future user has control over this.
    inf = -5;
    sup = 10;

    # For now, use one degree spacing; create linspace of same angles in radians
    alphas = np.linspace(-5,10,(sup-inf+1));
    alphas_rad = alphas*(np.pi/180);

    # Initialize coefficients' arrays
    c_l = np.zeros_like(alphas);
    c_mLE = np.zeros_like(alphas);
    c_mqc = np.zeros_like(alphas);

    # If NACA is in file-name, use definition for mean camberline equation for TAT
    if target_code in input_file_name:
        
        # Extract specific NACA code
        code = int(input_file_name.split("_")[1]);

        ## include if-else statement to distinguise 4-digit and 5-digit airfoils

        # NACA 4-digit codes with form 00XX are symmetric airfoils; Assume no impossible airfoil codes 
        symmetric = code // 100 == 0;

        if symmetric: 
            # Use known thin airfoil theory results for symmetric airfoils.
            # Lift coefficient, Leading Edge Moment coefficient, Quarter-Chord Moment coefficient
            c_l = 2 * np.pi * alphas_rad;
            c_mLE = (-1)*(c_l)/4;
            c_mqc = np.zeros_like(alphas);

        # If not symmetric, use mean camberline equations and derive Fourier coefficients via integration
        else:
            


            # Initialize Fourier coefficients 
            A_0 = np.zeros_like(alphas);
            A_1 = np.zeros();
            A_2 = np.zeros();
    

    
    # If not a NACA code, warn user that mean camberline would need to be approximated numerically from surface points
    else:
        print("This is not a NACA code. Non-NACA Airfoils are currently not handled by this solver.") # Placeholder
        sys.exit();
    




    # temp print statements to verify values w/ calculator
    print(f"alphas : {alphas}");
    print("---");
    print(f"lift coeff : {c_l}");
    print("---");
    print(f"LE moment coeff : {c_mLE}");
    print("---");
    print(f"QC moment coeff : {c_mqc}");

    # Temporary plotting to confirm visually
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(np.degrees(alphas), c_l,   label=r'$C_L$',       color='steelblue', linewidth=2)
    ax.plot(np.degrees(alphas), c_mLE, label=r'$C_{M,LE}$',  color='firebrick', linewidth=2)
    ax.plot(np.degrees(alphas), c_mqc, label=r'$C_{M,c/4}$', color='seagreen',  linewidth=2)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)'); ax.set_ylabel('Coefficient')
    ax.set_title('TAT — Quick Check'); ax.legend(); ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout(); plt.show()

