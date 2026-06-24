# This script will run Thin Airfoil Theory derivations and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# The chosen .csv file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.csv" for this script to execute successfully

# Currently, this thin airfoil theory generator can only handle NACA 4-digit series 

# For future implementation: if airfoil coordinates are provided and are not NACA series, 
# warn user that mean camber line will be approximated numerically which introduces uncertainty and noise to results

import sys;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.integrate as scpi;


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
            A_0 = np.zeros_like(alphas);
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
            c_mqc = c_mqc = np.full_like(alphas, (np.pi/4) * (A_2 - A_1), dtype=float);
               
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
    ax.plot(alphas, c_l,   label=r'$C_L$',       color='steelblue', linewidth=2)
    ax.plot(alphas, c_mLE, label=r'$C_{M,LE}$',  color='firebrick', linewidth=2)
    ax.plot(alphas, c_mqc, label=r'$C_{M,c/4}$', color='seagreen',  linewidth=2)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)'); ax.set_ylabel('Coefficient')
    ax.set_title('TAT — Quick Check'); ax.legend(); ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout(); plt.show()

