# This script will run Thin Airfoil Theory derivations and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# The chosen .csv file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.csv" for this script to execute successfully

# Currently, this thin airfoil theory generator can only handle NACA 4-digit series 

# For future implementation: if airfoil coordinates are provided and are not NACA series, 
# warn user that mean camber line will be approximated via spline which introduces uncertainty and noise to results

import sys;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.integrate as scpi;

def run_tat_solver(input_file_name, alphas):
    target_code = "NACA";

    # Since angle linespace creation is handled by the master script; extract angles in radians for computation.
    alphas_rad = alphas*(np.pi/180);

    # Initialize coefficients' arrays
    c_l = np.zeros_like(alphas);
    c_mLE = np.zeros_like(alphas);
    c_mqc = np.zeros_like(alphas);

    # If NACA is in file-name, use definition for mean camberline equation for TAT
    if target_code in input_file_name:     
        # Extract specific NACA code
        code_name = input_file_name.split("_")[1];
        # Check if NACA code has any letters/other characters
        try:
            code = int(code_name);
        except ValueError:
            # If error occurs (has non-numbers in code), 
            # warn user that mean camberline would need to be approximated via spline numerically from surface points
            result = handle_tat_fallback(input_file_name, alphas_rad);
            if result is None:
                return None;
            coeffs, angle_zero_lift = result;
        else:
            # If code is 4-digits, proceed with 4-digit algorithm.
            if len(code_name) == 4:
                # NACA 4-digit codes with form 00XX are symmetric airfoils; Assume for now no impossible airfoil codes 
                symmetric = code // 100 == 0;

                if symmetric: 
                    # Call SYMMETRIC solver function
                    coeffs, angle_zero_lift = sym_4digit_solver(code, alphas_rad)
                else:
                    # Call ASYMMETRIC solver function (integrates for Fourier coefficients and derives coefficients)
                    coeffs, angle_zero_lift = asym_4digit_solver(code, alphas_rad)
            
            # If a 5-digit airfoil, proceed with deriving MCL equation with 5-digit algorithm 
            elif len(code_name) == 5:
                # First search if name belongs to a standard series
                # If not, numerically derive parameters using function in NACA_geometric plotter script
                print("5-digit NACA airfoils not currently available"); # Placeholder

            # If NACA code doesn't abide by 4-digit or 5-digit series format, 
            # warn user that mean camberline would need to be approximated via spline numerically from surface points
            else: 
                result = handle_tat_fallback(input_file_name, alphas_rad);
                if result is None:
                    return None;
                coeffs, angle_zero_lift = result;
               
    # If not a NACA code, warn user that mean camberline would need to be approximated via spline numerically from surface points
    else:
        result = handle_tat_fallback(input_file_name, alphas_rad);
        if result is None:
            return None;
        coeffs, angle_zero_lift = result;
        # To be implemented

    # Extract coefficients from matrix
    c_l = coeffs[:,0];
    c_mLE = coeffs[:,1];
    c_mqc = coeffs[:,2];

    # print & plot values for verification (temp)
    printvals(alphas, c_l, c_mLE, c_mqc, angle_zero_lift);
    coeff_visualizer(alphas, c_l, c_mLE, c_mqc);

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

def sym_4digit_solver(code, alphas_rad):
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

def spline_TAT_solver(file, alphas_rad):
        # Spline solver on coordinates
        # placeholder
        print(file);
        c_l = np.zeros_like(alphas_rad);
        c_mLE = np.zeros_like(alphas_rad);
        c_mqc = np.zeros_like(alphas_rad);
        zero_lift_angle = 0;
        return np.column_stack((c_l, c_mLE, c_mqc)), zero_lift_angle;

# When non-NACA or improper NACA inputs are given in file name, ask if user wants to proceed/skip with approximated TAT results
def ask_user_tat_fallback():
    print("WARNING: Non-NACA or improper NACA airfoil name detected. TAT requires a mean camberline equation.")
    print("Options:")
    print("  [1] Proceed with numerical camberline approximation (results may have higher error)")
    print("  [2] Skip TAT solver and continue with VPM and XFOIL only")
    
    # Repeat until proper input is provided by user.
    while True:
        choice = input("Enter 1 or 2: ").strip()
        if choice == "1":
            return "numerical"
        elif choice == "2":
            return "skip"
        else:
            print("Invalid input, please enter 1 or 2.")

# Driver function for ask_user_tat_fallback() to decide what to do if non-NACA or improper NACA code 
# Kept separate from aforementioned function for manual testing reasons
def handle_tat_fallback(input_file_name, alphas_rad):
    # Call ask_user_tat_fallback() to give warning & get decision on spline solver or skipping
    decision = ask_user_tat_fallback()
    if decision == "skip":
        # Master script interprets "None" to skip TAT visualization 
        return None
    elif decision == "numerical":
        # Calls spline_TAT_solver returning coefficient matrix and zero-lift angle as a tuple 
        return spline_TAT_solver(input_file_name, alphas_rad)