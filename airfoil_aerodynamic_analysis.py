# This is the master script which the user will interact with. The user will never need to open TAT, VPM, or XFOIL scripts.
# The user will be able to select a saved airfoil coordinates file from the folder and choose to run an aerodynamic analysis
# The master script will call all three scripts for results and plot the c_l, c_m,LE, and c_m,c/4 values against AOA.
# The user can choose to save the results of all three scripts in the computation_results folder.

from thin_airfoil_theory import run_tat_solver
from constant_vpm import run_cvpm_solver
from xfoil_wrapper import run_xfoil_solver, run_xfoil_cp;
from panel_geometry import load_dat_coordinates
from post_processing import plot_coeffs, plot_cl_with_residuals, plot_cp_comparison;
import numpy as np;
import os;

def get_file_name():
    print();

# Function used to get user input for desired angle of attack mesh
def get_angle_params():
    print("-" * 10)
    print("Angle of Attack Range Configuration")
    print("-" * 3)
    print("Define the range over which all solvers will run.")
    print("Results will be computed and plotted across this range.")
    print("(All values in degrees)")
    print("-" * 10)

    # Upper and lower bounds of range is defined by user input
    try:
        inf = float(input("  Lower bound in whole degrees (e.g. -5): "))
    except ValueError:
        print("Invalid input. Using -5 deg. as default lower bound");
        inf = float(-5);

    try:
        sup = float(input("  Upper bound  (e.g. 15): "))
    except ValueError:
        print("Invalid input. Using 15 deg. as default upper bound");
        sup = float(15);


    # User determines step-size in degrees
    print()
    print("Step size determines how many points are computed per degree.")
    print("Finer step size = more points per degree, e.g.:")
    print("  1 → one point per degree  (fast, coarse)")
    print("  0.5 → two points per degree     (balanced)")
    print("  0.25 → four points per degree    (fine)")

    # Resolution of angular range expressed as step size
    try:
        step = float(input("  Step size in degrees (e.g. 1, 0.5, 0.25): "))
    except ValueError:
        print("Invalid input. Using default 1 degree step size.")
        step = 1.0;

    # Inform user of their desired range. Confirm with user. 
    while True:
    # Read user input, strip whitespace, and convert to uppercase
        choice = input(f"Solvers will run from {inf}° to {sup}° at {step}° intervals. Proceed? (Y/N)").strip().upper()
        if choice == 'Y':
            # If yes, return linspace from user-input angle parameters
            return np.arange(inf, sup + step, step)
        
        elif choice == 'N':
            # If no, ask user for angle parameters number via recursive function call;
            return get_angle_params();     
        else:
            print("Invalid choice. Please enter Y or N.");


angle_param = get_angle_params();

filename = ["NACA_0012_N100.dat",
            "NACA_2412_N100.dat"];

for k in filename:
    dat_path = os.path.join(os.path.dirname(__file__), "saved_airfoil_coords", k);

    code_name = k.split("_")[1];

    title_name = f"NACA {code_name}"

    # your solver
    geom_points = load_dat_coordinates(k);
    cL_KJ, cL_P, coeff_P_matrix = run_cvpm_solver(geom_points, angle_param, k);

    # XFOIL on the same file
    cl_xf, cmLE_xf, cmqc_xfv = run_xfoil_solver(dat_path, angle_param, len(geom_points));
    cl_xfV, cmLE_xfv, cmLE_xfV = run_xfoil_solver(dat_path, angle_param, len(geom_points), viscous=True);


    # Temporary sample file names to test cross-script function calls

    # Should print that input is a NACA airfoil
    #print(f"Using sample {sample_file_name1} : ")
    cl_tat, cmLE_tat, cmqc_tat = run_tat_solver(k, angle_param);
    #print("---");

    curves = {
    'Thin Airfoil':         cl_tat,
    'Kutta-Joukowski':      cL_KJ,
    'Pressure Integration': cL_P,
    'XFOIL, Inviscid':      cl_xf,
    }

    plot_coeffs(angle_param, cL_TAT=cl_tat, cL_KJ=cL_KJ, cL_P=cL_P, cL_xf=cl_xf, title=f'Lift Coefficient Results - {title_name}')
    plot_cl_with_residuals(angle_param, curves, ref_label='XFOIL, Inviscid', airfoil_name=code_name);

    xf_cp = run_xfoil_cp(dat_path, alpha=5.0, n_panels=len(geom_points));


# # Should print that input is not a NACA airfoil
# print(f"Using sample {sample_file_name2} : ")
# run_tat_solver(sample_file_name2);

# # Placeholder function calls
# run_vpm_solver(sample_file_name1);
# run_xfoil_solver(sample_file_name1);
