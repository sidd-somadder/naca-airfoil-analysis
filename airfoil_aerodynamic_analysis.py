# This is the master script which the user will interact with. The user will never need to open TAT, VPM, or XFOIL scripts.
# The user will be able to select a saved airfoil coordinates file from the folder and choose to run an aerodynamic analysis
# The master script will call all three scripts for results and plot the c_l, c_m,LE, and c_m,c/4 values against AOA.
# The user can choose to save the results of all three scripts in the computation_results folder.

from thin_airfoil_theory import run_tat_solver
from constant_vpm import run_cvpm_solver
from xfoil_wrapper import run_xfoil_solver;
from panel_geometry import load_dat_coordinates
from post_processing import plot_coeffs, plot_cl_with_residuals;
import numpy as np;
import os;

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

# Fill THIS array with filenames for the master script to loop through/
filenames = ["NACA_0012_N100.dat",
             "NACA_2412_N100.dat"];

# Initialize condition number array. 
cond_nums = np.zeros(len(filenames));

for k, file in enumerate(filenames):
    dat_path = os.path.join(os.path.dirname(__file__), "saved_airfoil_coords", file);

    code_name = file.split("_")[1];
    panel_count = 2  *int(((file.split("_")[2]).split(".")[0])[1:]) - 1;

    title_name = f"NACA {code_name}, {panel_count} panels"

    # Call from-scratch VPM solver
    geom_points = load_dat_coordinates(file);
    cL_KJ, cL_P, cmLE_VPM, cmqc_VPM, coeff_P_matrix, cond = run_cvpm_solver(geom_points, angle_param, file);

    cond_nums[k] = cond;

    # XFOIL on the same file
    cl_xf, cmLE_xf, cmqc_xf, a_xf = run_xfoil_solver(dat_path, angle_param, len(geom_points));
    cl_xfVisc, cmLE_xfVisc, cmqc_xfVisc, a_xfVisc = run_xfoil_solver(dat_path, angle_param, len(geom_points), viscous=True);
    cl_tat, cmLE_tat, cmqc_tat = run_tat_solver(file, angle_param);

    curves = {
    'Kutta-Joukowski':      cL_KJ,
    'Pressure Integration': cL_P,
    'XFOIL, Inviscid':      cl_xf,
    }

    plot_coeffs(angle_param, a_xfVisc=a_xfVisc, cL_TAT=cl_tat, cL_KJ=cL_KJ, cL_P=cL_P, cL_xf=cl_xf, cL_xfVisc=cl_xfVisc, title=fr'$C_L$ vs. $\alpha$ ({title_name})')
    plot_coeffs(angle_param, a_xfVisc=a_xfVisc, cmLE_TAT=cmLE_tat, cmLE_VPM=cmLE_VPM, cmLE_xf=cmLE_xf, title=rf'$C_{{M,LE}}$ vs. $\alpha$ ({title_name})', moment=True)
    plot_coeffs(angle_param, a_xfVisc=a_xfVisc, cmqc_TAT=cmqc_tat, cmqc_VPM=cmqc_VPM, cmqc_xf=cmqc_xf, title=rf'$C_{{M,c/4}}$  vs. $\alpha$ ({title_name})', moment=True)
    plot_cl_with_residuals(angle_param, curves, ref_label='XFOIL, Inviscid', airfoil_name=code_name, panel_count=panel_count);

print(f"Condition Number: {cond_nums}");