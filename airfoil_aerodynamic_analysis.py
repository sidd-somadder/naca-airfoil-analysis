# This is the master script which the user will interact with. The user will never need to open TAT, VPM, or XFOIL scripts.
# The user will be able to select a saved airfoil coordinates file from the folder and choose to run an aerodynamic analysis
# The master script will call all three scripts for results and plot the c_l, c_m,LE, and c_m,c/4 values against AOA.
# The user can choose to save the results of all three scripts in the computation_results folder.

from ThinAirfoilTheory import run_tat_solver
from VortexPanelMethod import run_vpm_solver
from xfoil_wrapper import run_xfoil_solver
import numpy as np;

def get_angle_params():
    print("-" * 10)
    print("Angle of Attack Range Configuration")
    print("-" * 3)
    print("Define the range over which all solvers will run.")
    print("Results will be computed and plotted across this range.")
    print("(All values in degrees)")
    print("-" * 10)

    inf = float(input("  Lower bound in whole degrees (e.g. -5): "))
    sup = float(input("  Upper bound  (e.g. 15): "))

    print()
    print("Step size determines how many points are computed per degree.")
    print("Finer step size = more points per degree, e.g.:")
    print("  1 → one point per degree  (fast, coarse)")
    print("  0.5 → two points per degree     (balanced)")
    print("  0.25 → four points per degree    (fine)")

    # Resolution of angular range expressed as step size
    step = float(input("  Step size in degrees (e.g. 1, 0.5, 0.25): "))

    print("-" * 10)
    print(f"Running solvers from {inf}° to {sup}° at {step}° intervals.")
    print("=" * 10)
    print()

    return np.arange(inf, sup + step, step)

angle_param = get_angle_params();

# Temporary sample file names to test cross-script function calls
sample_file_name1 = "ClarkY_N100.dat";
sample_file_name2 = "Clarky_N60.dat";

# Should print that input is a NACA airfoil
print(f"Using sample {sample_file_name1} : ")
run_tat_solver(sample_file_name1, angle_param);
print("---");

# # Should print that input is not a NACA airfoil
# print(f"Using sample {sample_file_name2} : ")
# run_tat_solver(sample_file_name2);

# # Placeholder function calls
# run_vpm_solver(sample_file_name1);
# run_xfoil_solver(sample_file_name1);
