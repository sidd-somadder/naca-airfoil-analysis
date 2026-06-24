# This is the master script which the user will interact with. The user will never need to open TAT, VPM, or XFOIL scripts.
# The user will be able to select a saved airfoil coordinates file from the folder and choose to run an aerodynamic analysis
# The master script will call all three scripts for results and plot the c_l, c_m,LE, and c_m,c/4 values against AOA.
# The user can choose to save the results of all three scripts in the computation_results folder.

from ThinAirfoilTheory import run_tat_solver
from VortexPanelMethod import run_vpm_solver
from xfoil_wrapper import run_xfoil_solver

# Temporary sample file names to test cross-script function calls
sample_file_name1 = "NACA_4412_N100.dat";
sample_file_name2 = "Clarky_N60.dat";

# Should print that input is a NACA airfoil
print(f"Using sample {sample_file_name1} : ")
run_tat_solver(sample_file_name1);
print("---");

# # Should print that input is not a NACA airfoil
# print(f"Using sample {sample_file_name2} : ")
# run_tat_solver(sample_file_name2);

# Placeholder function calls
run_vpm_solver(sample_file_name1);
run_xfoil_solver(sample_file_name1);
