# This script will run Thin Airfoil Theory derivations and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# The chosen .csv file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.csv" for this script to execute successfully

# Currently, this thin airfoil theory generator can only handle NACA 4-digit series 

# For future implementation: if airfoil coordinates are provided and are not NACA series, 
# warn user that mean camber line will be approximated numerically which introduces uncertainty and noise to results

# Master 
def run_tat_solver(input_file_name):
    target_code = "NACA";
    print(input_file_name) # Placeholder
    if target_code in input_file_name:
        print("This is a NACA code.");
    else:
        print("This is not a NACA code.")
