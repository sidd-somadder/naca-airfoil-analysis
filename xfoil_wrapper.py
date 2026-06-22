# This script will run a subprocess method on the XFOIL software to validate results from the VPM and TAT scripts
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# Calling this method is expected to yield the most accurate results since it combines VPM with integral boundary layer equations

def run_xfoil_solver(input_file_name):
    print("Placeholder");