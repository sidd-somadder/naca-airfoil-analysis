# This script will run a vortex panel computation method and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# Geometric surface points must be in the saved_airfoil_coords folder in order for successful execution
# Constant Strength Vortex Method is used here; Future implementations will include a linear-strength vortex method for comparison

import sys;
import numpy as np;
import matplotlib.pyplot as plt;
import scipy.integrate as scpi;
import pandas as pd;
import os;
import math;

def run_vpm_solver(input_file_name):
    print("Placeholder");

# Method that determines local induced velocity at collocation points (xi, zi) by vortex panels
def V2DC_induced_vel(loc_gam, x_i, z_i, x_j, z_j, x_jp1, z_jp1):
    # Note, Katz uses (x,z) for coordinates instead of (x,y), I use the same convention here.
    # Compute differences between collocation and jth panel node coordinates
    zdif1 = z_i - z_j;
    xdif1 = x_i - x_j;

    # Compute differences between collocation and j+1th panel node coordinates
    zdif2 = z_i - z_jp1
    xdif2 = x_i - x_jp1

    # Calculate geometric distance between collocation and each of jth and j+1th panels
    jdist = xdif1**2 + zdif1**2;
    jp1dist = xdif2**2 + zdif2**2;

    # Local induced velocity (Katz eqn. 11.44 & 11.45)
    u_i = (loc_gam/(2*np.pi))*(np.arctan2(zdif2/xdif2) - np.arctan2(zdif1/xdif1));
    w_i = (-loc_gam/(4*np.pi))*(np.log(jdist/jp1dist));
    
    # Returns local induced velocity in coordinate form
    return u_i, w_i;