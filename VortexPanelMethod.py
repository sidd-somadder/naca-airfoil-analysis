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

# Assume points are already processed into reverse Selig format
def run_vpm_solver(geom_points):

    print("Placeholder");



# Method that determines local induced velocity at collocation points (xi, zi) by panel between the jth and j+1th nodes
def V2DC_induced_vel(loc_gam, x_i, z_i, x_j, z_j, x_jp1, z_jp1, beta_j):
    # Note, Katz uses (x,z) for coordinates instead of (x,y), I use the same convention here. Note: jp1 = j+1.
    # First, derive translated coordinates of the ith collocation point
    x_i_prime = x_i - x_j;
    z_i_prime = z_i - z_j;

    # Using tangential angle beta of the j to j+1th panel, transform prime coordinates to panel coordinates (bar)
    x_i_bar = x_i_prime * np.cos(beta_j) + z_i_prime * np.sin(beta_j);
    z_i_bar = -x_i_prime * np.sin(beta_j) + z_i_prime * np.cos(beta_j);

    # Calculate panel length from panel node coordinates (x_j, z_j) and (x_jp1, z_jp1)
    p_L = np.sqrt((z_jp1 - z_j)**2 + (x_jp1 - x_j)**2);

    # Local induced velocity (Katz eqn. 11.44 & 11.45)
    ubar = loc_gam/(2*np.pi) * (np.arctan2(z_i_bar,(x_i_bar-p_L)) - np.arctan2(z_i_bar,x_i_bar)) 
    wbar = -loc_gam/(4*np.pi) * np.log(((x_i_bar)**2 + (z_i_bar)**2)/((x_i_bar - p_L)**2 + (z_i_bar)**2))
    
    # Use tangential angle again to rotate velocity vector to cartesian coordinates from panel coordinates
    u_i = ubar * np.cos(beta_j) - wbar * np.sin(beta_j);
    w_i = ubar * np.sin(beta_j) + wbar * np.cos(beta_j);
    
    # Returns local induced velocity in coordinate form
    return u_i, w_i;

# 
def V2DC_influence_matrix(geom_pts, midpoints, beta):
    N = len(geom_pts);

    A = np.zeros((N,N));

    for i in range(N-1):
        x_i = midpoints[i,0];
        z_i = midpoints[i,1];
        for j in range(N):
            # Define jth x,z coordinates required to calculate (i,j)th influence coefficient
            x_j = geom_pts[j,0]
            z_j = geom_pts[j,1]

            # If j = N-1, loop back to the first index;
            if j+1 == N:
               x_jp1 = geom_pts[0,0];
               z_jp1 = geom_pts[0,1];
            # Otherwise, continue down the array of points
            else:
               x_jp1 = geom_pts[j+1,0];
               z_jp1 = geom_pts[j+1,1];
            
            if i == j:
               # Using Katz eqn. 11.50
               A[i,j] = -0.5;
            else:
               # Calculate singularity vortex element influenced velocities
               u, w = V2DC_induced_vel(1, x_i, z_i, x_j, z_j, x_jp1, z_jp1, beta[j]);
               # Using Katz eqn. 11.49
               A[i,j] = u * np.cos(beta[i]) - w * np.sin(beta[i]);  
    
    # We replace the blunt TE panel influence equation with the Kutta Condiiton for the panels adjacent to TE points
    A[N-1,0] = 1;
    A[N-1,N-2] = 1;

    return A;

# Returns the RHS for the Ag = RHS matrix equation
# Takes the set of angle of attacks (alpha) and all local panel tangential angles (beta);
def V2DC_RHS_vec(alpha, beta):
    # Define Cartesian components of freestream velocity.
    # Since this solver returns coefficients instead of raw lift values, assume V = 1 
    U = np.cos(alpha);
    W = np.sin(alpha);

    # If there are N tangential angles, then there are N panels, and therefore the RHS has N elements
    RHS = np.zeros_like(beta);

    # (RHS)_i = - (U,W) dot (cos(beta[i] - sin(beta[i]))); Katz eqn. 11.51
    RHS = - U * np.cos(beta) + W * np.sin(beta);

    # To include the Kutta Condition for the Nth equation, we set RHS = 0 
    RHS[-1] = 0.0;

    return RHS;

def get_midpoints(geom_pts):

    # number of midpoints = number of panel nodes for a closed system of panels
    N = len(geom_pts);
    midpoints = np.zeros((N,2));
    
    # Use the face that the midpoint between two values is their sum divided by 2.
    for i in range(N):
        if i < (N-1):
            x_m = 0.5 *(geom_pts[i,0] + geom_pts[i+1,0]);
            z_m = 0.5 *(geom_pts[i,1] + geom_pts[i+1,1]);
            midpoints[i,:] = [x_m, z_m];
        # if i is on the N-1th index, we consider the 0 as the i+1 index to represent the first panel node
        else:
            x_m = 0.5 *(geom_pts[i,0] + geom_pts[0,0]);
            z_m = 0.5 *(geom_pts[i,1] + geom_pts[0,1]);
            midpoints[i,:] = [x_m, z_m];

    return midpoints;

# Using coordinate points, get the tangential angle and length of each panel.
def get_geom_params(geom_points):
    N = len(geom_points);

    # For N panels, there are N angles
    beta = np.zeros(N);
    p_lengths = np.zeros(N);

    for k in range(N):
        # Get kth panel node coordinates.
        x_k = geom_points[k,0];
        z_k = geom_points[k,1];

        # Get k+1th panel node coordinates; loop back to first node if the kth node is the Nth node.
        if (k+1) == N:
            x_kp1 = geom_points[0,0];
            z_kp1 = geom_points[0,1];        
        else:
            x_kp1 = geom_points[k+1,0];
            z_kp1 = geom_points[k+1,1];
        
        # Calculate vertical difference dz and horizontal difference dx.
        dx = x_kp1 - x_k;
        dz = z_kp1 - z_k;

        # Compute tangential angle using arctan.
        beta[k] = np.arctan2(dz,dx);

        # Calculate length by distance formula.
        p_lengths[k] = np.sqrt((dz)**2 + (dx)**2);
    
    return beta, p_lengths;



