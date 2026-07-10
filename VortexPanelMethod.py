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
def run_vpm_solver(geom_points, alphas):
    beta, panel_lengths, midpoints = get_geom_params(geom_points);

    # Each column (N) is a gamma distribution for a specific angle of attack
    # Each row (M) is a set of gamma values for each geometric tangential panel angle 
    N = len(alphas);
    M = len(beta);
    gamma_distribution = np.zeros((M,N));

    A = V2DC_influence_matrix(geom_points, midpoints, beta, panel_lengths);
    for k in range(N):
        # Define RHS vector for each alpha; has M elements
        RHS = V2DC_RHS_vec(alphas[k], beta);

        # resulting local vortex strength vector has M elements for the kth alpha
        gamma_distribution[:,k] = V2DC_solve_gamma_eqn(A, RHS);

    return gamma_distribution; # placeholder

# Method that determines local induced velocity at collocation points (xi, zi) by panel between the jth and j+1th nodes
def V2DC_induced_vel(loc_gam, x_i, z_i, x_j, z_j, beta_j, pL_j):
    # Note, Katz uses (x,z) for coordinates instead of (x,y), I use the same convention here. Note: jp1 = j+1.
    # First, derive translated coordinates of the ith collocation point
    x_i_prime = x_i - x_j;
    z_i_prime = z_i - z_j;

    # Using tangential angle beta of the j to j+1th panel, transform prime coordinates to panel coordinates (bar)
    x_i_bar = x_i_prime * np.cos(beta_j) + z_i_prime * np.sin(beta_j);
    z_i_bar = -x_i_prime * np.sin(beta_j) + z_i_prime * np.cos(beta_j);

    # Local induced velocity (Katz eqn. 11.44 & 11.45)
    # Note: x_jp1_bar = panel length (pL_j) 
    ubar = loc_gam/(2*np.pi) * (np.arctan2(z_i_bar,(x_i_bar-pL_j)) - np.arctan2(z_i_bar,x_i_bar)) 
    wbar = -loc_gam/(4*np.pi) * np.log(((x_i_bar)**2 + (z_i_bar)**2)/((x_i_bar - pL_j)**2 + (z_i_bar)**2))
    
    # Use tangential angle again to rotate velocity vector to cartesian coordinates from panel coordinates
    u_i = ubar * np.cos(beta_j) - wbar * np.sin(beta_j);
    w_i = ubar * np.sin(beta_j) + wbar * np.cos(beta_j);
    
    # Returns local induced velocity in coordinate form
    return u_i, w_i;

# Make influence coefficient using collocation points, panel nodes, panel lengths, and tangential angles,
def V2DC_influence_matrix(geom_pts, midpoints, beta, plengths):
    N = len(geom_pts);
    
    # Initialize square matrix with N equations and N local gammas
    A = np.zeros((N,N));

    for i in range(N-1):
        x_i = midpoints[i,0];
        z_i = midpoints[i,1];

        for j in range(N):
            # Define jth x,z coordinates required to calculate (i,j)th influence coefficient
            x_j = geom_pts[j,0]
            z_j = geom_pts[j,1]

            if i == j:
               # Using Katz eqn. 11.50
               A[i,j] = -0.5;
            else:
               # Calculate singularity vortex element influenced velocities
               u, w = V2DC_induced_vel(1, x_i, z_i, x_j, z_j, beta[j], plengths[j]);
               # Using Katz eqn. 11.49
               A[i,j] = u * np.cos(beta[i]) - w * np.sin(beta[i]);  
    
    # We replace the blunt TE panel influence equation with the Kutta Condition for the panels adjacent to TE points
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

# Using coordinate points, get the tangential angle and length of each panel.
def get_geom_params(geom_points):
    N = len(geom_points);

    # For N panels, there are N angles
    beta = np.zeros(N);
    p_lengths = np.zeros(N);
    midpoints = np.zeros((N,2));

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
    
        x_m = 0.5 *(x_k + x_kp1);
        z_m = 0.5 *(z_k + z_kp1);
        midpoints[k,:] = [x_m, z_m];
    
    return beta, p_lengths, midpoints;

# From other methods compute the coefficient matrix A and the RHS boundary conditions vector.
# Looped for each RHS vector for each alpha; main method saves local vortex strengths as matrix for each angle of attack
def V2DC_solve_gamma_eqn(A, RHS):
    #Solve the linear system of equations for all gamma values for each panel. 
    loc_gamma = np.linalg.solve(A, RHS);
    print(loc_gamma);
    return loc_gamma;