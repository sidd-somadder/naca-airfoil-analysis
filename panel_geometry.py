""" 
This script is used to retrieve important panel geometric parameters (angles in radians):
tangential angle (phi), outwards normal angle (beta), panel lengths (S), panel midpoints

Computes K and L normal and tangential geometric influence matrices based on panel geometry.
Math derivation found in accompanying report.

Handles loading coordinates into a N x 2 matrix from saved_coordinates folder with specified name.
Does not format from provided .dat file: only converts raw struture in file into coordinate array.
"""

import numpy as np
import os

# Using coordinate points, get the tangential angle and length of each panel.
def get_geom_params(geom_points):
    N = len(geom_points)

    # For N panels, there are N angles.
    phi = np.zeros(N)
    p_lengths = np.zeros(N)
    midpoints = np.zeros((N,2))

    for k in range(N):
        # Get kth panel node coordinates.
        x_k = geom_points[k,0]
        z_k = geom_points[k,1]

        # Get k+1th panel node coordinates; loop back to first node if the kth node is the Nth node.
        if (k+1) == N:
            x_kp1 = geom_points[0,0]
            z_kp1 = geom_points[0,1]        
        else:
            x_kp1 = geom_points[k+1,0]
            z_kp1 = geom_points[k+1,1]
        
        # Calculate vertical difference dz and horizontal difference dx.
        dx = x_kp1 - x_k
        dz = z_kp1 - z_k

        # Compute tangential angle using arctan.
        phi[k] = np.arctan2(dz,dx)

        # Calculate length by distance formula.
        p_lengths[k] = np.sqrt((dz)**2 + (dx)**2)

        # Calculate midpoint coordinates by averaging adjacent node coordinates.
        x_m = 0.5 *(x_k + x_kp1)
        z_m = 0.5 *(z_k + z_kp1)
        midpoints[k,:] = [x_m, z_m]
    
    # Define outwards normal angle from tangential angle
    # The VPM Neumann boundary condition RHS is built from beta, not phi.
    beta = phi + np.pi/2

    return phi, beta, p_lengths, midpoints

# Make influence coefficient using collocation points, panel nodes, panel lengths, and tangential angles.
def compute_KL_inf_matrices(geom_pts, midpoints, S, phi):
    """
    Assemble the normal (K) and tangential (L) influence coefficient matrices.

    Implements the geometric integral expressions derived in Section II-B of
    the accompanying report. Variable names A, B, C, D, E correspond to the
    intermediate terms defined there.

    Returns K, L: each (N, N), where N is the number of panels.
    """
    N = len(phi)
    
    # Initialize square matrix with N equations and N local gammas
    # K is the normal influence matrix used for the boundary condition
    K = np.zeros((N,N))
    # L is the tangential influence matrix used to find the tangential velocity per panel
    L = np.zeros((N,N))

    for i in range(N):
        # Retrieve ith collocation point coordinates
        x_i = midpoints[i, 0]
        z_i = midpoints[i, 1]
        for j in range(N):

            # Retrieve jth collocation point coordinates
            x_j = geom_pts[j,0]
            z_j = geom_pts[j,1]

            # Intermediate geometric integration terms 
            A = -(x_i - x_j)*np.cos(phi[j]) - (z_i - z_j)*np.sin(phi[j])
            B = (x_i - x_j)**2 + (z_i - z_j)**2
            C_k = - np.cos(phi[i] - phi[j])
            C_l = np.sin(phi[j] - phi[i])
            D_k = (x_i - x_j)*np.cos(phi[i]) + (z_i - z_j)*np.sin(phi[i])
            D_l = (x_i - x_j)*np.sin(phi[i]) - (z_i - z_j)*np.cos(phi[i])
            E = 0
            if B > A**2:
                E = np.sqrt(B - A**2)
            
            # Retrieve jth panel length
            s_j = S[j]
            
            # Diagonal terms = 0
            if not (i == j):        
                # Compute the logarithm terms for each resulting matrix entry 
                log_term = np.log((s_j**2 + 2*A*s_j + B)/B)

                K[i,j] = (C_k/2) * log_term
                L[i,j] = (C_l/2) * log_term

                if E != 0:
                    # Add arctan term only if E is nonzero.
                    atan_term = np.arctan2((s_j + A), E) - np.arctan2(A, E)
                    K[i,j] += ((D_k - A*C_k)/E)*atan_term
                    L[i,j] += ((D_l - A*C_l)/E)*atan_term

            # Zero out any remaining problem values if any (complex, NaN, inf)
            if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):
                K[i,j] = 0
            if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):
                L[i,j] = 0
    return K, L

# Reads a raw .dat coordinate file from the saved_airfoil_coords folder and returns an Nx2 numpy array of (x, y) points
# This is a raw parse only, assume in reverse Selig format.
def load_dat_coordinates(filename):
    input_dir = os.path.join(os.path.dirname(__file__), "saved_airfoil_coords")
    filepath = os.path.join(input_dir, filename)

    raw_points = []

    with open(filepath, "r") as f:
        for line in f:
            tokens = line.split()
            # Skip blank lines or anything that isn't exactly an (x, y) pair
            if len(tokens) != 2:
                continue
            try:
                x_val = float(tokens[0])
                y_val = float(tokens[1])
            except ValueError:
                # Catches the airfoil-name header line most .dat files start with
                continue
            raw_points.append((x_val, y_val))

    coords = np.array(raw_points)

    return coords