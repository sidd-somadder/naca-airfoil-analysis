""" 
This script is used to retrieve important panel geometric parameters (angles in radians):
tangential angle (phi), outwards normal angle (beta), panel lengths (S), panel midpoints

Computes K and L normal and tangential geometric influence matrices based on panel geometry.
Math derivation found in accompanying report.

Handles loading coordinates into a M x 2 matrix from saved_coordinates folder with specified name.
Does not format from provided .dat file: only converts raw structure in file into coordinate array.
"""

import numpy as np
import os

def get_geom_params(geom_points):
    """
    Computes panel geometric parameters from a set of airfoil surface coordinates.

    Inputs:
    Set of coordinates retrieved from .dat file (geom_points, size M x 2) in reverse-selig panel order.
    Panels are formed between consecutive nodes, with the final panel wrapping from the last node back to the first.

    Outputs:
    Tangential angle (phi) of each panel measured from the x-axis, in RADIANS (size M)
    Outwards normal angle (beta) of each panel, defined as phi + pi/2, in RADIANS (size M)
    Panel lengths (p_lengths) computed by distance formula between adjacent nodes (size M)
    Panel midpoint coordinates (midpoints) used as collocation points for the solver (size M x 2)
    """

    M = len(geom_points)

    # For M panels, there are M angles.
    phi = np.zeros(M)
    p_lengths = np.zeros(M)
    midpoints = np.zeros((M,2))

    for k in range(M):
        # Get kth panel node coordinates.
        x_k = geom_points[k,0]
        z_k = geom_points[k,1]

        # Get k+1th panel node coordinates; loop back to first node if the kth node is the Mth node.
        if (k+1) == M:
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

def compute_KL_inf_matrices(geom_pts, midpoints, S, phi):
    """
    Assemble the normal (K) and tangential (L) influence coefficient matrices.

    Implements the geometric integral expressions derived in Section II-B of
    the accompanying report. Variable names A, B, C, D, E correspond to the
    intermediate terms defined there.

    Returns K, L: each (M, M), where M is the number of panels.
    """
    M = len(phi)
    
    # Initialize square matrix with M equations and M local gammas
    # K is the normal influence matrix used for the boundary condition
    K = np.zeros((M,M))
    # L is the tangential influence matrix used to find the tangential velocity per panel
    L = np.zeros((M,M))

    for i in range(M):
        # Retrieve ith collocation point coordinates
        x_i = midpoints[i, 0]
        z_i = midpoints[i, 1]
        for j in range(M):

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

def load_dat_coordinates(filename):
    """
    Reads a .dat coordinate file from the saved_airfoil_coords folder into a coordinate array.

    Inputs:
    Filename (filename) as found in the saved_airfoil_coords folder, including the .dat extension.

    Outputs:
    Coordinate array (coords) of surface points in the order written in the file (size M x 2).
    Note: This is a raw parse only. Header lines and any row that is not exactly two floats are skipped.
    Coordinates are assumed to already be in reverse-selig order; no reordering is performed.
    """

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