# This script is used to retrieve important panel geometric parameters:
# tangential angle (phi)
# outwards normal angle (beta)
# panel lengths (S)
# panel midpoints

import numpy as np;

# Using coordinate points, get the tangential angle and length of each panel.
def get_geom_params(geom_points):
    N = len(geom_points);

    # For N panels, there are N angles
    phi = np.zeros(N);
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
        phi[k] = np.arctan2(dz,dx);

        # Define outwards normal angle from tangential angle
        beta = phi + np.pi/2;

        # Calculate length by distance formula.
        p_lengths[k] = np.sqrt((dz)**2 + (dx)**2);
    
        x_m = 0.5 *(x_k + x_kp1);
        z_m = 0.5 *(z_k + z_kp1);
        midpoints[k,:] = [x_m, z_m];
    
    return phi, beta, p_lengths, midpoints;

# Make influence coefficient using collocation points, panel nodes, panel lengths, and tangential angles,
def compute_KL_inf_matrices(geom_pts, midpoints, S, phi):
    N = len(geom_pts);
    
    # Initialize square matrix with N equations and N local gammas
    # K is the normal influence matrix used for the boundary condition
    K = np.zeros((N,N));
    # L is the tangential influence matrix used to find the tangential velocity per panel
    L = np.zeros((N,N));

    for i in range(N):
        for j in range(N):
            # Retrieve ith collocation point coordinates
            x_i = midpoints[i, 0];
            z_i = midpoints[i, 1];

            # Retrieve jth collocation point coordinates
            x_j = geom_pts[j,0];
            z_j = geom_pts[j,1];

            # Intermediate geometric integration terms 
            A = -(x_i - x_j)*np.cos(phi[j]) - (z_i - z_j)*np.sin(phi[j]);
            B = (x_i - x_j)**2 + (z_i - z_j)**2;
            C_k = - np.cos(phi[i] - phi[j]);
            C_l = np.sin(phi[j] - phi[i]);
            D_k = (x_i - x_j)*np.cos(phi[i]) + (z_i - z_j)*np.sin(phi[i]);
            D_l = (x_i - x_j)*np.sin(phi[i]) - (z_i - z_j)*np.cos(phi[i]);
            E = 0;
            if B > A**2:
                E = np.sqrt(B - A**2);
            
            # Retrieve jth panel length
            s_j = S[j];

            # Compute the logarithm and arctan terms for each resulting matrix entry 
            log_term = np.log((s_j**2 + 2*A*s_j + B)/B)
            atan_term = np.arctan2((s_j + A), E) - np.arctan2(A, E)
            
            # Diagonal terms = 0
            if not (i == j):
                K[i,j] = (C_k/2) * log_term;
                L[i,j] = (C_l/2) * log_term;

                K[i,j] += ((D_k - A*C_k)/E)*atan_term;
                L[i,j] += ((D_l - A*C_l)/E)*atan_term;

            # Zero out any remaining problem values
            if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):      # If K term is complex or a NAN or an INF
                K[i,j] = 0                                                          # Set K value equal to zero
            if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):      # If L term is complex or a NAN or an INF
                L[i,j] = 0                                                          # Set L value equal to zero
    return K, L;