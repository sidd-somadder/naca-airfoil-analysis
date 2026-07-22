# This script will run a vortex panel computation method and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# Geometric surface points must be in the saved_airfoil_coords folder in order for successful execution
# Constant Strength Vortex Method is used here; Future implementations will include a linear-strength vortex method for comparison

import numpy as np;
from panel_geometry import get_geom_params, compute_KL_inf_matrices;
from post_processing import plot_coeffs, plot_pressure, export_VPM_pressure;

# Assume points are already processed into reverse Selig format
def run_cvpm_solver(geom_points, alphas, input_file_name):

    # Retrieve tangential angles (phi), panel lengths, and panel midpoints (collocation points)
    phi, beta, panel_lengths, midpoints = get_geom_params(geom_points);

    # Convert angle of attack to radians
    alphas_rad = alphas * np.pi/180;

    N = len(alphas);
    M = len(phi);
    gamma_distribution = np.zeros((M,N));

    K,L = compute_KL_inf_matrices(geom_points, midpoints, panel_lengths, phi);

    K_solve = K.copy();
    K_solve[-1,:] = 0;
    K_solve[-1,0] = 1;
    K_solve[-1,M-2] = 1;

    for k in range(N):
        RHS = V2DC_RHS_vec(alphas_rad[k], beta);

        # The kth column of gamma_distribution corresponds to the kth angle in alphas_rad
        gamma_distribution[:,k] = V2DC_solve_gamma_eqn(K_solve, RHS);

    coeff_P_matrix = V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L);
    export_VPM_pressure(coeff_P_matrix, midpoints, alphas, input_file_name, method='CVPM');

    invP = invalid_panel_indices(M);
    c_l_KJ, c_l = VPM_get_coeffs(gamma_distribution, panel_lengths, coeff_P_matrix, beta, alphas_rad, invP);

    plot_coeffs(alphas, c_l_KJ=c_l_KJ, c_l_P=c_l, title='Const. VPM Lift Coefficient Comparison (KJ vs. Pressure)');
    plot_pressure(coeff_P_matrix, midpoints, alphas);

# Returns the RHS for the Kg = RHS matrix equation
# Takes specified angle of attack (alpha) and all local panel tangential angles (beta);
def V2DC_RHS_vec(alpha, beta):
    # If there are N outward normal angles, then there are N panels, and therefore the RHS has N elements
    RHS = np.zeros_like(beta);

    # (RHS)_i = 2pi * V_inf * cos(beta - alpha)
    RHS = 2 * np.pi * np.cos(beta - alpha)

    # If open, last row is overwritten by the Kutta condition; last row RHS is rewritten to be 0
    RHS[-1] = 0.0;
    
    return RHS;                          # length M

# From other methods compute the coefficient matrix A and the RHS boundary conditions vector.
# Looped for each RHS vector for each alpha; main method saves local vortex strengths as matrix for each angle of attack
def V2DC_solve_gamma_eqn(K_solve, RHS):
    # Solve the linear system of equations for all gamma values for each panel. 
    return np.linalg.solve(K_solve, RHS);

def V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L_matrix): 
    # initialize coefficient of pressure matrix using solved gammas matrix  
    M, N = gamma_distribution.shape;

    V_tang = np.zeros((M,N));

    # for each outwards normal angle in range
    for i in range(M-1):
        # for each angle of attack in range
        for k in range(N):
            sum_influence = 0;
            for j in range(M):
                sum_influence -= gamma_distribution[j, k]/(2*np.pi) * L_matrix[i, j]
            
            V_tang[i, k] = np.sin(beta[i] - alphas_rad[k]) + sum_influence + gamma_distribution[i,k]/2

    coeff_p = 1 - V_tang**2
    
    # Mask rather than delete, so every downstream array stays length M.
    coeff_p[invalid_panel_indices(M), :] = np.nan;
    return coeff_p;

def invalid_panel_indices(M):
    return [0, M-2, M-1];

# This method is called by the linear VPM script as well
# Since local vortex strengths, geometric parameters, and coefficient of pressures are solved for independently
def VPM_get_coeffs(gamma_distribution, p_lengths, coeff_p, beta, alphas_rad, invalid_panels):
    M, N = gamma_distribution.shape;
    
    # Derive lift coefficient in two ways; 
    # compute by summing circulation and using Kutta-Joukowski 
    c_l_KJ = np.zeros(N);
    # and integrate via pressure distribution
    c_l = np.zeros(N);

    # For high TE error, use masking on high outlier values on extremely small panels to make math smoother
    cp_mask = np.ones(M, dtype=bool);
    cp_mask[invalid_panels] = False;

    for k in range(N):
        # Circulation-based lift: every panel except the spurious closing one carries real circulation.
        c_l_KJ[k] = 2 * np.sum(gamma_distribution[:-1, k] * p_lengths[:-1]);

        # Pressure-based lift: only panels with a physically meaningful Cp.
        cp = coeff_p[cp_mask, k];
        s  = p_lengths[cp_mask];
        b  = beta[cp_mask];

        # Compute axial and normal force coefficients 
        c_n = -(cp * s * np.sin(b)).sum();
        c_a = -(cp * s * np.cos(b)).sum();

        c_l[k] = c_n * np.cos(alphas_rad[k]) - c_a * np.sin(alphas_rad[k]);

    return c_l_KJ, c_l;
