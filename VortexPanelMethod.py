# This script will run a vortex panel computation method and is called automatically by the master script
# .csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
# Geometric surface points must be in the saved_airfoil_coords folder in order for successful execution
# Constant Strength Vortex Method is used here; Future implementations will include a linear-strength vortex method for comparison

import numpy as np;
import matplotlib.pyplot as plt;
import pandas as pd;
import os;

# Assume points are already processed into reverse Selig format
def run_vpm_solver(geom_points, alphas, input_file_name):
    # Retrieve tangential angles (phi), panel lengths, and panel midpoints (collocation points)
    phi, panel_lengths, midpoints = get_geom_params(geom_points);

    # Define outwards normal angle from tangential angle
    beta = phi + np.pi/2;

    # Convert angle of attack to radians
    alphas_rad = alphas * np.pi/180;

    N = len(alphas);
    M = len(phi);
    gamma_distribution = np.zeros((M,N));

    K,L = V2DC_influence_matrices(geom_points, midpoints, panel_lengths, phi);
    for k in range(N):
        RHS = V2DC_RHS_vec(alphas_rad[k], beta);

        # The kth column of gamma_distribution corresponds to the kth angle in alphas_rad
        gamma_distribution[:,k] = V2DC_solve_gamma_eqn(K, RHS);

    coeff_P_matrix = V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L);
    export_pressure_results(coeff_P_matrix, midpoints, alphas, input_file_name);

    c_l = get_coeffs(gamma_distribution, panel_lengths);
    print(c_l);
    plot_coeffs(c_l, alphas);

    plot_pressure(coeff_P_matrix, midpoints, alphas);

    #return gamma_distribution, coeff_P_matrix, midpoints, c_l;

# Make influence coefficient using collocation points, panel nodes, panel lengths, and tangential angles,
def V2DC_influence_matrices(geom_pts, midpoints, plengths, phi):
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
            s_j = plengths[j];

            # Compute the logarithm and arctan terms for each resulting matrix entry 
            log_term = np.log((s_j**2 + 2*A*s_j + B)/B)
            atan_term = np.arctan2((s_j + A), E) - np.arctan2(A, E)
            
            # Diagonal terms = 0
            if not (i == j):
                K[i,j] = (C_k/2) * log_term;
                K[i,j] += ((D_k - A*C_k)/E)*atan_term;
    
                L[i,j] = (C_l/2) * log_term
                L[i,j] += ((D_l - A*C_l)/E)*atan_term;

    # We replace the blunt TE panel influence equation with the Kutta Condition for the panels adjacent to TE points
    K[N-1,:] = 0;
    K[N-1,0] = 1;
    K[N-1,N-2] = 1;

    return K, L;

# Returns the RHS for the Kg = RHS matrix equation
# Takes the set of angle of attacks (alpha) and all local panel tangential angles (beta);
def V2DC_RHS_vec(alpha, beta):

    # If there are N outward normal angles, then there are N panels, and therefore the RHS has N elements
    RHS = np.zeros_like(beta);

    # (RHS)_i = 2pi * V_inf * cos(beta - alpha)
    RHS = np.cos(beta - alpha)

    # To include the Kutta Condition for the Nth equation, we set RHS = 0 
    RHS[-1] = 0.0;

    # Factor in 2pi to simplify matrix equation
    return 2 * np.pi * RHS;

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

        # Calculate length by distance formula.
        p_lengths[k] = np.sqrt((dz)**2 + (dx)**2);
    
        x_m = 0.5 *(x_k + x_kp1);
        z_m = 0.5 *(z_k + z_kp1);
        midpoints[k,:] = [x_m, z_m];
    
    return phi, p_lengths, midpoints;

# From other methods compute the coefficient matrix A and the RHS boundary conditions vector.
# Looped for each RHS vector for each alpha; main method saves local vortex strengths as matrix for each angle of attack
def V2DC_solve_gamma_eqn(K, RHS):
    #Solve the linear system of equations for all gamma values for each panel. 
    loc_gamma = np.linalg.solve(K, RHS);
    return loc_gamma;

def V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L_matrix): 

    # initialize coefficient of pressure matrix using solved gammas matrix
    M, N = gamma_distribution.shape;

    V_tang = np.zeros((M,N));

    # for each outwards normal angle in range
    for i in range(M):
        # for each angle of attack in range
        for k in range(N):
            sum_influence = 0;
            for j in range(M):
                # for the gammas corresponding to the kth angle of attack
                sum_influence += gamma_distribution[j,k]/(2*np.pi) * L_matrix[i,j];
            V_tang[i,k] = np.sin(beta[i] - alphas_rad[k]) - sum_influence;

    coeff_p = 1 - V_tang**2
    return coeff_p;

# Interactively plots Cp vs x/c for a user-selected angle of attack, looping until the user quits.
# Takes inputs coeff_P_matrix (M X N), from V2DC_pressure_distribution
# and angle-of-attack array in DEGREES (alphas)
def plot_pressure(coeff_P_matrix, midpoints, alphas):
    x_coords = midpoints[:, 0];
    chord = np.max(x_coords) - np.min(x_coords);
    norm_x = x_coords / chord;

    # Reverse Selig order is TE -> lower -> LE -> upper -> TE, one continuous loop.
    # The leading edge is the point of minimum x/c; split the array there to color
    # the two surfaces independently.
    le_idx = np.argmin(norm_x);

    print(f"Available angles of attack: {alphas}");

    while True:
        user_input = input("Enter an angle of attack to plot Cp for (or 'q' to quit): ").strip();
        if user_input.lower() == 'q':
            break;

        try:
            requested_aoa = float(user_input);
        except ValueError:
            print("Invalid input -- please enter a numeric angle or 'q'.");
            continue;

        k = np.argmin(np.abs(alphas - requested_aoa));
        matched_aoa = alphas[k];

        if not np.isclose(matched_aoa, requested_aoa, atol=1e-6):
            print(f"No exact match for {requested_aoa}°; showing closest available angle: {matched_aoa:.2f}°.");

        fig, ax = plt.subplots(figsize=(9, 6));

        # Lower surface: index 0 through le_idx (inclusive, so the two segments share
        # the LE point and the line doesn't visibly break there).
        ax.plot(norm_x[:le_idx+1], coeff_P_matrix[:le_idx+1, k],
                color='blue', linewidth=2, marker='o', markersize=3, label='Lower surface');

        # Upper surface: le_idx through the end.
        ax.plot(norm_x[le_idx:], coeff_P_matrix[le_idx:, k],
                color='red', linewidth=2, marker='o', markersize=3, label='Upper surface');

        ax.axhline(0, color='black', linewidth=0.8, linestyle='--');
        ax.invert_yaxis();
        ax.set_xlabel(r'$x/c$');
        ax.set_ylabel(r'$C_p$');
        ax.set_title(rf'VPM — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)');
        ax.legend();
        ax.grid(True, linestyle=':', alpha=0.6);
        plt.tight_layout();
        plt.show();

# Exports the full Cp matrix (all panels x all angles of attack) as a single .csv.
def export_pressure_results(coeff_P_matrix, midpoints, alphas, input_file_name):
    output_dir = os.path.join(os.path.dirname(__file__), "computation_results");
    os.makedirs(output_dir, exist_ok=True);
    
    # Process file name and output path
    identifier = input_file_name.replace(".dat", "");
    output_path = os.path.join(output_dir, f"VPM_Cp_{identifier}_results.csv");

    # normalize x-axis values from 0-1 to be exactly x/c
    x_coords = midpoints[:, 0];
    chord = np.max(x_coords) - np.min(x_coords);
    norm_x = x_coords / chord;

    # .csv formatting and column titling
    column_names = [f"Cp_alpha_{a:.2f}" for a in alphas];
    df = pd.DataFrame(coeff_P_matrix, columns=column_names);
    df.insert(0, "x_over_c", norm_x);

    df.to_csv(output_path, index=False, float_format="%.6f");
    print(f"VPM pressure distribution results exported: {output_path}");

def get_coeffs(gamma_distribution, p_lengths):
    M, N = gamma_distribution.shape;
    
    c_l = np.zeros(N);

    for alpha in range(N):
        for k in range(M):
            c_l[alpha] += 2*gamma_distribution[k, alpha] * p_lengths[k];

    return c_l;

def plot_coeffs(c_l, a):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(a, c_l,   label=r'$C_L$',       color='steelblue', linewidth=2)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)'); ax.set_ylabel('c_l')
    ax.set_title('VPM — Quick Check'); ax.legend(); ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout(); plt.show()