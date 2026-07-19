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

    c_l_KJ, c_l = get_coeffs(gamma_distribution, panel_lengths, coeff_P_matrix, beta, alphas_rad);
    print(c_l_KJ);
    print(c_l);
    plot_coeffs(c_l_KJ, c_l, alphas);

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
                L[i,j] = (C_l/2) * log_term
                if E != 0:
                    K[i,j] += ((D_k - A*C_k)/E)*atan_term;
                    L[i,j] += ((D_l - A*C_l)/E)*atan_term;

            # Zero out any remaining problem values
            if (np.iscomplex(K[i,j]) or np.isnan(K[i,j]) or np.isinf(K[i,j])):      # If K term is complex or a NAN or an INF
                K[i,j] = 0                                                          # Set K value equal to zero
            if (np.iscomplex(L[i,j]) or np.isnan(L[i,j]) or np.isinf(L[i,j])):      # If L term is complex or a NAN or an INF
                L[i,j] = 0                                                          # Set L value equal to zero
    
    # We replace the blunt TE panel influence equation with the Kutta Condition for the panels adjacent to TE points
    K[-1,:] = 0;
    K[-1,0] = 1;
    K[-1,N-2] = 1;

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

def invalid_panel_indices(M):
    return [0, M-2, M-1];

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

def plot_pressure(coeff_P_matrix, midpoints, alphas):
    x_coords = midpoints[:, 0];

    # Reverse Selig order is TE -> lower -> LE -> upper -> TE.
    # Split at the leading edge (minimum x/c) to colour the surfaces separately.
    le_idx = np.argmin(x_coords);

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

        cp = coeff_P_matrix[:, k];

        fig, ax = plt.subplots(figsize=(9, 6));

        ax.plot(x_coords[:le_idx+1], cp[:le_idx+1],
                color='blue', linewidth=2, marker='o', markersize=3, label='Lower surface');
        ax.plot(x_coords[le_idx:], cp[le_idx:],
                color='red', linewidth=2, marker='o', markersize=3, label='Upper surface');

        # Set limits from the finite data only, so a stray value can never rescale the axes.
        lo, hi = np.nanmin(cp), np.nanmax(cp);
        pad = 0.1 * (hi - lo);
        ax.set_ylim(hi + pad, lo - pad);          # already inverted, so no invert_yaxis()

        ax.axhline(0, color='black', linewidth=0.8, linestyle='--');
        ax.set_xlabel(r'$x/c$');
        ax.set_ylabel(r'$C_p$');
        ax.set_title(rf'VPM — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)');
        ax.legend();
        ax.grid(True, linestyle=':', alpha=0.6);
        plt.tight_layout();
        plt.show();

def export_pressure_results(coeff_P_matrix, midpoints, alphas, input_file_name):
    output_dir = os.path.join(os.path.dirname(__file__), "computation_results");
    os.makedirs(output_dir, exist_ok=True);

    identifier  = input_file_name.replace(".dat", "");
    output_path = os.path.join(output_dir, f"VPM_Cp_{identifier}_results.csv");

    x_coords = midpoints[:, 0];
    chord    = np.max(x_coords) - np.min(x_coords);
    norm_x   = x_coords / chord;

    column_names = [f"Cp_alpha_{a:.2f}" for a in alphas];
    df = pd.DataFrame(coeff_P_matrix, columns=column_names);
    df.insert(0, "x_over_c", norm_x);

    # Masked panels export as empty cells; drop them so the CSV has no gaps.
    df = df.dropna();

    df.to_csv(output_path, index=False, float_format="%.6f");
    print(f"VPM pressure distribution results exported: {output_path}");

def get_coeffs(gamma_distribution, p_lengths, coeff_p, beta, alphas_rad):
    M, N = gamma_distribution.shape;
    
    c_l_KJ = np.zeros(N);
    c_l = np.zeros(N);

    cp_mask = np.ones(M, dtype=bool);
    cp_mask[[0, M-2, M-1]] = False;

    for k in range(N):
        # Circulation-based lift: every panel except the spurious closing one
        # carries real circulation, even the two Kutta-adjacent ones.
        c_l_KJ[k] = 2 * np.sum(gamma_distribution[:-1, k] * p_lengths[:-1]);

        # Pressure-based lift: only panels with a physically meaningful Cp.
        cp = coeff_p[cp_mask, k];
        s  = p_lengths[cp_mask];
        b  = beta[cp_mask];

        c_n = -(cp * s * np.sin(b)).sum();
        c_a = -(cp * s * np.cos(b)).sum();

        c_l[k] = c_n * np.cos(alphas_rad[k]) - c_a * np.sin(alphas_rad[k]);

    return c_l_KJ, c_l;

def plot_coeffs(c_l_KJ, c_l, a):
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(a, c_l_KJ,   label=r'$C_L (KJ)$',       color='steelblue', linewidth=2)
    ax.plot(a, c_l,   label=r'$C_L (P)$',       color='firebrick', linewidth=2)
    ax.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax.axvline(0, color='black', linewidth=0.8, linestyle='--')
    ax.set_xlabel(r'$\alpha$ (°)'); ax.set_ylabel('c_l')
    ax.set_title('Constant Strength VPM: Kutta Joukowski vs. Pressure Coefficient derived C_l'); ax.legend(); ax.grid(True, linestyle=':', alpha=0.6)
    plt.tight_layout(); plt.show()