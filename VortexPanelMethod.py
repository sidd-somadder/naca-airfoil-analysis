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
    beta, panel_lengths, midpoints = get_geom_params(geom_points);
    alphas_rad = alphas * np.pi/180;

    N = len(alphas);
    M = len(beta);
    gamma_distribution = np.zeros((M,N));

    A = V2DC_influence_matrix(geom_points, midpoints, beta, panel_lengths);
    for k in range(N):
        RHS = V2DC_RHS_vec(alphas_rad[k], beta);
        gamma_distribution[:,k] = V2DC_solve_gamma_eqn(A, RHS);

    coeff_P_matrix = V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution);
    export_pressure_results(coeff_P_matrix, midpoints, alphas, input_file_name);

    print(gamma_distribution);
    print(coeff_P_matrix);

    #return gamma_distribution, coeff_P_matrix, midpoints;

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
    return loc_gamma;

# Calculates coefficient of pressure for every panel, across every angle of attack, in one call.
# alphas_rad: array of angles of attack, in radians (matches the rest of this script's unit convention)
# beta: array of panel tangential angles (M,)
# gamma_distribution: solved vortex strengths, shape (M, N) -- M panels, N angles of attack
# Returns coeff_P: shape (M, N), same row/column convention as gamma_distribution
def V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution):
    # Broadcast beta (M,) against alphas_rad (N,) into an (M,N) matrix where
    # entry [j,k] = beta_j + alpha_k -- the (alpha + alpha_j) term in Katz eqn. 11.53.
    angle_matrix = beta[:, np.newaxis] + alphas_rad[np.newaxis, :];

    coeff_P = 1 - (np.cos(angle_matrix) + 0.5 * gamma_distribution)**2;

    return coeff_P;

# Interactively plots Cp vs x/c for a user-selected angle of attack, looping until the user quits.
# Takes inputs coeff_P_matrix (M X N), from V2DC_pressure_distribution
# and angle-of-attack array in DEGREES (alphas)
def plot_pressure(coeff_P_matrix, midpoints, alphas):
    x_coords = midpoints[:, 0];
    chord = np.max(x_coords) - np.min(x_coords);
    norm_x = x_coords / chord;

    print(f"Available angles of attack: {list(np.round(alphas, 2))}");

    # Keep asking user for angle of attack for pressure coefficient plot, loop until user quets.
    while True:
        user_input = input("Enter an angle of attack to plot Cp for (or 'q' to quit): ").strip();
        if user_input.lower() == 'q':
            break;

        # Check if user input is valid.
        try:
            requested_aoa = float(user_input);
        except ValueError:
            print("Invalid input -- please enter a numeric angle or 'q'.");
            continue;

        # Since user input likely doesn't match stored float, find nearest within small tolerance.
        k = np.argmin(np.abs(alphas - requested_aoa));
        matched_aoa = alphas[k];

        # If no angle is found, plot for the next nearest angle.
        if not np.isclose(matched_aoa, requested_aoa, atol=1e-6):
            print(f"No exact match for {requested_aoa}°; showing closest available angle: {matched_aoa:.2f}°.");

        fig, ax = plt.subplots(figsize=(9, 6));
        ax.plot(norm_x, coeff_P_matrix[:, k], color='steelblue', linewidth=2, marker='o', markersize=3);
        ax.axhline(0, color='black', linewidth=0.8, linestyle='--');
        ax.invert_yaxis();
        ax.set_xlabel(r'$x/c$');
        ax.set_ylabel(r'$C_p$');
        ax.set_title(rf'VPM — Pressure Distribution ($\alpha$ = {matched_aoa:.1f}°)');
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



    