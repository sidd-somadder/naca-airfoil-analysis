"""
This script will run a vortex panel computation method and is called automatically by the master script
.csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
Geometric surface points must be in the saved_airfoil_coords folder in order for successful execution
Constant Strength Vortex Method is used here; Future implementations will include a linear-strength vortex method for comparison
"""

import numpy as np
from panel_geometry import get_geom_params, compute_KL_inf_matrices
from post_processing import export_VPM_pressure

def run_cvpm_solver(geom_points, alphas, input_file_name):
    """
    Run the constant-strength vortex panel method solver on a given airfoil over angle of attack sweep. 

    Inputs: 
    Set of coordinates retrieved from .dat file (geom_points, size M x 2), Angle of Attack (alphas, size N) in DEGREES, 
    and .dat filename as found in the saved_airfoil_coords folder.
    Filename used specifically for labelling and titling exported files & plots. 

    Outputs:
    Coefficient of Lift derived from Kutta Joukowski theorem (cL_KJ) and from integrating Coefficient of Pressure (cL_P), (size N)
    Coefficient of Moment about Leading Edge (cmLE) and about Quarter Chord (cmqc) (size N)
    Coefficient of Pressure matrix (coeff_P_matrix) with values across each panel across each alpha (M x N)
    Condition number (condition_num) from running SVD on the K influence matrix with the Kutta Condition
    Midpoints of each panel (midpoints)
    """

    # Retrieve tangential angles (phi), panel lengths, and panel midpoints (collocation points)
    phi, beta, panel_lengths, midpoints = get_geom_params(geom_points)

    # Convert angle of attack to radians
    alphas_rad = alphas * np.pi/180

    # Save lengths of angle of attack array and tangential angle array
    N = len(alphas)
    M = len(phi)

    # Initialize local vortex strength matrix, organized with M rows and N columns.
    gamma_distribution = np.zeros((M,N))

    # Call imported method from panel_geometry.py for K and L matrices
    K,L = compute_KL_inf_matrices(geom_points, midpoints, panel_lengths, phi)

    # Modify the K infleunce matrix by replacing the last (open TE) panel influence equation with the Kutta condition. 
    # Note that the panel of index M-1 is the spurious open TE panel so we use the adjacent panel with index M-2. 
    K_solve = K.copy()
    K_solve[-1,:] = 0
    K_solve[-1,0] = 1
    K_solve[-1,M-2] = 1

    # Calculate the condition number of the K influence matrix after including the Kutta condition
    condition_num = np.linalg.cond(K_solve)

    for k in range(N):
        # For each fixed angle of attack, fill the RHS vector where each entry corresponds to each panel.
        # This loop is ran for all angles of attack, resulting in k sets of local vortex strengths.
        RHS = V2DC_RHS_vec(alphas_rad[k], beta)

        # The kth column of gamma_distribution corresponds to the kth angle of attack in alphas_rad
        gamma_distribution[:,k] = np.linalg.solve(K_solve, RHS)

    # Calculate & export pressure distribution
    coeff_P_matrix = V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L)
    export_VPM_pressure(coeff_P_matrix, midpoints, alphas, input_file_name, method='CVPM')

    invP = invalid_panel_indices(M)

    # Use pressure and local gamma strength distribution to derive key coefficients.
    cL_KJ, cL_P, cmLE, cmqc = VPM_get_coeffs(gamma_distribution, panel_lengths, coeff_P_matrix, beta, alphas_rad, midpoints, invP)

    # Return coefficients of lift (Kutta-Joukowski & Pressure derived), moments, pressure, and condition number
    return cL_KJ, cL_P, cmLE, cmqc, coeff_P_matrix, condition_num, midpoints

def V2DC_RHS_vec(alpha, beta):
    """
    Returns the RHS for the Kg = RHS matrix equation used to calculate local vortex strength. Helper method for main solver.

    Inputs:
    Specified Angle of Attack in RADIANS (alpha)
    Outwards normal angle (beta) for every panel in reverse-selig panel order (size M)

    Outputs:
    RHS vector for the Kg = RHS matrix equation with the last entry replaced for the Kutta condition (size M). 
    """

    # If there are M outward normal angles, then there are M panels, and therefore the RHS has M elements
    RHS = np.zeros_like(beta)

    # (RHS)_i = 2pi * V_inf * cos(beta - alpha)
    RHS = 2 * np.pi * np.cos(beta - alpha)

    # Last row is overwritten by the Kutta condition; last row RHS is rewritten to be 0
    RHS[-1] = 0.0

    # length M
    return RHS                         

def V2DC_pressure_distribution(alphas_rad, beta, gamma_distribution, L_matrix):
    """
    Returns coefficient of pressure distribution across each angle of attack and each panel in reverse-selig panel order.

    Inputs:
    Angle of Attack sweep (alphas_rad) in RADIANS (size N).
    Outward normal angle (beta) for every panel in reverse-selig panel order (size M)
    Local vortex strength distribution (gamma distribution) across each panel across each angle of attack (size M x N)
    L tangential influence matrix (L_matrix) calculated from compute_KL_inf_matrices from panel_geometry.py (size M x M)

    Output:
    Coefficient of Pressure distribution (coeff_p) across each angle of attack across each panel.
    Note: Trailing Edge panel and adjacent panels are masked due to spurious values.
    """
     
    # initialize coefficient of pressure matrix using solved gammas matrix  
    M, N = gamma_distribution.shape

    V_tang = np.zeros((M,N))

    # for each outwards normal angle in range
    for i in range(M-1):
        # for each angle of attack in range
        for k in range(N):
            sum_influence = 0
            for j in range(M):
                sum_influence -= gamma_distribution[j, k]/(2*np.pi) * L_matrix[i, j]
            
            V_tang[i, k] = np.sin(beta[i] - alphas_rad[k]) + sum_influence + gamma_distribution[i,k]/2

    coeff_p = 1 - V_tang**2
    
    # Mask rather than delete, so every downstream array stays length M.
    coeff_p[invalid_panel_indices(M), :] = np.nan
    return coeff_p

def invalid_panel_indices(M):
    return [0, M-2, M-1]

def VPM_get_coeffs(gamma_distribution, p_lengths, coeff_p, beta, alphas_rad, midpoints, invalid_panels=None):
    """
    Helper method. Uses local vortex strength, pressure, and geometric parameters to return lift and moment coefficients.

    Inputs:
    Local vortex strength distribution (gamma distribution) across each panel across each angle of attack (size M x N)
    Set of panel lengths (p_lengths) in reverse-selig panel order (size M)
    Coefficient of Pressure distribution (coeff_p) across each panel across each angle of attack (size M x N)
    Outwards normal angle (beta) for every panel in reverse-selig panel order (size M)
    Set of angle of attacks (alphas_rad) in RADIANS (size N)
    Set of panel midpoint coordinates in reverse-selig panel order (size M x 2)
    Array of masked indices (invalid_panels) to match masked panels used in pressure distribution matrix

    Outputs:
    Coefficient of Lift derived from Kutta Joukowski theorem (cL_KJ) and from integrating Coefficient of Pressure (cL_P), (size N)
    Coefficient of Moment about Leading Edge (cmLE) and about Quarter Chord (cmqc) (size N)
    """
    if invalid_panels is None:
        invalid_panels = []

    M, N = gamma_distribution.shape
    
    # Derive lift coefficient in two ways; 
    # compute by summing circulation and using Kutta-Joukowski 
    c_l_KJ = np.zeros(N)
    c_mLE = np.zeros(N)
    # and integrate via pressure distribution
    c_l_P = np.zeros(N)

    # For high TE error, use masking on high outlier values on extremely small panels to make math smoother
    cp_mask = np.ones(M, dtype=bool)
    cp_mask[invalid_panels] = False

    for k in range(N):
        # Circulation-based lift: every panel except the spurious closing one carries real circulation.
        c_l_KJ[k] = 2 * np.sum(gamma_distribution[:-1, k] * p_lengths[:-1])

        arm = (midpoints[:-1,0] * np.cos(alphas_rad[k]) + midpoints[:-1,1] * np.sin(alphas_rad[k]))
        c_mLE[k] = -2 * np.sum(gamma_distribution[:-1,k] * p_lengths[:-1] * arm)
        # Pressure-based lift: only panels with a physically meaningful Cp.
        cp = coeff_p[cp_mask, k]
        s  = p_lengths[cp_mask]
        b  = beta[cp_mask]

        # Compute axial and normal force coefficients 
        c_n = -(cp * s * np.sin(b)).sum()
        c_a = -(cp * s * np.cos(b)).sum()

        c_l_P[k] = c_n * np.cos(alphas_rad[k]) - c_a * np.sin(alphas_rad[k])

    cmqc = c_mLE + 0.25 * c_l_KJ

    return c_l_KJ, c_l_P, c_mLE, cmqc
