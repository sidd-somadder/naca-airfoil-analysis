"""
This script will run Thin Airfoil Theory derivations and is called automatically by the master script
.csv data for the lift, leading edge moment, and quarter-chord moment coefficients against angle of attacks is outputted
The chosen .dat file in saved_airfoil_coords folder must have proper naming: "NACA_XXXX_N#.dat" for this script to execute successfully
Currently, this thin airfoil theory solver's scope is to only handle NACA 4-digit series 
"""

import numpy as np
import scipy.integrate as scpi
from post_processing import export_tat_results

def run_tat_solver(input_file_name, alphas):
    """
    Runs the Thin Airfoil Theory solver on a NACA 4-digit airfoil over an angle of attack sweep.

    Inputs:
    Filename (input_file_name) as found in the saved_airfoil_coords folder, following the
    "NACA_XXXX_N#.dat" naming convention. The 4-digit code is parsed directly from the filename;
    the coordinate data itself is not read, since TAT operates on the analytical camber line.
    Angle of Attack sweep (alphas) in DEGREES (size N).

    Outputs:
    Coefficient of Lift (c_l), Coefficient of Moment about the Leading Edge (c_mLE),
    and about the Quarter Chord (c_mqc), each of size N.
    Results are also exported to .csv in the computation_results folder.
    Raises ValueError if the filename does not contain a valid 4-digit NACA code.
    """

    alphas_rad = alphas * (np.pi / 180)

    # Extract the numeric code from the required "NACA_XXXX_N#.dat" naming convention.
    code_name = input_file_name.split("_")[1]

    try:
        code = int(code_name)
    except ValueError:
        raise ValueError(f"'{code_name}' is not a valid NACA numeric code.")

    if len(code_name) != 4:
        raise ValueError(f"'{code_name}' is not a valid 4-digit NACA code.")

    symmetric = code // 100 == 0
    if symmetric:
        coeffs, angle_zero_lift = sym_4digit_solver(alphas_rad)
    else:
        coeffs, angle_zero_lift = asym_4digit_solver(code, alphas_rad)

    c_l = coeffs[:,0]
    c_mLE = coeffs[:,1]
    c_mqc = coeffs[:,2]

    export_tat_results(alphas, coeffs, input_file_name, angle_zero_lift)

    return c_l, c_mLE, c_mqc

def sym_4digit_solver(alphas_rad):
    """
    Returns thin airfoil theory coefficients for a symmetric (NACA 00XX) airfoil.

    Inputs:
    Angle of Attack sweep (alphas_rad) in RADIANS (size N).

    Outputs:
    Coefficient matrix (size N x 3) with columns ordered as [c_l, c_mLE, c_mqc],
    and the zero-lift angle of attack in DEGREES.
    Uses the closed-form results for zero camber: 
    c_l = 2*pi*alpha, c_mLE = -c_l/4, c_mqc = 0, and a zero-lift angle of 0 degrees.
    """

    # Use known thin airfoil theory results for symmetric airfoils.
    c_l = 2 * np.pi * alphas_rad
    c_mLE = (-1)*(c_l)/4
    c_mqc = np.zeros_like(alphas_rad)
    zero_lift_angle = 0.0

    return np.column_stack((c_l, c_mLE, c_mqc)), zero_lift_angle

def asym_4digit_solver(code, alphas_rad):
    """
    Returns thin airfoil theory coefficients for a cambered NACA 4-digit airfoil.

    Inputs:
    NACA 4-digit code as an integer (code), used to extract max camber and camber position.
    Angle of Attack sweep (alphas_rad) in RADIANS (size N).

    Outputs:
    Coefficient matrix (size N x 3) with columns ordered as [c_l, c_mLE, c_mqc],
    and the zero-lift angle of attack in DEGREES.

    The mean camber line slope is split at the max camber position and the Glauert
    transformation x/c = 0.5*(1 - cos(theta)) is applied. Fourier coefficients A_0, A_1,
    and A_2 are obtained by numerical quadrature over each segment, then used to derive
    the coefficients. Note that c_mqc = (pi/4)*(A_2 - A_1) is independent of angle of
    attack, so it is returned as a constant array.
    """

    # Extract digits from 4-digit series code
    m_dig = code // 1000
    p_dig = ((code % 1000) // 100)

    # Define mean camberline slope equations forward and aft of max camber location
    # Apply Glauert transformation (x/c = 0.5(1-cos(theta_0))) 
    def dz1(t0,m,p):
        return (2*m)/(p**2) * (p - 0.5*(1-np.cos(t0)))
            
    def dz2(t0,m,p):
        return (2*m)/((1-p)**2) * (p - 0.5*(1-np.cos(t0)))
                
    # Convert NACA digits into chordwise % values
    m = m_dig / 100
    p = p_dig / 10

    # Use Glauert transformation to find angle of max camber position
    t_c = np.arccos(1-2*p)

    # Integrate for Fourier coefficients
    A_0 = alphas_rad - (1/np.pi) * (scpi.quad(dz1, 0, t_c, args=(m,p))[0] + scpi.quad(dz2, t_c, np.pi, args=(m,p))[0])

    A_1 = (2/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * np.cos(t0), 0, t_c)[0] 
                        + scpi.quad(lambda t0: dz2(t0, m, p) * np.cos(t0), t_c, np.pi)[0])
    
    A_2 = (2/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * np.cos(2*t0), 0, t_c)[0] 
                        + scpi.quad(lambda t0: dz2(t0, m, p) * np.cos(2*t0), t_c, np.pi)[0])

    # Use Fourier coefficients to derive lift and moments coefficients of asymmetric airfoils
    c_l = np.pi * (2*A_0 + A_1)
    c_mLE = (-1)*(c_l/4) - (np.pi/4)*(A_1 - A_2)
    # Note, quarter-chord moment coefficient is theoretically constant
    c_mqc = np.full_like(alphas_rad, (np.pi/4) * (A_2 - A_1), dtype=float)

    # Compute zero lift angle of attack by expanding c_l and setting to 0
    zero_lift_angle = (1/np.pi) * (scpi.quad(lambda t0: dz1(t0, m, p) * (1-np.cos(t0)), 0, t_c)[0] 
                                + scpi.quad(lambda t0: dz2(t0, m, p) * (1-np.cos(t0)), t_c, np.pi)[0])

    # express zero lift angle in degrees
    zero_lift_angle = (180/np.pi)*zero_lift_angle

    # Return coefficients in matrix form and zero lift angle of attack as a tuple
    return np.column_stack((c_l, c_mLE, c_mqc)), zero_lift_angle
