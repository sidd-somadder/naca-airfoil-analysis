# This is the script for NACA 4-series airfoil plotting. 
# The geometric data points of the airfoil are exported in (x,y) coordinates in a .dat file

# Libraries for computation, numerical solving, plotting, and data exporting
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from pathlib import Path
import sys 

def naca_gen_script():
    # Collect NACA series airfoil from user input
    NACA_dig = get_naca4_input();

    # Get chordwise point count and total surface point count.
    N = get_mesh_param();

    # Use half-cosine spacing to make surface point distribution finer near trailing edge and leading edge
    beta = np.linspace(0, np.pi, N);
    x_axis = 0.5 * (1 - np.cos(beta));
    print("---");

    # At this point, assume that incorrect inputs have been filtered out by earlier input verification code block
    # Call function to get coordinate point matrix
    XY_coords = naca4gen(NACA_dig, x_axis);    

    print("---");

    # Call plotter method so user can see their airfoil visually
    naca_plot(NACA_dig, XY_coords);

    # Ask user if they want to save airfoil coordinates as .dat file
    save_coords_decision(NACA_dig, XY_coords);

# Plotter helper method called automatically in main script. 
def naca_plot(NACA_dig, XY_coords):
    print("Plotting airfoil...");
    # Extract X and Y values from output coordinate matrix
    x_plot = XY_coords[:, 0]
    y_plot = XY_coords[:, 1]

    code = "";
    if NACA_dig // 100 == 0:
        code = "00"
    code+= str(NACA_dig);

    # Plot airfoil geometry:
    plt.figure(figsize=(10, 4))

    # Including both line plot and scatterplot to show overall shape and individual points for visualization
    plt.plot(x_plot, y_plot, 'b-', linewidth=1.5)
    plt.scatter(x_plot, y_plot, s=5, color='blue', zorder=5) 

    # Format plot to be readable and prevent squashed look 
    plt.axis('equal') 
    plt.xlabel('x/c')
    plt.ylabel('y/c')
    plt.title(f'NACA {code} Airfoil Geometry')
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    print("---");

# Method that generates upper and lower coordinates using mean camberline and thickness distribution equations
def naca4gen(NACA_dig, x_axis):
    # For 4-digit airfoils, compute max camber (m), p (max camber location), and t (thickness) values in percentage chord
    m_dig = NACA_dig // 1000;
    p_dig = ((NACA_dig % 1000) // 100); 
    t_dig = NACA_dig % 100;

    # Convert digits into % chord values for each parameter 
    m = m_dig/100;
    p = p_dig/10;
    t = t_dig/100;

    # Print % chord values
    print(f"NACA {NACA_dig} details: ")
    print(f"Max Camber % : {m}");
    print(f"Max Camber Position % : {p}");
    print(f"Max Thickness % : {t}");

    # Handle symmetric vs asymmetric camber line calculations explicitly
    if m == 0 or p == 0:
        # Symmetric case: no camber line or slope variations
        y_mc = np.zeros_like(x_axis);
        dy_mc = np.zeros_like(x_axis);
    else:
        # Asymmetric case: define mean camberline & slope using standard equations
        y_mc = np.where((0 <= x_axis) & (x_axis < p), 
                        m/(p**2) * (2*p*x_axis - x_axis**2), 
                        m/((1-p)**2) * ((1-2*p) + 2*p*x_axis - x_axis**2));
        
        dy_mc = np.where((0 <= x_axis) & (x_axis < p), 
                         2*m/(p**2) * (p - x_axis), 
                         2*m/((1-p)**2) * (p - x_axis));

    # Calculate thickness distribution of NACA airfoil using known equation
    y_t = 5*t*(0.2969 * np.sqrt(x_axis) - 0.1260*x_axis - 0.3516*(x_axis**2) + 0.2843*(x_axis**3) - 0.1015*(x_axis**4));

    # Calculate tangential angle of the mean camberline slope
    theta = np.arctan(dy_mc);

    # Define upper and lower airfoil camber line coordinates
    x_upper = x_axis - y_t * np.sin(theta);
    x_lower = x_axis + y_t * np.sin(theta);
    y_upper = y_mc + y_t * np.cos(theta);
    y_lower = y_mc - y_t * np.cos(theta);

    # Process upper/lower x,y coordinates into a combined X,Y linspaces, respectively
    # Note that the points go from TE to LE via upper surface and then LE to TE via lower surface (Selig format).
    # Clip the first index of each upper coordinates arrays to prevent doubling the LE point before flipping.
    X = np.concatenate([np.flip(x_upper[1:]), x_lower]);
    Y = np.concatenate([np.flip(y_upper[1:]), y_lower]);

    # Normalize X to be 0 to 1
    chord = np.max(X) - np.min(X);
    X = X/chord;

    # Reverse the fully assembled Selig-ordered array to get reverse Selig:
    # TE -> lower surface -> LE -> upper surface -> TE (clockwise), matching the
    # traversal direction established for the VPM sign convention.
    X = np.flip(X);
    Y = np.flip(Y);

    # Combine combined X,Y coordinate linspaces into 2N-1 x 2 matrix
    XY_coords = np.column_stack((X,Y));
    return XY_coords;

def save_coords_decision(NACA_dig, XY_coords):
    # Before writing coordinate points to .dat file, check if dedicated folder exists to avoid conflicts
    Path("saved_airfoil_coords").mkdir(exist_ok=True)

    while True:
            # Read user input, strip whitespace, and convert to uppercase
            choice = input("Would you like to save your airfoil surface points to a .dat file? (Y/N): ").strip().upper()
        
            if choice == 'Y':
                # Build filename and path. Example format: "NACA_4412_N100.dat" would be written to saved_airfoils_coords
                
                code = "";
                if NACA_dig // 100 == 0:
                    code = "00"
                code+= str(NACA_dig);
                
                N = int(0.5*(len(XY_coords) + 1))
                filename = f"NACA_{code}_N{N}.dat"
                filepath = Path("saved_airfoil_coords") / filename

                with open(filepath, 'w') as f:
                    # Header line: airfoil designation (standard Selig format)
                    f.write(f"NACA {code}\n")
        
                    # Write x y coordinates, space-separated, 6 decimal places
                    for row in XY_coords:
                        f.write(f"  {row[0]:.12e}  {row[1]:.12e}\n")

                print(f"Coordinates saved to: {filepath}")
                break
            
            elif choice == 'N':
                # Do nothing
                print("Coordinates not saved.")
                break
                
            else:
                print("Invalid choice. Please enter Y or N.")

def get_mesh_param():
    while True:
        # Read user input for desired N chordwise points
        try:
            N = int(input("Enter desired number of chordwise points along x-axis (N): "))
        except ValueError:
            print("Invalid input. Please enter a whole number.")
            print("---")
            continue

        # Verify N is a positive integer greater than a practical minimum
        if N < 2:
            print(f"Invalid input: N must be at least 2. You entered {N}.")
            print("---")
            continue
        
        # Give user warning for small point count for practical reasons, but don't stop them until they confirm so 
        if N < 10:
            print(f"Warning: N = {N} is very low and may produce inaccurate geometry. Recommended minimum is N = 50.")

        print("---")

        # Inform user of expected number of surface points
        exp_ptnum = 2 * N - 1
        print(f"User Note: {N} chordwise points will be used for computation; "
              f"this generates {exp_ptnum} total surface points.")

        # Get user confirmation
        while True:
            choice = input("Proceed? (Y/N): ").strip().upper()

            if choice == 'Y':
                # Return chordwise point count and total surface point count
                # Inform user of their chosen number of chordwise points and how many coordinate points generated by plotter
                print(f"Proceeding to airfoil plotting with {N} chordwise points and {exp_ptnum} total surface points...");
                return N;

            elif choice == 'N':
                print("---");
                break;

            else:
                print("Invalid choice. Please enter Y or N.");

def get_naca4_input():
    while True:
        # Read as string to preserve leading zeros (e.g. 0012)
        raw = input("Please enter your NACA 4-digit airfoil code: ").strip()

        # Check exactly 4 digits, no letters or symbols
        if not (len(raw) == 4 and raw.isdigit()):
            print(f"Invalid input: '{raw}' is not a 4-digit code. Please enter exactly 4 digits (e.g. 0012, 2412).")
            print("---")
            continue

        # Parse individual digit groups
        m_dig = int(raw[0])      # max camber
        p_dig = int(raw[1])      # max camber location
        t_dig = int(raw[2:])     # max thickness (last two digits as a number)

        # first two digits must be both zero or both nonzero
        if (m_dig == 0) != (p_dig == 0):
            print(f"Invalid input: first two digits must be both zero (symmetric) or both nonzero (cambered).")
            print(f"  You entered m={m_dig}, p={p_dig}, one is zero and the other is not.")
            print("---")
            continue

        # thickness must be greater than 0
        if t_dig == 0:
            print(f"Invalid input: thickness digits '{raw[2:]}' give t=0, which is not a valid airfoil.")
            print("---")
            continue

        # Valid input — convert and return
        NACA_dig = int(raw)
        print(f"Your airfoil is: NACA {raw}")
        print("---")
        return NACA_dig;

naca_gen_script();
