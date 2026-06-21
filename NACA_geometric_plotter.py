# This is the script for NACA 4-series and 5-series airfoil plotting. 
# The geometric data points of the airfoil are exported in (x,y) coordinates in a .csv file

# Libraries for computation, numerical solving, plotting, and data exporting
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from pathlib import Path
import sys 

# Collect NACA series airfoil from user input
NACA_dig = int(input("Please enter your NACA airfoil 4-digit or 5-digit code: "))
#
# This needs verification step to handle unacceptable inputs (not 4 or 5 digits)
# Check if 5 digit input, that Q digit is valid (0 or 1)
# Check if 4 digit input has valid m and p values (either both or neither must be 0)
#
print(f"Your airfoil is: NACA {NACA_dig}")
print("---")

# Currently uniform spacing, will change to half-cosine spacing to be computationally efficient near LE/TE
can_proceed_disc = False;
while can_proceed_disc is False:
    # Read user input for desired N chordwise points
    N = int(input("Enter desired number of chordwise points (N): "));
    print("---");
    #
    # This needs verification step to handle unacceptable inputs
    #
    # Inform user of expected number of surface points as opposed to N chordwise points used for computation
    exp_ptnum = 2*N - 1;
    print(f"User Note: {N} points along x-axis will be used for computation; running this script generates expected {exp_ptnum} upper/lower surface points.");

    # Get user confirmation for input & expected number of points
    while True:
    # Read user input, strip whitespace, and convert to uppercase
        choice = input("Proceed? (Y/N): ").strip().upper()
    
        if choice == 'Y':
            can_proceed_disc = True;
            break
        
        elif choice == 'N':
            # Do nothing
            break
            
        else:
            print("Invalid choice. Please enter Y or N.")
    print("---");

print(f"Proceeding to airfoil plotting with {N} chordwise points and {exp_ptnum} total surface points...");
# Proceed with creating x-axis linspace with user-confirmed N points
x_axis = np.linspace(0, 1, N); 
print("---");

# At this point, assume that incorrect inputs have been filtered out by earlier input verification code block
# If input has 4 digits
if (NACA_dig // 10000) == 0:
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

    # Define mean camberline & slope using known equations for NACA 4-digit series airfoils
    y_mc = np.where((0 <= x_axis) & (x_axis < p), m/(p**2) * (2*p*x_axis - x_axis**2), m/((1-p)**2) * ((1-2*p) + 2*p*x_axis - x_axis**2));
    dy_mc = np.where((0 <= x_axis) & (x_axis < p), 2*m/(p**2) * (p - x_axis), 2*m/((1-p)**2) * (p - x_axis));

    # Define thickness across airfoil chord using max thickness % (t)
    y_t = 5*t*(0.2969 * np.sqrt(x_axis) - 0.1260*x_axis - 0.3516*(x_axis**2) + 0.2843*(x_axis**3) - 0.1036*(x_axis**4));
    theta = np.arctan(dy_mc);

    # Define upper and lower airfoil camber line coordinates
    x_upper = x_axis - y_t * np.sin(theta);
    x_lower = x_axis + y_t * np.sin(theta);
    y_upper = y_mc + y_t * np.cos(theta);
    y_lower = y_mc - y_t * np.cos(theta);

    # Process upper/lower x,y coordinates into a combined X,Y linspaces, respectively
    # Note that the points go from TE to LE via upper surface and then LE to TE via lower surface (Selig format).
    # Clip the first index of each upper coordinates arrays to prevent doubling the LE point.
    X = np.concatenate([np.flip(x_upper[1:]), x_lower]);
    Y = np.concatenate([np.flip(y_upper[1:]), y_lower]);

    # Combine combined X,Y coordinate linspaces into 2N-1 x 2 matrix
    XY_coords = np.column_stack((X,Y));

    # Close the trailing edge at (1,0) since it current surface point calculation leaves floating point arithmetic error
    # in addition Kutta condition is satisfied and VPM application is significantly simplified 
    XY_coords[0, :]  = [1, 0]   # upper trailing edge
    XY_coords[-1, :] = [1, 0]   # lower trailing edge

# Assume otherwise 5 digits, follow procedure for handling 5-digit series
else:
    # For NACA 5-digit airfoils, first extract L, P, Q, TT values via integer division & modulo
    L_dig  = NACA_dig // 10000
    P_dig  = (NACA_dig % 10000) // 1000
    Q_dig  = (NACA_dig % 1000) // 100
    t_dig  = NACA_dig % 100

    # Derive necessary parameters with digits and known convention
    cl_design = 3 * L_dig / 20
    x_mc = P_dig / 20
    t = t_dig / 100  
    # Boolean thats easier to interpret than just reading the value of Q
    reflexed  = (Q_dig == 1)    

    # Print % chord values
    print(f"NACA {NACA_dig} details: ")
    print(f"Design Lift Coefficient : {cl_design}");
    print(f"Max Camber Position % : {x_mc}");
    print(f"Max Thickness % : {t}");
    if reflexed:
        print("Reflex: Yes");
    else:
        print("Reflex: No");
    print("---");
    print("WARNING: This program has yet to implement a 5-digit plotter until computational scripts have been developed for the 4-digit series");
    
    # Putting pause on developing a 5-series plotter with user warning
    sys.exit()

    # Define equation to derive parameter r needed for mean camberline equation
    def r_equation(r, x_mc):
        return r * (1 - np.sqrt(r / 3)) - x_mc
    # Search in [0, 1] to get the smaller (physically relevant) root
    r_sol = brentq(r_equation, 0, 1, args=(x_mc,))
    
    print("...") # Placeholder

print("---");
print("Plotting airfoil...");
# Extract X and Y values from output coordinate matrix
x_plot = XY_coords[:, 0]
y_plot = XY_coords[:, 1]

# Plot airfoil geometry:
plt.figure(figsize=(10, 4))

# Including both line plot and scatterplot to show overall shape and individual points for visualization
plt.plot(x_plot, y_plot, 'b-', linewidth=1.5)
plt.scatter(x_plot, y_plot, s=5, color='blue', zorder=5) 

# Format plot to be readable and prevent squashed look 
plt.axis('equal') 
plt.xlabel('x/c')
plt.ylabel('y/c')
plt.title(f'NACA {NACA_dig} Airfoil Geometry')
plt.grid(True)
plt.tight_layout()
plt.show()
print("---");

# Before writing coordinate points to .csv file, check if dedicated folder exists to avoid conflicts
Path("saved_airfoil_coords").mkdir(exist_ok=True)

while True:
        # Read user input, strip whitespace, and convert to uppercase
        choice = input("Would you like to save your airfoil surface points to a .csv file? (Y/N): ").strip().upper()
    
        if choice == 'Y':
            # Build filename and path. Example format: "NACA_4412_N100.csv" would be saved to saved_airfoils_coords
            filename = f"NACA_{NACA_dig}_N{N}.csv"
            filepath = Path("saved_airfoil_coords") / filename

            # Create DataFrame from saved coordinate matrix and write to .csv
            df = pd.DataFrame(XY_coords, columns=["x", "y"])
            df.to_csv(filepath, index=False)

            print(f"Coordinates saved to: {filepath}")
            break
        
        elif choice == 'N':
            # Do nothing
            print("Coordinates not saved.")
            break
            
        else:
            print("Invalid choice. Please enter Y or N.")
