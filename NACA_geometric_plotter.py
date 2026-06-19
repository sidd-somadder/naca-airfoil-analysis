# This is the script for NACA 4-series and 5-series airfoil plotting. 
# The geometric data points of the airfoil are exported in (x,y) coordinates in a .csv file

# Libraries for computation, plotting, and data collection
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Collect NACA series airfoil from user input
NACA_dig = int(input("Please enter your NACA airfoil 4-digit or 5-digit code: "))
print(f"Your airfoil is: NACA {NACA_dig}")

# Currently uniform spacing, will change to half-cosine spacing to be computationally efficient near LE/TE
can_proceed_disc = False;
while can_proceed_disc is False:
    # Discretize the x axis from 0 to 1 with N points (user input)
    N = int(input("Enter desired number of chordwise points (N): "));
    x_axis = np.linspace(0, 1, N); 
    
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
            break
            
        else:
            print("Invalid choice. Please enter Y or N.")

print(f"Proceeding to airfoil plotting with {N} chordwise points and {exp_ptnum} total surface points...");

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
    print(f"Max Camber % : {m}");
    print(f"Max Camber Position % : {p}");
    print(f"Max Thickness % : {t}");

    # Define mean camberline & slope using known equations for NACA 4-digit series airfoils
    y_mc = np.where((0 <= x_axis) & (x_axis < p), m/(p**2) * (2*p*x_axis - x_axis**2), m/((1-p)**2) * ((1-2*p) + 2*p*x_axis - x_axis**2));
    dy_mc = np.where((0 <= x_axis) & (x_axis < p), 2*m/(p**2) * (p - x_axis), 2*m/((1-p)**2) * (p - x_axis));

    # Define thickness across airfoil chord using max thickness % (t)
    y_t = 5*t*(0.2969 * np.sqrt(x_axis) - 0.1260*x_axis - 0.3516*(x_axis**2) + 0.2843*(x_axis**3) - 0.1015*(x_axis**4));
    theta = np.arctan(dy_mc);

    # Define upper and lower airfoil camber line coordinates
    x_upper = x_axis - y_t * np.sin(theta);
    x_lower = x_axis + y_t * np.sin(theta);
    y_upper = y_mc + y_t * np.cos(theta);
    y_lower = y_mc - y_t * np.cos(theta);

    # Process upper/lower x,y coordinates into a combined X,Y linspaces, respectively
    X = np.concatenate([np.flip(x_lower[1:]), x_upper]);
    Y = np.concatenate([np.flip(y_lower[1:]), y_upper]);

    # Combine combined X,Y coordinate linspaces into 2N-1 x 2 matrix
    XY_coords = np.column_stack((X,Y));


# Assume otherwise user input has 5 digits, will handle improper inputs in future
else:
    # Handle NACA 5-digit airfoils (placeholder)
    # For now, notify that 5-digit airfoils are not available
    print("5-digit NACA airfoil handling is not implemented in this script yet.")

    # Placeholder: mean camberline is zeroes
    N = int(input("Enter desired number of points (N): "))
    y_mc = np.zeros_like(x_axis)
    print(x_axis)
    print(y_mc)