# This is the script for NACA 4-series and 5-series airfoil plotting. 
# The geometric data points of the airfoil are exported in (x,y) coordinates in a .csv file

# Libraries for computation, plotting, and data collection
import numpy as np
import pandas as pd
import matplotlib as plt

# Collect NACA series airfoil from user input
NACA_dig = int(input("Please enter your NACA airfoil 4-digit or 5-digit code: "))
print(f"Your airfoil is: NACA {NACA_dig}")

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
print(f"Thickness % : {t}");

# Discretize the x axis from 0 to 1 with N points (user input)
N = int(input("Enter desired number of points (N): "));
x_axis = np.linspace(0, 1, N);

# Define mean camberline using known equations for NACA 4-digit series airfoils
y_mc = np.where((0 <= x_axis) & (x_axis < p), m/(p**2) * (2*p*x_axis - x_axis**2), m/((1-p)**2) * ((1-2*p) + 2*p*x_axis - x_axis**2));

print(x_axis);
print(y_mc);