# Numerical Simulation of Unsteady One-Dimensional Heat Conduction in a Plate Using Python 

## Team Members:
Fan Huang(fhuang8@stevens.edu) and  Liye Guo(lguo15@stevens.edu)

## Problem Description
The physical system in this simulation is a steel plate with a thickness of 2L, its symmetric about the mid-plane. Internal heat generation will occur uniformly throughout the plate, and both surfaces exchange heat with the surrounding fluid by convection. 

The governing one-dimensional transient heat conduction equation is: 

∂T/∂t = α ∂²T/∂x² + q'''/(ρc) -> where α = k/(ρc) is the thermal diffusivity. 
          
The boundary conditions are:
- Symmetry at the centerline: ∂T/∂x = 0 at x = 0
- Convection at the surface: -k ∂T/∂x = h (T - T∞) at x = L

The given parameters are thermal conductivity k, density ρ, specific heat c, convective coefficient h, ambient temperature T∞, and uniform heat generation q'''. The goal is to compute the transient temperature distribution and verify the steady state profile and the time to steady state, and to compare the numerical steady-state solution with the analytical solution and lumped-capacitance estimation.
 
The implicit finite difference method is used for numerical stability and accuracy.

## Program structure
**material.py**:  Material class, that stores the physical properties: thermal conductivity, density, specific heat,
 convective coefficient, internal heat generation, ambient temperature, and plate half thickness used in the 1D plate conduction problem.
 The program also computes derived properties such as thermal diffusivity and implements operator overloading and formatted output.

**heat_plate.py**: Defines the HeatPlate class, which performs the numerical simulation. This class assembles the coefficient matrix for
  the implicit finite difference scheme, updates the temperature field in time, and implements lambda-based boundary conditions. A method 
  provides a generator for temperature output, and the class includes plotting utilities and an analytical steady-state solution for comparison.

**io_utils.py**: Contains utility functions to read simulation parameters from a CSV input file and construct Material and HeatPlate objects. 

**main.py:** This file is the  main executable script, which reads input parameters, runs the implicit solver, estimates the time to steady state
  using lumped-capacitance theory, and generates all required plots.
  
**input_params.csv:** The Input file, which contains material properties, numerical parameters, and simulation time settings.

## To run program 
 ### Option 1: Run with the provided input parameters
 
1. **Download the project**
   - Click **Code → Download ZIP** on GitHub
   - Extract the files to your local machine
     
2. **Open a terminal**
   - Navigate to the project directory:
     ```bash
     cd heat_conduction
     ```
3. **Run the program**
   ```bash
   python main.py

### Option 2: Run with your own input parameters (optional)
1. **After downloading the project, open the file input_params.csv**
  
2. **Modify the  simulation parameters**
   - Update the values for the material and numerical settings:
     - `k` — thermal conductivity  
     - `rho` — density  
     - `c` — specific heat  
     - `h` — convective heat transfer coefficient  
     - `q_gen` — internal heat generation  
     - `T_inf` — ambient temperature  
     - `L` — plate half-thickness  
     - `n_nodes` — number of spatial nodes  
     - `dt` — time step size  
     - `t_final` — total simulation time  
     - `T_init` — initial temperature (optional)
       
3. **Save the file**
   
4. **Open a terminal**
   - Navigate to the project directory:
     ```bash
     cd heat_conduction
     ```
5. **Run the program**
   ```bash
   python main.py

## Notes 
- Python 3.x is required
- Required libraries: numpy, matplotlib

## Team Contributions
Fan Huang: heat_plate.py(Numerical model simulation solver) , material.py,input_params.csv

Liye Guo: heat_plate.py (Plots for the model simulation), main.py, io_utils.py
