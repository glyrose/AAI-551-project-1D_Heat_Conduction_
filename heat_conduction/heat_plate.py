"""Defines the HeatPlate class, which performs the numerical simulation. 
This class assembles the coefficient matrix for the implicit finite difference scheme 
updates the temperature field in time, implements lambda-based boundary conditions,
provides a generator for temperature output, and includes plotting utilities and an analytical steady-state solution for comparison. """

import math
import numpy as np
import matplotlib.pyplot as plt
from typing import List, Dict, Tuple
from material import Material


class HeatPlate:
    

    def __init__(self,material: Material,n_nodes,dt,t_final,T_init):
        
        #check parameters and ensure values are within range
        if n_nodes < 2:
            raise ValueError("n_nodes must be at least 2.")
        if dt <= 0 or t_final <= 0:
            raise ValueError("dt and t_final must be positive.")

        self.material = material
        self.n_nodes = n_nodes  
        self.dt = dt
        self.t_final = t_final
        self.dx = material.L / n_nodes

        #fourier number
        self.r = material.alpha() * dt / self.dx**2

        #create the spatial grid for the plate 
        self.x: List[float] = [i * self.dx for i in range(n_nodes + 1)]

        
        self.metadata: Dict[str, object] = {
            "dx": self.dx,
            "r": self.r,
            "x_range": (0.0, material.L),
            "description": "1D plate with internal heat generation and convection",
        }

        # intitialize initial temperature field  
        if T_init is None:
            T_init = material.T_inf
            
        self.T_current = np.full(n_nodes + 1, float(T_init))
        self.time_points: List[float] = [0.0]
        self.T_history: List[np.ndarray] = [self.T_current.copy()]

        # Lambda functions for the  boundary condition coefficients  
        self.center_bc = lambda r: (1.0 + 2.0 * r, -2.0 * r)
        self.surface_bc = lambda r, h, dt, dx, rho, c: (
            1.0 + 2.0 * r + 2.0 * h * dt / (dx * rho * c),
            -2.0 * r,
        )
        
        self.A = self._assemble_matrix()

         

 
    #Assemble the (n+1)x(n+1) coefficient matrix A for the implicit finite-difference scheme.
    def _assemble_matrix(self):
    
        n = self.n_nodes
        r = self.r
        k = self.material.k
        rho = self.material.rho
        c = self.material.c
        h = self.material.h
        dt = self.dt
        dx = self.dx

        A = np.zeros((n + 1, n + 1), dtype=float)

        #define center node  x = 0, symmetry: dT/dx = 0
        a00, a01 = self.center_bc(r)
        A[0, 0] = a00
        A[0, 1] = a01

        #define interior nodes: 1 .. n-1
        for i in range(1, n):
            A[i, i] = 1.0 + 2.0 * r
            A[i, i - 1] = -r
            A[i, i + 1] = -r

        #define surface node x = L, convection
        ann, annm1 = self.surface_bc(r, h, dt, dx, rho, c)
        A[n, n] = ann
        A[n, n - 1] = annm1

        return A

 
    #Allow HeatPlate object to be called like a function to run the simulation.
    def __call__(self, steady_tol):
        
        #calls the solver and  returns its result 
        return self.run(steady_tol)
  
    #perform one implicit time step given the temperature at the previous time level.
    def step_implicit(self, T_old) :
     
        n = self.n_nodes
        rho = self.material.rho
        c = self.material.c
        h = self.material.h
        q_gen = self.material.q_gen
        T_inf = self.material.T_inf
        dt = self.dt
        dx = self.dx

        #right hand side vector b
        b = np.zeros(n + 1, dtype=float)

        #volumetric generation term that appears through the formula
        q_term = q_gen * dt / (rho * c)

        #define center node
        b[0] = T_old[0] + q_term

        #define interior nodes
        for i in range(1, n):
            b[i] = T_old[i] + q_term

        #define surface node (i=n) with convection contribution
        b[n] = T_old[n] + 2.0 * h * dt / (dx * rho * c) * T_inf + q_term

        #Solve the linear system A T_new = b -> T_new = b\A
        T_new = np.linalg.solve(self.A, b)
        return T_new

    #run the simulation until t_final or until the maximum nodal change is below steady_tol.
    def run(self, steady_tol) :
   
        n_steps_max = int(math.ceil(self.t_final / self.dt))
        T_old = self.T_current.copy()

 
        for step in range(1, n_steps_max + 1):
            t = step * self.dt
            T_new = self.step_implicit(T_old)

            self.time_points.append(t)
            self.T_history.append(T_new.copy())

            #convergence check to stop early if nearly steady
            max_diff = np.max(np.abs(T_new - T_old))
            if max_diff < steady_tol:
                break

            T_old = T_new

        times = np.array(self.time_points)
        T_all = np.vstack(self.T_history)

        #update current state
        self.T_current = T_all[-1, :]
        return times, T_all

    
  
    #yield temperature profiles one by one as time and T_profile pairs. 
    def temperature_generator(self):
        for t, T in zip(self.time_points, self.T_history):
            yield t, T

 
    #analytical steady state solution for comparison.
    def analytical_steady_profile(self):
       
        q = self.material.q_gen
        k = self.material.k
        h = self.material.h
        T_inf = self.material.T_inf
        L = self.material.L

        #T(x) = T_inf + (q/(2k)) (L² - x²) + (q L / h)
        coeffs: Tuple[float, float] = (q / (2.0 * k), q * L / h)
        a, b = coeffs

        x_arr = np.array(self.x)
        return T_inf + a * (L**2 - x_arr**2) + b

  
    