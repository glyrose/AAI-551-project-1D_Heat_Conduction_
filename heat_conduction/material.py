"""
 Material class,that stores the physical properties : thermal conductivity, density, specific heat,
 convective coefficient, internal heat generation, ambient temperature, and plate half-thickness used in the 1D plate conduction problem.
 Also computes derived properties  such as thermal diffusivity and implements operator overloading and formatted output.
"""
 
#Material class to  stores physical parameters and provides utility methods.
class Material:
    def __init__(self,k,rho,c,h,q_gen,T_inf,L):
        
        self.k = k 
        self.rho = rho 
        self.c = c
        self.h = h
        self.q_gen = q_gen
        self.T_inf = T_inf
        self.L = L
      
        # Basic validation of parameters
        if self.k <= 0 or self.rho <= 0 or self.c <= 0:
            raise ValueError("k, rho and c must all be positive.")
        if self.L <= 0:
            raise ValueError("Plate half-thickness L must be positive.")
        if self.h < 0:
                raise ValueError("Convective coefficient h cannot be negative.")
 
    
    #compute the thermal diffusivity α = k/(ρc).
    def alpha(self):
        return self.k / (self.rho * self.c)
    
    #Scale internal heat generation by a scalar for  parametric testing 
    def __mul__(self, factor):
        
        if not isinstance(factor, (int, float)):
            raise TypeError("Material can only be multiplied by a scalar.")
        
        return Material(
            k = self.k,
            rho = self.rho,
            c = self.c,
            h = self.h,
            q_gen = self.q_gen * factor,
            T_inf = self.T_inf,
            L = self.L,
        )
    #Format description of all properties to display
    def __str__(self):
        lines = [
            "Material properties:",
            f"  k      = {self.k:.3g} W/(m·K)",
            f"  rho    = {self.rho:.3g} kg/m³",
            f"  c      = {self.c:.3g} J/(kg·K)",
            f"  h      = {self.h:.3g} W/(m²·K)",
            f"  q_gen  = {self.q_gen:.3g} W/m³",
            f"  T_inf  = {self.T_inf:.3g} °C",
            f"  L      = {self.L:.3g} m (half thickness)",
            f"  alpha  = {self.alpha():.3e} m²/s",
        ]
        return "\n".join(lines)
