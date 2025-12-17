"""
Main executable script for the 1D unsteady heat conduction project.
"""
 
from io_utils import load_data_from_csv
 


def main():
    
    #create material and plate objects with data from input_params.csv
    csv_path = "input_params.csv"
    material, plate = load_data_from_csv(csv_path)

    #Print material properties with material __str__ method
    print(material)  
 
    # Run implicit solver  
    solve_implicit(plate, steady_tol=1e-4)


    # estimate lumped-capacitance time to reach 95% of steady state
    rho = material.rho
    c = material.c
    h = material.h
    # Use math.log to estimate the lumped capacitance time scale
    # theta(t)/theta_end = 0.95 = 1 - exp(-h t / (rho c)) ->  t = -(rho c / h) * ln(0.05)
    t_lumped = -(rho * c / h) * log(0.05)

    print(f"Estimated steady state time (lumped capacitance) ≈ {t_lumped:.2f} s")

    # Generate all figures and print the numerical and analytical comparison
    plot_results(plate)

    


if __name__ == "__main__":
    main()
