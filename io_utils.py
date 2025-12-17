 
"""
This file contains the utility functions to read simulation parameters from CSV, 
building Material object, and Heatplate object.
"""
import os
import pandas as pd
from material import Material
from heat_plate import HeatPlate

#read the simulation data from CSV to create Material and HeatPlate objects
def load_data_from_csv(csv_path) :

    #check if file exist 
    if not os.path.exists(csv_path):
        #rasie error if file does not exist 
         raise FileNotFoundError(f"Input CSV file not found: {csv_path}")
     
     
    #if file exist read into a pandas dataframe
    df = pd.read_csv(csv_path)
    #if dataframe is empty raise an error
    if df.empty:
        raise ValueError("Input CSV file is empty.")

    #read the first row of CSV file
    row = df.iloc[0]

    #create material object
    material = Material(
        k = float(row["k"]),
        rho = float(row["rho"]),
        c = float(row["c"]),
        h = float(row["h"]),
        q_gen = float(row["q_gen"]),
        T_inf = float(row["T_inf"]),
        L = float(row["L"]),
        
    )

    # Read simulation parameters from CSV and convert to the appropriate data type 
    n_nodes = int(row["n_nodes"])
    dt = float(row["dt"])
    t_final = float(row["t_final"])
    T_init = float(row["T_init"]) if "T_init" in row and not pd.isna(row["T_init"]) else None

    #create HeatPlate object
    plate = HeatPlate(
        material = material,
        n_nodes = n_nodes,
        dt = dt,
        t_final = t_final,
        T_init = T_init,
    )
    return material, plate
