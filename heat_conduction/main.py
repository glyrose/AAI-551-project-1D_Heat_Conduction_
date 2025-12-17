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

    


if __name__ == "__main__":
    main()
