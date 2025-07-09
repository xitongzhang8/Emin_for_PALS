*****************************************************************************************************************************
Free Energy Minimization Code for the Manuscript:
**Tunable Wetting of Droplets on Patterned Liquid Surfaces Supplementary Materials**

This code is designed to simulate the static wetting behavior of a droplet on the Patterned Liquid Surfaces (PALS).
*****************************************************************************************************************************

**Features:**
- Compute the equilibrium state of a four-phase fluid system using free energy minimization.
- Apply constraints (volume/pressure) to different phases.
- Utilize a confinement potential to stabilize the lubricant layer and form the desired pattern.
- Support parameter configuration through the `input.txt` file.

-----------------------------------------------------------------------------------------------------------------------------
**Requirements:**
To use this script, ensure you have the following environment:
- Python 3.4 or higher
- Fortran 90 compiler

-----------------------------------------------------------------------------------------------------------------------------
**Governing file:**
The Python script energy_minimization.py is the core driver for generating the necessary data files and managing the 
computational workflow. 

This script is responsible for:
- Computational domain and parameters
- Fluid properties
- Phase constraints
- Phase initialization
- Simulation logic flow

-----------------------------------------------------------------------------------------------------------------------------
**Running the Script:**

1. Configure parameters in the `input.txt` file and the `energy_minimization.py` script.
2. Execute the Python script `energy_minimization.py`.  
    - The script reads simulation parameters from `input.txt`, creates a `data` folder, and generates the `data.in` file.
    - It then calls the executables `./GMIN` and `./gmin` to perform the computations.

-----------------------------------------------------------------------------------------------------------------------------
**Key User Inputs:**
- **NPOSTX, NPOSTY**: Number of patterns in the x and y directions.
- **WIDTHX, WIDTHY**: Dimensions of the patterns in the x and y directions.
- **HEIGHT1**: Depth of the lubricant layer.
- **GRIDX, GRIDY, GRIDZ**: Domain size in the x, y, and z directions.
- **N_PHASE**: Number of phases to simulate.
- **alpha**: Interface width.
- **Interfacial Tensions**: Tension values between different fluid phases.
- **Phase Constraints**: Volume or pressure constraints applied to different phases.

-----------------------------------------------------------------------------------------------------------------------------
**Output Files:**

- **data_log**: Logs the simulation parameters and energy values for each run.
- **coords_before**: Stores the system's coordinates prior to droplet initialization.
- **coords_after** : Stores the system's coordinates after droplet initialization.
- **coords.step**  : Stores the system's coordinates during the process of free energy minimization.
- **lowest_test_1**: A binary file containing the final minimized configuration.  
  - Includes the distribution of order parameters for the three phases (excluding the bulk phase) in the format: `[x, y, z, c1, c2, c3]`.  
  - The order parameter for the bulk phase (`c4`) can be calculated as:  
    `c4 = 1.0 - c1 - c2 - c3`.

-----------------------------------------------------------------------------------------------------------------------------
**License:**
This code is provided as-is. You may modify it as needed, but no warranty or liability is provided.
