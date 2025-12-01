# 🥊 Force IBI Project

This project implements Iterative Boltzmann Inversion (IBI) and related coarse-grained force field optimization techniques. The code has been updated to use Julia as the main computation engine for better performance, while maintaining Python for pre/post-processing.

> ⚠️ **Note:** This repository is for personal use. Dependencies and environment setup are assumed to be manually managed.  
---
# ⚙️ Usage Guide
## Quick Start for Main Inversion

If you're already familiar with the setup process, you can jump directly to running the main inversion routine:

1. Ensure your parameter file `./inputs/params.json` is properly configured.
2. Execute the following command from the project root:
    ```bash
    julia -t auto forceIBI/grinter_parallel.jl
    ```
3. The results will be saved in the `./outputs/` directory for further analysis.

For detailed steps, refer to the sections below.
## Full Workflow
### 1. Prepare Input Configurations
Place your configuration files in the `./inputs/configs/` directory. These files will serve as the starting point for the inversion process.

### 2. Generate Targets for Inversion
Run the following scripts from the project root to prepare the target data for the inversion:

1. Minimize the input data:
    ```bash
    python3 forceIBI/min_d.py
    ```

2. Create histograms for the radial distribution function (RDF):
    ```bash
    python3 forceIBI/gr_histo.py <dr>
    ```
    Replace `<dr>` with the desired bin width for the RDF.

3. Generate weighted target data:
    ```bash
    python3 forceIBI/weighted_gen_target.py
    ```

4. Create a dummy `.json` parameter file:
    ```bash
    python3 init_dummy_json.py
    ```

### 3. Run the Main Computation
Edit the parameter file `./inputs/paramters.json`. 
Execute the main iterative Boltzmann inversion routine using Julia:
```bash
julia -t auto forceIBI/grinter_parallel.jl
```
This will perform the coarse-grained force field optimization.

### 4. Analyze and Visualize Results
To plot and analyze the results, run:
```bash
python3 forceIBI/plot_results.py
```
The generated plots and data will be saved automatically in the `./outputs/` directory.

---

## 📁 Folder Structure
```
project-root/
├── forceIBI/                # Core scripts and modules
│   ├── gr_borgis.jl         # Contains the main function for Borgis Formula implementation
│   ├── grinter_parallel.jl  # Contains the main potential convergence loop
│   ├── utils.jl             # Auxiliary functions for the inversion loop
|   ├── ...
|   ├── *.py                 # All python files used for pre- and post-processing
│   └── README.md            # This file
|
├── inputs/                  
│   ├── configs/             # Where configurations should be stored
│   ├── configs_bin/         # Auto-generated directory with binary-converted input files
│   └── params.json          # Params for the inversion
|
└── outputs/                 # Output directory (auto-generated)
```