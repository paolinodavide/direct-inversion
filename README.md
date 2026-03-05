# 🥊 Force IBI Project

This project implements Iterative Boltzmann Inversion (IBI) by exploiting the force formula proposed by Borgis et al. The code has been updated to use Julia as the main computation engine for better performance, while maintaining Python for pre/post-processing.

> ⚠️ **Note:** This repository is for personal use. Dependencies and environment setup are assumed to be manually managed.  
---
---
# ⚙️ Usage Guide
## Quick Start for Main Inversion

1. Ensure your parameter file is configured at `<YOUR_DIR>/inputs/params.json`.
2. Execute the following command from the project root:
    ```bash
    julia -t auto forceIBI/grinter_parallel.jl --directory <YOUR_DIR>
    ```

## Full Workflow
### 1. Prepare Input Configurations
Place your configuration files in the `<YOUR_DIR>/inputs/configs/` directory.

### 2. Generate Targets for Inversion
Run the following scripts from the project root, providing your data directory with the `-d` flag:

1. Minimize the input data:
    ```bash
    python3 forceIBI/min_d.py -d <YOUR_DIR>
    ```

2. Create histograms for the radial distribution function (RDF):
    ```bash
    python3 forceIBI/gr_histo.py -d <YOUR_DIR>
    ```

3. Generate weighted target data:
    ```bash
    python3 forceIBI/weighted_gen_tgt.py -d <YOUR_DIR>
    ```

4. Create a dummy `.json` parameter file:
    ```bash
    python3 init_dummy_json.py -d <YOUR_DIR>
    ```

### 3. Run the Main Computation
The inversion routine now looks for `inputs/params.json` relative to your provided directory:
```bash
julia --project=forceIBI -t auto forceIBI/grinter_parallel.jl --directory <YOUR_DIR>

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
│   ├── grinter_parallel.jl  # Contains the main potential inversion loop
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