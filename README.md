# Direct Boltzmann inversion method from particle configurations at arbitrary state points
This is the code repository for the article ["Direct Boltzmann inversion method from particle configurations at arbitrary state points"](https://doi.org/10.48550/arXiv.2603.12081) by Coquand, Paolino, Berthier.
This is the code utilised to produce data and figures for the article, while a more recent version of the code is available at [https://github.com/paolinodavide/direct-inversion.git](https://github.com/paolinodavide/direct-inversion.git).

This projects builds a tool to reconstruct the interaction between particles starting from observations of the system. The main computation loop is written in Julia, while there are some Python scripts for pre- and post-processing. 
 
---
---
# ⚙️ Usage Guide
## Quick Start for Main Inversion
If your target directory is already fully prepared with configurations and a `params.json` file at `<YOUR_DIR>/inputs/params.json`, you can launch the inversion directly from the project root:

```bash
julia -t <NUM_THREADS or auto> forceIBI/grinter_parallel.jl -d <YOUR_DIR>
```

## Full Workflow
### 1. Prepare Input Configurations
Place your configuration files in the `<YOUR_DIR>/inputs/configs/` directory.

If you don't have any configurations ready, you can use the one provided in `direct-inversion/examples/` as a template. 
Run the following command to copy the example configuration to your working directory:
```bash
python3 format_data.py -i examples/lj_92_2.dat -d <YOUR_DIR>
```
This will create a `inputs/` folder in your specified directory with the necessary input structure.

### 2. Generate Targets for Inversion
Run the following scripts from the project root, providing your data directory with the `-d` flag:

1. Minimize the input data:
    ```bash
    python3 direct-inversion/min_d.py -d <YOUR_DIR>
    ```

2. Create histograms for the radial distribution function (RDF):
    ```bash
    python3 direct-inversion/gr_histo.py -d <YOUR_DIR> --dr <BIN_WIDTH>
    ```

3. Generate weighted target data:
    ```bash
    python3 direct-inversion/weighted_gen_tgt.py -d <YOUR_DIR>
    ```

4. Create a dummy `.json` parameter file:
    ```bash
    python3 direct-inversion/init_dummy_json.py -d <YOUR_DIR>
    ```
5. Edit the generated `params.json` file to set the correct parameters for your inversion. If you don't know how, you can read the original article https://doi.org/10.48550/arXiv.2603.12081. 

### 3. Run the Main Computation
The inversion routine now looks for `inputs/params.json` relative to your provided directory:
```bash
julia --project=direct-inversion -t <NUM_THREADS or auto> direct-inversion/grinter_parallel.jl -d <YOUR_DIR>
```

### 4. Analyze and Visualize Results
To plot and analyze the results, run:
```bash
python3 direct-inversion/plot_results.py -d <YOUR_DIR>
```
or 
```bash
python3 direct-inversion/convergence_plot.py -d <YOUR_DIR>
```

The generated plots and data will be saved automatically in the `./outputs/` directory.

---

## 📁 Directory Architecture

The pipeline is designed to be highly flexible. You can keep the repository codebase in one location and run the inversion on **any external working directory** on your machine using the `-d` or `--directory` flags. 

### 1. Source Code Repository
This is the structure of the provided codebase:
```text
project-root/
├── gr_borgis.jl         # Main function for Borgis Formula implementation
├── grinter_parallel.jl  # Main potential inversion loop
├── utils.jl             # Auxiliary functions for the inversion loop
├── ...                  # Additional Julia modules
├── *.py                 # Python files for pre- and post-processing
├── examples/*           # Example input files for testing and demonstration
│
└── README.md            # This documentation
```

### 2. User Working Directory (```YOUR_DIR```)
You must prepare a specific working directory for your data. When running the scripts, point them to this directory. The structure should look like this:
```text
<YOUR_DIR>/
├── inputs/                  
│   ├── configs/             # ⚠️ REQUIRED: Place your initial configuration files here
│   ├── configs_bin/         # (Auto-generated) Binary-converted input files
│   └── params.json          # (Auto-generated) Parameters for the inversion, to be edited by you
└── outputs/                 # (Auto-generated) Results, output data, and plots
```
