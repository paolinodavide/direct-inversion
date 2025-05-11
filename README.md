# 🥊 Force IBI Project

This project implements Iterative Boltzmann Inversion (IBI) and related coarse-grained force field optimization techniques. The code is modular and optimized for Python with optional Cython acceleration.

> ⚠️ **Note:** This repository is for personal use. Dependencies and environment setup are assumed to be manually managed.  
---
## ⚙️ Usage Guide
> 🧵 All steps assume execution from within the `forceIBI/` directory.
### 1. Compile Cython Extensions
```bash
python3 setup.py build_ext --inplace
```
### 2. Set Number of Threads
Edit the file `forceIBI/gr_borgis.py` and set the number of cores (recommended: 8): 
```python
NUM_THREADS = 8  # Example
```
### 3. Prepare Configuration Module
Ensure that the system configuration files are correctly compiled in `.npy` format.
Edit your parameters file and generate the dictionary:
```bash
python3 make_config_iter.py
```
### 4. Run Main Computation
Execute the main routine:
```bash
python3 grinter_parallel.py
```
### 5. Generate and Parse Results
To organize and analyze the output:
```bash
bash Gen_res.sh
```
Results will be saved automatically in the `Results/` directory. 

---

## 📁 Folder Structure
```
project-root/
├── forceIBI/                # Core scripts and modules
│   ├── gr_borgis.py         # Contains the main function for Borgis Formula implementation
│   ├── grinter_parallel.py  # Contains the main potential convergence loop
│   ├── *.py                 # Python scripts
│   ├── *.pyx                # Cython code
│   ├── *.cpp                # (Deprecated) C++ scripts
│   ├── *.sh                 # Utility shell scripts
│   ├── README.md            # This file
│   └── utils/               # Utility modules
│       └── *.py
├── configs/                 # Input configuration folders (e.g., config_iter_0)
├── configs_npy/             # Serialized config data (NumPy)
└── Results/                 # Output directory (auto-generated)
```
---
## 📝 Notes
- Do **not** place configuration folders inside subdirectories — they must be located directly in the project root.
- The `gen_pos.py` script is now deprecated and no longer required.
- The `format_data.py` script now also converts files from `.dat` to `.npy`. 
- The `Results/` folder will be created automatically if it doesn't exist.
- Target RDFs and waiting time lists should be saved in the `project-root` directory.
- Cache is stored to improve performance; however, outdated caches might lead to ineffective runs. To address this, `grinter_parallel` includes a routine to reset the `forceIBI/__pycache__` directory at each run, clearing the cache. However, this may not always suffice. An additional run should resolve the issue in such cases.   