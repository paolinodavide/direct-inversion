# Force IBI Project

This project contains code for performing iterative Boltzmann inversion (IBI) or similar force field optimization techniques.
The code is organized into modular Python and Cython scripts, and is intended to be executed from within the `src/` directory, where a detailed `README.md` guides execution.

> **Note:** This repository is for personal use. No installation or dependency details are provided here.
> _You can edit this section to specify the scientific purpose or include environment setup instructions._

---

## 📁 Folder Structure
```
project-root/
├── src/ # Core scripts and modules
│ ├── *.py / *.pyx # Python & Cython code
│ ├── *.sh # Shell scripts for compiling/running
│ ├── README.txt # Usage instructions
│ └── utils/ # Utility scripts
│     └── *.py
├── Results/ # Automatically generated results
├── config_iter_/ # Input configuration folders
└── README.md # This file
```

---

## 📌 Notes

- The `Results/` folder will be automatically created in the root directory when the code is run.
- Any configuration folders like `config_iter_0`, `config_iter_1`, etc., **must be placed directly inside the master directory** (i.e., next to `src/`).
- For detailed instructions on running the code, see `src/README.txt`.

---
