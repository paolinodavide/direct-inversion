How to use:

0) Adapt the number of Threads in the cpp source for the parallel computation
 (recommended 8).

1) Compile the cpp file to create the library.

$ bash Compile_gr_pair_parallel.sh

2) Compile the Cython libraries.

$ python3 setup.py build_ext --inplace

3) Edit the config file and create the module.

$ python3 make_config_iter.py

4) Choose the appropriate subset of snapshots in list_wt.dat.
Once done, create the json file for positions.

$ python3 gen_pos.py

5) The code is ready to run.

$ python3 grinter_parallel.py

6) Parse the result files.

$ python3 convert_json.py
