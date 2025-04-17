g++ -o gr_pair_parallel.o -c gr_pair_parallel.cpp -fopenmp -std=gnu++11 -Wall -fPIC `python -m pybind11 --includes` -ljsoncpp -O3
g++ -o gr_pair_parallel.so -shared gr_pair_parallel.o -ljsoncpp -fopenmp -O3
