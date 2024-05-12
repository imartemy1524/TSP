# Different TSP solvers


This project contains different TSP solvers, including:
- Any colony optimization (ACO) [ant_colony.py](src%2Falgorithms%2Fant_colony.py)
- Simulated annealing (SA) [sim_anneling.py](src%2Falgorithms%2Fsim_anneling.py)
- Lin-Kernighan algorithm(LK) [lin_kernighan.cpp](src%2Falgorithms%2Fopt%2Flin_kernighan.cpp)
- OPT-1 algorithm (otimize the graph, tryting moving each node in `O(n^2)`) [opt1.py](src%2Falgorithms%2Fopt%2Fopt1.py)


For some other algorithms, like OPT-1 and Ant Colony, the code is written in Python, using numba NJit to speed up the code. For the Lin-Kernighan algorithm, the code is written in C++, ported to python using `extern C` methods.

