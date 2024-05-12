# Different TSP solvers

This project contains different TSP solvers, including:

- Any colony optimization (ACO) [ant_colony.py](src%2Falgorithms%2Fant_colony.py)
- Simulated annealing (SA) [sim_anneling.py](src%2Falgorithms%2Fsim_anneling.py)
- Lin-Kernighan algorithm(LK) [lin_kernighan.cpp](src%2Falgorithms%2Fopt%2Flin_kernighan.cpp)
- OPT-1 algorithm (otimize the graph, tryting moving each node in `O(n^2)`) [opt1.py](src%2Falgorithms%2Fopt%2Fopt1.py)

For some other algorithms, like OPT-1 and Ant Colony, the code is written in Python, using numba NJit to speed up the
code. For the Lin-Kernighan algorithm, the code is written in C++, ported to python using `extern C` methods.

## Preview comparation

As you can see in [out/output.txt](out/output.txt) file statistics, the Lin-Kernighan algorithm is the fastest one and
usually gives the best results, BUT, its result highly depends on inital input.
The best results it gives, when passing as an input result of Ant colony Algorithm.

### TODO

- Add superants to the Ant Colony Algorithm.
- Try to improve `SA` on big amounts of data (it is super uneffictive on over 100+ nodes)
- Compare with existing libraries, like `ortools` or `networkx` (already in progress) and branch and bound method.



## How to run
To run the test (tested only on linux) do the following:

1. Create python venv:
```bash
python3 -m venv venv
```
2. Install the libriries:
```bash
pip install -r requirements.txt
```
3. Run the test:
```bash
python3 ./checker.py
```
