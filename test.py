import pathlib
import time

from z3 import *

from src.reader import read_data

data, sol = read_data(pathlib.Path(__file__).parent.parent / "datasets/bays29.txt")


def exactly_one_of(*args):
    g = [Or(*args)]
    for i in range(len(args)):
        for j in range(i + 1, len(args)):
            g.append(Or(Not(args[i]), Not(args[j])))
    return g


solver = Solver()
variables = []
variables_bools = []
for n in range(len(data)):
    from_n = Int(f"from_{n}")
    bools = []
    for m in range(len(data)):
        if n == m: bools.append(None); continue
        chosen = Bool(f"chosen_{n}_{m}")
        bools.append(chosen)
        solver.add(Implies(chosen, from_n == data[n, m]))
    variables.append(from_n)
    variables_bools.append(bools)
    print(n)
for i in range(len(data)):
    solver.add(exactly_one_of(*(variables_bools[j][i] for j in range(len(data)) if i != j)))
    solver.add(exactly_one_of(*(variables_bools[i][j] for j in range(len(data)) if i != j)))

for n in range(len(data)):
    for m in range(len(data)):
        if n == m: continue
        solver.add(Implies(variables_bools[n][m], Not(variables_bools[m][n])))
solver.add(Sum(variables) <= sol+1)
t = time.time()
print("Started")
print(
    solver.check()
)
answer = (
    solver.model()
)
print(time.time() - t)
ans = 0
for variable in answer:

    if is_true(answer[variable]):
        v = variable.name()
        from_, to_ = map(int, v.split('_')[1:])
        print("From", from_, "to", to_, "cost", data[from_, to_])
        ans += data[from_, to_]
print(ans, sol)

