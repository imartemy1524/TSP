import os
import pathlib
from ctypes import *
from datetime import datetime

import networkx as nx
import numpy as np

from src.reader import read_data

out_path = str(pathlib.Path(__file__).parent / 'dll.so')
try:
    os.system(f" gcc ./dll.c -shared -o {out_path}")

    lib_c = cdll.LoadLibrary(out_path)
    lib_c.tsp_branch.argtypes = [c_int, c_void_p]
    lib_c.tsp_branch.restype = c_int


    def tsp_branch(n, py_arr, lib):
        if n < 2:
            return {}
        flatten_arr = list(np.concatenate(py_arr).flat)
        l = [-1] * (n * n - len(flatten_arr))
        flatten_arr = (flatten_arr + l)[:n * n:]
        int_arr = (c_int * (n * n))(*flatten_arr)
        res = lib.tsp_branch(n, byref(int_arr))
        if res > 0:
            l = list(int_arr)[:res:]
            return {'len': l.pop(0), 'steps': l.pop(0), 'path': l}
        else:
            return {}


    INF = -1

    input_matrix, best_solution = read_data(pathlib.Path(__file__).parent.parent / 'datasets/bays29.txt')
    for i in range(len(input_matrix)):
        input_matrix[i][i] = -1
    n = input_matrix.shape[0]

    start_time = datetime.now()
    res1 = tsp_branch(n, input_matrix, lib_c)
    print(datetime.now() - start_time)
    print(f"{best_solution=}")
    print('min_len =', res1['len'], ', steps =', res1['steps'], res1['path'])

    if 'path' in res1:
        d1 = []
        for i, v in enumerate(res1['path']):
            d1.append([res1['path'][i - 1], res1['path'][i]])
finally:
    os.unlink(out_path)