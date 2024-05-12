import sys

import numpy
import numpy as np
from numba import njit, typed

from src.algorithms._abs import AbsSolver
from src.algorithms.utils import nearest_neighbour_solution
from src.element import Element


class Opt1Improver(AbsSolver):
    """
    1-Opt algorithm improver.

    Goes through each edge and tries relinking it somewhere else.

    Complexity: O(n^2)
    """

    def __init__(self, data: np.ndarray, start_point: int, start_data: list[int] = None):
        super().__init__(data, start_point)
        self.optimized = start_data or nearest_neighbour_solution(data)

    def solve(self):
        self.optimized = optimize(self.data, typed.List(self.optimized))
        return [
            Element(self.optimized[i], self[self.optimized[i - 1], self.optimized[i]]) for i in
            range(len(self.optimized))
        ]

@njit
def current_size(data, optimized):
    sum_ = 0
    for i in range(len(optimized)):
        sum_ += data[optimized[i - 1], optimized[i]]
    return sum_


@njit
def optimize(data, optimized):
    for i in range(999999999999):
        # yield False, " ".join("Optimizing solution:", current_size(data, optimized), "#", i, '\r')
        for move_node_index in range(len(optimized)):
            optimized2 = can_improve(optimized, data, move_node_index)
            if optimized2 is not None:
                optimized = optimized2
                break
        else:
            break
    return optimized


@njit
def at(optimized, index: int):
    return optimized[index % len(optimized)]


@njit
def can_improve(optimized, data, move_node_index: int):
    element = optimized[move_node_index]
    # How much improvement would we get, removing removing node with index `move_node_index`
    gain = (data[at(optimized, move_node_index - 1), element]
            + data[element, at(optimized, move_node_index + 1)]
            - data[at(optimized, move_node_index - 1), at(optimized, move_node_index + 1)])
    for i in range(len(optimized)):
        a_index, b_index = ((i - 1) % len(optimized)), i
        if a_index == move_node_index or b_index == move_node_index: continue
        a, b = optimized[a_index], optimized[b_index]
        gain_curr = (
                data[a, element]
                + data[element, b]
                - data[a, b]
        )
        # Check if replacing edge move_node_index between edges a - b would be beneficial
        if gain_curr < gain:
            # Great, replace it!
            if move_node_index < a_index:
                optimized.insert(b_index if b_index != 0 else len(optimized), element)
                optimized.remove(element)
            else:
                optimized.remove(element)
                optimized.insert(b_index if b_index != 0 else len(optimized), element)
            return optimized
    return None
