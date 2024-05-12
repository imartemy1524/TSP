import warnings

import numpy as np

from ._abs import AbsSolver
from .utils import nearest_neighbour_solution
from ..element import Element


class SimulatedAnnealing(AbsSolver):
    def __init__(self, data: np.ndarray, start_point: int):
        super().__init__(data, start_point)
        self.temp =self.data.sum() / len(self.data) / 5
        self.data_length = len(data[0])
        self.iterations = self.data_length * 1000
        self.multiplier = 0.999
        self._random = np.random.RandomState(start_point)
        self._current_way = nearest_neighbour_solution(data)
        self._current_length = self.weight(self._current_way)
        self.q = self.data_length

    def try_accept(self, data):
        """Checks if the state transition will execute."""
        candidate_weight = self.weight(data)
        if candidate_weight < self._current_length:
            self._current_length = candidate_weight
            self._current_way = data
            if candidate_weight < self.minimum[0]:
                self.minimum = (candidate_weight, data)

        else:
            if self.random > min(1, np.exp((self._current_length - candidate_weight) / self.temp)):
                self._current_length = candidate_weight
                self._current_way = data
            else:
                pass
    def weight(self, way):
        return sum(self.data[i, j] for i, j in zip(way, way[1:] + [way[0]]))

    def solve(self):

        self.minimum = (self._current_length, self._current_way)
        now_min = (self.minimum[0])
        for _ in range(self.iterations):
            candidate = self._current_way[::]
            l = self._random.randint(2, max(int(self.q - 1), 3))
            i = self._random.randint(0, self.data_length - l)

            candidate[i: (i + l)] = (candidate[i: (i + l)])[::-1]

            self.try_accept(candidate)
            self.temp /= self.multiplier
            # self.q *= 0.999
        if self.minimum[0] == now_min:
            warnings.warn("No better solution found")

        return [Element(i, self.data[i, j]) for i, j in
                zip(self.minimum[1], self.minimum[1][1:] + [self.minimum[1][0]])]

    @property
    def random(self):
        return self._random.random()

    @property
    def random_indexes(self):
        a = self._random.randint(0, len(self.data[0]))
        b = self._random.randint(0, len(self.data[0]) - 1)
        if b >= a:
            b += 1
        return a, b
