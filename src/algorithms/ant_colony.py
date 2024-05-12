import typing
from multiprocessing import Manager

import numpy
import numpy as np
from numba import njit, prange

from ._abs import AbsSolver
from ..element import Element

T = typing.TypeVar('T')


class AntColony(AbsSolver):
    def __init__(self, data: np.ndarray, start_point: int):
        super().__init__(data, start_point)
        self.alpha = 1.0  # Weight for pheromone level
        self.beta = 2.0  # Weight for heuristic information (distance)
        self.rho = 0.5  # Pheromone evaporation rate
        self.Q = 100  # Pheromone deposit constant
        self.iterations = 1000
        self.repeat = 35

        # Initialize pheromone levels to a small constant value
        self.pheromones = np.ones((len(data), len(data))) * 0.1

        self.weights = self.build_weights()
        self.randoms = [numpy.random.RandomState(i + (start_point + 5) ** 5) for i in range(self.iterations)]

    def build_weights(self):
        weights = np.zeros_like(self.data, dtype=float)
        for i in range(len(self.data)):
            for j in range(len(self.data)):
                if self.data[i, j] != 0:
                    # Combine pheromone and heuristic information
                    weights[i, j] = (self.pheromones[i, j] ** self.alpha) * \
                                    ((1 / self.data[i, j]) ** self.beta)
                else:
                    weights[i, j] = 0
        return weights

    def solve(self):
        best_path = ()
        manager = Manager()
        return_dict = [None] * self.iterations
        for iteration in range(self.repeat):
            def run_ant_for_ant(ant, return_dict, count):
                random_numbers = self.randoms[ant].random_sample(len(self.data) * count)
                data = run_ants(count, self.data, self.weights, random_numbers)
                data = [(list(path), path_length) for path, path_length in data]
                # for i in range(len(data)):
                #     intarray = data[i][0]
                #     data[i] = ([Element(to=intarray[j], count=self.data[intarray[j-1], intarray[j]]) for j in range(len(intarray))], data[i][1])
                return_dict[ant:ant + count] = data

            run_ant_for_ant(0, return_dict, self.iterations)

            # do some code here
            for number in range(len(return_dict)):
                if not best_path or return_dict[number][1] < best_path[1]:
                    best_path = return_dict[number]
            # for ant in range(1000):
            #     path, path_length = self.run_ant()
            #     all_paths[ant] = (path, path_length)
            #     if not best_path or path_length < best_path[1]:
            #         best_path = (path, path_length)
            # Update pheromone based on paths
            self.update_pheromone(return_dict)

            # Recalculate weights after updating pheromone
            self.weights = self.build_weights()
            print(
                f"Done {iteration + 1}/{self.repeat} ({(iteration + 1) / self.repeat:0.0%}) iteration: best {best_path[1]}",
                end='\r'
            )
        # Return the best path found
        return [Element(best_path[0][i], self.data[best_path[0][i - 1], best_path[0][i]]) for i in
                range(len(best_path[0]))]

    def update_pheromone(self, all_paths):
        # Evaporate pheromone
        self.pheromones *= (1 - self.rho)

        # Deposit pheromone based on paths
        for path, path_length in all_paths:
            pheromone_to_add = self.Q / path_length
            for index, element in enumerate(path):
                from_node = path[index - 1]
                to_node = element
                self.pheromones[from_node, to_node] += pheromone_to_add

    def _get_max(self, from_: int, visited: list[bool]) -> int:
        max_index = -1
        max_value = -float('inf')
        for i in range(len(self.pheromones)):
            if not visited[i] and self.pheromones[from_, i] > max_value:
                max_value = self.pheromones[from_, i]
                max_index = i
        return max_index


@njit
def run_ant(data: numpy.ndarray, weights: numpy.ndarray, random_numbers):
    path = numpy.array([-1] * len(data), dtype='int32')
    visited_nodes = np.array([False] * len(data), dtype='bool')
    start_point = int(random_numbers[0] * (len(path)))  # [0, len(path)-1]  #self.start_point
    current_node = start_point
    path_length = 0

    visited_nodes[current_node] = True

    for i in range(len(data) - 1):
        next_node = _next_to_visit(weights, current_node, visited_nodes, random_numbers[i + 1])
        path[i] = next_node
        path_length += data[current_node, next_node]
        visited_nodes[next_node] = True
        current_node = next_node
    path_length += data[current_node, start_point]
    path[-1] = start_point
    return path, path_length


@njit
def _next_to_visit(wights: np.ndarray, from_: int, visited_nodes: list[bool], random_number) -> int:
    probabilities = wights[from_].copy()
    c = 0
    for i in range(len(visited_nodes)):
        if visited_nodes[i]:
            probabilities[i] = 0
            c += 1
    if c == len(visited_nodes):
        raise ValueError("All nodes are visited")
    total_probability = sum(probabilities)
    choice = random_choice(probabilities, total_probability, random_number)

    return choice


@njit
def random_choice(probabilities: 'list[T]', total_probability: float, random_number: float) -> T:
    r = random_number * total_probability
    cumulative = 0
    for i, probability in enumerate(probabilities):
        cumulative += probability
        if cumulative >= r:
            return i
    raise ValueError("Random value was incorrect")


@njit(parallel=True)
def run_ants(count: int, data: numpy.ndarray, weights: numpy.ndarray, random_numbers: numpy.ndarray):
    paths_lengths = numpy.array([0] * count, dtype='int64')
    paths = numpy.array([[0] * len(data)] * count, dtype='int32')
    for i in prange(count):
        path, path_length = run_ant(data, weights, random_numbers[i * len(data):])
        paths[i] = path
        paths_lengths[i] = path_length
        # paths.append((path, path_length))
    return [(path, path_length) for path, path_length in zip(paths, paths_lengths)]
