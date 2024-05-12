import os
import pathlib
import warnings
from copy import deepcopy
from ctypes import CDLL, c_int, c_longlong

import numpy as np
from numpy.ctypeslib import ndpointer

from .._abs import AbsSolver
from ...element import Element


class Path:
    """
    Class, describes a path
    """

    def __init__(self, path: 'list[int]', cost: float = 0):
        self.path = path
        self.cost = cost

    def __hash__(self):
        return hash(tuple(self.path))

    def __eq__(self, other):
        return self.path == other.path


class Tour:
    """
    Class to represent a tour in LKH.
    """

    def __init__(self, tour: 'list[int]', data: 'np.ndarray'):
        self.tour_array = tour
        self.tour_indexes = {node: i for i, node in enumerate(tour)}
        self.length = self._cost(data)
        self.edges = set(
            (self.tour_array[i - 1], self.tour_array[i]) if self.tour_array[i - 1] < self.tour_array[i] else (
                self.tour_array[i], self.tour_array[i - 1])
            for i in range(self.size)
        )

    @property
    def size(self):
        return len(self.tour_array)

    def _cost(self, data: 'np.ndarray'):
        cost = sum(data[self.tour_array[i - 1], self.tour_array[i]] for i in range(self.size))
        return cost

    def __getitem__(self, item) -> int:
        return self.tour_array[item]

    def __contains__(self, edge):
        return edge in self.edges or edge[:-1] in self.edges

    def index(self, item):
        """
        Return the index of a node in a tour.
        """
        return self.tour_indexes[item]

    def around(self, node):
        """
        Return the predecessor and successor of the current node, given by index.
        Return: tuple of (predecessor, successor)
        """
        index = self.index(node)

        prev_child = index - 1
        next_child = index + 1 if index + 1 < self.size else 0
        return self[prev_child], self[next_child]

    def generate(self, broken: 'set[(int,int)]', joined: 'set[(int,int)]'):
        """
        Generate a tour with the exclusions (broken) and inclusions (joined).
        """
        edges = (self.edges - broken) | joined

        # If we do not have enough edges, we cannot form a tour -- should not happen within LKH
        if len(edges) < self.size:
            warnings.warn("Not enough edges to form a tour.")
            return False, []
        successors = self._build_successors(edges)
        if len(successors) < self.size:
            warnings.warn("Not enough successors to form a tour.")
            return False, []

        succ = successors[0]
        new_tour = [0]
        visited = set(new_tour)

        # If we already encountered a node it means we have a loop
        while succ not in visited:
            visited.add(succ)
            new_tour.append(succ)
            succ = successors[succ]

        # If all nodes are visited without a loop - onw have a tour
        return len(new_tour) == self.size, new_tour

    def _build_successors(self, edges: 'set[(int,int)]'):
        """
        Build a dictionary of successors for each node in the tour.
        :param edges: set of edges in the tour
        :returns: dictionary of successors
        """
        successors = {}
        node = 0

        # Build the list of successors
        while edges:
            for from_, to_ in edges:
                if from_ == node:
                    successors[node] = to_
                    node = to_
                    break
                elif to_ == node:
                    successors[node] = from_
                    node = from_
                    break
            edges.remove((from_, to_))
        return successors


class LinKernighanOld(AbsSolver):
    """
    Lin-Kernighan algorithm implementation.
    Works only for symmetric matrices.
    """

    def __init__(self, data: np.ndarray, start_point: int, start_path: 'list[int]' = None):
        super().__init__(data, start_point)

        self.paths: set[Path] = set()
        self.cost = 9999999999
        self.heuristic_path: 'list[int]' = start_path if start_path is not None else [7, 26, 23, 15, 18, 14, 3, 9, 19,
                                                                                      1, 20, 0, 27, 5, 11, 8, 4, 25, 28,
                                                                                      2, 12, 24, 6, 22, 10, 21, 13, 17,
                                                                                      16]  # nearest_neighbour_solution(data)
        self.neighbours = {}
        for i in self.heuristic_path:
            el = self[i]
            self.neighbours[i] = [
                j for j, dist in enumerate(el) if dist > 0
            ]

    def solve(self):
        improve = True
        while improve:
            improve = self.improve()
            # Paths always begin at 0 so this should manage to find duplicate
            # solutions
            self.paths.add(Path(self.heuristic_path, self.cost))
        answer = min(self.paths, key=lambda x: x.cost).path
        return [Element(i, self[i, j]) for i, j in zip(answer, answer[1:] + [answer[0]])]

    def closest(self, t2i, tour, gain, broken, joined) -> 'list[tuple[int, tuple[int, int]]]':
        """
        Find the closest neighbours of a node ordered by potential gain. As a side effect, also compute the partial improvement of joining a node.

        :param t2i: node to relink from
        :param tour: current tour to optimise
        :param gain: current gain
        :param broken: set of edges to remove (X)
        :param joined: set of edges to join (Y)
        :return: sorted list of neighbours based on potential improvement with next omission
        """
        neighbours: 'dict[int, tuple[int, int]]' = {}

        # Create the neighbours of t_2i
        for node in self.neighbours[t2i]:
            yi = (t2i, node) if t2i < node else (node, t2i)
            Gi: int = gain - self[t2i, node]

            # Any new edge has to have a positive running sum, not be a broken
            # edge and not belong to the tour.
            if Gi <= 0 or yi in broken or yi in tour:
                continue

            for succ in tour.around(node):
                xi = (node, succ) if node < succ else (succ, node)

                # valid first thing in `x` so this should be sufficient
                #
                # Check that "x_i+1 exists"
                if xi not in broken and xi not in joined:
                    diff = self[node, succ] - self[t2i, node]
                    if node in neighbours and diff > neighbours[node][0]:
                        neighbours[node] = (diff, neighbours[node][1])
                    else:
                        neighbours[node] = diff, Gi

        # Sort the neighbours by potential gain
        return sorted(neighbours.items(), key=lambda x: x[1][0], reverse=True)

    def improve(self):
        """
        Start the LKH algorithm with the current tour.
        """
        tour = Tour(self.heuristic_path, self.data)

        # Find all valid 2-opt moves and try them
        for t1 in self.heuristic_path:
            around = tour.around(t1)

            for t2 in around:
                broken = {(t1, t2) if t1 < t2 else (t2, t1)}
                # Initial savings
                gain = self[t1, t2]

                close = self.closest(t2, tour, gain, broken, set())

                # Number of neighbours to try
                tries = 5

                for t3, (_, Gi) in close:
                    # Make sure that the new node is none of t_1's neighbours, so it does not belong to the tour.
                    if t3 in around:
                        continue

                    joined = {(t2, t3) if t2 < t3 else (t3, t2)}

                    # The positive Gi is taken care of by `closest()`

                    if self.choice_remove(tour, t1, t3, Gi, broken, joined):
                        # Return to Step 2, that is the initial loop
                        return True
                    # Else try the other options

                    tries -= 1
                    # Explored enough nodes, change t_2
                    if tries <= 0: break

        return False

    def choice_remove(self, tour: 'Tour', t1: int, last: int, gain: int, broken: 'set[tuple[int,int]]',
                      joined: 'set[tuple[int,int]]'):
        """
        Choose an edge to omit from the tour.

        :param tour: current tour to optimise
        :param t1: starting node for the current k-opt
        :param last: tail of the last edge added (t_2i-1)
        :param gain: current gain (Gi)
        :param broken: potential edges to remove (X)
        :param joined: potential edges to add (Y)
        :return: whether we found an improved tour
        """
        if len(broken) == 4:
            pred, succ = tour.around(last)

            # Give priority to the longest edge for x_4
            if self[pred, last] > self[succ, last]:
                around = [pred]
            else:
                around = [succ]
        else:
            around = tour.around(last)

        for t2i in around:
            xi = (last, t2i) if last < t2i else (t2i, last)
            # Gain at current iteration
            Gi = gain + self[last, t2i]

            # Verify that X and Y are disjoint, though I also need to check
            # that we are not including an x_i again for some reason.
            if xi not in joined and xi not in broken:
                added = deepcopy(joined)
                removed = deepcopy(broken)

                removed.add(xi)
                added.add((t2i, t1) if t2i < t1 else (t1, t2i))  # Try to relink the tour

                relink = Gi - self[t2i, t1]
                is_tour, new_tour = tour.generate(removed, added)

                # The current solution does not form a valid tour
                if not is_tour and len(added) > 2:
                    continue

                # Stop the search if we come back to the same solution
                if Path(new_tour) in self.paths:
                    return False

                # Save the current solution if the tour is better, we need
                # `is_tour` again in the case where we have a non-sequential
                # exchange with i = 2
                if is_tour and relink > 0:
                    self.heuristic_path = new_tour
                    self.cost -= relink

                    return True
                else:
                    # Pass on the newly "removed" edge but not the relink
                    choice = self.choice_add(tour, t1, t2i, Gi, removed, joined)

                    if len(broken) == 2 and choice:
                        return True
                    else:
                        # Single iteration for i > 2
                        return choice

        return False

    def choice_add(self, tour: Tour, t1: int, t2i: int, gain: int, broken: 'set[tuple[int,int]]',
                   joined: set[tuple[int, int]]):
        """
        Choose an edge to add to the new tour.
        :param tour: current tour to optimise
        :param t1: starting node for the current k-opt
        :param t2i: tail of the last edge removed (t_2i)
        :param gain: current gain (Gi)
        :param broken: potential edges to remove (X)
        :param joined: potential edges to add (Y)
        :returns: whether we found an improved tour
        """
        ordered = self.closest(t2i, tour, gain, broken, joined)

        if len(broken) == 2:
            # Check the five nearest neighbours when i = 2
            top = 5
        else:
            # Otherwise the closest only
            top = 1

        for node, (_, Gi) in ordered:
            yi = (t2i, node) if t2i < node else (node, t2i)
            added = deepcopy(joined)
            added.add(yi)

            # Stop at the first improving tour
            if self.choice_remove(tour, t1, node, Gi, broken, added):
                return True

            top -= 1
            # Tried enough options
            if top == 0:
                return False

        return False


fle = pathlib.Path(__file__).parent / "lin_kernighan_so.so"
fle_cpp = pathlib.Path(__file__).parent / "lin_kernighan.cpp"
if not fle.is_file():
    os.system("g++ -shared -fPIC -o " + str(fle) + " " + str(fle_cpp))
lib = CDLL(str(fle))
lib.solve.restype = c_int
lib.solve.argtypes = [
    c_int,
    ndpointer(c_longlong, flags="C_CONTIGUOUS"),
    ndpointer(c_longlong, flags="C_CONTIGUOUS"),
    ndpointer(c_longlong, flags="C_CONTIGUOUS")
]


class LinKernighan(AbsSolver):

    def __init__(self, data: np.ndarray, start_point: int, start_path: 'list[int]' = None):
        super().__init__(data, start_point)
        self.start_path = start_path

    def solve(self) -> list[Element]:
        n = self.data.shape[0]
        data = self.data.flatten()
        start_path = np.array(self.start_path or [-1], dtype=np.int64)
        result = np.empty(n, dtype=np.int64)
        answer = lib.solve(n, data, start_path, result)
        if answer == -1:
            raise ValueError("Error in c++ library")
        answer_list = result.tolist()
        return [Element(i, int(self.data[i, j])) for i, j in zip(answer_list, answer_list[1:] + [answer_list[0]])]
