import random

import numpy as np


def nearest_neighbour_solution(dist_matrix: 'np.ndarray') -> list:
    """
    Computes the initial solution (nearest neighbour strategy)
    """
    assert len(dist_matrix.shape) == 2, "Matrix should be 2D"
    assert dist_matrix.shape[0] == dist_matrix.shape[1], "Matrix should be square"
    node = random.randrange(dist_matrix.shape[0])
    result = [node]

    nodes_to_visit = list(range(dist_matrix.shape[0]))
    nodes_to_visit.remove(node)

    while nodes_to_visit:
        nearest_node = min([(dist_matrix[node, j], j) for j in nodes_to_visit], key=lambda x: x[0])
        node = nearest_node[1]
        nodes_to_visit.remove(node)
        result.append(node)

    return result
