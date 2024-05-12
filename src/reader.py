
import numpy as np


# data2 - http://comopt.ifi.uni-heidelberg.de/software/TSPLIB95/atsp/
# data - https://people.sc.fsu.edu/~jburkardt/datasets/tsp/tsp.html
# minimal length is 937
def format_data(data: list):
    if not all(len(i) == 3 for i in data): return data
    answer = [[0] * len(data) for _ in range(len(data))]
    for i in range(len(data)):
        for j in range(i, len(data)):
            answer[i][j] = answer[j][i] = round(
                np.sqrt((data[i][1] - data[j][1]) ** 2 + (data[i][2] - data[j][2]) ** 2)
            )
    return answer


def read_data(file_path: str):
    with open(file_path, 'r') as file:
        data = file.readlines()
    solution = int(data[-1].split(': ')[-1])
    data = [list(map(int, filter(None, line.strip().split(' ')))) for line in data[:-2]]
    data = format_data(data)
    for i in range(len(data[0])):
        data[i][i] = 0
    return np.array(data), solution
