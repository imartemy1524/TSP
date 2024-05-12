import pathlib
import time
import unittest
import warnings

import numpy as np

from src.algorithms import LinKernighan, Opt1Improver, SimulatedAnnealing, AntColony
from src.element import Element
from src.reader import read_data

FILES = ['bays29.txt', 'data5.txt', 'pr124.txt', 'data7.txt', 'data.txt', 'data1.txt', 'data2.txt', 'data3.txt',
         'data6.txt']

COUNT = 10


def verify_path(path: list, data: np.ndarray):
    assert len(path) == len(set(path)), f"Path contains duplicates: {path}"
    all_nodes = set(range(data.shape[0]))
    sum_ = 0
    for i in range(len(path)):
        sum_ += data[path[i - 1], path[i]]
        all_nodes.remove(path[i])
    assert not all_nodes, f"Not all nodes were visited: {all_nodes}"
    return sum_


class MyTestCase(unittest.TestCase):
    @property
    def _each_file(self):
        for i in FILES:
            data, solution = read_data(str(pathlib.Path(__file__).parent / "datasets" / i))
            yield data, solution, i

    def test_any_colony(self):
        self._test_it(AntColony, "Ant colony")

    def test_sumulated_annealing(self):
        self._test_it(SimulatedAnnealing, "Simulated Annealing")

    def test_lin_kernighan(self):
        self._test_it(LinKernighan, "Lin-Kernighan", False, True)

    def test_1opt(self):
        self._test_it(Opt1Improver, "1-opt", False)

    def _test_it(self, cls, name: str, try_to_improve: bool = True, only_symmetric: bool = False):
        print("Running ", name + "...", end="\n\n")
        bests = []
        for data, solution, fname in self._each_file:
            if only_symmetric and not (data == data.transpose()).all():
                bests.append((0, 0, fname))
                warnings.warn(f"Data is not symmetric, thus, it wouldn't work with {name}")
                continue
            answers = []
            for i in range(COUNT):
                start_time = time.time()
                method = cls(data, i)
                answer = method.solve()
                assert len(answer) == len(set(i.to for i in answer)), f"Invalid answer: {answer}"
                best = sum(i.count for i in answer)
                assert verify_path([i.to for i in answer], data) == best, "Method returned invalid path"
                if try_to_improve:
                    answer2 = self._improve_answer(data, answer)
                    best2 = sum(i.count for i in answer2)
                    assert verify_path([i.to for i in answer2], data) == best2, "Method returned invalid path"
                else:
                    best2 = best
                if best2 < solution:
                    print(f"\n\n\n\n\nWE FOUND BETTER SOLITION!!!!!\n\n\n\n {[i.to for i in answer2]}\n\n\n\n")
                delta_time = time.time() - start_time
                print(
                    f"[{fname:11s}] Best: {str(solution):7s}| my: {str(best2):7s} ({str(best):7s} original) | {best2 / solution:0.3%} | {delta_time:0.2f}s\t\t\t")
                answers.append(best2)
                # if best2 == solution: break
            bests.append((min(answers), min(answers) / solution, fname))
        print(f"\n\n {name} result\n\tAverage: {sum(i[1] for i in bests) / len(bests)}")
        print(f'\n', *(i[1] for i in bests), sep='\n')
        with open('/shared/Projects/C++/HafmanAlgorithm/python/out/output.txt', 'a', encoding='utf-8') as file:
            file.write(name)
            file.write('\n')
            for answer, percent, fname in bests:
                file.write(
                    f"[{fname:11s}] {f'{percent:0.3%}':10s} {str(answer):10s} (best: {round(answer / percent) if percent else 0})\n")

    @staticmethod
    def _improve_answer(data: np.ndarray, answer: list[Element]):
        if (data == data.transpose()).all():
            improver = LinKernighan(data, 0, [i.to for i in answer])
        else:
            # if data is not symmetric, we can't use LinKernighan (opt-k improver) - so we use opt-1
            improver = Opt1Improver(data, 0, [i.to for i in answer])
        ans = improver.solve()
        assert len(ans) == len(answer), f"Invalid answer: {answer}, {ans}"
        return ans


if __name__ == '__main__':
    unittest.main()
