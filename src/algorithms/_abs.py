from abc import ABC, abstractmethod
from typing import overload, TypeVar, Annotated, Literal

import numpy as np
import numpy.typing as npt

from ..element import Element

DType = TypeVar("DType", bound=np.generic)

Array2D = Annotated[npt.NDArray[DType], Literal["N", "N"]]
Array1D = Annotated[npt.NDArray[DType], Literal["N"]]


class AbsSolver(ABC):

    def __init__(self, data: Array2D[np.int64], start_point: int):
        self.data = data
        self.start_point = start_point

    @abstractmethod
    def solve(self) -> list[Element]:
        raise NotImplementedError('Method solve must be implemented in subclass')

    def __len__(self):
        return len(self.data)

    @overload
    def __getitem__(self, item: int) -> Array1D[np.int64]: ...
    @overload
    def __getitem__(self, item: (int, int)) -> np.int64: ...

    def __getitem__(self, item):
        return self.data[item]
