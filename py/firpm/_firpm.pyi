# The private extension behind the firpm package; see firpm/__init__.py.
from typing import Callable, ClassVar, Sequence, Tuple, Union

import numpy as np
import numpy.typing as npt

__version__: str

class Strategy:
    UNIFORM: ClassVar[Strategy]
    SCALING: ClassVar[Strategy]
    AFP: ClassVar[Strategy]
    __members__: ClassVar[dict[str, Strategy]]
    def __init__(self, value: int) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __hash__(self) -> int: ...
    def __index__(self) -> int: ...
    def __int__(self) -> int: ...
    @property
    def name(self) -> str: ...
    @property
    def value(self) -> int: ...

_BandValue = Union[float, Callable[[float], float]]
_BandSpec = Tuple[float, float, _BandValue, _BandValue]

def firpm(n: int, f: Sequence[float], a: Sequence[float], w: Sequence[float],
          type: str, eps: float, nmax: int, strategy: Strategy, depth: int,
          rstrategy: Strategy, prec: int) -> npt.NDArray[np.float64]: ...
def firpm_bands(n: int, bands: Sequence[_BandSpec], nyquist: float,
                eps: float, nmax: int, strategy: Strategy, depth: int,
                rstrategy: Strategy, prec: int) -> npt.NDArray[np.float64]: ...

# Only in builds with MPFR.
def firpm_mp(n: int, f: Sequence[float], a: Sequence[float],
             w: Sequence[float], type: str, eps: float, nmax: int,
             strategy: Strategy, depth: int, rstrategy: Strategy,
             prec: int) -> npt.NDArray[np.float64]: ...
def firpm_bands_mp(n: int, bands: Sequence[_BandSpec], nyquist: float,
                   eps: float, nmax: int, strategy: Strategy, depth: int,
                   rstrategy: Strategy, prec: int) -> npt.NDArray[np.float64]: ...
