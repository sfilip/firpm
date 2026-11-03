"""Robust Parks-McClellan FIR filter design.

There are two design routines, which differ only in how a specification is
written: :func:`remez` takes the arguments of ``scipy.signal.remez``, and
:func:`firpm` those of MATLAB's ``firpm``. Both also accept a sequence of
:class:`Band`, whose amplitude and weight may vary across the band.
"""

from __future__ import annotations

import math
from typing import (Any, Callable, List, Literal, NamedTuple, Optional,
                    Sequence, Tuple, Union, cast)

import numpy as np
import numpy.typing as npt

from . import _firpm
from ._firpm import Strategy, __version__

__all__ = ["Band", "Strategy", "firpm", "remez", "__version__"]

BandValue = Union[float, Callable[[float], float]]
BandLike = Union["Band", Tuple[Any, ...]]


class Band(NamedTuple):
    """One band of a :func:`remez` or :func:`firpm` specification.

    ``f0`` and ``f1`` are the band edges, in the frequency units of the
    routine: units of ``fs`` for :func:`remez`, and normalized so that 1 is
    the Nyquist frequency for :func:`firpm`. ``amplitude`` (the desired
    amplitude) and ``weight`` are each a number or a callable. Callables
    are called with one frequency in ``[f0, f1]``, in the same units, and
    return a number.
    """
    f0: float
    f1: float
    amplitude: BandValue
    weight: BandValue = 1.0


_OPTIONS = """\
    eps : float
        Convergence threshold: roughly, the relative accuracy of the final
        minimax error.
    nmax : int
        Degree of the local approximations used to find the error extrema.
    strategy : Strategy
        How the initial reference is chosen: ``UNIFORM``, ``SCALING``
        (reference scaling, from a shorter filter) or ``AFP``.
    depth : int
        With ``SCALING``, the number of scaling levels. It must be at least
        1 for ``SCALING`` to have any effect.
    rstrategy : Strategy
        With ``SCALING``, how the shortest filter's reference is chosen:
        ``UNIFORM`` or ``AFP``.
    precision : {"double", "mpfr"}
        The arithmetic used. ``"mpfr"`` needs a build with MPFR; the
        published wheels are built without it.
    prec : int
        With ``precision="mpfr"``, the precision in bits.
"""

_RETURNS = """\
    Raises
    ------
    ValueError
        If the specification is invalid (e.g. amplitudes not matching the
        band edges, or band edges out of order), the design failed, or
        ``precision`` is not available; as ``scipy.signal.remez`` does.

    Warns
    -----
    RuntimeWarning
        If the exchange algorithm did not converge.

    Notes
    -----
    Bands with only numeric amplitudes and weights give exactly the design
    of the equivalent band-edge specification. Callables are called a few
    thousand times per design, one frequency at a time, so keep them cheap.
    An exception raised inside a callable is re-raised.
"""


class _Options(NamedTuple):
    eps: float
    nmax: int
    strategy: Strategy
    depth: int
    rstrategy: Strategy
    prec: int


def _routines(precision: str) -> Tuple[Any, Any]:
    if precision == "double":
        return _firpm.firpm, _firpm.firpm_bands
    if precision == "mpfr":
        if not hasattr(_firpm, "firpm_mp"):
            raise ValueError("precision='mpfr' needs firpm built with MPFR; "
                             "the published wheels are built without it")
        return _firpm.firpm_mp, _firpm.firpm_bands_mp
    raise ValueError(
        f"precision must be 'double' or 'mpfr', not {precision!r}")


def _design_bands(n: int, bands: Sequence[BandLike], nyquist: float,
                  precision: str,
                  options: _Options) -> npt.NDArray[np.float64]:
    """Designs from Bands with edges in [0, nyquist]."""
    flat, banded = _routines(precision)
    specs = [Band(*band) for band in bands]
    if any(callable(b.amplitude) or callable(b.weight) for b in specs):
        return cast(npt.NDArray[np.float64], banded(
            n, [(b.f0, b.f1, b.amplitude, b.weight) for b in specs],
            nyquist, *options))
    f: List[float] = []
    a: List[float] = []
    for b in specs:
        f += [b.f0 / nyquist, b.f1 / nyquist]
        a += [cast(float, b.amplitude)] * 2
    w = [cast(float, b.weight) for b in specs]
    return cast(npt.NDArray[np.float64],
                flat(n, f, a, w, "bandpass", *options))


def remez(numtaps: int,
          bands: Union[Sequence[float], Sequence[BandLike]],
          desired: Optional[Sequence[float]] = None, *,
          weight: Optional[Sequence[float]] = None,
          type: Literal["bandpass", "differentiator", "hilbert"] = "bandpass",
          fs: Optional[float] = None,
          eps: float = 0.01,
          nmax: int = 4,
          strategy: Strategy = Strategy.UNIFORM,
          depth: int = 0,
          rstrategy: Strategy = Strategy.UNIFORM,
          precision: Literal["double", "mpfr"] = "double",
          prec: int = 165) -> npt.NDArray[np.float64]:
    """Design a linear-phase FIR filter, as ``scipy.signal.remez`` does.

    ``remez(numtaps, bands, desired, weight=..., type=..., fs=...)`` takes
    the specification of ``scipy.signal.remez``. ``remez(numtaps, bands)``,
    with a sequence of :class:`Band` (edges in units of ``fs``), allows
    amplitudes and weights that vary across a band.

    Parameters
    ----------
    numtaps : int
        The number of taps. With ``type="bandpass"``, odd ``numtaps`` gives
        a type I filter and even ``numtaps`` a type II filter; an even
        ``numtaps`` whose last band has nonzero amplitude at Nyquist (which
        a type II filter cannot have) gives the type I filter with
        ``numtaps + 1`` taps. Otherwise, odd ``numtaps`` gives type III and
        even ``numtaps`` type IV.
    bands : sequence of float, or sequence of Band
        Either band edges, a monotonic sequence of pairs from 0 to
        ``fs / 2``, or the bands (see :class:`Band`; plain ``(f0, f1,
        amplitude[, weight])`` tuples also work).
    desired : sequence of float, optional
        With band edges: the desired amplitude in each band.
    weight : sequence of float, optional
        With band edges: the weight of each band (1 by default).
    type : {"bandpass", "differentiator", "hilbert"}
        ``"bandpass"`` gives a symmetric filter (types I and II) with a flat
        amplitude in each band. The others give antisymmetric filters
        (types III and IV): ``"differentiator"`` has an amplitude of
        ``desired * f / fs`` in each band, and ``"hilbert"`` a flat
        amplitude. They need band edges, not :class:`Band`.
    fs : float, optional
        The sampling frequency (1 by default, so that Nyquist is 0.5).
OPTIONS
    Returns
    -------
    numpy.ndarray
        The ``numtaps`` taps (``numtaps + 1`` after promotion to type I).

RETURNS
    For differentiators, ``scipy.signal.remez`` divides the weight by the
    frequency in every band with a nonzero slope; here, as in MATLAB, only
    the first band's weight is. The two agree for the usual differentiator,
    whose only band with a nonzero slope is the first.
    """
    if type not in ("bandpass", "differentiator", "hilbert"):
        raise ValueError("type must be 'bandpass', 'differentiator' or "
                         f"'hilbert', not {type!r}")
    if numtaps < 1:
        raise ValueError(f"numtaps must be at least 1, not {numtaps}")
    if fs is None:
        fs = 1.0
    elif np.ndim(fs) != 0 or not fs > 0:
        # remez's wording, so that code checking it keeps working.
        raise ValueError("Sampling frequency fs must be a single scalar.")
    nyquist = 0.5 * float(fs)
    options = _Options(eps, nmax, strategy, depth, rstrategy, prec)

    if desired is None:
        if weight is not None:
            raise TypeError("weight is given by each Band; it can't also "
                            "be passed to remez()")
        if type != "bandpass":
            raise ValueError(f"type={type!r} needs band edges and desired "
                             "amplitudes, not Bands")
        return _design_bands(numtaps - 1, cast(Sequence[BandLike], bands),
                             nyquist, precision, options)

    edges = [float(e) for e in cast(Sequence[float], bands)]
    nbands = len(edges) // 2
    if weight is None:
        weight = [1.0] * nbands
    if len(edges) % 2 or len(desired) != nbands or len(weight) != nbands:
        raise ValueError(
            f"{len(edges)} band edges need {nbands} values each of desired "
            f"and weight, not {len(desired)} and {len(weight)}")
    if type == "bandpass":
        return _design_bands(
            numtaps - 1, [Band(edges[2*i], edges[2*i + 1], desired[i],
                               weight[i]) for i in range(nbands)],
            nyquist, precision, options)

    flat, _ = _routines(precision)
    f = [e / nyquist for e in edges]
    if type == "differentiator":
        # remez's slope is per unit of normalized frequency, f / fs.
        a = [desired[i // 2] * edges[i] / float(fs) for i in range(len(f))]
    else:
        a = [desired[i // 2] for i in range(len(f))]
    w = [float(x) for x in weight]
    if type == "differentiator":
        # The library divides the first band's weight by omega; remez
        # divides it by f / fs, which is omega / (2 pi).
        w[0] *= 2 * math.pi
    return cast(npt.NDArray[np.float64],
                flat(numtaps - 1, f, a, w, type, *options))


def firpm(n: int,
          f: Union[Sequence[float], Sequence[BandLike]],
          a: Optional[Sequence[float]] = None,
          w: Union[Sequence[float], str, None] = None,
          ftype: Optional[Literal["hilbert", "differentiator"]] = None, *,
          eps: float = 0.01,
          nmax: int = 4,
          strategy: Strategy = Strategy.UNIFORM,
          depth: int = 0,
          rstrategy: Strategy = Strategy.UNIFORM,
          precision: Literal["double", "mpfr"] = "double",
          prec: int = 165) -> npt.NDArray[np.float64]:
    """Design a linear-phase FIR filter, as MATLAB's ``firpm`` does.

    ``firpm(n, f, a)``, ``firpm(n, f, a, w)``, ``firpm(n, f, a, ftype)``
    and ``firpm(n, f, a, w, ftype)`` take the specifications of MATLAB's
    ``firpm``. ``firpm(n, bands)``, with a sequence of :class:`Band` (edges
    normalized so that 1 is Nyquist), allows amplitudes and weights that
    vary across a band.

    Parameters
    ----------
    n : int
        The filter order; ``n + 1`` taps are returned. Without ``ftype``,
        even ``n`` gives a type I filter and odd ``n`` a type II filter; an
        odd ``n`` whose last band has nonzero amplitude at Nyquist (which a
        type II filter cannot have) gives the type I filter of order
        ``n + 1``. With ``ftype``, even ``n`` gives type III and odd ``n``
        type IV.
    f : sequence of float, or sequence of Band
        Either band edges, a monotonic sequence of pairs from 0 to 1 (the
        Nyquist frequency), or the bands (see :class:`Band`; plain ``(f0,
        f1, amplitude[, weight])`` tuples also work).
    a : sequence of float, optional
        With band edges: the desired amplitude at each edge. Within a band,
        the amplitude is linear between its edges' values.
    w : sequence of float, optional
        With band edges: the weight of each band (1 by default). A string
        here is taken as ``ftype``, as in MATLAB.
    ftype : {"hilbert", "differentiator"}, optional
        An antisymmetric filter (types III and IV). For differentiators, the
        first band's weight is divided by the frequency. They need band
        edges, not :class:`Band`.
OPTIONS
    Returns
    -------
    numpy.ndarray
        The ``n + 1`` taps (``n + 2`` after promotion to type I).

RETURNS"""
    if isinstance(w, str):
        if ftype is not None:
            raise TypeError("firpm() takes ftype either in place of w, or "
                            "after it, not both")
        w, ftype = None, cast(Literal["hilbert", "differentiator"], w)
    if ftype not in (None, "hilbert", "differentiator"):
        raise ValueError("ftype must be 'hilbert' or 'differentiator', not "
                         f"{ftype!r}")
    if n < 0:
        raise ValueError(f"n must be at least 0, not {n}")
    options = _Options(eps, nmax, strategy, depth, rstrategy, prec)

    if a is None:
        if w is not None:
            raise TypeError("w is given by each Band; it can't also be "
                            "passed to firpm()")
        if ftype is not None:
            raise ValueError(f"ftype={ftype!r} needs band edges and "
                             "amplitudes, not Bands")
        return _design_bands(n, cast(Sequence[BandLike], f), 1.0, precision,
                             options)

    edges = [float(e) for e in cast(Sequence[float], f)]
    weights = [1.0] * (len(edges) // 2) if w is None else list(w)
    if len(edges) % 2 or len(a) != len(edges) or \
            len(weights) != len(edges) // 2:
        raise ValueError(
            f"{len(edges)} band edges need as many amplitudes and "
            f"{len(edges) // 2} weights, not {len(a)} and {len(weights)}")
    flat, _ = _routines(precision)
    return cast(npt.NDArray[np.float64],
                flat(n, edges, list(a), weights, ftype or "bandpass",
                     *options))


for _routine in (remez, firpm):
    _routine.__doc__ = cast(str, _routine.__doc__).replace(
        "OPTIONS\n", _OPTIONS).replace("RETURNS", _RETURNS)
del _routine
