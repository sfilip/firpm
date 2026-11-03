# Changelog

## 0.2.2 (Nov. 3, 2026)

The first release published to PyPI, as `firpm`.

### Added

- Python bindings: the `firpm` package, built with pybind11, with binary
  wheels for Linux (x86_64, aarch64), macOS (arm64) and Windows (x86_64).
  `firpm.remez()` takes the arguments of `scipy.signal.remez`, and
  `firpm.firpm()` those of MATLAB's `firpm`; both return NumPy arrays and
  support differentiators and Hilbert transformers. Either also takes one
  `Band` per band, whose amplitude and weight may be callables of
  frequency. Invalid specifications and failed designs raise
  `ValueError`, as in `scipy.signal.remez`, and convergence warnings are
  reported as `RuntimeWarning`. With MPFR available at build time,
  `precision="mpfr"` computes in arbitrary precision; the wheels are built
  without it. The package is typed, and has a `__version__`.
- A `firpm()` overload taking `band_t` bands with arbitrary amplitude and
  weight callables, for non-flat magnitude specifications such as CIC
  droop correction. Exceptions raised in Python callbacks are re-raised on
  the calling thread.
- A numerics test executable, `firpmlib_numerics_test`, which checks the
  numerical kernels against 256-bit MPFR.

### Changed

- This is the first release under the BSD license, which replaced GPLv3
  upstream in April 2024.
- C++20 is required.
- Eigen 5.x is supported as well as 3.x.
- The tests build and run without MPFR.

### Fixed

- A design could come back wrong by as much as its ripple while reporting
  success: the final filter was levelled with the error of the previous
  iteration's reference, whose sign flips when the last iteration shifts
  the reference by one node. Symmetric Hilbert transformers hit this most
  easily (e.g. `remez(11, [0.1, 0.4], [1], type="hilbert")`). Every design
  also loses the small spurious component this left in its taps.
- With the default `nmax=4`, the extremum search placed the error's
  extrema slightly off their true peaks, so designs converged short of the
  minimax optimum (by 2e-4 relative on scipy's `test_compare` designs). The
  located extrema are now refined in the final iterations, and designs
  reach the optimum at the default `nmax`, for about 10% more work in those
  iterations.
- `uniform()` could give the last band a negative reference count and
  write out of bounds for three or more bands of uneven width. It now
  always allocates exactly `n` references, and rejects layouts that cannot
  be filled with `std::domain_error`.
- The leveled error `delta` and the error function are now computed
  without cancellation against the desired response. Designs whose error
  approaches double's rounding floor (e.g. comb filters) used to converge
  only by luck, or report success on a wrong reference.
- The barycentric weights can no longer overflow or underflow at any
  filter length, and the error is evaluated stably for ill-conditioned
  reference sets.
- Extremum search: nodes rounded just outside a band are pinned to the
  band edge, and references are counted correctly when a band has none.
- Concurrent designs in separate threads (e.g. from Python, which releases
  the GIL) no longer share exchange state.
- MSVC builds: OpenMP loops use signed indices, and pi no longer relies on
  `M_PI`. `long double` now uses a full-precision pi.

## 0.2.1 and earlier

Not recorded here; see the git history.
