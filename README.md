firpm library
==================================

## Python

The library is also available from Python, as the `firpm` package:

        pip install firpm

Binary wheels are published for Linux (x86_64, aarch64), macOS (arm64) and
Windows (x86_64). They compute in double precision only; see *Building from
source* below for MPFR support.

The package has two design routines, which differ only in how a
specification is written:

|                    | `firpm.remez`                  | `firpm.firpm`                    |
|--------------------|--------------------------------|----------------------------------|
| follows            | `scipy.signal.remez`           | MATLAB's `firpm`                 |
| length             | number of taps, `numtaps`      | filter order, `n` (`n + 1` taps) |
| band edges         | from 0 to `fs / 2` (`fs=1`)    | from 0 to 1 (Nyquist)            |
| desired amplitude  | one per band                   | one per edge (linear in a band)  |
| antisymmetric      | `type="hilbert"` etc.          | `ftype` (`"hilbert"` etc.)       |

### Designing a filter

```python
import numpy as np
import firpm

# 101-tap type I lowpass: passband [0, 0.2], stopband [0.25, 0.5] (in
# cycles/sample, so Nyquist is 0.5), with the stopband weighted 10x.
h = firpm.remez(101, [0.0, 0.2, 0.25, 0.5], [1.0, 0.0], weight=[1.0, 10.0])

# The same filter, in MATLAB's terms.
h = firpm.firpm(100, [0.0, 0.4, 0.5, 1.0], [1.0, 1.0, 0.0, 0.0], [1.0, 10.0])

H = np.fft.rfft(h, 4096)
```

`remez` takes the arguments of `scipy.signal.remez`, and `firpm` those of
MATLAB's `firpm` (`firpm(n, f, a)`, `firpm(n, f, a, w)`, and either with an
`ftype`). Odd `numtaps` (even `n`) gives a type I filter, and even
`numtaps` (odd `n`) a type II filter. `weight` and `w` default to 1 in each
band.

### Arbitrary amplitude and weight functions

Either routine also takes the specification band by band, as a sequence of
`firpm.Band(f0, f1, amplitude, weight=1.0)` (plain tuples work too). The
amplitude and the weight are each a number or a callable, so the target can
vary across a band; for example, to correct the passband droop of a CIC
decimator:

```python
import numpy as np
import firpm

R = 8   # CIC decimation ratio (third order)

def cic_inverse(f):
    if f == 0:
        return 1.0
    w = 2 * np.pi * f
    return (np.sin(w / 2) / (R * np.sin(w / (2 * R)))) ** -3

h = firpm.remez(201, [firpm.Band(0.0, 0.04, cic_inverse),
                      firpm.Band(0.06, 0.5, 0.0)])
```

Band edges, and the frequencies passed to callables, are in each routine's
own units: units of `fs` for `remez`, and normalized to Nyquist for `firpm`.
Callables are called with one frequency at a time; a design takes a few
thousand calls, so keep them cheap. Bands with only numeric amplitudes and
weights make no calls into Python at all. An exception raised inside a
callable is re-raised from the design routine.

### Differentiators and Hilbert transformers

`remez(..., type="differentiator")` or `type="hilbert"`, and
`firpm(..., ftype)` likewise, give an antisymmetric filter (type III for
odd `numtaps`, type IV for even):

```python
import firpm

h = firpm.remez(32, [0.0, 0.45], [1.0], type="differentiator")
h = firpm.remez(31, [0.05, 0.45], [1.0], type="hilbert")
h = firpm.firpm(30, [0.1, 0.9], [1.0, 1.0], "hilbert")
```

These take band edges only, not `Band`s. A `remez` differentiator's
amplitude is `desired * f / fs` in each band; a `firpm` one's is given at
the band edges, as in MATLAB. For differentiators, `scipy.signal.remez`
divides the weight by the frequency in every band with a nonzero slope,
but this library (like MATLAB) divides only the first band's. The two
agree for the usual differentiator, whose first band is the only one with
a slope. `scipy.signal.remez`'s `maxiter` and `grid_density` have no
counterpart here, since the library does not use a frequency grid.

### Options and errors

Both routines take these keyword options:

* `strategy` selects how the initial reference is chosen:
  `firpm.Strategy.UNIFORM` (the default), `firpm.Strategy.SCALING` or
  `firpm.Strategy.AFP`. Reference scaling (`strategy=firpm.Strategy.SCALING,
  depth=1`) usually needs far fewer iterations for long filters, and is more
  robust for them. `depth` is the number of scaling levels, and must be at
  least 1 for `SCALING` to have any effect.
* `eps` is the convergence threshold and `nmax` the degree of the local
  approximations used to locate the error extrema; the defaults rarely need
  changing.
* `precision="mpfr"` computes in MPFR arbitrary precision (`prec` bits, 165
  by default), when `firpm` was built with MPFR.

Both routines return the taps as a NumPy array. As in `scipy.signal.remez`,
an invalid specification or a failed design raises `ValueError`; a design
that did not fully converge emits a `RuntimeWarning` (use `warnings.simplefilter("error")` to
turn it into an exception). `help(firpm.remez)` and `help(firpm.firpm)`
have the details.

### Building from source

From a checkout, `pip install .` builds and installs the package. The
library's own dependencies (a C++20 compiler with OpenMP, Eigen) must be
available; see the installation instructions below. The published wheels
leave out MPFR, as `pip install .` does; for `precision="mpfr"`, use a CMake
build with MPFR, GMP and mpreal present. A CMake build also builds the Python
package whenever pybind11 is found, under `python/` in the build directory,
and `make test` then runs the Python tests too (these need `pytest` and
`numpy`; with `scipy` installed, they also compare designs against
`scipy.signal.remez`).

The package's version is `firpm.__version__`, and it ships type
annotations.

## Folder organization

Prior to November 2019, there were three different folders containing
implementations with different numerical precision used for the computations:
* double (64-bit)
* long double (80-bit on x86 architectures)
* MPFR-based custom extendend precision

They have now been merged into one template version.

## Installation instructions

This code has been tested on recent Mac OS X (>= 10.14) and Linux installations.
In order to compile and use it, a recent C++ compiler with C++20 and OpenMP support is necessary (recent g++ and clang compilers have been shown to work). An active internet
connection and some external utilities and libraries must also be installed and
available on your system search paths:
* CMake version 3.16 or newer
* Eigen version 3.3 or newer
* (optional) GMP version 6.0 or newer for multiple precision support
* (optional) MPFR version 4.0 or newer for multiple precision support
* (optional) [mpreal](https://github.com/advanpix/mpreal) wrapper for MPFR
* (optional) doxygen to generate the accompanying code documentation
* Google gtest framework for generating the test executables will be downloaded during CMake configuration

Assuming these prerequisites are taken care of and you are located at the base
folder, building the library can be done using the following commands
(with the appropriate user privileges):

        mkdir build
        cd build
        cmake ..
        make all
        make install

This series of steps will install the library on your system. Usually the default
location is /usr/local, but this behavior can be changed when calling CMake by
specifying an explicit install prefix:

        cmake -DCMAKE_INSTALL_PREFIX=/your/install/path/here ..


With gtest downloaded in the build directory, the *make all* command should
have generated two test executables:
* firpmlib_scaling_test : should contain code to generate all the example filters used in [1] and tests for subsequent bugs
* firpmlib_extensive_test : contains a more extensive set of over 50 different filters, giving an iteration count comparison for uniform initialization, reference scaling and AFP initialization, when applicable

The *make test* command will launch these two test executables on your system.

The *make all* target also generates the documentation if Doxygen was found on
your system when running CMake. It can also be generated individually by running
the command

        make doc

after CMake was called.


Running the test executables will print out information regarding the final reference error
(i.e., final delta value) obtained when executing the Parks-McClellan exchange algorithm,
along with iteration count information. Depending on your machine, these tests can take
some time to finish.

An example output for one test case from firpm_scaling_test (using the long double
instantiation) of the routines looks like this:

        [ RUN      ] firpm_scaling_test/1.combfir
        START Parks-McClellan with uniform initialization
        Final Delta     = 1.6039185593030203561e-07
        Iteration count = 13
        FINISH Parks-McClellan with uniform initialization
        START Parks-McClellan with reference scaling
        Final Delta     = 1.6066986792983147255e-07
        Iteration count = 3
        FINISH Parks-McClellan with reference scaling
        START Parks-McClellan with AFP
        Final Delta     = 1.606565353078586782e-07
        Iteration count = 4
        FINISH Parks-McClellan with AFP
        Iteration count reduction for final filter  RS: 0.76923076923076916245
        Iteration count reduction for final filter AFP: 0.69230769230769229061
        [       OK ] firpm_scaling_test/1.combfir (378 ms)

## Use

Examples of how to use the library can be found in the **test** folder.

## Licensing

Up until April 7, 2024 the code was GPLv3+ licensed. It has since switched to a more permissive 3-Clause BSD license.

## References
[1] S.-I. Filip, A robust and scalable implementation of the Parks-McClellan
algorithm for designing FIR filters, ACM Trans. Math. Softw., vol. 43,
no. 1, Aug. 2016, Art. no. 7.

[2] S.-I. Filip, Robust tools for weighted Chebyshev approximation and
applications to digital filter design, Ph.D. dissertation, ENS de Lyon, France, 2016.

## Acknowledgements
* Graeme Smecher [@gsmecher]( https://github.com/gsmecher ) - setting up the template version of the code base
* [@jlwehle]( https://github.com/jlwehle ) - various bug reports and fixes
