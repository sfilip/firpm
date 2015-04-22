firpm library
==================================

## Folder organization and installation instructions

This directory contains three different folders corresponding to the
three implementations that are mentioned in the submitted article:

* firpm_d   : corresponds to the double (64-bit) precision version of the routines
* firpm_ld  : corresponds to the long double (80-bit on x86 architectures) precision version
* firpm_mp  : corresponds to the multiple precision (MPFR-based) version

Inside each folder you will find a README.md file describing the external dependencies of that
particular version along with instructions on how to compile and install it, generate the
corresponding Doxygen documentation or the test executables, which are designed to
use Google's gtest unit testing framework.

The code has been tested *only* on machines with Linux-based installations. For example, a recent
version of Ubuntu Linux (like versions 14.4 and 14.10) should work without problem.

We opted to provide the end user with the choice of which version to use. This is mostly due to the
fact that the MPFR multiple precision version of the routines require more external dependencies (GMP and MPFR)
than the double and long double code, which the end user might not necessarily be willing to use.

## Licensing

The provided code is primarily GPLv3+ licensed.
