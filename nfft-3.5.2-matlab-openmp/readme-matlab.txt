This archive contains the Matlab interface of NFFT 3.5.2
compiled for 64-bit Windows using GCC 9.3.0 x86_64-w64-mingw32
with -march=haswell and Matlab R2018b,
x86_64 Linux using GCC 8.3.0 with -march=haswell
and Matlab R2017b Update 9,
x86_64 macOS using GCC 9.3.0 with -march=haswell
and Matlab R2017b.
In addition, the binaries for Octave 5.2.0 on Windows are included.

Please note that since the binaries were compiled with gcc flag -march=haswell,
they may not work on older CPUs (below Intel i3/i5/i7-4xxx or
AMD Excavator/4th gen Bulldozer) as well as on some Intel Atom/Pentium CPUs.


NFFT - Nonequispaced FFT
=========================

Overview
--------
NFFT is a software library, written in C, for computing non-equispaced fast
Fourier transforms and related variations. It implements the following
transforms:

1. Non-equispaced fast Fourier transform (NFFT)
    - forward transform *(NFFT)*, i.e. frequency to time/space domain
    - adjoint transform *(adjoint NFFT)*, i.e. time/space to frequency domain

2. Generalisations
    - to arbitrary nodes in time *and* frequency domain *(NNFFT)*
    - to real-valued data, i.e. (co)sine transforms, *(NFCT, NFST)*
    - to the sphere S^2 *(NFSFT)*
    - to the rotation group *(NFSOFT)*
    - to the hyperbolic cross *(NSFFT)*

3. Generalised inverse transformations based on iterative methods, e.g. CGNR/CGNE

Some examples for application of these transforms are provided:

1. Medical imaging
    - magnetic resonance imaging (mri)
    - computerised tomography (radon)

2. Summation schemes
    - fast summation (fastsum)
    - fast Gauss transform (FGT)
    - singular kernels
    - zonal kernels

3. polar FFT, discrete Radon transform, ridgelet transform

Detailed API documentation in HTML format can be found in
`doc/html/index.html`, if you are working from a release tarball.
When working from a source repository, the documentation can be
generated with Doxygen.
```
make doc
```

Building
--------
The NFFT depends on the [FFTW](https://fftw.org) library, which is available for many Linux distros, Homebrew on macOS and MSYS2 on Windows. If you compile the FFTW yourself, it should be configured `--enable-shared`.

When working from a source repository, you need to run libtoolize and autoreconf first. A bash script to do this is provided.
```
./bootstrap.sh
```

The rest of the build process is standard.
```
./configure --enable-all --enable-openmp [add options as necessary, see below]
```

Alternatively, you might run the configure script for Matlab.
```
./configure --enable-all --enable-openmp --with-matlab=/path/to/matlab
```

Here are some useful optional flags for `./configure`:
* `--enable-all` specifies that all modules should be compiled,
* `--enable-openmp` enables the multicore support and
* `--enable-julia` specifies that the julia interface will be compiled.
* `--with-matlab=path/to/matlab` specifies a path of Matlab, and
* `--with-octave=path/to/octave` does the same for GNU Octave.
* For a list of all available options, run `./configure --help`.

Build the software.
```
make
```

Optionally, unit tests may be run.
```
make check
```

Optionally, install NFFT on your system.
```
make install
```

Citing
------
The most current general paper, the one that we recommend if you wish to cite NFFT, is *Keiner, J., Kunis, S., and Potts, D.
''Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms''
ACM Trans. Math. Software,36, Article 19, 1-30, 2009*.

Feedback
--------
Your comments are welcome! This is the third version of the library and may
not be as robust or well documented as it should be. Please keep track of bugs
or missing/confusing instructions and report them to
[Daniel Potts](mailto:potts@mathematik.tu-chemnitz.de).
The postal address is

```
  Prof. Dr. Daniel Potts
  TU Chemnitz, Fakultaet fuer Mathematik
  Reichenhainer Str. 39
  09107 Chemnitz
  GERMANY
```

Alternatively, you might contact
[Stefan Kunis](mailto:stefan.kunis@math.uos.de)
or
[Jens Keiner](mailto:jens@nfft.org).

If you find NFFT useful, we would be delighted to hear about what application
you are using NFFT for!

Legal Information & Credits
---------------------------
Copyright (c) 2002, 2017 Jens Keiner, Stefan Kunis, Daniel Potts

This software was written by Jens Keiner, Stefan Kunis and Daniel Potts.
It was developed at the Mathematical Institute, University of
Luebeck, and at the Faculty of Mathematics, Chemnitz University of Technology.

NFFT3 is free software. You can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation; either version 2 of the License, or (at your option) any later
version. If not stated otherwise, this applies to all files contained in this
package and its sub-directories.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA 

FFTW
----
The compiled NFFT files contain parts of the FFTW library (http://www.fftw.org)
Copyright (c) 2003, 2007-14 Matteo Frigo
Copyright (c) 2003, 2007-14 Massachusetts Institute of Technology
