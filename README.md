# FaVeST: Fast Vector Spherical Harmonic Transforms
This is our Matlab implementation for FaVeST in the paper

>Q. T. Le Gia, M. Li, Y. G. Wang. [FaVeST: Fast Vector Spherical Harmonic Transforms](https://arxiv.org/abs/1908.00041). arXiv preprint arXiv:1908.00041, 2019.

Author: Quoc Thong Le Gia (qleqia@unsw.edu.au), Ming Li (ming.li.ltu@gmail.com;mingli@zjnu.edu.cn), Yu Guang Wang (yuguang.wang@unsw.edu.au).

## Abstract
Vector spherical harmonics on $\mathbb{S}^{2}\subset \mathbb{R}^3$ have wide applications in geophysics, quantum mechanics and astrophysics. In the representation of a tangent field, one needs to evaluate the expansion and the Fourier coefficients of vector spherical harmonics. In this paper, we develop fast algorithms (FaVeST) for vector spherical harmonic transforms for these evaluations. The forward FaVeST which evaluates the Fourier coefficients has computational steps proportional to $N\log \sqrt{N}$ for $N$ number of evaluation points. The adjoint FaVeST which evaluates a linear combination of vector spherical harmonics with degree up to $\sqrt{M}$ for $M$ evaluation points is proportional to $M\log\sqrt{M}$. Numerical examples illustrate the accuracy and efficiency of FaVeST.

## Citation 
If you use find our package useful, please cite our paper:
```
@article{FaVeST,
  title={FaVeST: Fast Vector Spherical Harmonic Transforms},
  author={Le Gia, Quoc T. and Li, Ming and Wang, Yu Guang},
  journal={arXiv preprint arXiv:1908.00041},
  year={2019}
}
```
## Environment Requirement
The code has been tested in Matlab environment in Windows, Linux and MacOS. 

## Functions
* **utils**: the folder contains some basic tools/resources/auxiliary functions used for implementing our main functions, including
   1. SD: it contains six examples of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html) used in our numerical experiments.
   2. tangent_field: the folder that contains several functions for generating three vector fields and their visualization used in our paper. These functions come from E. J. Fuselier and G. B. Wright [*Stability and error estimates for vector field interpolation and decomposition on the sphere with RBFs. SIAM Journal on Numerical Analysis, 47(5):3213-39*](https://epubs.siam.org/doi/abs/10.1137/080730901).
   3. m_map: [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this package for visualization of tangent fields. 
   4. QpS2.m: the function is used for computing the weights and quadrature nodes (for a given degree and a specific type of quadrature rule) in either Cartesian coordinates or spherical coordinates. 

* **nfft-3.5.2-matlab-openmp**: The pre-compiled Matlab interfaces of NFFT 3.5.2 with AVX2 and OpenMP support, downloaded from [NFFT library](https://www-user.tu-chemnitz.de/~potts/nfft/): Keiner, J., Kunis, S., and Potts, D. [*Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms, ACM Trans. Math. Software,36, Article 19, 1-30, 2009*](https://dl.acm.org/citation.cfm?id=1555388). As stated in [NFFT Downloads](https://www-user.tu-chemnitz.de/~potts/nfft/download.php), this version can be compiled on 64-bit Windows using GCC 9.3.0 x86_64-w64-mingw32 with -march=haswell and Matlab R2018b, x86_64 Linux using GCC 8.3.0 with -march=haswell and Matlab R2017b Update 9, or x86_64 macOS using GCC 9.3.0 with -march=haswell and Matlab R2017b. The package also includes the binaries for Octave 5.2.0 on Windows.

* **Demo.m**: the function tests **FaVeST_fwd.m** and **FaVeST_adj.m** on a tangent field. It is used to test whether users have successfully configured NFFT packages by **Setup.m**. 

* **FaVeST_fwd.m**: forward FFTs computing Fourier coefficients associated with a quadrature rule (e.g., Algorithm 1 in the paper):
```
T - tangent field samples
L - degree for vector spherical harmonic
X,w - quadrature rule used for evaluating FFT
```

* **FaVeST_adj.m**: adjoint FFTs for vector spherical harmonic expansion with given inputs (e.g., Algorithm 2 in the paper): 
```
alm -  Fourier coefficients for divergent-free part
blm - Fourier coefficients of curl-free part
X - evaluation points on the sphere
```


* **Fig2a,2b,2c.m**, **Fig3a,3b,3c.m**, **Table1.m**, **Table2_Fig4.m**: these routines are used to reproduce the numerical results of the corresponding figures and tables of the paper.

## Demo
We provide a demo **Demo.m** for **FaVeST**, which can be run in Windows, Linus or MacOS. Other numerical tests in our paper are also availabel. Users can try, for example, **Fig3a.m**. 

<img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_rec_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_err_gl.png" width="250">


## Acknowledgements
We would thank E. J. Fuselier and G. B. Wright for providing their MATLAB program which generates simulated tangent fields. The NFFT package is used for the **FaVeST** package.  M. Li acknowledges support from the National Natural Science Foundation of China under Grant 61802132, and the Australian Research Council under Discovery Project
DP160101366 when he worked with P. Broadbridge and A. Olenko at La Trobe University. Q. T. Le Gia and Y. G. Wang acknowledge support from the Australian Research Council under Discovery Project DP180100506.

## Contributing
Copyright (c) <2020> <NeurIPS>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
