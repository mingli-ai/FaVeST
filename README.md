# FaVeST: Fast Vector Spherical Harmonic Transforms
This is our Matlab implementation for FaVeST in the paper

>Q. T. Le Gia, M. Li, Y. G. Wang. [FaVeST: Fast Vector Spherical Harmonic Transforms](https://arxiv.org/abs/1908.00041). arXiv preprint arXiv:1908.00041, 2019.

Author: Dr. Quoc Thong Le Gia (qleqia@unsw.edu.au), Dr. Ming Li (ming.li.ltu@gmail.com;mingli@zjnu.edu.cn), Dr. Yu Guang Wang (yuguang.wang@unsw.edu.au).

## Abstract
Vector spherical harmonics on $\mathbb{S}^{2}\subset \mathbb{R}^3$ have wide applications in geophysics, quantum mechanics and astrophysics. In the representation of a tangent field, one needs to evaluate the expansion and the Fourier coefficients of vector spherical harmonics. In this paper, we develop fast algorithms (FaVeST) for vector spherical harmonic transforms for these evaluations. The forward FaVeST which evaluates the Fourier coefficients has computational steps proportional to $N\log \sqrt{N}$ for $N$ number of evaluation points. The adjoint FaVeST which evaluates a linear combination of vector spherical harmonics with degree up to $\sqrt{M}$ for $M$ evaluation points is proportional to $M\log\sqrt{M}$. Numerical examples illustrate the accuracy and efficiency of FaVeST.

## Citation 
If you use our codes and datasets, please cite:
```
@article{FaVeST,
  title={FaVeST: Fast Vector Spherical Harmonic Transforms},
  author={Le Gia, Quoc T. and Li, Ming and Wang, Yu Guang},
  journal={arXiv preprint arXiv:1908.00041},
  year={2019}
}
```
## Environment Requirement
The code has been tested in Matlab environment. 

## Functions and Folders
* **utils**: This folder contains some basic tools/resources/auxiliary functions used for implementing our main functions, including
   1. SD: A folder that contains six examples of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html), which are used in our paper. 
   2. tangent_field: A folder that contains several functions for generating three vector fields and their visualization used in our paper. These functions come from E. J. Fuselier and G. B. Wright as they use in "[*Stability and error estimates for vector field interpolation and decomposition on the sphere with RBFs. SIAM Journal on Numerical Analysis, 47(5):3213-39*](https://epubs.siam.org/doi/abs/10.1137/080730901)".
   3. m_map: A [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this package for visualization of tangent fields. 
   4. QpS2.m: A function used for computing the weights and quadrature nodes (for a given degree and a specific type of quadrature rule) in both Cartesian and spherical coordinates. 

* **nfft-3.5.2-matlab-openmp**: The pre-compiled Matlab interfaces of NFFT 3.5.2 with AVX2 and OpenMP support, downloaded from [NFFT library](https://www-user.tu-chemnitz.de/~potts/nfft/): Keiner, J., Kunis, S., and Potts, D. "[*Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms, ACM Trans. Math. Software,36, Article 19, 1-30, 2009*](https://dl.acm.org/citation.cfm?id=1555388)". This version was compiled on 64-bit Windows using GCC 9.3.0 x86_64-w64-mingw32 with -march=haswell and Matlab R2018b, x86_64 Linux using GCC 8.3.0 with -march=haswell and Matlab R2017b Update 9, x86_64 macOS using GCC 9.3.0 with -march=haswell and Matlab R2017b. In addition, the binaries for Octave 5.2.0 on Windows were included.

* **Demo.m**: It tests **FaVeST_fwd.m** and **FaVeST_adj.m** on a tangent field. It is used to test whether users have successfully configured NFFT packages by **Setup.m**. 

* **FaVeST_fwd.m**: Main function for implementing forward FFTs computing Fourier coefficients associated with a quadrature rule: T - tangent field samples; L - degree for vector spherical harmonic; X,w - quadrature rule used for evaluating FFT. See Algorithm 1 in our paper.

* **FaVeST_adj.m**: Main function for implementing adjoint FFTs for vector spherical harmonic expansion with given inputs: alm -  Fourier coefficients for divergent-free part; blm - Fourier coefficients of curl-free part; X - evaluation points on the sphere.


* **Fig2a,2b,2c.m**, **Fig3a,3b,3c.m**, **Table1.m**, **Table2_Fig4.m**: These routines are used to reproduce the numerical results of the corresponding figures and tables of our paper.

## Demo
Users can try different settings by running **Demo.m**, or test our simulation programs by running other m-scripts (e.g. **Fig3a.m**). 

<img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_rec_gl.png" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/vf_1_err_gl.png" width="250">


## Acknowledgements
We would thank E. J. Fuselier and G. B. Wright for providing their MATLAB program which generates simulated tangent fields. The NFFT package is used for the **FaVeST** package. M. Li acknowledges support from the Australian Research Council under Discovery Project DP160101366. Q. T. Le Gia and Y. G. Wang acknowledge support from the Australian Research Council under Discovery Project DP180100506.

## Notes
The package **FaVeST** may be used for any research purposes under the following conditions:
* The user must acknowledge the use of **FaVeST** in publications resulting from the use of the functions/tools.
* The user may not redistribute **FaVeST** without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from us.
