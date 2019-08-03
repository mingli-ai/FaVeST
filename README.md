# FaVeST: Fast Vector Spherical Harmonic Transforms
This is our Matlab implementation for the paper:

>Q. T. Le Gia, M. Li, Y. G. Wang. [FaVeST: Fast Vector Spherical Harmonic Transforms](https://arxiv.org/abs/1908.00041). arXiv preprint arXiv:1908.00041, 2019.

Author: Dr. Quoc Thong Le Gia (qleqia@unsw.edu.au), Dr. Ming Li (ming.li.ltu@gmail.com), Dr. Yu Guang Wang (yuguang.wang@unsw.edu.au).

## Abstract
Vector spherical harmonics on $\mathbb{S}^{2}\subset \mathbb{R}^3$ have wide applications in geophysics, quantum mechanics and astrophysics. In the representation of a tangent field, one needs to evaluate the expansion and the Fourier coefficients of vector spherical harmonics. In this paper, we develop fast algorithms (FaVeST) for vector spherical harmonic transforms for these evaluations. The forward FaVeST which evaluates the Fourier coefficients has computational steps proportional to $N\log \sqrt{N}$ for $N$ number of evaluation points. The adjoint FaVeST which evaluates a linear combination of vector spherical harmonics with degree up to $\sqrt{M}$ for $M$ evaluation points is proportional to $M\log\sqrt{M}$. Numerical examples illustrate the accuracy and efficiency of FaVeST.

## Citation 
If you want to use our codes and datasets in your research, please cite:
```
@article{FaVeST,
  title={FaVeST: Fast Vector Spherical Harmonic Transforms},
  author={Le Gia, Quoc T. and Li, Ming and Wang, Yu Guang},
  journal={arXiv preprint arXiv:1908.00041},
  year={2019}
}
```
## Environment Requirement
The code has been tested running in Matlab. To run the FaVeST functions and demos, users need to download the NFFT package:
* [NFFT library](https://www-user.tu-chemnitz.de/~potts/nfft/): Keiner, J., Kunis, S., and Potts, D. "[*Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms, ACM Trans. Math. Software,36, Article 19, 1-30, 2009*](https://dl.acm.org/citation.cfm?id=1555388)".  [Download Website](https://www-user.tu-chemnitz.de/~potts/nfft/download.php)
* We provide a script **Setup.m** to help download the recent version: 	*nfft-3.5.0-mexw64-openmp.zip* or *nfft-3.5.0-mexw32-openmp.zip* for **Windows binaries**, or *nfft-3.5.0-mexmaci64-openmp.zip* for **MacOS binaries**. In this manner, users can just run **Setup.m** in Matlab when using for the first time. 

## Functions and Folders
* **utils**: This folder contains some basic tools/resources/auxiliary functions used for implementing our main functions, such as:
   1. SD: A folder that saves some examples of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html), i.e., corresponding to the six cases used in our paper. 
   2. tangent_field: A folder that contains several functions for generating three vector fields and the visualization used in our paper. These functions come from E. J. Fuselier and G. B. Wright who have simulated in "[*Stability and error estimates for vector field interpolation and decomposition on the sphere with RBFs. SIAM Journal on Numerical Analysis, 47(5):3213-39*](https://epubs.siam.org/doi/abs/10.1137/080730901)".
   3. m_map: A [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this tool in the visualization of the vector fields. 
   4. QpS2.m: A function used for computing the weights and quadrature nodes (for a given degree and type of quadrature points) in both Cartesian and spherical coordinates. 

* **Setup.m**: The script used for downloading NFFT package (matching the current operating system and computing environment) and unzipping the package and putting in the current folder with FaVeST. If successful, users can use FaVeST by referring to the demos and examples. 

* **Demo.m**: This is used to test whether users have successfully configured NFFT packages by running **Setup.m**, as an immediate trial  for how to run **FaVeST_fwd.m** and **FaVeST_adj.m** with a given vector field. 

* **FaVeST_fwd.m**: Main function for implementing forward FFTs computing Fourier coefficients with given inputs: T-the vector field samples; L-the degree for vector spherical harmonic; X,w -the quadrature rule used for evaluating FFT. See Algorithm 1 in our paper.

* **FaVeST_adj.m**: Main function for implementing adjoint FFTs for vector spherical harmonic expansion with given inputs: alm- Fourier coefficients for divergent-free part; blm-Fourier coefficients of curl-free part; X -given quadrature rule points on the sphere.


* **Fig.2*.m**, **Fig.3*.m**, **Table1.m**, **Table2_Fig4.m**: These functions are used to reproduce the simulation results demonstrated in our paper.

## Demo
We provide a simple demonstration by running `Demo.m`.
* Both GL and SD points satisfying quadrature rule that is exact for spherical polynomials of certain degree can be tried by specifying `QN = 'GL'` or `QN = 'SD'`
* For a toy illustration, we only use `L=20` in the simulation, corresponding to 882 GL points and 864 SD points.

Running `Demo.m` can produce the following basic results, showing the effectiveness of FaVeST:

#### GL points
<img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Tar_VF_GL.jpg" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Appro_VF_GL.jpg" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Error_VF_GL.jpg" width="250">

#### SD points
<img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Tar_VF_SD.jpg" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Appro_VF_SD.jpg" width="250"><img src="https://github.com/mingli-ai/FaVeST/blob/master/images/Error_VF_SD.jpg" width="250">



## Acknowledgement
We would thank E. J. Fuselier and G. B. Wright for providing their MATLAB codes implementing the artificial tangent fields and Matlab routines for the visualization of vector fields. We also thank J. Keiner, S. Kunis, and D. Potts, for their NFFT library implementing FFTs for scalar spherical harmonics.

## Notes
The package **FaVeST** may be used for any research purposes under the following conditions:
* The user must acknowledge the use of **FaVeST** in publications resulting from the use of the functions/tools.
* The user may not redistribute **FaVeST** without separate permission.
* The user may not use this information for any commercial or revenue-bearing purposes without first obtaining permission from us.
