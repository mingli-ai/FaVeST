# FaVeST: Fast Vector Spherical Harmonic Transforms
This is our Matlab implementation for the paper:

>Q. T. Le Gia, M. Li, Y. G. Wang. [FaVeST: Fast Vector Spherical Harmonic Transforms](https://arxiv.org/abs/xxxxxx). arXiv....., 2019.

Author: Dr. Quoc Thong Le Gia (qleqia@unsw.edu.au), Dr. Ming Li (ming.li.ltu@gmail.com), Dr. Yu Guang Wang (yuguang.wang@unsw.edu.au).

## Abstract
Vector spherical harmonics on $\mathbb{S}^{2}\subset \mathbb{R}^3$ have wide applications in geophysics, quantum mechanics and astrophysics. In the representation of a tangent field, one needs to evaluate the expansion and the Fourier coefficients of vector spherical harmonics. In this paper, we develop fast algorithms (FaVeST) for vector spherical harmonic transforms for these evaluations. The forward FaVeST which evaluates the Fourier coefficients has computational steps proportional to $N\log \sqrt{N}$ for $N$ number of evaluation points. The adjoint FaVeST which evaluates a linear combination of vector spherical harmonics with degree up to $\sqrt{M}$ for $M$ evaluation points is proportional to $M\log\sqrt{M}$. 
Numerical examples illustrate the accuracy and efficiency of FaVeST.

## Citation 
If you want to use our codes and datasets in your research, please cite:
```
@article{FaVeST,
  title={Fast Vector Spherical Harmonic Transforms},
  author={},
  journal={arXiv PrePrint},
  volume={xx,
  number={xx},
  pages={xx
  year={2019}
}
```
## Environment Requirement
The code has been tested running in Matlab. To run the codes, users need to install the required packages are as follows:
* [NFFT library](https://www-user.tu-chemnitz.de/~potts/nfft/): Keiner, J., Kunis, S., and Potts, D. ''Using NFFT 3 - a software library for various nonequispaced fast Fourier transforms'' ACM Trans. Math. Software,36, Article 19, 1-30, 2009.
* After a successful installation of NFFT library fitting your OS, users can add the NFFT package folder 'xx' into the current working directory (i.e., ../FaVeST/xx).

## Functions and Folders
* SD points: This folder saves an example of [symmetric spherical design points](https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/ss.html), i.e., ss041.mat, containing 864 points satisfying quadrature rule that is exact for spherical polynomials of degree 41. 

* m_map: A [mapping package](https://www.eoas.ubc.ca/~rich/map.html#ack) for Matlab. We have used some functions of this tool in visualizing vector fields. 

* *QpS2.m*: function used for computing the weights and quadrature nodes (for a given degree and type of quadrature points) in both Cartesian and spherical coordinates. 

* *Vecfield_sim.m*: function used for generating simulated vector fileds with vector spherical harmonics with random coefficients. We use this toy example to display how to run the main FaVeST functions.

* *FaVeST_fwd.m*: Main function for implementing forward FFTs computing Fourier coefficients with given inputs: T-the vector field samples; L-the degree for vector spherical harmonic; X,w -the quadrature rule used for evaluating FFT. See Algorithm 1 in our paper.

* *FaVeST_adj.m*: Main function for implementing adjoint FFTs for vector spherical harmonic expansion with given inputs: alm- Fourier coefficients for divergent-free part; blm-Fourier coefficients of curl-free part; X -given points on the sphere.

* *Demo.m*: This is used to demonstrate how to run **FaVeST_fwd.m** and **FaVeST_adj.m** with a given vector field. 

* *Visualization.m*: This is used to illstrate the simulation results with three plots including the target vector field, the approximated vector field and the error.

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
