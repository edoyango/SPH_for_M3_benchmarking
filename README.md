# SPH_basic_para
basic WCSPH with RK4 time-integration. Written in FORTRAN.
Parallelised using MPI and domain-partitioning parallelisation strategy.
Continual improvement project of my master's work (Yang et al., 2020)

Todo:
- Improve ORB partitioning (remove need to copy entire grid)
- Translate code style (line length, all subroutines in modules, comments etc.) to serial code.
- Try to obtain better low-core scaling performance
- GPU code

Reference:
Yang E et al. (2020) A scalable parallel computing SPH frame-work for predictions of geophysical granular flows. Comput Geotech 121:103474
