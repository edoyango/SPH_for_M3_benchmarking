!---------------------------------------------------------
!     Including file for parameters and constants used 
!     in the entire SPH software packages.
!---------------------------------------------------------

module param

	use constants, only: g
	
	! double or single precision (chance f to match)
	integer,parameter,public:: df = kind(1.d0)
	integer,parameter,public:: sf = kind(1.)
	integer,parameter,public:: f = df
	
	!dim : Dimension of the problem (1, 2 or 3)
	integer,parameter,public:: dim = 2
	
	!Smoothing kernel function 
	!skf = 1, cubic spline kernel by W4 - Spline (Monaghan 1985)
	!    = 2, Gauss kernel   (Gingold and Monaghan 1981) 
	integer,parameter,public:: skf = 1
	
	!spacing and kernel radii parameters
	real(f),parameter,public:: dxo = 0.1d0, kappa = 1.2d0, v_max = 22.15d0
	
	!material density (per particle)
	real(f),parameter,public:: irho = 1000d0
	
	!derived parameters. c: speed of sound, hsml: smoothing length, dt: time-step size, mass: mass per particle
	real(f),parameter,public:: c = 10d0*v_max,hsml = kappa*dxo, dt = 1.5d0*hsml/c,mass=irho*dxo**dim
	
	integer,parameter,public:: mp = 250, np = 250, op = 750, pp = 400
	
	! state equation parameter,public::s
	real(f),parameter,public:: rh0 = irho
	integer,parameter,public:: gamma = 7
	
	! artificial viscosity parameter,public::s
	real(f),parameter,public:: alpha = 0.01d0, beta = 0.d0, etq = 0.1d0
	
	! repulsive force parameter,public::s
	real(f),parameter,public:: rr0 = dxo,dd = 5d0*g*25d0
	integer,parameter,public:: p1=4,p2=2
	
	character(len=200),parameter,public:: output_directory = "outputdata"
	
	character(len=3),parameter,public:: output_flt_type='dbl'
	logical,parameter,public:: output_phys(2) = (/.true.,.true./)
	logical,parameter,public:: output_halo(2) = (/.true.,.false./)
	logical,parameter,public:: output_virt(2) = (/.true.,.false./)
	logical,parameter,public:: output_CPUtime = .false.
	logical,parameter,public:: output_boundary = .false.

end
