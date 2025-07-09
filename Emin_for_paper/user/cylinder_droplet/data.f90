!PF_NPHASE_SURFPOT_SQUARES
subroutine set_defaults()
  use potential
  
  
  !#################!
  ! User Parameters !
  !#################!	
  
  OPEN(1,FILE='data.in')
  READ(1,*) NPOSTX 
  READ(1,*) NPOSTY
  READ(1,*) WIDTHX
  READ(1,*) WIDTHY 
  READ(1,*) TOPEXX 
  READ(1,*) TOPEXY 
  READ(1,*) HEIGHT1 
  READ(1,*) HEIGHT2
  READ(1,*) LIPX 
  READ(1,*) LIPY
  READ(1,*) LIPZ 

  READ(1,*) GRIDX
  READ(1,*) GRIDY
  READ(1,*) GRIDZ 
  READ(1,*) N_PHASE
  READ(1,*) LAMBDA  
  CLOSE(1)
 
  
  ALPHA=1.0
  
  K44=0.0
  
  CONF_STRENGTH=-1.0D0
  
  XMIRROR=.TRUE. !If .FALSE., simulate periodic x-boundaries, otherwise simulate mirror x-boundaries
  YMIRROR=.TRUE. !If .FALSE., simulate periodic y-boundaries, otherwise simulate mirror y-boundaries
  ZMIRROR=.TRUE. !If .FALSE., simulate periodic z-boundaries, otherwise simulate mirror z-boundaries
 
  GRIDSIZE=1 !Size of lattice spacing

  N = GRIDX*GRIDY*GRIDZ*(N_PHASE-1) !Number of degrees of freedom

  !###################!
  ! L-BFGS Parameters !
  !###################!
  M = 10
  max_iterations  = 1000000
  convergence_rms = 1d-6
  max_step_size   = 0.25d0
  H0init          = 1d1
  relative_energy_check = .true.
  dE_max          = 1d-1
  ! Print convergence information to GMIN_OUT
  debug = .true.

  !##########################!
  ! Basin Hopping Parameters !
  !##########################!

  ! Number of basin hopping steps to take
  max_runs          = 1
  max_similarity    = 1d-4
  temperature       = 1000d8
  binary_io         = .false.
  ! Save only the lowest i minima. 0 Saves all
  minima_save_limit = 0
end subroutine set_defaults
