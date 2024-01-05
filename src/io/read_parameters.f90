subroutine read_parameters(Lx,Ly,boundary_condition,beta,epsilon, &
  therm_cycles,n_bins,mc_sweeps,S,J_perp,J_par,magnetic_field,crystal_field) &
  bind(c,name='read_parameters')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: Lx
  integer(kind=c_int),    intent(inout) :: Ly
  character(kind=c_char), intent(inout) :: boundary_condition

  real(kind=c_double),    intent(inout) :: beta
  real(kind=c_double),    intent(inout) :: epsilon
  integer(kind=c_long),   intent(inout) :: therm_cycles
  integer(kind=c_int),    intent(inout) :: n_bins
  integer(kind=c_long),   intent(inout) :: mc_sweeps

  real(kind=c_double),    intent(inout) :: S
  real(kind=c_double),    intent(inout) :: J_perp
  real(kind=c_double),    intent(inout) :: J_par
  real(kind=c_double),    intent(inout) :: magnetic_field
  real(kind=c_double),    intent(inout) :: crystal_field

  namelist /Lattice/ Lx,Ly,boundary_condition
  namelist /Simulation/ beta,epsilon,therm_cycles,n_bins,mc_sweeps
  namelist /Hamiltonian/ S,J_perp,J_par,magnetic_field,crystal_field

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Lattice)
  read(unit=100,nml=Simulation)
  read(unit=100,nml=Hamiltonian)
  close(unit=100)
  
endsubroutine read_parameters

