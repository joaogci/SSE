subroutine read_parameters_analysis(n_rebin) &
  bind(c,name='read_parameters_analysis')

  use,intrinsic :: iso_c_binding,only:c_char,c_int,c_double,c_long
  implicit none

  integer(kind=c_int),    intent(inout) :: n_rebin

  namelist /Analysis/ n_rebin

  open(unit=100,file='parameters',status='old')
  read(unit=100,nml=Analysis)
  close(unit=100)
  
endsubroutine read_parameters_analysis

