Module Variables

  ! Critical temperature for a 3D Heisenberg model is Tc = 0.945d0.

  ! lx=10, ly=1, lz=1 means a 1D chain lattice with size 10.
  ! lx=10, ly=10, lz=1 means a 2D square lattice.
  ! lx=10, ly=10, lz=10 means a 3D cubic lattice.

  Integer*4, parameter :: lx = 2                     ! Number of spins in x.
  Integer*4, parameter :: ly = 1                     ! Number of spins in y.
  Integer*4, parameter :: lz = 1                     ! Number of spins in z.
  Integer*4, parameter :: Nbins = 20                  ! Averages written to file after mcsteps.
  Integer*4, parameter :: mcsteps = 1e4               ! Monte Carlo steps.
  Integer*4, parameter :: term_steps = 1e4            ! Termalization steps.
  Integer*4, parameter :: N = lx * ly * lz            ! Total number of spins.
  Real*8,    parameter :: temp_ini = 1.05d0            ! Initial temperature.
  Integer*4, parameter :: t_steps = 10                ! Number of temperature steps.
  Real*8,    parameter :: dt = -0.025d0                 ! Size of temperature steps.

  Integer*4 :: NH                   ! Number of H-operator.
  Integer*4 :: L                    ! Maximum string size. (Don't change it)
  Integer*4 :: d                    ! Dimension. (Selected by the program)
  Integer*4 :: Nb                   ! Number o bonds between spins. (N = d*lx*ly*lz)
  Integer*4 :: Nn

  Integer*4, Allocatable, Dimension (:)   :: spin                  ! Spin value.
  Integer*4, Allocatable, Dimension (:)   :: opstring              ! Operator string ID.
  Integer*4, Allocatable, Dimension (:,:) :: bound_spin_ID         ! List of the spin sites connect between the bound_spin_ID(1, bound_idx) and bound_spin_ID(2, bound_idx)
  Integer*4, Allocatable, Dimension (:)   :: vertex_link           ! List of vertex links.
  Integer*4, Allocatable, Dimension (:)   :: first_vertex_visitted ! First operator on each site in linked vertex list.
  Integer*4, Allocatable, Dimension (:)   :: last_vertex_visitted  ! Last operator on each site in linked vertex list.

  Integer*4 :: NH_max
  Character*60 :: arq

  Real*8    :: beta
  Real*8    :: rem_prob
  Real*8    :: add_prob
  Real*8    :: n_opH = 0.0d0, n_opH2 = 0.0d0
  Real*8    :: ususc = 0.0d0, staggered = 0.0d0

End Module Variables

Program Main

  Use Variables
  Implicit None

  Integer*4 t  

  write(arq, '("results_",I0,"x",I0,"x",I0,".dat")') lx, ly, lz
  Open(40, file=arq, position='append')

  write(40, '(4x,"Beta",4x,"Energy",4x,"Specific Heat",4x,"Uniform Susceptibility", &
           &  4x,"Magnetization^2",4x,"Number of H-operators",    &
           &  4x,"(Number of H-operators)^2" )')

  Do Nn = 0, 100

  Do t = 0, t_steps - 1

    beta = 1.0d0 / (temp_ini + dt * t)

    NH = 0; NH_max = 0;
    L = max(4, N / 4)

    Allocate(opstring(0:L - 1))
    Allocate(first_vertex_visitted(N), last_vertex_visitted(N))
    Allocate(vertex_link(0:4 * L - 1))

    opstring = 0

    Call lattice
    Call init

    add_prob = 0.5d0 * beta * Nb
    rem_prob = 1.0d0 / (0.5d0 * beta * Nb)

    Call termalization

    Call qmc_steps

    Call results

    Call make_dist

    Call free_memory

    Call flush()

  end do

  end do

  Close(40)

End Program Main

Subroutine termalization

  Use Variables
  Implicit None

  Integer*4 i

  Do i = 1, term_steps

    Call diagonalupdate
    Call linkvertices
    Call loopupdate
    Call adjustcutoff(i)

  end do

  write(*,'(1x,"Series cut-off: ",I0)') L
  write(*,*)'----------------------------'
  write(*,*)' '

End Subroutine termalization

Subroutine qmc_steps

  Use Variables
  Implicit None

  Integer*4 i, j

  Do j = 1, Nbins

    Do i = 1, mcsteps

      Call diagonalupdate
      Call linkvertices
      Call loopupdate
      Call measure

    end do

    Call write_results

  end do 

  Close(10)
  Close(20)
  Close(30)

End Subroutine qmc_steps

Subroutine lattice

  Use Variables
  Implicit None

  Integer*4 :: i, j, k, bound_idx, spin_idx

  ! mod(x, L), in this case, returns 0 if x == Lx or returns x if x != Lx.

  write(*,*) " "
  write(*,*)"heisenberg_sse.f90"
  write(*,*) " "

  if ((lx > 1 .and. ly > 1) .and. (lz == 1)) then

    print*, '2D with Periodic boundary '

    d = 2                     
    Nb = d * lx*ly            
    Allocate(bound_spin_ID(2, Nb))

    Do i = 0, lx - 1
      Do j = 0, ly - 1
        spin_idx = 1 + i + j * lx

        bound_idx = 1 + i + j * lx

        bound_spin_ID(1, bound_idx) = spin_idx
        bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx) + j * lx

        bound_spin_ID(1, bound_idx + lx * ly) = spin_idx
        bound_spin_ID(2, bound_idx + lx * ly) = 1 + i + mod(j + 1, ly) * lx
      end do
    end do

  elseif ( (lx > 1) .and. (ly == 1 .and. lz == 1) ) then

    print*, '1D with Periodic boundary '

    d = 1                     
    Nb = d * lx               
    Allocate(bound_spin_ID(2, Nb))

    Do i = 0, lx - 1
       spin_idx = i + 1

       bound_idx = i + 1

       bound_spin_ID(1, bound_idx) = spin_idx
       bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx)
    end do

  elseif (lx > 1 .and. ly > 1 .and. lz > 1) then

    print*, '3D with Periodic boundary '

    d = 3                     
    Nb = d * lx * ly * lz
    Allocate(bound_spin_ID(2, Nb))

    Do i = 0, lx - 1
      Do j = 0, ly - 1
        Do k = 0, lz - 1
          
          spin_idx = 1 + i + j * lx + k * lx * ly

          bound_idx = 1 + i + j * lx + k * lx * ly

          bound_spin_ID(1, bound_idx) = spin_idx
          bound_spin_ID(2, bound_idx) = 1 + mod(i + 1, lx) + j * lx + k * lx * ly

          bound_spin_ID(1, bound_idx + lz * lx * ly) = spin_idx
          bound_spin_ID(2, bound_idx + lz * lx * ly) = 1 + i + mod(j + 1, ly) * lx + k * lx * ly

          bound_spin_ID(1, bound_idx + 2 * lz * lx * ly) = spin_idx
          bound_spin_ID(2, bound_idx + 2 * lz * lx * ly) = 1 + i + j * lx + mod(k + 1, lz) * lx * ly

        end do
      end do
    end do

  else

    print*, 'For a 1D lattice: lx > 0, ly==1 and lz==1'
    print*, 'For a 2D lattice: lx > 0, ly > 0 and lz==1'
    print*, 'For a 3D lattice: lx > 0, ly > 0 and lz > 0'
    stop
  end if

End Subroutine lattice

Subroutine init

  Use Variables
  Implicit None

  Integer*4 :: i

  Allocate(spin(N))

  ! Spins with random values of -1 or 1.

  Do i = 1, N
    spin(i) = (-1)**(mod(i - 1, lx) + (i - 1) / lx)
  End do

  print*,' '
  write(*,'(1x,"Lattice: ", 2x,"lx=",I0,2x," ly=",I0,2x," lz=",I0)') lx, ly, lz
  print*,' '
  write(*,'(1x,"Number of spins: ",I0)') N
  write(*,'(1x,"Number of bins: ",I0)') Nbins
  write(*,'(1x,"Monte Carlo steps: ",ES7.1)') dble(mcsteps)
  write(*,'(1x,"Termalization steps: ",ES7.1)') dble(term_steps)
  write(*,'(1x,"Temperatures: ",F6.4," -> ",F6.4)') 1.0d0 / beta, temp_ini + (t_steps - 1) * dt
  write(*,'(1x,"N: ",I0)') Nn

  write(arq, '("bin_",I0,"x",I0,"x",I0,"_T=",F6.4,".dat")') lx, ly, lz, 1.0d0 / beta
  Open(10, file=arq)

  write(10, '(4x,"Energy",4x,"Specific Heat",4x,"Uniform Susceptibility", &
           &  4x,"Magnetization^2",4x,"Number of H-operators",    &
           &  4x,"(Number of H-operators)^2" )')
  

  write(arq, '("raw_",I0,"x",I0,"x",I0,"_T=",F6.4,".dat")') lx, ly, lz, 1.0d0 / beta
  Open(20, file=arq)

  write(20, '(4x,"Number of H-operators")')

  write(arq, '("cut-off_",I0,"x",I0,"x",I0,"_T=",F6.4,".dat")') lx, ly, lz, 1.0d0 / beta
  Open(30, file=arq)

  write(30, '(4x,"Step",4x,"Cut-off",4x,"Number of H-operators", &
           &  4x,"Maximum Number of H-operators" )')

End Subroutine init

Subroutine diagonalupdate

  Use Variables
  Implicit None

  Integer*4 :: p, bound_idx, op
  Real*8 :: ran

  ! b(p) is the bound index. 
  ! a(p) is the operator type, here 1 or 2.

  ! opstring = 2 * b(p) + a(p) - 1

  ! opstring == 0    => Identity operator.  
  ! opstring == even => Diagonal operator.
  ! opstring == odd  => Off-diagonal operator.

  ! Since 2 * b(p) is always integer even, (a(p) - 1)
  ! is the term which defines the operator type.
  ! If a(p) = 1, opstring is even.
  ! If a(p) = 2, opstring is odd.

  ! From the opstring value is possible to get
  ! bouth a(p) or b(p) using the relations,
  ! a(p) = MOD(opstring[p], 2) + 1, and
  ! b(p) = opstring[p] / 2.

  ! print*, step

  Do p = 0, L - 1

    op = opstring(p)

    if (op == 0) then ! 0 == Identity operator.

      ! Add a diagonal operator if the spins connected
      ! by a bound, chosen at random, are anti-parallel.

      Call random_number(ran)

      bound_idx = int(ran * Nb) + 1

      if (spin(bound_spin_ID(1, bound_idx)) /= spin(bound_spin_ID(2, bound_idx))) then

        Call random_number(ran)

        if (ran * (NH + 1 - Nn) <= (add_prob / (L - NH)) * (NH + 1) ) then

          opstring(p) = 2 * bound_idx
          NH = NH + 1
        end if

      end if

    elseif (mod(op, 2) == 0) then ! Even == Diagonal operator.

      ! Remove a diagonal operator.
      
      Call random_number(ran)

      if (ran * NH <= (L - NH + 1) * rem_prob * (NH - Nn) ) then

        opstring(p) = 0
        NH = NH - 1
      end if

    else !Odd == Off-diagonal operator.

      ! Change off-diagonal operators is more complex. Here,
      ! we only make a basic flip of spins every time a off-diagonal
      ! operator is found. This is allowed because a valid configuration
      ! in this system always have with a zero or a odd number of 
      ! off-diagonal operators. An even number of off-diagonal operators
      ! breaks the boundary condition imposed by the trace \sum <\sigma|H|\sigma>. 

      bound_idx = op / 2

      spin(bound_spin_ID(1, bound_idx)) = -spin(bound_spin_ID(1, bound_idx))
      spin(bound_spin_ID(2, bound_idx)) = -spin(bound_spin_ID(2, bound_idx))

    end if

  end do

End Subroutine diagonalupdate

Subroutine linkvertices

  Use Variables
  Implicit None

  Integer*4 :: v0, p, op, bound_idx, first_vertex, last_vertex
  Integer*4 :: last_vertex1, last_vertex2
  Integer*4 :: spin_idx, spin_idx1, spin_idx2

  first_vertex_visitted(:) = -1
  last_vertex_visitted(:) = -1

  Do v0 = 0, 4 * L - 1, 4
    
    p = v0 / 4

    op = opstring(p)
    
    if (op /= 0) then
       
       bound_idx = op / 2
       
       spin_idx1 = bound_spin_ID(1, bound_idx)
       spin_idx2 = bound_spin_ID(2, bound_idx)
       
       last_vertex1 = last_vertex_visitted(spin_idx1)
       last_vertex2 = last_vertex_visitted(spin_idx2)
       
       if (last_vertex1 /= -1) then

          vertex_link(last_vertex1) = v0
          vertex_link(v0) = last_vertex1

       else

          first_vertex_visitted(spin_idx1) = v0

       endif

       if (last_vertex2 /= -1) then

          vertex_link(last_vertex2) = v0 + 1
          vertex_link(v0 + 1) = last_vertex2

       else

          first_vertex_visitted(spin_idx2) = v0 + 1

       endif

       last_vertex_visitted(spin_idx1) = v0 + 2
       last_vertex_visitted(spin_idx2) = v0 + 3

    else

       vertex_link(v0:v0 + 3) = -1

    endif
  enddo

  Do spin_idx = 1, N
    
    first_vertex = first_vertex_visitted(spin_idx)
    
    if (first_vertex /= -1) then

        last_vertex = last_vertex_visitted(spin_idx)

        vertex_link(last_vertex) = first_vertex
        vertex_link(first_vertex) = last_vertex

    endif
  enddo
  
End Subroutine linkvertices

Subroutine loopupdate

  Use Variables
  Implicit None

  Real*8 :: ran
  Integer*4 :: spin_idx, v0, vertex_in, vertex_out

  ! -1 => Visited vertex but not flipped.
  ! -2 => Visited and flipped vertex.

  Do v0 = 0, 4 * L - 1, 2

    ! If the vertex v0 was visited, cycle the loop.
    if (vertex_link(v0) < 0) cycle

    vertex_in = v0

    Call random_number(ran)

    if (ran < 0.5d0) then

      do 

        ! Change the operator type.
        opstring(vertex_in / 4) = ieor(opstring(vertex_in / 4), 1)
        
        ! Mark as visitted
        vertex_link(vertex_in) = -2
        
        ! Next vertex
        vertex_out = ieor(vertex_in, 1)

        vertex_in = vertex_link(vertex_out)

        vertex_link(vertex_out) = -2
        
        ! If the loop is closed, exit.
        if (vertex_in == v0) exit

      end do

    else

      do 

        vertex_link(vertex_in) = -1

        vertex_out = ieor(vertex_in, 1)

        vertex_in = vertex_link(vertex_out)

        vertex_link(vertex_out) = -1

        if (vertex_in == v0) exit

      end do

    end if

  end do

  ! Flip free spin lines. (Lines without operators)

  Do spin_idx = 1, N

    if (first_vertex_visitted(spin_idx) /= -1) then

      if (vertex_link(first_vertex_visitted(spin_idx)) == -2) then
        spin(spin_idx) = -spin(spin_idx)
      end if

    else
      
      Call random_number(ran)

      if (ran < 0.5d0) spin(spin_idx) = -spin(spin_idx)

    end if

  end do 

End Subroutine loopupdate

Subroutine adjustcutoff(step)

  Use Variables
  Implicit None

  Integer*4, Allocatable, Dimension (:) :: copy_opstring
  Integer*4 :: L_new, step

  L_new = NH + NH / 3

  ! Uncomment for save data about the cut-off adjust.

  ! Call cut_off_analisys(L_new, step)

  if (L_new <= L) return

  Allocate(copy_opstring(0:L - 1))
  
  copy_opstring(:) = opstring(:)

  Deallocate(opstring)
  Allocate(opstring(0:L_new - 1))
  
  opstring(0:L - 1) = copy_opstring(:)
  opstring(L:L_new - 1) = 0
  
  Deallocate(copy_opstring)

  L = L_new
  
  Deallocate (vertex_link)
  Allocate(vertex_link(0:4 * L - 1))

End Subroutine adjustcutoff

Subroutine measure

  Use Variables
  Implicit None


  integer :: i, bound_idx, op, spin_idx1, spin_idx2, stag_mag
  real(8) :: stag_mag2

  stag_mag = 0

  Do i = 1, N
     stag_mag = stag_mag + spin(i) * (-1)**(mod(i - 1, lx) + (i - 1) / lx)
  end do      

  stag_mag = stag_mag / 2
  stag_mag2 = 0.d0

  Do i= 0, L - 1

     op = opstring(i)

     if (mod(op,2) == 1) then  

        bound_idx = op / 2

        spin_idx1 = bound_spin_ID(1, bound_idx)
        spin_idx2 = bound_spin_ID(2, bound_idx)

        spin(spin_idx1) = -spin(spin_idx1)
        spin(spin_idx2) = -spin(spin_idx2)

        stag_mag = stag_mag + 2 * spin(spin_idx1) * (-1)**(mod(spin_idx1 - 1, lx) + (spin_idx1 - 1) / lx)

     end if

     stag_mag2 = stag_mag2 + dfloat(stag_mag)**2

  end do

  stag_mag2 = stag_mag2 / dble(L)

  staggered = staggered + stag_mag2
  ususc = ususc + dble(sum(spin) / 2)**2

  n_opH = n_opH + dble(NH)
  n_opH2 = n_opH2 + dble(NH)**2

  write(20, *) NH, - ( NH / (beta * N) - 0.25d0 * dble(Nb) / dble(N))

End Subroutine measure

Subroutine write_results

  Use Variables
  Implicit None

  Real*8 :: Cv, ener

  ! print*, n_opH, n_opH2, mcsteps

  n_opH = n_opH / dble(mcsteps)
  n_opH2 = n_opH2 / dble(mcsteps)
  ususc = ususc / dble(mcsteps)
  staggered = staggered / dble(mcsteps)

  ener = - ( n_opH / (beta * N) - 0.25d0 * dble(Nb) / dble(N))

  Cv = (n_opH2 - n_opH**2 - n_opH) / dble(N)

  ! Cv = (n_opH**2 - n_opH2 - n_opH) / dble(N)

  staggered = 3.d0 * staggered / dble(N)**2

  ususc = beta * ususc / dble(N)

  n_opH = n_opH / dble(N)
  n_opH2 = n_opH2 / dble(N)

  write(10, *) ener, Cv, ususc, staggered, n_opH, n_opH2

  n_opH = 0.0d0
  n_opH2 = 0.0d0
  ususc = 0.0d0
  staggered = 0.0d0

End Subroutine write_results

Subroutine results

  Use Variables
  Implicit None

  Real*8 :: e(Nbins), c(Nbins), n_op(Nbins), n_op2(Nbins), x(Nbins), stag(Nbins)
  Integer*4 :: i

  write(arq, '("bin_",I0,"x",I0,"x",I0,"_T=",F6.4,".dat")') lx, ly, lz, 1.0d0 / beta
  Open(10, file=arq)

  read(10, '(A)')

  Do i = 1, Nbins

    read(10, *) e(i), c(i), x(i), stag(i), n_op(i), n_op2(i)

  end do 

  Close(10)

  e = e / dble(Nbins)
  n_op = n_op / dble(Nbins)
  n_op2 = n_op2 / dble(Nbins)
  c = c / dble(Nbins)
  x = x / dble(Nbins)
  stag = stag / dble(Nbins)

  write(40, *) beta, sum(e), sum(c), sum(x), sum(stag), sum(n_op), sum(n_op2)

End Subroutine results

Subroutine make_dist

  Use Variables
  Implicit None

  Integer*4 opn(Nbins * mcsteps), i, idx
  Real*8  energy(Nbins * mcsteps)
  Integer*4, Allocatable, Dimension(:) :: n_dist, e_dist
  Real*8 energy_bin, e_min, e_max

  energy_bin = 1.0d0 / (beta * N)

  write(arq, '("raw_",I0,"x",I0,"x",I0,"_T=",F6.4,".dat")') lx, ly, lz, 1.0d0 / beta
  Open(20, file=arq)

  read(20, '(A)')

  Do i = 1, Nbins * mcsteps

    read(20, *) opn(i), energy(i)

  end do

  Close(20)

  Allocate(n_dist(minval(opn):maxval(opn)))
  n_dist = 0

  Do i = 1, Nbins * mcsteps

    n_dist(opn(i)) = n_dist(opn(i)) + 1

  end do

  write(arq, '("dist_n_",I0,"x",I0,"x",I0,"_T=",F6.4,"Nn=",I0,".dat")') lx, ly, lz, 1.0d0 / beta, Nn
  Open(50, file=arq)

  Do i = minval(opn), maxval(opn)

    write(50, *) i, n_dist(i)

  End do

  Close(50)

  e_min = minval(energy)
  energy = energy - e_min

  e_max = maxval(energy)

  Allocate(e_dist(0: int(e_max / energy_bin)))

  e_dist = 0 

  Do i = 1, Nbins * mcsteps

    idx = int(energy(i) / energy_bin)

    e_dist(idx) = e_dist(idx) + 1

  end do

  write(arq, '("dist_e_",I0,"x",I0,"x",I0,"_T=",F6.4,"Nn=",I0,".dat")') lx, ly, lz, 1.0d0 / beta, Nn
  Open(50, file=arq)

  Do i = 0, int(e_max / energy_bin)

    write(50, *) i * energy_bin + e_min, e_dist(i)

  End do

  Close(50)  


End Subroutine make_dist

Subroutine cut_off_analisys(L_new, step)

  Use Variables
  Implicit None

  Integer*4, intent(in) :: L_new, step

  if ( NH > NH_max) NH_max = NH

  if (L_new <= L) then

    write(30, *) step, L, NH, NH_max
    return

  else 

    write(30, *) step, L_new, NH, NH_max

  end if

End Subroutine cut_off_analisys

Subroutine free_memory

  Use Variables
  Implicit None

  Deallocate(spin)
  Deallocate(opstring)
  Deallocate(bound_spin_ID)
  Deallocate(first_vertex_visitted)
  Deallocate(last_vertex_visitted)
  Deallocate(vertex_link)

End Subroutine free_memory

