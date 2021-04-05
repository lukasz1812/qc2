!         __ _       __
!       / /| |     / /
!     / / | | /| / /   this program was written by Łukasz Wantoch
!   / /__| |/ |/ /     During the laboratory course of WP 4 in summerterm 2021
! /_____/__/|__/       at the Rheinrich-Wilhelm-University of Bonn

module scf_main
    !> Include standard Fortran environment for IO
    use iso_fortran_env, only : output_unit, error_unit

    ! ------------------------------------------------------------------------
    !> library functions provided by your lab assistents:

    !> interface to LAPACK's double precision symmetric eigenvalue solver (dspev)
    !  examples:
    !  call solve_spev(mat, eigval, eigvec)
    use linear_algebra, only : solve_spev

    !> expansion of slater-functions into contracted gaussians,
    !  coefficients and primitive exponents are taken from R.F. Stewart, JCP, 1970
    !  example:
    !  call expand_slater(zeta, alpha, coeff)
    use slater, only : expand_slater

    !> calculates one-electron integrals and two-electron integrals over
    !  spherical gaussians (s-functions). One-electron quanities supported
    !  are overlap, kinetic energy and nuclear attraction integrals.
    !  Two-electron integrals are provided in chemist notation.
    !  examples:
    !  call oneint(xyz, chrg, r_a, r_b, alp, bet, ca, ca, s, t, v)
    !  call twoint(r_a, r_b, r_c, r_d, alp, bet, gam, del, ca, cb, cc, cd, g)
    use integrals, only : oneint, twoint

    !> prints a matrix quantity to screen
    !  examples:
    !  call write_vector(vec, name='vector')
    !  call write_matrix(mat, name='matrix')
    !  call write_matrix(mat, name='packed matrix')
    use print_matrix, only : write_vector, write_matrix

    !> other tools that may help you jump ahead with I/O-heavy tasks
    !  example:
    !  call read_line(input, line)
    use io_tools, only : read_line

    !> Always declare everything explicitly
    implicit none

    !> All subroutines within this module are not exported, except for scf_prog
    !  which is the entry point to your program
    private
    public :: scf_prog

    !> Selecting double precision real number
    integer, parameter :: wp = selected_real_kind(15)


contains


!> This is the entry point to your program, do not modify the dummy arguments
!  without adjusting the call in lib/prog.f90
subroutine scf_prog(input)

    !> Always declare everything explicitly
    implicit none

    !> IO unit bound to the input file
    integer, intent(in) :: input

    !> System specific data
    !> Number of atoms
    integer :: nat

    !> Number of electrons
    integer :: nel

    !> Atom coordinates of the system, all distances in bohr
    real(wp), allocatable :: xyz(:,:)

    !> Nuclear charges
    real(wp), allocatable :: chrg(:)

    !> Number of basis functions
    integer :: nbf

    !>Number of Gausiian functions for Slater expansion
    integer :: ng

    !> Slater exponents of basis functions
    real(wp),allocatable :: zeta(:)

    !> Nuclear repulsion energy
    real(wp) :: Erep

    !> Pointer for result file
    integer :: io2

    !> Hartree-Fock energy
    real(wp) :: escf

    !  declarations may not be complete, so you have to add your own soon.
    !  Create a program that reads the input and prints out final results.
    !  And, please, indent your code.

    !  Write the self-consistent field procedure in a subroutine.


    open(file="results.txt", newunit=io2)
    call banner(io2)

    call input_reader (nat,nel,nbf,xyz,chrg,zeta,io2)

    call NucRep(nat,xyz,chrg,Erep,io2)



    write(output_unit, '(a)') 'Here could start a Hartree-Fock calculation'

end subroutine scf_prog

!> Below all subroutines used in the scf_prog subroutine are declared, the subroutines are listed in the same order as apearance in the scf_prog
!In order to make the program structured the beginning of each subroutine is defined by a line of ł the full name of subroutine inbetween
!the end of each subroutine is followed by a lines od $ with an "-End-" inbetweeen


!łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł   PROGRAM TITLE     łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł
subroutine banner(io2)
  integer :: io2

  write (*,*) "======================================================================================================="
  write (*,*) "                      ___                    ___                           ___   "
  write (*,*) "                     /  /\                  /  /\                         /  /\  "
  write (*,*) "                    /  /:/_                /  /:/                        /  /:/_ "
  write (*,*) "                   /  /:/ /\              /  /:/                        /  /:/ /\"
  write (*,*) "                  /  /:/ /::\            /  /:/  ___                   /  /:/ /:/"
  write (*,*) "                 /__/:/ /:/\:\          /__/:/  /  /\                 /__/:/ /:/ "
  write (*,*) "                 \  \:\/:/~/:/          \  \:\ /  /:/                 \  \:\/:/  "
  write (*,*) "                  \  \::/ /:/            \  \:\  /:/                   \  \::/   "
  write (*,*) "                   \__\/ /:/              \  \:\/:/                     \  \:\   "
  write (*,*) "                     /__/:/                \  \::/                       \  \:\  "
  write (*,*) "                     \__\/ E L F            \__\/ O N S I S T E N T       \__\/ I E L D "
  write(*,*) ""
  write(*,*) "                         A Hartree-Fock program with use of Rothan-Haal equations"
  write (*,*) "======================================================================================================="
  write(*,*) ""
  write (io2,*) "======================================================================================================="
  write (io2,*) "                      ___                    ___                           ___   "
  write (io2,*) "                     /  /\                  /  /\                         /  /\  "
  write (io2,*) "                    /  /:/_                /  /:/                        /  /:/_ "
  write (io2,*) "                   /  /:/ /\              /  /:/                        /  /:/ /\"
  write (io2,*) "                  /  /:/ /::\            /  /:/  ___                   /  /:/ /:/"
  write (io2,*) "                 /__/:/ /:/\:\          /__/:/  /  /\                 /__/:/ /:/ "
  write (io2,*) "                 \  \:\/:/~/:/          \  \:\ /  /:/                 \  \:\/:/  "
  write (io2,*) "                  \  \::/ /:/            \  \:\  /:/                   \  \::/   "
  write (io2,*) "                   \__\/ /:/              \  \:\/:/                     \  \:\   "
  write (io2,*) "                     /__/:/                \  \::/                       \  \:\  "
  write (io2,*) "                     \__\/ E L F            \__\/ O N S I S T E N T       \__\/ I E L D "
  write (io2,*) ""
  write (io2,*) "                         A Hartree-Fock program with use of Rothan-Haal equations"
  write (io2,*) "======================================================================================================="
  write(io2,*) ""
end subroutine banner

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END PROGRAM TITLE $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! -----------------------------------------------------------------------------------------------------------------------

!łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł   INPUT READER     łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł
subroutine input_reader (nat,nel,nbf,xyz,chrg,zeta,io2)

  !>declaration of local variables
  character(len = 100):: input
  integer :: nat, nel, nbf, io, i, j,k,  dim,io2
  real(wp), allocatable :: xyz(:,:),chrg(:), zeta(:), basis(:)


  dim=3
  k=1
  j=1
  write(*,*) "Give the input file name :"
  read(*,*) input

  open(file="molecules/"//input, newunit=io)
  read(io,*) nat,nel,nbf
  allocate (xyz(dim,nat), chrg(nat), basis(nat), zeta(nbf))


    do i=1,nat

      read (io,*) xyz(j:dim,i),chrg(i), basis(i)

      do while (k<=sum(basis))

        read(io,*) zeta(k)
        k=k+1

      end do

    end do

  close(io)
!>variable check

write(*,*)
write(*,*) "                      ==================================================="
write(*,*) "                      ||               System parameters               ||"
write(*,*) "                      |–––––––––––––––––––––––––––––––––––––––––––––––––|"
write(*,*) "                      |Number of atoms          |", nat,"          |"
write(*,*) "                      |Number of electrons      |", nel,"          |"
write(*,*) "                      |Number of basis functions|", nbf,"          |"
write(*,*) "                      ==================================================="
write(*,*)
write(*,*)
write(io2,*) "                      ==================================================="
write(io2,*) "                      ||               System parameters               ||"
write(io2,*) "                      |–––––––––––––––––––––––––––––––––––––––––––––––––|"
write(io2,*) "                      |Number of atoms          |", nat,"          |"
write(io2,*) "                      |Number of electrons      |", nel,"          |"
write(io2,*) "                      |Number of basis functions|", nbf,"          |"
write(io2,*) "                      ==================================================="
write(io2,*)
write(io2,*)
Write(io2,*)
call write_matrix(xyz, "       ===== Atom positions/[Bohr] =====", io2)
write(io2,*)"        ================================="
write(io2,*)
write(io2,*)
call write_vector(chrg,  "       == Charge/[q] ==", io2)
write(io2,*)"        ================"
write(io2,*)
write(io2,*)
call write_vector(basis,"        == Basis fct. ==", io2)
write(io2,*)"        ================"
write(io2,*)
write(io2,*)
call write_vector(zeta,"       ==== ζ exp. ====",io2)
write(io2,*)"        ================"

end subroutine input_reader

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END INPUT READER $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! -----------------------------------------------------------------------------------------------------------------------

!łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł   NUCLEAR REPULSION ENERGY    łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł
subroutine NucRep(nat,xyz,chrg,Erep,io2)

  !>Declaration of global variables
  integer :: nat, io2
  real(wp) :: Erep
  real(wp), allocatable :: xyz(:,:), chrg(:)

  !>Declaration of local variables
  integer :: i, j
  real(wp) :: distance,partRep
  real(wp) :: diff(3)

  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                 Calculating nuclear repulsion energy"
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                 Calculating nuclear repulsion energy"
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)

  j=1
  i=j+1

  if (nat==1) then
    Erep=0
  else
    do while (i<=nat)
      diff=xyz(1:3,j)-xyz(1:3,i)
      distance=sqrt(sum(diff**2))
      write(*,*) distance
      partRep=(chrg(j)*chrg(i))/distance**2
      Erep=Erep+partRep
      j=j+1
      i=i+1
    end do
  end if
  write(*,*) "                      ==================================================="
  write(*,*) "                      ||            Nuclear repulsion energy           ||"
  write(*,*) "                      |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(*,*) "                      |", Erep, "H                     |"
  write(*,*) "                      ==================================================="
  write(io2,*) "                      ==================================================="
  write(io2,*) "                      ||            Nuclear repulsion energy           ||"
  write(io2,*) "                      |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(io2,*) "                      |", Erep, "H                     |"
  write(io2,*) "                      ==================================================="


end subroutine NucRep

!$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$ END NUCLEAR REPULSION ENERGY $$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

! -----------------------------------------------------------------------------------------------------------------------

!łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł   SLATER EXPANSION   łłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłłł
subroutine expansion(ng, nbf, zeta, exponents, coefficients)

  !>Declaration of global variables
  integer :: ng, nbf
  real(wp), allocatable:: zeta(:), exponents(:), coefficients(:)
end module scf_main
