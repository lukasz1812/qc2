!        __  __       __
!      / /  | |     / /
!     / /  | | /| / /   this program was written by Łukasz Wantoch
!   / /_ _| |/ |/ /     During the laboratory course of WP 4 in summerterm 2021
!  /_____/__/|__/       at the Rheinrich-Wilhelm-University of Bonn

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

    !> Gaussian exponents and coeefficients of slater expansion
    real(wp),allocatable :: exponents(:)
    real(wp),allocatable :: coefficients(:)

    !> Nuclear repulsion energy
    real(wp) :: Erep

    !> One electron overlap integral
    real(wp),allocatable :: sab(:,:), packsab(:)


    !> One electron overlap integral
    real(wp),allocatable :: tab(:,:), packtab(:)

    !> One electron overlap integral
    real(wp),allocatable :: vab(:,:), packvab(:)

    !> Two electron integrals
    real(wp),allocatable :: twointeg(:)

    !> Fock matrix
    real(wp),allocatable :: Fock(:,:)

    !> Coefficients matrix
    real(wp), allocatable :: cab(:,:)

    !> Density matrix
    real(wp), allocatable :: pab(:,:)

    !> Eigenvalues and eigenvectors for LAPACK eigenvalue solver
    real(wp),allocatable :: eigval(:)
    real(wp),allocatable :: eigvec(:,:)

    !>Symmetric orhonormalizer
    real(wp),allocatable :: xab(:,:)

    !> Pointer for result file
    integer :: io2


    !> Hartree-Fock energy
    real(wp) :: escf

    !  declarations may not be complete, so you have to add your own soon.
    !  Create a program that reads the input and prints out final results.
    !  And, please, indent your code.

    !  Write the self-consistent field procedure in a subroutine.
    integer:: i



    open(file="results.txt", newunit=io2)
    call banner(io2)


    call input_reader(nat,nel,nbf,xyz,chrg,zeta,io2)



    call NucRep(nat,xyz,chrg,Erep,io2)

    ng=3
    allocate (exponents(ng*nbf), coefficients(ng*nbf), sab(nbf,nbf), tab(nbf,nbf), vab(nbf,nbf), packsab(nbf*(1+nbf)/2),packtab(nbf*(1+nbf)/2),packvab(nbf*(1+nbf)/2))
    allocate(eigval(nbf),eigvec(nbf,nbf),xab(nbf,nbf),Fock(nbf,nbf),cab(nbf,nbf))
    allocate(twointeg(((nbf*(nbf-1)/2+nbf)*(nbf*(nbf-1)/2+nbf)/2+(nbf*(nbf-1)/2+nbf)-1)))

    call expansion(ng, nbf, zeta, exponents, coefficients,io2)

    call oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2)

    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                           Packing Matrices                                           "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"<<<<<<<<<<<<<<<<<     Packing order: (I) S matrix, (II) T Matrix, (III) V Matrix    >>>>>>>>>>>>>>>>>>>"
    write(io2,*)"======================================================================================================="

    call packer(sab,packsab, nbf,io2)
    call packer(tab,packtab, nbf,io2)
    call packer(vab,packvab, nbf,io2)

  call orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf)
  write(*,*)
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                     Initial guess of Fock Matrix                                      "
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                             Initial guess                                      "
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

  Fock=tab+vab
  write(io2,*)
  call write_matrix(Fock,"      ============================       Fock Matrix      ============================",io2)
  write(io2,*)"      =================================================================================="

  call coeff(cab,Fock,xab,nbf,io2)

  call twoIntegrals(xyz, exponents, coefficients, twointeg, ng,nbf)


  deallocate(coefficients, exponents)

end subroutine scf_prog

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!> Below all subroutines used in the scf_prog subroutine are declared, the subroutines are listed in the same order as apearance in the scf_prog
!In order to make the program structured the beginning of each subroutine is defined by a line of ł the full name of subroutine inbetween
!the end of each subroutine is followed by a lines od $ with an "-End-" inbetweeen


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    PROGRAM TITLE     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END PROGRAM TITLE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   INPUT READER   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END INPUT READER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   NUCLEAR REPULSION ENERGY    >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
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


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NUCLEAR REPULSION ENERGY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    SLATER EXPANSION   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine expansion(ng, nbf, zeta, exponents, coefficients,io2)

  !>Declaration of global variables
  integer :: ng, nbf,io2
  real(wp) ::zeta(:)
  real(wp) :: coefficients(ng*nbf), exponents(ng*nbf)

  !>Declaration of local variables
  integer :: i


  do i=0,nbf
    call expand_slater(ng, zeta(i+1), exponents(i*ng+1:(i+1)*ng),coefficients(i*ng+1:(i+1)*ng))
  end do
  write(*,*)
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                            Slater Expansion                                           "
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)
  Write(*,*) nbf, "Slater functions transformed into ", ng*nbf ,"Gaussian functions"
  write(io2,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                            Slater Expansion                                           "
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)
  Write(io2,*) nbf, "Slater functions transformed into ", ng*nbf ,"Gaussian functions"
  write(io2,*)
end subroutine expansion

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SLATER EXPANSION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    ONE ELECTRON INTEGRALS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2)

  !>Declaration of global variables
  real(wp) :: xyz(:,:), chrg(:), coefficients(:), exponents(:), sab(:,:), tab(:,:), vab(:,:)
  integer:: nbf, ng

  !>Declaration of local variables
  integer :: i, j,io2

  write(*,*)
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                 Calculation of one electron integrals                                 "
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                 Calculation of one electron integrals                                 "
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

  do i=1, nbf

    do j=1, nbf

      call oneint(xyz,chrg,xyz(1:3,i), xyz(1:3,j),exponents(ng*(i-1)+1:ng*i),exponents(ng*(j-1)+1:ng*j),coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), sab(i,j), tab(i,j), vab(i,j))

    end do

  end do

  write(*,*)
  write(*,*)"                     –––––––   Done, all matrices saved in results.txt   –––––––"
  write(*,*)
  write(io2,*)
  write (io2,*) "======================================================================================================="
  call write_matrix(sab, "  Overlap Matrix S   ",io2)
  write (io2,*) "======================================================================================================="
  write(io2,*)
  write (io2,*) "======================================================================================================="
  call write_matrix(tab, "Kinetic Matrix T", io2)
  write (io2,*) "======================================================================================================="
  write(io2,*)
  write (io2,*) "======================================================================================================="
  call write_matrix(vab, "Nuclear Attraction Matrix V", io2)
  write (io2,*) "======================================================================================================="

end subroutine oneelint

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SLATER EXPANSION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    PACKING OF SYMMETRIC MATRICES   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine packer(matrix,packmatrix,nbf,io2)

    !>Declaration of local variables
    integer :: nbf, i,j,k,io2
    real(wp) :: matrix(:,:),packmatrix(:)
    real(wp), allocatable:: singlem(:,:)

    allocate(singlem(nbf,nbf))
    singlem=sngl(matrix)
    k=1
    j=1
    if(all(transpose(singlem)==singlem)) then
     do i=1,nbf
       j=1
       do while (j<=i)
         packmatrix(k)=matrix(i,j)
         k=k+1
         j=j+1
       end do
     end do
     write(io2,*)
     write(io2,*)"                                       Matrix succsesfully packed."
    else
     write(io2,*) "                                        Matrix is not symmetric."
    end if


    deallocate (singlem)
  end subroutine packer
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END PACKING OF SYMMETRIC MATRICES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    SYMMETRIC ORTHONORMALIZER  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf)

    !>Declaration of local variables
    integer:: stat, io2, nbf,i,j
    real(wp) :: packsab(:),sab(:,:),eigval(:),eigvec(:,:), xab(nbf,nbf)
    real(wp), allocatable :: matrix(:), halfs(:,:), check(:,:)

    allocate(matrix(nbf*(nbf+1)/2), halfs(nbf,nbf),check(nbf,nbf))

    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                                Calculation of symmetric orthonormalizer                               "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                Calculation of symmetric orthonormalizer                               "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

    ! Taking care of no overwriting of the matrix
    matrix=packsab

    ! Calculating eigenvalues and eigenfun ctions of the overlap matrix
    call solve_spev(matrix,eigval,eigvec,stat)
    !call write_vector(eigval, "Eigenvalue of SAB")




    do i=1,nbf
      do j=1,nbf

        if(i==j) then
          halfs(i,j)=1_wp/sqrt(eigval(i))
        else
          halfs(i,j)=0
        end if
      end do
    end do

  ! caLCULATING SYMMETRIC ORTHONORMALIZER
  xab=matmul(halfs,transpose(eigvec))
  xab=matmul(eigvec,xab)

      write(*,*)
       write(*,*)"                         –––––––   Symmetric orthonormalizer obtained   –––––––"
       write(io2,*)
       call write_matrix(xab,"      ====================    Symmetric orthonormalizer obtained    ====================",io2)
       write(io2,*)"      =================================================================================="

  !>Check of orthonormalizer
  !check=matmul(sab,xab)
  !check=matmul(transpose(xab),check)
  !call write_matrix(check, "Check of orthonormalizer")

  deallocate(halfs,matrix,check)
  end subroutine orthonormalizer

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SYMMETRIC ORTHONORMALIZER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW COEFFICIENTS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine coeff(cab,Fock,xab,nbf,io2)

    !>Declaration of local variables
    integer :: io2,nbf,stat,i,j,k,l, ij,kl,ijkl, m, n, final
    real(wp) :: cab(:,:),Fock(:,:),xab(:,:)
    real(wp), allocatable :: Fockprim(:,:),packFockprim(:),eigval(:),eigvec(:,:)

    allocate(Fockprim(nbf,nbf), packFockprim(nbf*(nbf+1)/2), eigval(nbf), eigvec(nbf,nbf))

    Fockprim=matmul(Fock,xab)
    Fockprim=matmul(transpose(xab),Fockprim)

    call packer(Fockprim,packFockprim,nbf,io2)



    call solve_spev(packFockprim,eigval,eigvec,stat)

    cab=matmul(xab,eigvec)
    write(*,*)
    write(*,*)
    write(*,*)"                            ––––––    Orbital coefficients obtained    –––––––"
    write(io2,*)
    write(io2,*)
    call write_matrix(cab,"         ========================    Orbital coefficients    =========================",io2)
    write(io2,*)"      =================================================================================="





  end subroutine coeff

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW COEFFICIENTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  subroutine twoIntegrals(xyz, exponents, coefficients, twointeg, ng,nbf)

    integer :: i, j, k, l, ij, kl, ijkl, nbf,final,ng
    real(wp) :: xyz(:,:), exponents(:), coefficients(:), twointeg(:)

    i=1
    j=1
    k=1
    l=1

    do i=1, nbf

       do j=1, i

         !> First indexing
          if(i>j)then
            ij=i*(i-1)/2+j
          else
            ij=j*(j-1)/2+i
          endif

           do k=1, i

             if (i==k) then
               final=j
             else
               final=k
             end if

                do l=1, final

                  !> Second indexing
                  if(l>k) then
                    kl=k*(k-1)/2+l
                  else
                    kl=l*(l-1)/2+k
                  endif

                  !> FinaL index
                  if (ij>kl) then
                    ijkl=ij*(ij-1)/2+kl
                  else
                    ijkl=kl*(kl-1)/2+ij
                  end if


                  call twoint(xyz(1:3,i), xyz(1:3,j), xyz(1:3,k), xyz(1:3,l), exponents(ng*(i-1)+1:ng*i), exponents(ng*(j-1)+1:ng*j), exponents(ng*(k-1)+1:ng*k), exponents(ng*(l-1)+1:ng*l), coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), coefficients(ng*(k-1)+1:ng*k), coefficients(ng*(l-1)+1:ng*l), twointeg(ijkl))
                  write(*,*)ijkl!, twointeg(ijkl)


              end do

            end do

         end do

      end do

      call write_matrix(twointeg, "Two electron integrals")
  end subroutine twoIntegrals

end module scf_main
