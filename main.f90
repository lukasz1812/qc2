!      ____       __
!     / /| |     / /
!    /   | | /| / /   this program was written by Łukasz Wantoch
!   / /__| |/ |/ /     During the laboratory course of WP 4 in summerterm 2021
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

    !>New Fock matrix
    real(wp),allocatable :: Fock_new(:,:)

    !>G Tensor
    real(wp), allocatable :: gabcd(:,:)

    !> Coefficients matrix
    real(wp), allocatable :: cab(:,:)

    !> Density matrix
    real(wp), allocatable :: pab(:,:)

    !>Core Hamiltonian matrix
    real(wp), allocatable :: hab(:,:)


    !> Eigenvalues and eigenvectors for LAPACK eigenvalue solver
    real(wp),allocatable :: eigval(:)
    real(wp),allocatable :: eigvec(:,:)

    !>Symmetric orhonormalizer
    real(wp),allocatable :: xab(:,:)

    !> Pointer for result file
    integer :: io2,io3

    !>time variables
    real :: start, finish, starttei,finishtei, startscf,finishscf
    !> Hartree-Fock energy
    real(wp) :: escf

    !> Hartree-Fock energy in SCF cycle
    real(wp) :: newescf

    !>End condition of SCF cycle
    real(wp):: delta

    !> Counter
    integer:: i

    !  declarations may not be complete, so you have to add your own soon.
    !  Create a program that reads the input and prints out final results.
    !  And, please, indent your code.

    !  Write the self-consistent field procedure in a subroutine.


    !>generating file with results
    open(file="results.txt", newunit=io2)

    !> Printing Banner
    call banner(io2)

    !> Reading input file
    call input_reader(nat,nel,nbf,xyz,chrg,zeta,io2)
    call cpu_time(start)


    !>Calculatin nuclear repulsion
    call NucRep(nat,xyz,chrg,Erep,io2)

    !>Parameter for choose od the Basis Sets Type STO-(ng)G
    ng=6

    !>allocate memory for used arrays
    allocate (exponents(ng*nbf), coefficients(ng*nbf), sab(nbf,nbf), tab(nbf,nbf), vab(nbf,nbf), packsab(nbf*(1+nbf)/2),packtab(nbf*(1+nbf)/2),packvab(nbf*(1+nbf)/2))
    allocate(eigval(nbf),eigvec(nbf,nbf),xab(nbf,nbf),Fock(nbf,nbf),cab(nbf,nbf),pab(nbf,nbf), hab(nbf,nbf), Fock_new(nbf,nbf),gabcd(nbf,nbf))
    allocate(twointeg(((nbf*(nbf-1)/2+nbf)*(nbf*(nbf-1)/2+nbf)/2+(nbf*(nbf-1)/2+nbf)-1)))

    !>Slater Expansion into Gaussians
    call expansion(ng, nbf, zeta, exponents, coefficients,io2)

    !>Calculation of one electron integrals
    call oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2)


    !> Text for result File
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                           Packing Matrices                                           "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"<<<<<<<<<<<<<<<<<     Packing order: (I) S matrix, (II) T Matrix, (III) V Matrix    >>>>>>>>>>>>>>>>>>>"
    write(io2,*)"======================================================================================================="

    !>Packin one electron matrices
    call packer(sab,packsab, nbf,io2)
    call packer(tab,packtab, nbf,io2)
    call packer(vab,packvab, nbf,io2)


    !> Calculating orthonormalizer with use of symmetric procedure
    call orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf)


    !>Initial Guess Annoucing (stdout+file)
    write(*,*)
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                                             Initial guess                                       "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                             Initial guess                                      "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"


    !>Setting initial Fock Matrix as a core Hamiltonian
    Fock=tab+vab
    write(*,*)
    write(*,*)"                                         ✓ Fock matrix obtained        "
    write(io2,*)
    call write_matrix(Fock,"      ============================       Fock Matrix      ============================",io2)
    write(io2,*)"      =================================================================================="

    !>calculating coefficients from the initial Fock matrix
    call coeff(cab,Fock,xab,nbf,io2)

    !>Calculating two electron integrals
    call twoIntegrals(xyz, exponents, coefficients, twointeg, ng,nbf,starttei,finishtei)

    !Calculating new density matrix
    call new_density(nel, nbf,cab,pab,io2)
    write(*,*)
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                 <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"               <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

    !HCore matrix as sum of kinetic and attraction matrices
    hab=tab+vab

    !Calculating Hartree Fock energy
    call HFenergy(nbf,escf,hab,Fock,pab)

    !Printing the initial Fock energy
    write(*,*)
    write(*,*) "                           *************************************************"
    write(*,*) "                           *  ", "initial E_{HF}=", escf, "H  *"
    write(*,*) "                           *************************************************"
    write(*,*)
    write(io2,*)
    write(io2,*) "                           *************************************************"
    write(io2,*) "                               ", "initial E_{HF}=", escf, "H"


    !Calculating the G Tensor and new Fock Matrix
    call newFock(nbf,pab,hab,Fock_new,twointeg,gabcd,io2)

    !Printing results and initial HF energy /stdout+file
    call write_matrix(gabcd, "G Tensor",io2)
    write(io2,*)
    call write_matrix(Fock_new, "Fock Matrix",io2)
    write(io2,*)
    call HFenergy(nbf,escf,hab,Fock_new,pab)
    write(io2,*) escf

    !>Starting SCF Procedure
    write(*,*)
    write(*,*)"    ==============================================================================================="
    write(*,*)"    ||                                        SCF ITERATIONS                                     ||"
    write(*,*)"    ==============================================================================================="
    write(*,*)
    write(io2,*)
    write(io2,*)
    write(io2,*)"    ==============================================================================================="
    write(io2,*)"    ||                                       SCF ITERATIONS                                      ||"
    write(io2,*)"    ==============================================================================================="
    write(io2,*)


    call cpu_time(startscf)

    !Setting stop criterion for at least one SCF Run
    delta=1

    !Setting run counter
    i=1

    open(file="SCF-results.txt", newunit=io3)
    !Definition of Convergence
    do while (delta>0.0000000000001.or.i==50)

      !Calling the steps
      call iterations(io2,cab,Fock_new,xab,nbf,nel,pab,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3)

      !Calculating energy change
      delta=escf-newescf

      !Overwriting the HF Energy
      escf=newescf

      !Counter inrease
      i=i+1

    end do

    !Printing information about end of SCF Procedure
    if (i<50) then

      write(*,*)"    ==============================================================================================="
      write(*,*)"    ||                                      SCF SUCCSESFULL                                      ||"
      write(*,*)"    ==============================================================================================="
      write(*,*)
      write(*,*)
      write(*,*)
      write(*,*) "                          ==================================================="
      write(*,*) "                          ||                 Energy Values                 ||"
      write(*,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
      write(*,*) "                          |Nuc. rep.=", Erep, "H           |"
      write(*,*) "                          |E HF     =", Escf, "H           |"
      write(*,*) "                          |.................................................|"
      write(*,*) "                          |E Tot.   =", Erep+escf, "H           |"
      write(*,*) "                          ==================================================="
      write(io2,*)"    ==============================================================================================="
      write(io2,*)"    ||                                      SCF SUCCSESFULL                                      ||"
      write(io2,*)"    ==============================================================================================="
      write(io2,*)
      write(io2,*)
      write(io2,*)
      write(io2,*) "                          ==================================================="
      write(io2,*) "                          ||                 Energy Values                 ||"
      write(io2,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
      write(io2,*) "                          |Nuc. rep.=", Erep, "H           |"
      write(io2,*) "                          |E HF     =", Escf, "H           |"
      write(io2,*) "                          |.................................................|"
      write(io2,*) "                          |E Tot.   =", Erep+escf, "H           |"
      write(io2,*) "                          ==================================================="
      write(io3,*)"    ==============================================================================================="
      write(io3,*)"    ||                                      SCF SUCCSESFULL                                      ||"
      write(io3,*)"    ==============================================================================================="
      write(io3,*)

    else
      write(*,*)"    ==============================================================================================="
      write(*,*)"    ||                                        SCF FAILED                                         ||"
      write(*,*)"    ==============================================================================================="

    end if

    call cpu_time(finishscf)
    write(io3,*) "======================================================================================================="
    write(io3,*) "SCF iterations time:", finishscf-startscf, "s"

    call cpu_time(finish)
    write(*,*)
    write(*,*)
    write(*,*)
    write(*,*)"Total computational time:",finish-start, "s"
    write(*,*)"Thereof: for tei         ", finishtei-starttei, "s"
    write(*,*)"         for SCF         ", finishscf-startscf, "s"
    write(io2,*)
    write(io2,*)
    write(io2,*)
    write(io2,*)"Total computational time:",finish-start, "s"
    write(io2,*)"Thereof: for tei         ", finishtei-starttei, "s"
    write(io2,*)"         for SCF         ", finishscf-startscf, "s"
 !deallocate(coefficients, exponents)
 deallocate (exponents, coefficients, sab, tab, vab, packsab,packtab,packvab)
 deallocate(eigval,eigvec,xab,Fock,pab, hab, Fock_new,gabcd,twointeg)
end subroutine scf_prog

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!> Below all subroutines used in the scf_prog subroutine are declared, the subroutines are listed in the same order as apearance in the scf_prog
!In order to make the program structured the beginning of each subroutine is defined by a line of ł the full name of subroutine inbetween
!the end of each subroutine is followed by a lines od $ with an "-End-" inbetweeen


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    PROGRAM TITLE     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine banner(io2)

  !>Declaration of local variables
  integer :: io2

  !>Print Banner /stdout+file
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
  write(*,*) "                       A Hartree-Fock program with use of Roothaan-Haal equations"
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
  write (io2,*) "                       A Hartree-Fock program with use of Roothaan-Haal equations"
  write (io2,*) "======================================================================================================="
  write(io2,*) ""
end subroutine banner

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END PROGRAM TITLE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   INPUT READER   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine input_reader (nat,nel,nbf,xyz,chrg,zeta,io2)      !EXERCISE 2

  !>declaration of local variables
  character(len = 100):: input
  integer :: nat, nel, nbf, io, i, j,k,  dim,io2
  real(wp), allocatable :: xyz(:,:),chrg(:), zeta(:), basis(:)

  !>Set of XYZ coordinates
  dim=3
  k=1
  j=1

  !>Input file name over stdin
  write(*,*) "Give the input file name :"
  read(*,*) input

  !>opening input file
  open(file="molecules/"//input, newunit=io)

  !Reading # of atoms, elecons and basis functions
  read(io,*) nat,nel,nbf

  !Allocate needed memory
  allocate (xyz(dim,nat), chrg(nat), basis(nat), zeta(nbf))

    !Run ober all atoms
    do i=1,nat

      !Read positions, nuclear charge and number of basis functions for each atom
      read (io,*) xyz(j:dim,i),chrg(i), basis(i)

      !Run over all basis functions
      do while (k<=sum(basis))

        !Read all slater exponens
        read(io,*) zeta(k)
        k=k+1

      end do

    end do

  close(io)


  !>variable check stdout+file
  write(*,*)
  write(*,*) "                          ==================================================="
  write(*,*) "                          ||               System parameters               ||"
  write(*,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(*,*) "                          |Number of atoms          |", nat,"          |"
  write(*,*) "                          |Number of electrons      |", nel,"          |"
  write(*,*) "                          |Number of basis functions|", nbf,"          |"
  write(*,*) "                          ==================================================="
  write(*,*)
  write(*,*)
  write(io2,*) "                          ==================================================="
  write(io2,*) "                          ||               System parameters               ||"
  write(io2,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(io2,*) "                          |Number of atoms          |", nat,"          |"
  write(io2,*) "                          |Number of electrons      |", nel,"          |"
  write(io2,*) "                          |Number of basis functions|", nbf,"          |"
  write(io2,*) "                          ==================================================="
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
subroutine NucRep(nat,xyz,chrg,Erep,io2)     !EXERCISE 3

  !>Declaration of global variables
  integer :: nat, io2,i,j
  real(wp) :: Erep, distance, partrep, diff(3)
  real(wp), allocatable :: xyz(:,:), chrg(:)

  !>Printing some info text /stdout+file
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                 Calculating nuclear repulsion energy"
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                 Calculating nuclear repulsion energy"
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)

  !Setting counters and taking care about selfcounting
  j=1
  i=j+1

  !treating a single atom case
  if (nat==1) then
    Erep=0
  else
    !looping over all atom pairs
    do while (i<=nat)

      !>Calculating nuclear repulsion energy between two attoms
      diff=xyz(1:3,j)-xyz(1:3,i)
      distance=sqrt(sum(diff**2))
      partRep=(chrg(j)*chrg(i))/distance
      !>Summing energies to get the repulsion of the system
      Erep=Erep+partRep

      j=j+1
      i=i+1

    end do

  end if


  !>Printing subroutine results /stdout+file
  write(*,*) "                         ==================================================="
  write(*,*) "                         ||            Nuclear repulsion energy           ||"
  write(*,*) "                         |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(*,*) "                         |", Erep, "H                     |"
  write(*,*) "                         ==================================================="
  write(io2,*) "                          ==================================================="
  write(io2,*) "                          ||            Nuclear repulsion energy           ||"
  write(io2,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
  write(io2,*) "                          |", Erep, "H                     |"
  write(io2,*) "                          ==================================================="


end subroutine NucRep


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NUCLEAR REPULSION ENERGY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    SLATER EXPANSION   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine expansion(ng, nbf, zeta, exponents, coefficients,io2)      !EXERCISE 4

  !>Declaration of variables
  integer :: ng, nbf,io2,i
  real(wp) ::zeta(:)
  real(wp) :: coefficients(ng*nbf), exponents(ng*nbf)

  !>Loop over all basis functions to calculate one electron integrals
  do i=0,nbf

    call expand_slater(ng, zeta(i+1), exponents(i*ng+1:(i+1)*ng),coefficients(i*ng+1:(i+1)*ng))

  end do
  !>Print the results of Slater expansion
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
subroutine oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2)  !EXERCISE 5

  !>Declaration of variables
  real(wp) :: xyz(:,:), chrg(:), coefficients(:), exponents(:), sab(:,:), tab(:,:), vab(:,:)
  integer:: nbf, ng
  integer :: i, j,io2


  !Print some procedure title /Stdout+file
  write(*,*)
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                 Calculation of one electron integrals                                 "
  write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                 Calculation of one electron integrals                                 "
  write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"


  !Loop over all electron pairs
  do i=1, nbf

    do j=1, nbf

      !CAlculating one electron integrals (overlap, kinetic enegry, e-nuc attraction)
      call oneint(xyz,chrg,xyz(1:3,i), xyz(1:3,j),exponents(ng*(i-1)+1:ng*i),exponents(ng*(j-1)+1:ng*j),coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), sab(i,j), tab(i,j), vab(i,j))

    end do

  end do

  !Printing some results /file
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

  subroutine packer(matrix,packmatrix,nbf,io2)      !EXERCISE 6

    !>Declaration of variables
    integer :: nbf, i,j,k,io2
    real(wp) :: matrix(:,:),packmatrix(:)
    real(wp), allocatable:: singlem(:,:)

    !Allocatte needed memory
    allocate(singlem(nbf,nbf))

    !>creating a matrix equivalent with single precission for symmetry check
    singlem=sngl(matrix)

    k=1
    j=1

    !>Symmetry check
    if(all(transpose(singlem)==singlem)) then

      !>Packing matrices into a vector
     do i=1,nbf
       j=1
       do while (j<=i)
         packmatrix(k)=matrix(i,j)
         k=k+1
         j=j+1
       end do
     end do

     !>Text output into file
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
 subroutine orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf)     !EXERCISE 7

    !>Declaration of variables
    integer:: stat, io2, nbf,i,j
    real(wp) :: packsab(:),sab(:,:),eigval(:),eigvec(:,:), xab(nbf,nbf)
    real(wp), allocatable :: matrix(:), halfs(:,:), check(:,:)

    !Allocation of needed memory
    allocate(matrix(nbf*(nbf+1)/2), halfs(nbf,nbf),check(nbf,nbf))

    !Print some procedure title /Stdout+file
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                                Calculation of symmetric orthonormalizer                               "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                Calculation of symmetric orthonormalizer                               "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

    ! Taking care of no overwriting of the matrix
    matrix=packsab

    ! Calculating eigenvalues and eigenfunctions of the overlap matrix
    call solve_spev(matrix,eigval,eigvec,stat)
    !call write_vector(eigval, "Eigenvalue of SAB")



    !Calculating s^(-1/2)
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
    xab=matmul(eigvec,halfs)
  xab=matmul(xab,transpose(eigvec))


      !Prinrt some success text
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
  subroutine coeff(cab,Fock,xab,nbf,io2) !Exercise 8.3

    !>Declaration of local variables
    integer :: io2,nbf,stat,i,j,k,l, ij,kl,ijkl, m, n, final
    real(wp) :: cab(:,:),Fock(:,:),xab(:,:)
    real(wp), allocatable :: Fockprim(:,:),packFockprim(:),eigval(:),eigvec(:,:)

    !Allocation of needed memory
    allocate(Fockprim(nbf,nbf), packFockprim(nbf*(nbf+1)/2), eigval(nbf), eigvec(nbf,nbf))

    !Calculation of F' matrix
    Fockprim=matmul(Fock,xab)
    Fockprim=matmul(transpose(xab),Fockprim)

    !Print the F' matrix into result File
    write(io2,*)
    call write_matrix(Fockprim,"           ====================      F'ock matrix obtained     ====================     ",io2)
    write(io2,*)"      =================================================================================="

    !Packing F' into vector
    call packer(Fockprim,packFockprim,nbf,io2)


    !Solving the eigenvalue problem of the F'ock matrix
    call solve_spev(packFockprim,eigval,eigvec,stat)

    !Calculating coefficients
    cab=matmul(xab,eigvec)

    !Printing results stdout+file
    write(*,*)
    write(*,*)"                                    ✓ Orbital coefficients obtained    "
    write(io2,*)
    write(io2,*)
    call write_matrix(cab,"         ========================    Orbital coefficients    =========================",io2)
    write(io2,*)"      =================================================================================="



    deallocate(Fockprim, packFockprim)

  end subroutine coeff

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW COEFFICIENTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW DENISTY MATRIX   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine new_density(nel, nbf,cab,pab,io2) !Exercise 8.4

    !Declaration of variables
    integer:: i, j, nel,nbf,io2
    real(wp)::cab(nbf,nbf), pab(nbf,nbf)
    real(wp), allocatable :: nocc(:,:)

    !Allocation of needed memory
    allocate(nocc(nbf,nbf))

    !Preparring the occupation matrix
    nocc=0
    do i=1,nel/2
      do j=1,nel/2
        if(i==j) then
          nocc(i,j)=2
        end if
      end do
    end do

    !CAlculation of the new density matrix
    pab=matmul(nocc,transpose(cab))
    pab=matmul(cab,pab)
    write(*,*)

    !Printing succes text /stdout+file
    write(*,*)"                                       ✓ Density matrix obtained       "
    write(*,*)
    write(io2,*)
    write(io2,*)
    call write_matrix(pab,"            ========================    Density Matrix    =========================",io2)
    write(io2,*)"      =================================================================================="

    deallocate(nocc)

  end subroutine new_density


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW DENSITY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  HF ENERGY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine HFenergy(nbf,escf,hab,Fock,pab)

    !>Declaration of local variables
    integer:: i, j,nbf
    real(wp):: escf
    real(wp):: hab(nbf,nbf),Fock(nbf,nbf),pab(nbf,nbf)
    real(wp), allocatable :: Energy(:,:)

    !Setting energy to 0
    escf=0

    !>Allocate array space
    allocate(Energy(nbf,nbf))

    !>Calculating matrix for energy calculation
    Energy=matmul((hab+Fock), pab)


    !>Calculate the E_{HF} as a half of the matrix trace
    do i=1,nbf

      do j=1,nbf

        if(i==j) then

          escf=escf+0.5_wp*(Energy(i,i))

        end if

      end do

    end do

    !>Free array space
    deallocate(Energy)

  end subroutine HFenergy


  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END HF ENERGY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
    !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  TWO ELECTRON INTEGRALS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine twoIntegrals(xyz, exponents, coefficients, twointeg, ng,nbf,starttei,finishtei)

    integer :: i, j, k, l, ij, kl, ijkl, nbf,final,ng
    real:: starttei, finishtei
    real(wp) :: xyz(:,:), exponents(:), coefficients(:), twointeg(:)

    call cpu_time(starttei)

    i=1
    j=1
    k=1
    l=1


    !Run over all different integrals
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
                  if(k>l) then
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

                  !Calculating two electron integrals
                  call twoint(xyz(1:3,i), xyz(1:3,j), xyz(1:3,k), xyz(1:3,l), exponents(ng*(i-1)+1:ng*i), exponents(ng*(j-1)+1:ng*j), exponents(ng*(k-1)+1:ng*k), exponents(ng*(l-1)+1:ng*l), coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), coefficients(ng*(k-1)+1:ng*k), coefficients(ng*(l-1)+1:ng*l), twointeg(ijkl))


              end do

            end do

         end do

      end do

      call cpu_time(finishtei)

  end subroutine twoIntegrals
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ ENSD TWO ELECTRON INTEGRALS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

      subroutine newFock(nbf,pab,hab,Fock_new,twointeg,gabcd,io2) !Exercise12.1


        !Declaration of variables
        integer:: i,j,k,l, ij, kl, ijkl, il, kj, ilkj, nbf,io2
        real(wp):: pab(:,:), hab(:,:), Fock_new(:,:),twointeg(:), gabcd(:,:)

        !Setiing new Fovck equal 0
        Fock_new=0

        !>Run over all basis functions
        do i=1, nbf
          do j=1, nbf
            do k=1, nbf
              do l=1, nbf

            !Indexing for the additive part
            !> First indexing
             if(i>j)then
               ij=i*(i-1)/2+j
             else
               ij=j*(j-1)/2+i
             endif

             !> Second indexing
             if(k>l) then
               kl=k*(k-1)/2+l
             else
               kl=l*(l-1)/2+k
             endif

             !> FinaL indexing
             if (ij>kl) then
               ijkl=ij*(ij-1)/2+kl
             else
               ijkl=kl*(kl-1)/2+ij
             end if

             !Indexing for the negative part
             !> First indexing
              if(i>l)then
                il=i*(i-1)/2+l
              else
                il=l*(l-1)/2+i
              endif
              !> Second indexing
              if(k>j) then
                kj=k*(k-1)/2+j
              else
                kj=j*(j-1)/2+k
              endif
              !> FinaL index
              if (il>kj) then
                ilkj=il*(il-1)/2+kj
              else
                ilkj=kj*(kj-1)/2+il
              end if


              !Calculating the G Tensor
             Fock_new(i,j)=Fock_new(i,j)+pab(l,k)*(twointeg(ijkl)-0.5*twointeg(ilkj))


           end do

         end do

       end do

     end do



    !Setiing G tensor equal the new Fock
    gabcd=Fock_new

    !>Adding the core Hamiltonian part
    Fock_new=Fock_new+hab

end subroutine newFock

subroutine iterations(io2,cab,Fock_new,xab,nbf,nel,pab,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3)

  integer :: i,io2,io3,nbf,nel,ng
  real(wp) :: xyz(:,:), exponents(:), coefficients(:), twointeg(:)
  real(wp) :: Erep,newescf
  real(wp) :: cab(:,:),Fock_new(:,:),xab(:,:), pab(:,:),hab(:,:), gabcd(:,:)

  !Printing iteration counter
  write(*,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)"                                         ITERATION", i
  write(*,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(*,*)
  write(io2,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)"                                         ITERATION", i
  write(io2,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io2,*)
  write(io3,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io3,*)"                                         ITERATION", i
  write(io3,*)"               –––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  write(io3,*)


  !Calculation of the transformed F matrix and new coefficients
  call coeff(cab,Fock_new,xab,nbf,io3)
  write(io2,*)"                                 ✓ New orbital coefficients obtained    "

  !Formation of new denisty matrix from obtained coefficients
  call new_density(nel, nbf,cab,pab,io3)
  write(io2,*)"                                   ✓ New density matrix obtained       "

  !Calculation of new G tensors and Fock Matrix
  call newFock(nbf,pab,hab,Fock_new,twointeg,gabcd,io3)

  !Printing some results /std+file
  write(*,*)"                                          ✓ G tensor obtained       "
  write(io2,*)"                                      ✓ New G tensor obtained       "
  call write_matrix(gabcd,"                =========================      G Tensor      ==========================",io3)
  write(io3,*)"      =================================================================================="
  write(io3,*)
  write(*,*)
  write(*,*)"                                      ✓ New Fock Matrix obtained       "
  write(io2,*)"                                      ✓ New Fock Matrix obtained       "
  call write_matrix(Fock_new,"               ========================      Fock Matrix      =========================",io3)
  write(io3,*)"      =================================================================================="
  write(io3,*)

  !Calculating new HF Energy
  call HFenergy(nbf,newescf,hab,Fock_new,pab)

  !Printing new HF energy
  write(*,*)
  write(*,*) "                           *************************************************"
  write(*,*) "                                   ", "E_{HF}=", newescf, "H                     "
  write(*,*)
  write(*,*)
  write(io2,*)
  write(io2,*) "                           *************************************************"
  write(io2,*) "                                   ", "E_{HF}=", newescf, "H                     "
  write(io2,*)
  write(io2,*)
  write(io3,*)
  write(io3,*) "                           *************************************************"
  write(io3,*) "                                   ", "E_{HF}=", newescf, "H                     "
  write(io3,*)
  write(io3,*)

end subroutine iterations



end module scf_main
