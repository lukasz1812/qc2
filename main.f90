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

    !> Input name
    integer, intent(in) :: input

    !> Case (RHF or UHF)
    character(len = 100):: input_name
    character(:), allocatable :: output_name

    !> Case (RHF or UHF)
    character(len = 3):: version

    !> System specific data
    !> Number of atoms
    integer :: nat

    !> Number of electrons for RHF
    integer :: nel

    !> Number of electrons for UHF
    integer :: nelalpha, nelbeta

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

    !>basis functions per atom
    real(wp), allocatable:: basis(:)

    !> Gaussian exponents and coeefficients of slater expansion
    real(wp),allocatable :: exponents(:)
    real(wp),allocatable :: coefficients(:)

    !> Nuclear repulsion energy
    real(wp) :: Erep

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
    integer :: io,io2,io3

    !>Control over print on screen
    integer:: ausgabe_erfolgt

    !>time variables
    real :: start, finish, starttei,finishtei, startscf,finishscf
    !> Hartree-Fock energy
    real(wp) :: escf

    !> Hartree-Fock energy in SCF cycle
    real(wp) :: newescf

    !>End condition of SCF cycle
    real(wp):: delta

    !> Counter
    integer:: i, j

    !> strings length
    integer::length

    !> Test of directory
    logical:: dirExists

    !  declarations may not be complete, so you have to add your own soon.
    !  Create a program that reads the input and prints out final results.
    !  And, please, indent your code.

    !  Write the self-consistent field procedure in a subroutine.






    !> Printing Banner
    call banner

    !>Input file name over stdin
    write(*,*) "Give the input file name :"
    read(*,*) input_name

    !Setting prefix of the output file
    length=len(trim(input_name))
    allocate(character(length) :: output_name)  ! Note the correct form
    output_name=input_name(1:length-3)

    ! Check if the directory exists first
    inquire( file=output_name//"-results"//'/.', exist=dirExists )
    write(*,*)dirExists, "dirExist"

    if (dirExists) then

    else
    CALL EXECUTE_COMMAND_LINE("mkdir "//output_name//"-results")
  end if
    !>generating file with results
    open(file=output_name//"-results/"//output_name//"-results.out", newunit=io2)

    !> Print banner in output file
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

    open(file=output_name//"-results/"//output_name//"-SCF-results.out", newunit=io3)

    !>opening input file
    open(file="molecules/"//input_name, newunit=io)

    !>Reading whether the restricted or Unrestricted version should proceed
    read(io,*) version

    close(io)


    !Starting the case distinction
    select case (version)

    case ("RHF")

      !Printing what we're gonna do /stdout+file
      write(*,*)
      write(*,*)"#######################################################################################################"
      write(*,*)"#########################                     Restricted                     ##########################"
      write(*,*)"#######################                     Hartree - Fock                     ########################"
      write(*,*)"#######################################################################################################"
      write(*,*)
      write(io2,*)
      write(io2,*)"#######################################################################################################"
      write(io2,*)"#########################                     Restricted                     ##########################"
      write(io2,*)"#######################                     Hartree - Fock                     ########################"
      write(io2,*)"#######################################################################################################"
      write(io2,*)

      !> Reading input file
      call input_reader(version,input_name,io,nat,nel,nbf,xyz,chrg,zeta,io2,basis)
      call cpu_time(start)

      ausgabe_erfolgt=0

      !>Calculating nuclear repulsion
      call NucRep(nat,xyz,chrg,Erep,io2,ausgabe_erfolgt)

      !>Parameter for choose od the Basis Sets Type STO-(ng)G
      ng=6

      !>allocate memory for used arrays
      allocate (exponents(ng*nbf), coefficients(ng*nbf))
      allocate(eigval(nbf),eigvec(nbf,nbf),xab(nbf,nbf),Fock(nbf,nbf),cab(nbf,nbf),pab(nbf,nbf), hab(nbf,nbf), Fock_new(nbf,nbf),gabcd(nbf,nbf))

      !>Slater Expansion into Gaussians
      call expansion(ng, nbf, zeta, exponents, coefficients,io2,ausgabe_erfolgt)


      !>Starting the restricted Hartree Fock Routine
      call restrictedHF(nbf,ng,xyz,chrg,coefficients,exponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,escf, ausgabe_erfolgt)


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

      !call  mulliken(nbf,pab,sab)

      !call opt_geo(nbf,ng,xyz,chrg,coefficients,exponents,basis,nat,nel,erep, finishtei, starttei, finishscf, startscf,escf, ausgabe_erfolgt,io2)

      !call opt_expo(zeta,nbf,ng,io2, nat,nel,io3,xyz, basis,erep,chrg,escf)

      !deallocate (exponents, coefficients)
      deallocate(eigval,eigvec,xab,Fock,pab, hab, Fock_new,gabcd)

    case ("UHF")

      !Printing what we're gonna do /stdout+file
      write(*,*)
      write(*,*)"#######################################################################################################"
      write(*,*)"#########################                    Unrestricted                    ##########################"
      write(*,*)"#######################                     Hartree - Fock                     ########################"
      write(*,*)"#######################################################################################################"
      write(*,*)
      write(io2,*)
      write(io2,*)"#######################################################################################################"
      write(io2,*)"#########################                    Unrestricted                    ##########################"
      write(io2,*)"#######################                     Hartree - Fock                     ########################"
      write(io2,*)"#######################################################################################################"
      write(io2,*)


      !> Reading input file
      call unrestricted_input_reader(version,input_name,io,nat,nelalpha,nelbeta,nbf,xyz,chrg,zeta,io2,basis)
      call cpu_time(start)

      ausgabe_erfolgt=0

      !>Calculating nuclear repulsion
      call NucRep(nat,xyz,chrg,Erep,io2,ausgabe_erfolgt)

      !>Parameter for choose od the Basis Sets Type STO-(ng)G
      ng=6

      !>allocate memory for used arrays
      allocate (exponents(ng*nbf), coefficients(ng*nbf))
      allocate(eigval(nbf),eigvec(nbf,nbf),xab(nbf,nbf),Fock(nbf,nbf),cab(nbf,nbf),pab(nbf,nbf), hab(nbf,nbf), Fock_new(nbf,nbf),gabcd(nbf,nbf))

      !>Slater Expansion into Gaussians
      call expansion(ng, nbf, zeta, exponents, coefficients,io2,ausgabe_erfolgt)

      !>Starting the restricted Hartree Fock Routine
      call unrestrictedHF(nbf,ng,xyz,chrg,coefficients,exponents,io2,basis,nat, nelalpha, nelbeta,erep, finishtei, starttei, finishscf, startscf,io3,escf, ausgabe_erfolgt)

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




    case default
      write(*,*)
      write(*,*) "ERROR!!"
      write(*,*)"No version chosen, check your input file and start again."

    end select

end subroutine scf_prog

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!> Below all subroutines used in the scf_prog subroutine are declared, the subroutines are listed in the same order as apearance in the scf_prog
!In order to make the program structured the beginning of each subroutine is defined by a line of ł the full name of subroutine inbetween
!the end of each subroutine is followed by a lines od $ with an "-End-" inbetweeen


!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    PROGRAM TITLE     >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine banner



  !>Print Banner /stdout
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

end subroutine banner

!!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END PROGRAM TITLE @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   INPUT READER   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine input_reader(version,input,io,nat,nel,nbf,xyz,chrg,zeta,io2,basis)     !EXERCISE 2

  !>declaration of local variables
  character(len = 100):: input
  character(len = 3):: version
  integer :: nat, nel, nbf, io, i, j,k,  dim,io2
  real(wp), allocatable :: xyz(:,:),chrg(:), zeta(:), basis(:)

  !>Set of XYZ coordinates
  dim=3
  k=1
  j=1


    !>opening input file
    open(file="molecules/"//input, newunit=io)

    !>Reading whether the restricted or Unrestricted version should proceed
    read(io,*) version

      !Reading # of atoms, elecons and basis functions
      read(io,*) nat,nel,nbf

      !Allocate needed memory
      allocate (xyz(dim,nat), chrg(nat), basis(nat), zeta(nbf))

      !Run over all atoms
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
subroutine NucRep(nat,xyz,chrg,Erep,io2, ausgabe_erfolgt)     !EXERCISE 3

  !>Declaration of global variables
  integer :: nat, io2,i,j, ausgabe_erfolgt
  real(wp) :: Erep, distance, partrep, diff(3)
  real(wp), allocatable :: xyz(:,:), chrg(:)

  Erep=0
  if(ausgabe_erfolgt==0) then
    !>Printing some info text /stdout+file
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                                 Calculating nuclear repulsion energy"
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                 Calculating nuclear repulsion energy"
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
  end if

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
      distance=sum(sqrt(diff**2))
      partRep=(chrg(j)*chrg(i))/distance
      !>Summing energies to get the repulsion of the system
      Erep=Erep+partRep


      j=j+1
      i=i+1

    end do

  end if

if(ausgabe_erfolgt==0) then
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
else
end if


end subroutine NucRep


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NUCLEAR REPULSION ENERGY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    SLATER EXPANSION   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine expansion(ng, nbf, zeta, exponents, coefficients,io2,ausgabe_erfolgt)      !EXERCISE 4

  !>Declaration of variables
  integer :: ng, nbf,io2,i, ausgabe_erfolgt
  real(wp) ::zeta(:)
  real(wp) :: coefficients(:), exponents(:)


coefficients=0
exponents=0

  !>Loop over all basis functions to calculate one electron integrals
  do i=0,nbf-1

    call expand_slater(ng, zeta(i+1), exponents(i*ng+1:(i+1)*ng),coefficients(i*ng+1:(i+1)*ng))

  end do


  if(ausgabe_erfolgt==0) then
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

end if

end subroutine expansion

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SLATER EXPANSION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    RESTRICTED FOCK >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine restrictedHF(nbf,ng,xyz,chrg,coefficients,exponents,io2,basis,nat,nel,erep, finishtei, starttei, finishscf, startscf,io3,escf, ausgabe_Erfolgt)


   integer :: nbf, nel,nat, io2, ng, io3, i, ausgabe_erfolgt
   real(wp) :: delta, Erep, escf, newescf, evenergy
   real :: finishtei, starttei,finishscf, startscf
   real (wp):: xyz(:,:), chrg(:), coefficients(:), exponents(:)
   real(wp), allocatable:: packsab(:), packtab(:), packvab(:), sab(:,:), tab(:,:), vab(:,:), xab(:,:),hab(:,:), Fock(:,:), Fock_new(:,:)
   real(wp),allocatable :: eigval(:), cab(:,:), pab(:,:), gabcd(:,:), twointeg(:), aufpunkt(:,:)
   real(wp),allocatable :: eigvec(:,:), basis(:)

   allocate(packsab(nbf*(1+nbf)/2),packtab(nbf*(1+nbf)/2),packvab(nbf*(1+nbf)/2),xab(nbf,nbf), eigval(nbf), eigvec(nbf,nbf), Fock(nbf,nbf),Fock_new(nbf,nbf))
   allocate(hab(nbf,nbf), sab(nbf,nbf), tab(nbf,nbf), vab(nbf,nbf), cab(nbf,nbf), pab(nbf,nbf), gabcd(nbf,nbf))
   allocate(twointeg(((nbf*(nbf-1)/2+nbf)*(nbf*(nbf-1)/2+nbf)/2+(nbf*(nbf-1)/2+nbf)-1)))

    !>Calculation of one electron integrals
    call oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2,basis,nat,aufpunkt,ausgabe_Erfolgt)


    if(ausgabe_erfolgt==0) then
    !> Text for result File
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                           Packing Matrices                                           "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"<<<<<<<<<<<<<<<<<     Packing order: (I) S matrix, (II) T Matrix, (III) V Matrix    >>>>>>>>>>>>>>>>>>>"
    write(io2,*)"======================================================================================================="
  end if

    !>Packin one electron matrices
    call packer(sab,packsab, nbf,io2,ausgabe_erfolgt)
    call packer(tab,packtab, nbf,io2,ausgabe_erfolgt)
    call packer(vab,packvab, nbf,io2,ausgabe_erfolgt)


    !> Calculating orthonormalizer with use of symmetric procedure
    call orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf,ausgabe_Erfolgt)


    if(ausgabe_erfolgt==0) then
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
  end if

    !>Setting initial Fock Matrix as a core Hamiltonian
    Fock=tab+vab

    if(ausgabe_erfolgt==0) then
    write(*,*)
    write(*,*)"                                         ✓ Fock matrix obtained        "
    write(io2,*)
    call write_matrix(Fock,"      ============================       Fock Matrix      ============================",io2)
    write(io2,*)"      =================================================================================="
    end if

    !>calculating coefficients from the initial Fock matrix
    call coeff(cab,Fock,xab,nbf,io2,ausgabe_erfolgt)

    !>Calculating two electron integrals
    call twoIntegrals(aufpunkt, exponents, coefficients, twointeg, ng,nbf,starttei,finishtei)


    if(ausgabe_erfolgt==0) then
    !Calculating new density matrix
    call new_density(nel, nbf,cab,pab,io2,ausgabe_erfolgt)
    write(*,*)
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                 <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"               <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

    end if

    !HCore matrix as sum of kinetic and attraction matrices
    hab=tab+vab

    !Calculating Hartree Fock energy
    call HFenergy(nbf,escf,hab,Fock,pab)

    if(ausgabe_erfolgt==0) then
    !Printing the initial Fock energy
    write(*,*)
    write(*,*) "                           *************************************************"
    write(*,*) "                           *  ", "initial E_{HF}=", escf, "H  *"
    write(*,*) "                           *************************************************"
    write(*,*)
    write(io2,*)
    write(io2,*) "                           *************************************************"
    write(io2,*) "                               ", "initial E_{HF}=", escf, "H"

    end if

    !Calculating the G Tensor and new Fock Matrix
    call newFock(nbf,pab,hab,Fock_new,twointeg,gabcd,io2)

    if(ausgabe_erfolgt==0) then
    !Printing results and initial HF energy /stdout+file
    call write_matrix(gabcd,"           ===========================      G Tensor       ===========================     ",io2)
    write(io2,*)"      =================================================================================="
    write(io2,*)
    call write_matrix(Fock_new,"           ==========================       Fock Matrix      ==========================     ",io2)
    write(io2,*)"      =================================================================================="
    write(io2,*)
    end if


    call HFenergy(nbf,escf,hab,Fock_new,pab)

    if(ausgabe_erfolgt==0) then
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
    end if


    call cpu_time(startscf)

    !Setting stop criterion for at least one SCF Run
    delta=1

    !Setting run counter
    i=1


    !Definition of Convergence
    do while (delta>0.0000000000001.and.i<50)

      !Calling the steps
      call iterations(io2,cab,Fock_new,xab,nbf,nel,pab,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3,ausgabe_erfolgt)

      !Calculating energy change
      delta=escf-newescf

      !Overwriting the HF Energy
      escf=newescf
      evenergy=(Erep+escf)*1.042*10.**(-21)

      !Counter increase
      i=i+1

    end do

    !Printing information about end of SCF Procedure
    if (i<50) then
      if(ausgabe_Erfolgt==0) then

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
        write(*,*) "                          |         =", (Erep+escf)*27.2114, "eV          |"
        write(*,1)"                           |         =  ", evenergy,"     kcal        |"
        1 format (a,e20.14,a)
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
        write(io2,*) "                          |         =", (Erep+escf)*27.2114, "eV          |"
        write(io2,1)"                           |         =  ", evenergy,"     kcal        |"
        write(io2,*) "                          ==================================================="
        write(io3,*)"    ==============================================================================================="
        write(io3,*)"    ||                                      SCF SUCCSESFULL                                      ||"
        write(io3,*)"    ==============================================================================================="
        write(io3,*)
      end if
    else
      if(ausgabe_erfolgt==0) then
        write(*,*)"    ==============================================================================================="
        write(*,*)"    ||                                        SCF FAILED                                         ||"
        write(*,*)"    ==============================================================================================="
      end if
    end if

    if(ausgabe_erfolgt==0) then
      call cpu_time(finishscf)
      write(io3,*) "======================================================================================================="
      write(io3,*) "SCF iterations time:", finishscf-startscf, "s"
    end if
    escf=escf+erep
    write(*,*)"final SCF energy",escf
    call mulliken(nat, nel, nbf,pab,sab, chrg, basis)

  deallocate(packsab,packtab,packvab,xab,eigval,eigvec,Fock,Fock_new)
  deallocate(hab,sab,tab,vab,cab,pab,gabcd, twointeg)

  end subroutine restrictedHF
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END RESTRICTED FOCK @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    ONE ELECTRON INTEGRALS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
subroutine oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2, basis,nat,aufpunkt,ausgabe_erfolgt)  !EXERCISE 5

  !>Declaration of variables
  real(wp) :: xyz(:,:), chrg(:), coefficients(:), exponents(:), sab(:,:), tab(:,:), vab(:,:),basis(:)
  real(wp), allocatable::aufpunkt(:,:)
  integer:: nbf, ng,nat, ausgabe_erfolgt
  integer :: i, j,io2,k

!Allocate needed memory
  allocate(aufpunkt(3,nbf))
  k=0
  j=1

  sab=0
  tab=0
  vab=0
  !Run over all atoms
  do i=1,nat

    !Set counter of basis functions per atom
    k=k+basis(i)
    do while(j<=k)

      !Set aufpunkt for each orbital
      aufpunkt(1:3,j)=xyz(1:3,i)

      j=j+1
    end do
  end do

  if(ausgabe_erfolgt==0) then
    !Print some procedure title /Stdout+file
    write(*,*)
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(*,*)"                                 Calculation of one electron integrals                                 "
    write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    write(io2,*)"                                 Calculation of one electron integrals                                 "
    write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
  end if

  !Loop over all electron pairs
  do i=1, nbf

    do j=1, nbf

      !CAlculating one electron integrals (overlap, kinetic enegry, e-nuc attraction)
      call oneint(xyz,chrg,aufpunkt(1:3,i), aufpunkt(1:3,j),exponents(ng*(i-1)+1:ng*i),exponents(ng*(j-1)+1:ng*j),coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), sab(i,j), tab(i,j), vab(i,j))
    end do

  end do

  if(ausgabe_erfolgt==0) then
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
  end if
end subroutine oneelint

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SLATER EXPANSION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    PACKING OF SYMMETRIC MATRICES   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine packer(matrix,packmatrix,nbf,io2,ausgabe_erfolgt)      !EXERCISE 6

    !>Declaration of variables
    integer :: nbf, i,j,k,io2,ausgabe_erfolgt
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

     if(ausgabe_erfolgt==0) then
       !>Text output into file
       write(io2,*)
       write(io2,*)"                                       Matrix succsesfully packed."
     end if
   else
     if(ausgabe_erfolgt==0) then
       write(io2,*) "                                        Matrix is not symmetric."
     end if
   end if

    deallocate (singlem)
  end subroutine packer
  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END PACKING OF SYMMETRIC MATRICES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<    SYMMETRIC ORTHONORMALIZER  >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
 subroutine orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf, ausgabe_erfolgt)     !EXERCISE 7

    !>Declaration of variables
    integer:: stat, io2, nbf,i,j, ausgabe_Erfolgt
    real(wp) :: packsab(:),sab(:,:),eigval(:),eigvec(:,:), xab(nbf,nbf)
    real(wp), allocatable :: matrix(:), halfs(:,:), check(:,:)

    !Allocation of needed memory
    allocate(matrix(nbf*(nbf+1)/2), halfs(nbf,nbf),check(nbf,nbf))


    if(ausgabe_erfolgt==0) then
      !Print some procedure title /Stdout+file
      write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
      write(*,*)"                                Calculation of symmetric orthonormalizer                               "
      write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
      write(io2,*)
      write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
      write(io2,*)"                                Calculation of symmetric orthonormalizer                               "
      write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
    end if

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


    if(ausgabe_erfolgt==0) then
      !Print some success text
      write(*,*)
       write(*,*)"                         –––––––   Symmetric orthonormalizer obtained   –––––––"
       write(io2,*)
       call write_matrix(xab,"      ====================    Symmetric orthonormalizer obtained    ====================",io2)
       write(io2,*)"      =================================================================================="
    end if

  !>Check of orthonormalizer
  !check=matmul(sab,xab)
  !check=matmul(transpose(xab),check)
  !call write_matrix(check, "Check of orthonormalizer")

  deallocate(halfs,matrix,check)
  end subroutine orthonormalizer

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SYMMETRIC ORTHONORMALIZER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW COEFFICIENTS   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine coeff(cab,Fock,xab,nbf,io2,ausgabe_erfolgt) !Exercise 8.3

    !>Declaration of local variables
    integer :: io2,nbf,stat,i,j,k,l, ij,kl,ijkl, m, n, final,ausgabe_erfolgt
    real(wp) :: cab(:,:),Fock(:,:),xab(:,:)
    real(wp), allocatable :: Fockprim(:,:),packFockprim(:),eigval(:),eigvec(:,:)

    !Allocation of needed memory
    allocate(Fockprim(nbf,nbf), packFockprim(nbf*(nbf+1)/2), eigval(nbf), eigvec(nbf,nbf))

    !Calculation of F' matrix
    Fockprim=matmul(Fock,xab)
    Fockprim=matmul(transpose(xab),Fockprim)

  if(ausgabe_erfolgt==0) then
    !Print the F' matrix into result File
    write(io2,*)
    call write_matrix(Fockprim,"           ====================      F'ock matrix obtained     ====================     ",io2)
    write(io2,*)"      =================================================================================="
  end if

    !Packing F' into vector
    call packer(Fockprim,packFockprim,nbf,io2,ausgabe_erfolgt)


    !Solving the eigenvalue problem of the F'ock matrix
    call solve_spev(packFockprim,eigval,eigvec,stat)

    !Calculating coefficients
    cab=matmul(xab,eigvec)

    if(ausgabe_erfolgt==0) then
      !Printing results stdout+file
      write(*,*)
      write(*,*)"                                    ✓ Orbital coefficients obtained    "
      write(io2,*)
      write(io2,*)
      call write_matrix(cab,"         ========================    Orbital coefficients    =========================",io2)
      write(io2,*)"      =================================================================================="
    end if


    deallocate(Fockprim, packFockprim)

  end subroutine coeff

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW COEFFICIENTS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW DENISTY MATRIX   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine new_density(nel, nbf,cab,pab,io2,ausgabe_erfolgt) !Exercise 8.4

    !Declaration of variables
    integer:: i, j, nel,nbf,io2,ausgabe_erfolgt
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


    if(ausgabe_erfolgt==0) then
      !Printing succes text /stdout+file
      write(*,*)
      write(*,*)"                                       ✓ Density matrix obtained       "
      write(*,*)
      write(io2,*)
      write(io2,*)
      call write_matrix(pab,"            ========================    Density Matrix    =========================",io2)
      write(io2,*)"      =================================================================================="
    end if

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

  subroutine twoIntegrals(aufpunkt, exponents, coefficients, twointeg, ng,nbf,starttei,finishtei)

    integer :: i, j, k, l, ij, kl, ijkl, nbf,final,ng
    real:: starttei, finishtei
    real(wp) :: aufpunkt(:,:), exponents(:), coefficients(:), twointeg(:)

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
                  call twoint(aufpunkt(1:3,i), aufpunkt(1:3,j), aufpunkt(1:3,k), aufpunkt(1:3,l), exponents(ng*(i-1)+1:ng*i), exponents(ng*(j-1)+1:ng*j), exponents(ng*(k-1)+1:ng*k), exponents(ng*(l-1)+1:ng*l), coefficients(ng*(i-1)+1:ng*i),coefficients(ng*(j-1)+1:ng*j), coefficients(ng*(k-1)+1:ng*k), coefficients(ng*(l-1)+1:ng*l), twointeg(ijkl))


              end do

            end do

         end do

      end do


      call cpu_time(finishtei)

  end subroutine twoIntegrals
      !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END TWO ELECTRON INTEGRALS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
      !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
      !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW FOCK MATRIX >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

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
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW FOCK @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SCF ITERATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine iterations(io2,cab,Fock_new,xab,nbf,nel,pab,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3,ausgabe_erfolgt)

  integer :: i,io2,io3,nbf,nel,ng, ausgabe_erfolgt
  real(wp) :: xyz(:,:), exponents(:), coefficients(:), twointeg(:)
  real(wp) :: Erep,newescf
  real(wp) :: cab(:,:),Fock_new(:,:),xab(:,:), pab(:,:),hab(:,:), gabcd(:,:)

  if (ausgabe_erfolgt==0) then
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
  end if

  !Calculation of the transformed F matrix and new coefficients
  call coeff(cab,Fock_new,xab,nbf,io3,ausgabe_erfolgt)

  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                 ✓ New orbital coefficients obtained    "
  end if

  !Formation of new denisty matrix from obtained coefficients
  call new_density(nel, nbf,cab,pab,io3,ausgabe_erfolgt)

  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                   ✓ New density matrix obtained       "
  end if

  !Calculation of new G tensors and Fock Matrix
  call newFock(nbf,pab,hab,Fock_new,twointeg,gabcd,io3)

  if(ausgabe_erfolgt==0) then
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
  end if


  !Calculating new HF Energy
  call HFenergy(nbf,newescf,hab,Fock_new,pab)

  if(ausgabe_erfolgt==0) then
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
  end if

end subroutine iterations

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SCF ITERATIONS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MULLIKEN CHARGES >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine mulliken(nat,nel, nbf,pab,sab, chrg, basis)

    !>Declaration of local variables
    integer::i,j,nbf, begin, finish, nat, nel
    real(wp)::pab(:,:), sab(:,:), mlkn(nbf,nbf), chrg(:), basis(:), diag(nbf)
    real(wp), allocatable :: Mlkn_charge(:)

    !>allocate needed memory
    allocate(Mlkn_charge(nbf))

    !>Matrix multiplikation for Mulliken analysis
    mlkn=matmul(pab,sab)



    do i=1, nbf
      do j=1,nbf
        if(i==j) then
          diag(i)=mlkn(i,j)
        end if
   end do

 end do


  begin=1
  finish=0
  Mlkn_charge=chrg

 J=0
  !Mulliken Population Calculation
    do i=1,nat

       begin=finish+1

       finish=finish+(basis(i))


       j=begin
       do while(j<=finish)

         Mlkn_charge(i)=Mlkn_charge(i)-diag(j)


         j=j+1
       end do

     end do

call write_vector(Mlkn_charge)
end subroutine mulliken


!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END MULLIKEN CHARGES @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< START GEOMETRY OPTIMIZATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

   subroutine opt_geo(nbf,ng,xyz,chrg,coefficients,exponents,basis,nat,nel,erep, finishtei, starttei, finishscf, startscf,escf,ausgabe_erfolgt,io2)
     !>Exercise 15 && 16

     !>Declaration of local variables
     integer :: nbf, nel,nat, io2, ng, ausgabe_erfolgt, i
     real(wp):: Erep, plus, erepminus, der,escf, plusescf, minusescf, eplus, eminus, energy, opt_start, opt_end
     real :: finishtei, starttei,finishscf, startscf, delta
     real (wp):: xyz(:,:),  coefficients(:), exponents(:)
     real(wp), allocatable:: basis(:),xyzplus(:,:), xyzminus(:,:), geo_gradient(:,:),chrg(:)

     !>allocate needed memory
     allocate(xyzplus(3,nat), xyzminus(3,nat), geo_gradient(3,nat))

     call cpu_time(opt_start)

     !> Some stdout+file outpt here
     write (*,*)
     write (*,*) "#######################################################################################################"
     write(*,*)
     write(*,*)"    ==============================================================================================="
     write(*,*)"    ||                                    OPTIMIZATION RUN                                       ||"
     write(*,*)"    ==============================================================================================="
     write(*,*)
     write (io2,*)
     write (io2,*) "#######################################################################################################"
     write(io2,*)
     write(io2,*)"    ==============================================================================================="
     write(io2,*)"    ||                                   OPTIMIZATION RUN                                        ||"
     write(io2,*)"    ==============================================================================================="
     write(io2,*)
     write(*,*)"           ========================     GRADIENT PREPARATION        ========================="
     write(io2,*)"           ========================     GRADIENT PREPARATION        ========================="
     write(io2,*)

     !>chose of proper gradient vector
     !Loop over all atoms
      do i=1, nat
        !Set geometry change as 0.1 of the start positions
        geo_gradient(1:3,i)=0.1*xyz(1:3,i)
     end do

     !>Moving back and forth
     xyzplus=xyz+geo_gradient
     xyzminus=xyz-geo_gradient


     !>Reducing screen input
    ausgabe_erfolgt=1

    !>Calculating new values of repulsion and HF Energy in further step
     call NucRep(nat,xyzplus,chrg,plus,io2,ausgabe_erfolgt)
     call restrictedHF(nbf,ng,xyzplus,chrg,coefficients,exponents,io2,basis,nat,nel,plus, finishtei, starttei, finishscf, startscf,io2,plusescf,ausgabe_erfolgt)

    !>Calculating new values of repulsion and HF Energy in backward step
    call NucRep(nat,xyzminus,chrg,Erepminus,io2,ausgabe_erfolgt)
     call restrictedHF(nbf,ng,xyzminus,chrg,coefficients,exponents,io2,basis,nat,nel,erepminus, finishtei, starttei, finishscf, startscf,io2, minusescf,ausgabe_erfolgt)

     !> Calculating al total energies
     energy=erep+escf
     eplus=plus+plusescf
     eminus=erepminus+minusescf

       i=1
       !Chosing the optimal start gradient
       do while(energy<eplus.and.energy<eminus)

         !>Printing the actual step stdout+file
         write(*,*)
         write(*,*)"           -----------------------   PREPARATION RUN",i,"   -----------------------"
         write(io2,*)"           -----------------------   PREPARATION RUN",i,"   -----------------------"
         write(io2,*)

         !> Restart of energies
         plus=0
         erepminus=0
         plusescf=0
         minusescf=0

         !Setting new gradient value
         geo_gradient=0.5*geo_gradient

         !>moving back and forth
         xyzplus=xyz+geo_gradient
         xyzminus=xyz-geo_gradient

          !>Calculating new values of repulsion and HF Energy in further step
         call NucRep(nat,xyzplus,chrg,plus,io2,ausgabe_erfolgt)
         call restrictedHF(nbf,ng,xyzplus,chrg,coefficients,exponents,io2,basis,nat,nel,plus, finishtei, starttei, finishscf, startscf,io2,plusescf,ausgabe_erfolgt)

        !>Calculating new values of repulsion and HF Energy in further step
         call NucRep(nat,xyzminus,chrg,Erepminus,io2, ausgabe_erfolgt)
         call restrictedHF(nbf,ng,xyzminus,chrg,coefficients,exponents,io2,basis,nat,nel,erepminus, finishtei, starttei, finishscf, startscf,io2, minusescf,ausgabe_erfolgt)

         !>Counter increase
         i=i+1

         !> Calculating total energies
         eplus=plus+plusescf
         eminus=erepminus+minusescf

       end do

       !> success information
       write(*,*)
       write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
       write(*,*)"            <<<<<<<<<<<<<<      Gradient preparation ended succesfully      >>>>>>>>>>>>>>"
       write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
       write(*,*)
       call write_matrix(geo_gradient,"               ========================    Gradient    =========================")
       write(*,*)"      =================================================================================="
       write(io2,*)
       write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
       write(io2,*)"            <<<<<<<<<<<<<<      Gradient preparation ended succesfully      >>>>>>>>>>>>>>"
       write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
       write(io2,*)
       call write_matrix(geo_gradient,"               ========================    Gradient    =========================",io2)
       write(io2,*)"      =================================================================================="

       der=(plusescf+plus-minusescf-erepminus)/2*sum(sqrt(geo_gradient**2))

       !>Taking care the loop will run at least one time
       delta=1

       !>Counter reset
       i=1
       !>Starting the steepest descent procedure
          do while(delta>0.0001.and.i<100)

            !Setting new values for further and backward step
            xyzminus=xyz-(0.5**i)*geo_gradient
            xyzplus=xyz+(0.5**i)*geo_gradient

            !CAlculating new energies values
            call NucRep(nat,xyzplus,chrg,plus,io2,ausgabe_erfolgt)
            call restrictedHF(nbf,ng,xyzplus,chrg,coefficients,exponents,io2,basis,nat,nel,plus, finishtei, starttei, finishscf, startscf,io2,plusescf,ausgabe_erfolgt)
            call NucRep(nat,xyzminus,chrg,Erepminus,io2, ausgabe_erfolgt)
            call restrictedHF(nbf,ng,xyzminus,chrg,coefficients,exponents,io2,basis,nat,nel,erepminus, finishtei, starttei, finishscf, startscf,io2, minusescf,ausgabe_erfolgt)
            eminus=minusescf+erepminus
            eplus=plus+plusescf
            der=(eplus-eminus)/2*sum(sqrt(geo_gradient**2))


            !>Control movement on the PES
              if(ENERGY>eminus)then

                !Calculating convergence criterion
                delta=energy-eminus
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                   ", "E_{HF}=", Eminus, "H                     "
                write(*,*) "                                       ", "∆=", delta, "H                     "
                write(*,*)
                write(io2,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                   ", "E_{HF}=", Eminus, "H                     "
                write(io2,*) "                                       ", "∆=", delta, "H                     "
                write(io2,*)

                !Setting new Values
                Energy=eminus
                xyz=xyzminus
                i=i+1
              elseif (Energy>eplus) then

                !Calculating convergence criterion
                delta=energy-eplus
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                   ", "E_{HF}=", Eplus, "H   plus             "
                write(*,*) "                                       ", "∆=", delta, "H                     "
                write(*,*)
                write(io2,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                   ", "E_{HF}=", Eplus, "H                   "
                write(io2,*) "                                            ∆=", delta, "H                     "
                write(io2,*)


                Energy=eplus
                xyz=xyzplus
                i=i+1
              else
                write(io2,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                             run failed                     "
                write(io2,*)
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                             run failed                     "
                write(*,*)
                delta=1
                i=i+1
          end if


        end do

        call write_matrix(xyz,"            ========================  Final geometry    =========================")
        write(*,*)"      =================================================================================="

        call write_matrix(xyz,"            ========================  Final geometry    =========================",io2)
        write(io2,*)"      =================================================================================="

        call cpu_time(opt_end)

        !>Output over geometry optimization run
        write(*,*) "Geometry optimization done in:", opt_end-opt_start, "s"
        write(*,*)
        write(io2,*) "Geometry optimization done in:", opt_end-opt_start, "s"
        write(io2,*)


  end subroutine opt_geo

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END GEOMETRY OPTIMIZATION @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
!<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<< START EXPONENTS OPTIMIZATION >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

  subroutine opt_expo(zeta,nbf,ng,io2, nat,nel,io3,xyz, basis,erep,chrg,escf)

    !>Declaration of local variables
    integer ::  nbf, ng,io2, nat, nel,io3, ausgabe_erfolgt,i
    real(wp):: zeta(:),xyz(:,:), erep,chrg(:), minusescf, plusescf,escf, delta
    real::finishtei, starttei,finishscf, startscf
    real (wp), allocatable:: expo_gradient(:), zeta1(:), zeta2(:), plusexponents(:), pluscoefficients(:), minuexponents(:), minucoefficients(:),basis(:)


    !>allocate needed memory
    allocate(expo_gradient(nbf), zeta1(nbf),plusexponents(ng*nbf), pluscoefficients(ng*nbf))
                         allocate(zeta2(nbf),minuexponents(ng*nbf), minucoefficients(ng*nbf))

    ausgabe_erfolgt=1

    !>Setting starting value for exponent gradient
    expo_gradient=0.5*zeta

    !>Increasing and decreasing Slater exponents
    zeta1=zeta-expo_gradient
    zeta2=zeta+expo_gradient



      !>Expand Slater to Gaussians and calculate RHF Energy
      call expansion(ng, nbf, zeta1, minuexponents, minucoefficients,io2,ausgabe_erfolgt)
      call restrictedHF(nbf,ng,xyz,chrg,minucoefficients,minuexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,minusescf, ausgabe_erfolgt)

      zeta1=zeta+expo_gradient

      !>Expand Slater to Gaussians and calculate RHF Energy
      call expansion(ng, nbf, zeta2, plusexponents, pluscoefficients,io2,ausgabe_erfolgt)
      call restrictedHF(nbf,ng,xyz,chrg,pluscoefficients,plusexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,plusescf, ausgabe_erfolgt)


      i=1
      do while((escf<plusescf.and.escf<minusescf))

        !>Printing the actual step stdout+file
        write(*,*)
        write(*,*)"           -----------------------   PREPARATION RUN",i,"   -----------------------"
        write(io2,*)"           -----------------------   PREPARATION RUN",i,"   -----------------------"
        write(io2,*)

        !>Setting starting value for exponent gradient
        expo_gradient=0.5*expo_gradient

        minuexponents=0
        minucoefficients=0
        minusescf=0
        plusexponents=0
        pluscoefficients=0
        plusescf=0

        !>Increasing and decreasing Slater exponents
        zeta1=zeta-expo_gradient
        zeta2=zeta+expo_gradient


          !>Expand Slater to Gaussians and calculate RHF Energy
         call expansion(ng, nbf, zeta1, minuexponents, minucoefficients,io2,ausgabe_erfolgt)
          call restrictedHF(nbf,ng,xyz,chrg,minucoefficients,minuexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,minusescf, ausgabe_erfolgt)

          !>Expand Slater to Gaussians and calculate RHF Energy
          call expansion(ng, nbf, zeta2, plusexponents, pluscoefficients,io2,ausgabe_erfolgt)
         call restrictedHF(nbf,ng,xyz,chrg,pluscoefficients,plusexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,plusescf, ausgabe_erfolgt)
         i=i+1
        end do

        !>Taking care the loop will run at least one time
        delta=1

        !>Counter reset
        i=1
        !>Starting the steepest descent procedure
           do while(delta>0.0001)

             !Setting new values for further and backward step
             zeta1=zeta-(0.5**i)*expo_gradient
             zeta2=zeta+(0.5**i)*expo_gradient

             !>Expand Slater to Gaussians and calculate RHF Energy
            call expansion(ng, nbf, zeta1, minuexponents, minucoefficients,io2,ausgabe_erfolgt)
             call restrictedHF(nbf,ng,xyz,chrg,minucoefficients,minuexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,minusescf, ausgabe_erfolgt)

             !>Expand Slater to Gaussians and calculate RHF Energy
             call expansion(ng, nbf, zeta2, plusexponents, pluscoefficients,io2,ausgabe_erfolgt)
            call restrictedHF(nbf,ng,xyz,chrg,pluscoefficients,plusexponents,io2,basis,nat, nel,erep, finishtei, starttei, finishscf, startscf,io3,plusescf, ausgabe_erfolgt)

            !>Control movement on the PES
              if(escf>minusescf)then

                !Calculating convergence criterion
                delta=escf-minusescf
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                   ", "E_{HF}=", minusescf, "H                     "
                write(*,*) "                                       ", "∆=", delta, "H                     "
                write(io2,*)
                write(*,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                   ", "E_{HF}=", minusescf, "H                     "
                write(io2,*) "                                       ", "∆=", delta, "H                     "
                write(io2,*)

                !Setting new Values
                escf=minusescf
                zeta=zeta1
                i=i+1
              elseif (Escf>plusescf) then

                !Calculating convergence criterion
                delta=escf-plusescf
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                   ", "E_{HF}=", plusescf, "H                     "
                write(*,*) "                                       ", "∆=", delta, "H                     "
                write(io2,*)
                write(*,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                   ", "E_{HF}=", plusescf, "H                     "
                write(io2,*) "                                       ", "∆=", delta, "H                     "
                write(io2,*)

                !Setting new Values
                escf=plusescf
                zeta=zeta2
                i=i+1

              else
                write(io2,*)
                write(io2,*)"                                          Optimization run", i
                write(io2,*) "                           *************************************************"
                write(io2,*) "                                             run failed                     "
                write(io2,*)
                write(*,*)
                write(*,*)"                                          Optimization run", i
                write(*,*) "                           *************************************************"
                write(*,*) "                                             run failed                     "
                write(*,*)
                delta=1
                i=i+1
          end if
        end do
        write(io2,*)
        call write_vector(zeta,"           ========================  Final exponents     =========================",io2)
        write(io2,*)"      =================================================================================="
        call write_vector(zeta,"           ========================  Final exponents     =========================")
        write(*,*)"      =================================================================================="




  end subroutine opt_expo

  !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!              SUBROUTINES FOR UNRESTRICTED CASE              !!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   INPUT READER   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  subroutine unrestricted_input_reader(version,input,io,nat,nelalpha,nelbeta,nbf,xyz,chrg,zeta,io2,basis)

    !>declaration of local variables
    character(len = 100):: input
    character(len = 3):: version
    integer :: nat, nelalpha,nelbeta, nbf, io, i, j,k,  dim,io2, nelchange
    real(wp), allocatable :: xyz(:,:),chrg(:), zeta(:), basis(:)

    !>Set of XYZ coordinates
    dim=3
    k=1
    j=1


      !>opening input file
      open(file="molecules/"//input, newunit=io)

      !>Reading whether the  or Unrestricted version should proceed
      read(io,*) version

        !Reading # of atoms, elecons and basis functions
        read(io,*) nat,nelalpha,nelbeta,nbf

        !Allocate needed memory
        allocate (xyz(dim,nat), chrg(nat), basis(nat), zeta(nbf))

        !Run over all atoms
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

      !Check id number of electrons is correct
      if(nelbeta>nelalpha) then

        write(io2,*)
        write(io2,*)"ERROR! Number of α electrons is lower than β! The values have been switched."

        !Set values as default
        nelchange=nelalpha
        nelalpha=nelbeta
        nelbeta=nelchange

      end if


      !>variable check stdout+file
      write(*,*)
      write(*,*) "                          ==================================================="
      write(*,*) "                          ||               System parameters               ||"
      write(*,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
      write(*,*) "                          |Number of atoms          |", nat,"          |"
      write(*,*) "                          |Number of α electrons    |", nelalpha,"          |"
      write(*,*) "                          |Number of β electrons    |", nelbeta,"          |"
      write(*,*) "                          |Number of basis functions|", nbf,"          |"
      write(*,*) "                          ==================================================="
      write(*,*)
      write(*,*)
      write(io2,*) "                          ==================================================="
      write(io2,*) "                          ||               System parameters               ||"
      write(io2,*) "                          |–––––––––––––––––––––––––––––––––––––––––––––––––|"
      write(io2,*) "                          |Number of atoms          |", nat,"          |"
      write(io2,*) "                          |Number of α electrons    |", nelalpha,"          |"
      write(io2,*) "                          |Number of β electrons    |", nelbeta,"          |"
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


  end subroutine unrestricted_input_reader

  !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END INPUT READER @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  subroutine unrestrictedHF(nbf,ng,xyz,chrg,coefficients,exponents,io2,basis,nat,nelalpha, nelbeta,erep, finishtei, starttei, finishscf, startscf,io3,escf, ausgabe_Erfolgt)


    integer :: nbf, nelalpha, nelbeta,nat, io2, ng, io3, i, ausgabe_erfolgt
    real(wp) :: delta, Erep, escf, newescf
    real :: finishtei, starttei,finishscf, startscf, evenergy
    real (wp):: xyz(:,:), chrg(:), coefficients(:), exponents(:)
    real(wp), allocatable:: packsab(:), packtab(:), packvab(:), sab(:,:), tab(:,:), vab(:,:), xab(:,:),hab(:,:), Fockalpha(:,:), Fockbeta(:,:),Fock_newalpha(:,:),Fock_newbeta(:,:)
    real(wp),allocatable :: eigval(:), cabalpha(:,:), cabbeta(:,:), pabalpha(:,:), pabbeta(:,:), gabcd(:,:), twointeg(:), aufpunkt(:,:)
    real(wp),allocatable :: eigvec(:,:), basis(:)

    allocate(packsab(nbf*(1+nbf)/2),packtab(nbf*(1+nbf)/2),packvab(nbf*(1+nbf)/2),xab(nbf,nbf), eigval(nbf), eigvec(nbf,nbf), Fockalpha(nbf,nbf),Fock_newalpha(nbf,nbf), Fockbeta(nbf,nbf),Fock_newbeta(nbf,nbf))
    allocate(hab(nbf,nbf), sab(nbf,nbf), tab(nbf,nbf), vab(nbf,nbf), cabalpha(nbf,nbf), cabbeta(nbf,nbf), pabalpha(nbf,nbf), pabbeta(nbf,nbf), gabcd(nbf,nbf))
    allocate(twointeg(((nbf*(nbf-1)/2+nbf)*(nbf*(nbf-1)/2+nbf)/2+(nbf*(nbf-1)/2+nbf)-1)))

     !>Calculation of one electron integrals
     call oneelint(nbf, ng, xyz, chrg, coefficients, exponents, sab, tab, vab,io2,basis,nat,aufpunkt,ausgabe_Erfolgt)


     if(ausgabe_erfolgt==0) then
     !> Text for result File
     write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
     write(io2,*)"                                           Packing Matrices                                           "
     write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
     write(io2,*)
     write(io2,*)"<<<<<<<<<<<<<<<<<     Packing order: (I) S matrix, (II) T Matrix, (III) V Matrix    >>>>>>>>>>>>>>>>>>>"
     write(io2,*)"======================================================================================================="
   end if

     !>Packin one electron matrices
     call packer(sab,packsab, nbf,io2,ausgabe_erfolgt)
     call packer(tab,packtab, nbf,io2,ausgabe_erfolgt)
     call packer(vab,packvab, nbf,io2,ausgabe_erfolgt)


     !> Calculating orthonormalizer with use of symmetric procedure
     call orthonormalizer(packsab,sab,eigval,eigvec,xab,io2,nbf,ausgabe_Erfolgt)


     if(ausgabe_erfolgt==0) then
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
   end if

     !>Setting initial Fock Matrix as a core Hamiltonian
     Fockalpha=tab+vab
     Fockbeta=tab+vab

     if(ausgabe_erfolgt==0) then
     write(*,*)
     write(*,*)"                                         ✓ Fock matrices obtained        "
     write(io2,*)
     call write_matrix(Fockalpha,"      ============================       Fock Matrix      ============================",io2)
     write(io2,*)"      =================================================================================="
     end if

     !>calculating coefficients from the initial Fock matrix
     call coeff(cabalpha,Fockalpha,xab,nbf,io2,ausgabe_erfolgt)
     call coeff(cabbeta,Fockbeta,xab,nbf,io2,ausgabe_erfolgt)

     !>Calculating two electron integrals
     call twoIntegrals(aufpunkt, exponents, coefficients, twointeg, ng,nbf,starttei,finishtei)


     if(ausgabe_erfolgt==0) then

     !Calculating new density matrices for both spins
     call unres_new_density(nelalpha, nbf,cabalpha,pabalpha,io2,ausgabe_erfolgt)
     call unres_new_density(nelbeta, nbf,cabbeta,pabbeta,io2,ausgabe_erfolgt)
     write(*,*)
     write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
     write(*,*)"                 <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
     write(*,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
     write(io2,*)
     write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"
     write(io2,*)"               <<<<<<<<<<<<<<      Initial guess ended succesfully      >>>>>>>>>>>>>>                 "
     write(io2,*)"–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––"

     end if

     !HCore matrix as sum of kinetic and attraction matrices
     hab=tab+vab

     !Calculating Hartree Fock energy
     call UHFenergy(nbf,escf,hab,Fockalpha, Fockbeta,pabalpha, pabbeta)

     if(ausgabe_erfolgt==0) then
     !Printing the initial Fock energy
     write(*,*)
     write(*,*) "                           *************************************************"
     write(*,*) "                           *  ", "initial E_{HF}=", escf, "H  *"
     write(*,*) "                           *************************************************"
     write(*,*)
     write(io2,*)
     write(io2,*) "                           *************************************************"
     write(io2,*) "                               ", "initial E_{HF}=", escf, "H"

     end if

     !Calculating the G Tensor and new Fock Matrix
     call newuFock(nbf,pabalpha,pabbeta,hab,Fock_newalpha,twointeg,gabcd,io2)
    call newuFock(nbf,pabbeta,pabalpha,hab,Fock_newbeta,twointeg,gabcd,io2)

     if(ausgabe_erfolgt==0) then
       !Printing results and initial HF energy /stdout+file
       call write_matrix(gabcd,"           ===========================      G Tensor       ===========================     ",io2)
       write(io2,*)"      =================================================================================="
       write(io2,*)
       call write_matrix(Fock_newalpha,"         ==========================       Fock α Matrix      ==========================     ",io2)
       write(io2,*)"      =================================================================================="
       write(io2,*)
       write(io2,*)
       call write_matrix(Fock_newbeta,"         ==========================       Fock β Matrix      ==========================     ",io2)
       write(io2,*)"      =================================================================================="
       write(io2,*)
     end if


     call UHFenergy(nbf,escf,hab,Fockalpha, Fockbeta,pabalpha, pabbeta)

    if(ausgabe_erfolgt==0) then
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
    end if


     call cpu_time(startscf)

     !Setting stop criterion for at least one SCF Run
     delta=1

     !Setting run counter
     i=1


     !Definition of Convergence
     do while (delta>0.0000001.and.i<50)


       !Calling the steps
       call unrest_iterations(io2,cabalpha, cabbeta,Fock_newalpha,Fock_newbeta ,xab,nbf,nelalpha,nelbeta,pabalpha, pabbeta,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3,ausgabe_erfolgt)

       !Calculating energy change
       delta=escf-newescf
       delta=sqrt(delta**2)
       !Overwriting the HF Energy
       escf=newescf
       evenergy=(Erep+escf)*1.042*10.**(-21)

       !Counter increase
       i=i+1

     end do

     !Printing information about end of SCF Procedure
     if (i<50) then
       if(ausgabe_Erfolgt==0) then

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
         write(*,*) "                          |         =", (Erep+escf)*27.2114, "eV          |"
         write(*,1)"                           |         =  ", evenergy,"     kcal        |"
         1 format (a,e20.14,a)
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
         write(io2,*) "                          |         =", (Erep+escf)*27.2114, "eV          |"
         write(io2,1)"                           |         =  ", evenergy,"     kcal        |"
         write(io2,*) "                          ==================================================="
         write(io3,*)"    ==============================================================================================="
         write(io3,*)"    ||                                      SCF SUCCSESFULL                                      ||"
         write(io3,*)"    ==============================================================================================="
         write(io3,*)
       end if
     else
       if(ausgabe_erfolgt==0) then
         write(*,*)"    ==============================================================================================="
         write(*,*)"    ||                                        SCF FAILED                                         ||"
         write(*,*)"    ==============================================================================================="
       end if
     end if

     if(ausgabe_erfolgt==0) then
       call cpu_time(finishscf)
       write(io3,*) "======================================================================================================="
       write(io3,*) "SCF iterations time:", finishscf-startscf, "s"
      end if
      escf=escf+erep
      write(*,*)"final SCF energy",escf



   end subroutine unrestrictedHF
 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END UNRESTRICTED FOCK @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

 !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  NEW UNRES. DENISTY MATRIX   >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine unres_new_density(nel, nbf,cab,pab,io2,ausgabe_erfolgt) !Exercise 8.4

   !Declaration of variables
   integer:: i, j, nel,nbf,io2,ausgabe_erfolgt
   real(wp)::cab(nbf,nbf), pab(nbf,nbf)
   real(wp), allocatable :: nocc(:,:)

   !Allocation of needed memory
   allocate(nocc(nbf,nbf))

   !Preparring the occupation matrix
   nocc=0
   do i=1,nel
     do j=1,nel
       if(i==j) then
         nocc(i,j)=1
       end if
     end do
   end do

   !CAlculation of the new density matrix
   pab=matmul(nocc,transpose(cab))
   pab=matmul(cab,pab)


   if(ausgabe_erfolgt==0) then
     !Printing succes text /stdout+file
     write(*,*)
     write(*,*)"                                       ✓ Density matrix obtained       "
     write(*,*)
end if
   deallocate(nocc)

 end subroutine unres_new_density


 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END UNRES. NEW DENSITY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  UHF ENERGY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine UHFenergy(nbf,escf,hab,Fockalpha, Fockbeta,pabalpha, pabbeta)

   !>Declaration of local variables
   integer:: i, j,nbf
   real(wp):: escf
   real(wp):: hab(nbf,nbf),Fockalpha(nbf,nbf),pabalpha(nbf,nbf),Fockbeta(nbf,nbf),pabbeta(nbf,nbf)

   !Setting energy to 0
   escf=0


   !>Calculate the E_{HF} as a half of the matrix trace
   do i=1,nbf

     do j=1,nbf

         escf=escf+0.5_wp*(pabalpha(i,j)*(hab(j,i)+fockalpha(j,i)))+0.5_wp*(pabbeta(i,j)*(hab(j,i)+fockbeta(j,i)))

     end do

   end do


 end subroutine UHFenergy


 !@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END UHF ENERGY @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 !–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
 !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<   MATRIX >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

 subroutine newuFock(nbf,pabalpha,pabbeta,hab,Fock_newalpha,twointeg,gabcd,io2) !Exercise12.1


   !Declaration of variables
   integer:: i,j,k,l, ij, kl, ijkl, il, kj, ilkj, nbf,io2
   real(wp):: pabalpha(:,:), pabbeta(:,:), hab(:,:), Fock_newalpha(:,:),twointeg(:), gabcd(:,:)

   !Setiing new Fock equal 0
   Fock_newalpha=0

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
            Fock_newalpha(i,j)=Fock_newalpha(i,j)+(pabalpha(l,k)+pabbeta(l,k))*twointeg(ijkl)-pabalpha(l,k)*twointeg(ilkj)


          end do

        end do

      end do

    end do



    !Setiing G tensor equal the new Fock
    gabcd=Fock_newalpha

    !>Adding the core Hamiltonian part
    Fock_newalpha=Fock_newalpha+hab


  end subroutine newuFock
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END NEW FOCK @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
  !<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<  SCF ITERATIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

subroutine  unrest_iterations(io2,cabalpha, cabbeta,Fock_newalpha,Fock_newbeta ,xab,nbf,nelalpha,nelbeta,pabalpha, pabbeta,newescf,Erep,hab,gabcd,xyz, exponents,coefficients, twointeg,i,io3,ausgabe_erfolgt)

  integer :: i,io2,io3,nbf,nelalpha,nelbeta,ng, ausgabe_erfolgt
  real(wp) :: xyz(:,:), exponents(:), coefficients(:), twointeg(:)
  real(wp) :: Erep,newescf,escf
  real(wp) :: cabalpha(:,:), cabbeta(:,:),Fock_newalpha(:,:),Fock_newbeta(:,:),xab(:,:), pabalpha(:,:),pabbeta(:,:), hab(:,:), gabcd(:,:)

  if (ausgabe_erfolgt==0) then
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
  end if

  !Calculation of the transformed F matrix, coefficients and density for alpha spin
  write(*,*)"                                     .............................."
  write(*,*)"                                           α spin calculation"
  write(*,*)"                                     .............................."
  write(io2,*)"                                     .............................."
  write(io2,*)"                                           α spin calculation"
  write(io2,*)"                                     .............................."
  call coeff(cabalpha,Fock_newalpha,xab,nbf,io3,ausgabe_erfolgt)

  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                 ✓ New orbital coefficients obtained    "
  end if

  call unres_new_density(nelalpha, nbf,cabalpha,pabalpha,io2,ausgabe_erfolgt)
  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                   ✓ New density matrix obtained       "
  end if

  !Calculation of the transformed F matrix, coefficients and density for beta spin
  write(*,*)"                                     .............................."
  write(*,*)"                                           β spin calculation"
  write(*,*)"                                     .............................."
  write(io2,*)"                                     .............................."
  write(io2,*)"                                           β spin calculation"
  write(io2,*)"                                     .............................."
  call coeff(cabbeta,Fock_newbeta,xab,nbf,io3,ausgabe_erfolgt)

  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                 ✓ New orbital coefficients obtained    "
  end if

  call unres_new_density(nelbeta, nbf,cabbeta,pabbeta,io2,ausgabe_erfolgt)

  if(ausgabe_erfolgt==0) then
    write(io2,*)"                                   ✓ New density matrix obtained       "
  end if

  !Calculation of new G tensors and Fock Matrix
  call newuFock(nbf,pabalpha,pabbeta,hab,Fock_newalpha,twointeg,gabcd,io2)
  call newuFock(nbf,pabbeta,pabalpha,hab,Fock_newbeta,twointeg,gabcd,io2)

  if(ausgabe_erfolgt==0) then
    !Printing some results /std+file
    write(*,*)"                                     .............................."
    write(io2,*)"                                     .............................."
    write(*,*)
    write(*,*)"                                          ✓ G tensor obtained       "
    write(io2,*)"                                      ✓ New G tensor obtained       "
    call write_matrix(gabcd,"                =========================      G Tensor      ==========================",io3)
    write(io3,*)"      =================================================================================="
    write(io3,*)
    write(*,*)
    write(*,*)"                                      ✓ New Fock α Matrix obtained       "
    write(*,*)
    write(*,*)"                                      ✓ New Fock β Matrix obtained       "
    write(io2,*)"                                    ✓ New Fock α Matrix obtained       "
    call write_matrix(Fock_newalpha,"               ========================     Fock α Matrix     =========================",io3)
    write(io3,*)"      =================================================================================="
    write(io3,*)
    write(io2,*)"                                    ✓ New Fock β Matrix obtained       "
    call write_matrix(Fock_newbeta,"               ========================      Fock β Matrix      =========================",io3)
    write(io3,*)"      =================================================================================="
    write(io3,*)
  end if

  !Calculating new HF Energy
   call UHFenergy(nbf,newescf,hab,Fock_newalpha, Fock_newbeta,pabalpha, pabbeta)

  if(ausgabe_erfolgt==0) then
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
  end if

end subroutine unrest_iterations

!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ END SCF ITERATIONS @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––

end module scf_main
