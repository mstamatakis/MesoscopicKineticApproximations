module constants_module

implicit none

!maximum number of two body interaction terms
integer, parameter :: inttermsmax = 66

! for repetitive warnings: the maximum number of message that will be issued for each such warning
integer, parameter :: maxrepwarnings = 10

real(8), parameter :: kboltz = 1.380650321752056d-23*6.24150974d+18 ! eV/K

! commenting character
character(1), parameter :: remchar = '#'

! double the maximum number of words separated by a delimiter that can be read
integer, parameter :: maxwords2 = 6000

! record length (in number of characters) that can be read
integer, parameter :: lengthrecinp = 8192

! allowed length (in characters) of the species', site-types' and mechanism-steps' names
integer, parameter :: nnam0 = 64

! Record length buffer for formatted input/output.
! Used to fix recl to a default values across compilers
! to fix some tests with NAG compilers.
! The NAG compiler creates code that exits with buffer overflow if a text line
! is longer than recl. However this does not affect Intel which uses a much
! smaller recl. Apparently this is implemented differently.
! For formatted files the default value is 1024 using NAG,
! Intel 132 (https://software.intel.com/en-us/node/511275)
! 1073741824 using gfortran (https://gcc.gnu.org/onlinedocs/gcc-4.9.1/gfortran.pdf)
! This needs to store all the mechanism names in one line
integer, parameter :: defaultrecl = nnam0*1000

! *** File names and unit numbers

! The input file names are not constant (parameter) to allow passing them from the
! command line. These are the default values.

! simulation input file
character(100) :: csimfname = 'calculation_input.dat'
integer, parameter :: iread = 103

! restart file
character(100) :: crestartfname = 'restart.inf'
integer :: irestart = 201

! general output file
character(30) :: cgenoutfname = 'calculation_setup.txt'
integer :: iwrite = 202
end module constants_module
