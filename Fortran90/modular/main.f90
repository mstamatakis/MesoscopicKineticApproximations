	program k2nnc1
	use constants
	implicit none
	logical check
	integer i, j
	integer, dimension(nsi) :: state
	integer, allocatable, dimension(:) :: vec
	real*8 hamiltonian, part, corr
	real*8 cov
	real*8, dimension(4) :: x

	! Calculating the MF correction parameters
	check=.false.
	allocate(vec(1))
	state=0
	call confs(state,0)
	x(1)=0.1d0
	x(2)=0.1d0
	x(3)=0.1d0
	x(4)=0.1d0
	chemp=-1.4d0
	do i=1,240
	 chemp=chemp+0.01d0
     write(*,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~'
     write(*,*) 'Chemical potential', chemp
	 call solver(x,4,check) 
	 ! Calculating the coverage
	 cov=0.d0
	 do j=1,13
	  vec(1)=j
	  cov=cov+corr(vec,1,1,x)/part(4,x)
	 end do
	 cov=cov/13.d0
	 write(16,*) chemp, cov
	end do
	end program
