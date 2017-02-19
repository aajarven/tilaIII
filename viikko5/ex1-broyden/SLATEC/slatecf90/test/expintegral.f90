program expintegral

!
! SLATEC exponential integral test.
!
! Compilation: gfortran expintegral.f90 -lslatec -L../src
!

  implicit none
  integer,parameter :: rk=selected_real_kind(10,40),MAXBUF=200
  integer,parameter :: kode=1,m=1
  real(rk) :: x1,x2,dx,x,tol,en(m),d1mach
  integer :: n,nz,ierr,ip,i
  character(len=MAXBUF) :: argu

  if (command_argument_count()/=5) then
     call get_command_argument(0,argu)
     write(0,'(a,a,a)') 'usage: ',trim(argu),' n x1 x2 dx tol'
     stop
  endif
  
  call get_command_argument(1,argu) ;read(argu,*) n
  call get_command_argument(2,argu) ;read(argu,*) x1
  call get_command_argument(3,argu) ;read(argu,*) x2
  call get_command_argument(4,argu) ;read(argu,*) dx
  call get_command_argument(5,argu) ;read(argu,*) tol
  
  ip=(x2-x1)/dx

  do i=1,ip
     x=x1+dx*(i-1)
     call dexint(x,n,kode,m,tol,en,nz,ierr)
     write(6,'(2g20.10,5x,2i4)') x,en(1),nz,ierr
  end do

  stop
end program expintegral
