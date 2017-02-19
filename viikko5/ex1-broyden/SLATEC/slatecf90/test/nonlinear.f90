!---------------------------------------------------
!
! Demonstration of the SLATEC routine DNSQE.
!
! Find the root of nonlinear equations:
!   f1(x1,x2)=0, f2(x1,x2)=0
! where
!   f1(x1,x2)=exp(-x1**2-x2**2)-1/8
!   f2(x1,x2)=sin(x1)-cos(x2)
!
! Compilation 
!   gfortran nonlinear.f90  -lslatec -L../src
!
!         Antti Kuronen, 2011
!         2017: Corrected a bug: Jacobian was transposed
!
!---------------------------------------------------

module sizes
  integer,parameter :: rk=selected_real_kind(10,40),MAXBUF=200
  integer :: fcn_calls,jac_calls
end module sizes


! --- Calculate function values ---

subroutine fcn(n,x,fvec,iflag)
  use sizes
  implicit none
  integer :: n,iflag
  real(rk) :: x(n),fvec(n)
  fvec(1)=exp(-x(1)**2-x(2)**2)-0.125
  fvec(2)=sin(x(1))-cos(x(2))
  if (iflag==0) print '(4g20.10,2x,g20.10)',x,fvec,sum(fvec**2)
  fcn_calls=fcn_calls+1
  return
end subroutine fcn


! --- Calculate Jacobian ---

subroutine jac(n,x,fvec,fjac,ldfjac,iflag)
  use sizes
  implicit none
  integer :: n,ldfjac,iflag
  real(rk) :: x(n),fvec(n),fjac(ldfjac,n)
  fjac(1,1)=-2.0*x(1)*exp(-x(1)**2-x(2)**2)
  fjac(1,2)=-2.0*x(2)*exp(-x(1)**2-x(2)**2)
  fjac(2,1)=cos(x(1))
  fjac(2,2)=sin(x(2))
  !fjac=reshape([1.0,0.0,0.0,1.0],[2,2])
  jac_calls=jac_calls+1
  return
end subroutine jac


! --- Main program ---

program nonlinear

  use sizes
  implicit none
  integer,parameter :: n=2
  integer :: iopt,nprint,info,lwa
  real(rk) :: tol,x(n),fvec(n)
  real(rk),allocatable :: wa(:)
  external :: fcn,jac
  character(len=MAXBUF) :: arg

  lwa=(3*n**2+13*n)/2
  allocate(wa(lwa))
  iopt=1

  call get_command_argument(1,arg); read(arg,*) x(1)
  call get_command_argument(2,arg); read(arg,*) x(2)
  call get_command_argument(3,arg); read(arg,*) tol
  call get_command_argument(4,arg); read(arg,*) nprint
  
  call dnsqe(fcn,jac,iopt,n,x,fvec,tol,nprint,info,wa,lwa)
  
  print '(a,i0)','# finished, info=',info
  print '(4g20.10,2x,i0)',x,fvec
  print *,fcn_calls,jac_calls
  
  stop
end program nonlinear
