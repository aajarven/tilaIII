
!
! SLATEC Bessel function test.
!
! Compilation: gfortran besseltest.f90 -lslatec -L../src
!
!

program besseltest

  implicit none
  integer,parameter :: rk=selected_real_kind(kind(1.0d0))
  real(rk) :: x,dbesk1
  integer :: i
 
  do i=1,100
     x=0.01*i
     print *,x,dbesk1(x)
  end do


end program besseltest
