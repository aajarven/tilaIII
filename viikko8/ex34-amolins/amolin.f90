program amolin
! fit a sraight line to data using the L-infinity norm
  use amoeba
  implicit none
  interface
    real function funk(p, ndim)
      integer, intent(in) :: ndim
      real, dimension(ndim) :: p
    end function
  end interface

  integer, parameter :: ndim = 2     ! dimension of the space
  real p(ndim+1, ndim), &
       y(ndim+1)
  integer	:: Ncalls

  ! vertices of the initial simplex
  p(1,:) = (/ -1.0, -1.0 /)
  p(2,:) = (/ 0.0, 3.0 /)
  p(3,:) = (/ 3.0, 0.0 /)

  call simplex(funk, 2, p, y)

  write(*,*) p(1,:), y(1)
  write(*,*) p(2,:), y(2)
  write(*,*) p(3,:), y(3)


end program

real function funk(p, ndim)
  implicit none
  integer, intent(in) :: ndim
  real, dimension(ndim) :: p
  
  funk= exp(norm2(p))
end function
