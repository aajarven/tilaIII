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
  integer :: Ncalls

  ! vertices of the initial simplex
  p(1,:) = (/ -5.0, 0.0 /)
  p(2,:) = (/ 5.0, 0.0 /)
  p(3,:) = (/ 0.0, 5.0 /)

  call simplex(funk, 2, p, y)

  write(*,*) p(1,:), y(1)
  write(*,*) p(2,:), y(2)
  write(*,*) p(3,:), y(3)


end program

real function funk(p, ndim)
  implicit none
  integer, intent(in) :: ndim
  real, dimension(ndim) :: p
  real, dimension(6, 2) :: points
  integer :: i

  points(1,:) = (/ 0.0, 0.0 /)
  points(2,:) = (/ 1.0, 1.0 /)
  points(3,:) = (/ 2.0, 0.0 /)
  points(4,:) = (/ 3.0, 2.0 /)
  points(5,:) = (/ 4.0, 4.0 /)
  points(6,:) = (/ 5.0, 5.0 /)
  
  funk = 0
  do i=1,6
    funk = funk + abs(points(i, 2) - p(1)*points(i, 1) - p(2))
  end do
end function
