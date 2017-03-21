
module amoeba
  implicit none
  integer, parameter :: NMAX = 1000     ! max number of steps
  real, parameter :: FTOL = 0.0001, &   ! termination tolerance
                     ALPHA = 1.0, &     ! reflection factor
                     BETA = 0.5, &      ! reduction factor
                     GAMMA = 2.0        ! expansion factor
  integer :: nfunk                      ! number of function calls
             
contains

subroutine simplex(f, ndim, p, y)
! find the minimum of the function f by the Nelder-Mead
! polytope or simplex method
! ndim is the dimension of the space
! p(ndim+1, ndim) contains the values of the function
  interface
    real function f(p, ndim)
      integer, intent(in) :: ndim
      real, dimension(ndim) :: p
    end function
  end interface
  integer, intent(in) :: ndim
  real, dimension (ndim+1, ndim) :: p
  real, dimension (ndim+1) :: y
!
  integer :: mpts         ! number of points of the simples
  integer i,j, ilo,ihi,inhi
  real ytry, ysave, sum, rtol, psum(ndim)

  mpts = ndim+1
  
  ! initial values of the function
  do i=1,mpts
    y(i) = f(p(i,:), ndim)
  end do

  nfunk=0

  ! sums of coordinates; the centre of gravity 
  ! opposite to apex i is (psum-p(i))/ndim
  do j=1,ndim
    sum=0.0
    do i=1,mpts
      sum = sum+p(i,j) 
    end do
    psum(j) = sum
  end do

  main: do
    ! finfd the best (ilo), worst (ihi) and
    ! second worst (inhi) apex of the simplex
    ilo=1
    if (y(1) > y(2)) then
       inhi=2; ihi=1
    else
       inhi=1; ihi=2
    end if
    do i=1,mpts
      if (y(i) < y(ilo)) ilo=i
      if (y(i) > y(ihi)) then
         inhi=ihi; ihi=i
      else if (y(i) > y(inhi)) then
         if (i /= ihi) inhi=i
      end if
    end do

    ! stop if the best and worst values are almost equal 
    ! NB: relative error; the function must not be zero
    rtol=2.0*abs(y(ihi)-y(ilo)) / (abs(y(ihi))+abs(y(ilo)))
    if (rtol < FTOL) then
           write(6,'("rtol=",f10.5)') rtol
           write(6,'("Number of function calls =",i5)') nfunk
           exit main
    end if
    ! stop if too many iterations
    if (nfunk >= NMAX) then
           write(6,'("nfunk=",i5)') nfunk
           exit main
    end if

    ! reflection
    ytry=amotry(-ALPHA)

    ! if the result improves, move the new apex even further
    if (ytry <= y(ilo)) then
      ytry=amotry(GAMMA)
    else if (ytry >= y(inhi)) then
      ! the new point worse than the second worst;
      ! shrink the simplex
      ysave=y(ihi)
      ytry=amotry(BETA)
      if (ytry >= ysave) then
        ! still too big; shrink the simplex
        ! w.r.t the best point
        do i=1,mpts
          if (i /= ilo) then
            do j=1,ndim
              psum(j)=0.5*(p(i,j)+p(ilo,j))
              p(i,j)=psum(j)
            end do
            y(i)=f(psum, ndim)              
          end if
        end do
        nfunk = nfunk+ndim
        ! update sums of coordinates
        do j=1,ndim
          sum=0.0
          do i=1,mpts
            sum = sum+p(i,j) 
          end do
          psum(j) = sum
        end do
      end if
    end if 
  end do main

contains

real function amotry(fac)
  implicit none
  real, intent(in) :: fac
  integer j
  real :: fac1, fac2, ytry, ptry(ndim)
  fac1=(1.0-fac)/ndim
  fac2=fac-(1-fac)/ndim
  do j=1,ndim
     ptry(j)=psum(j)*fac1+p(ihi,j)*fac2
  end do
  ytry=f(ptry, ndim)
  nfunk = nfunk+1
  if (ytry < y(ihi)) then
    y(ihi)=ytry
    do j=1,ndim
      psum(j) = psum(j)+ptry(j)-p(ihi,j)
      p(ihi,j)=ptry(j)
    end do
  end if
  amotry=ytry
end function

end subroutine

end module










