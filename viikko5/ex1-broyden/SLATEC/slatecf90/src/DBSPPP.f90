subroutine DBSPPP (T, A, N, K, LDC, C, XI, LXI, WORK)
!
!! DBSPPP converts the B-representation of a B-spline to the piecewise ...
!  polynomial (PP) form.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BSPPP-S, DBSPPP-D)
!***KEYWORDS  B-SPLINE, PIECEWISE POLYNOMIAL
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!         DBSPPP is the BSPLPP routine of the reference.
!
!         DBSPPP converts the B-representation (T,A,N,K) to the
!         piecewise polynomial (PP) form (C,XI,LXI,K) for use with
!         DPPVAL.  Here XI(*), the break point array of length LXI, is
!         the knot array T(*) with multiplicities removed.  The columns
!         of the matrix C(I,J) contain the right Taylor derivatives
!         for the polynomial expansion about XI(J) for the intervals
!         XI(J)  <=  X  <=  XI(J+1), I=1,K, J=1,LXI.  Function DPPVAL
!         makes this evaluation at a specified point X in
!         XI(1)  <=  X  <=  XI(LXI+1)
!
!     Description of Arguments
!
!         Input      T,A are double precision
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K  >=  1
!          LDC     - leading dimension of C, LDC  >=  K
!
!         Output     C,XI,WORK are double precision
!          C       - matrix of dimension at least (K,LXI) containing
!                    right derivatives at break points
!          XI      - XI break point vector of length LXI+1
!          LXI     - number of break points, LXI  <=  N-K+1
!          WORK    - work vector of length K*(N+3)
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  DBSPDR, DBSPEV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSPPP
!
  INTEGER ILEFT, INEV, K, LDC, LXI, N, NK
  DOUBLE PRECISION A, C, T, WORK, XI
!     DIMENSION T(N+K),XI(LXI+1),C(LDC,*)
!     HERE, * = THE FINAL VALUE OF THE OUTPUT PARAMETER LXI.
  DIMENSION T(*), A(*), WORK(*), XI(*), C(LDC,*)
!***FIRST EXECUTABLE STATEMENT  DBSPPP
  if ( K < 1) go to 100
  if ( N < K) go to 105
  if ( LDC < K) go to 110
  call DBSPDR(T, A, N, K, K, WORK)
  LXI = 0
  XI(1) = T(K)
  INEV = 1
  NK = N*K + 1
  DO 10 ILEFT=K,N
    if (T(ILEFT+1) == T(ILEFT)) go to 10
    LXI = LXI + 1
    XI(LXI+1) = T(ILEFT+1)
    call DBSPEV(T,WORK(1),N,K, K,XI(LXI),INEV,C(1,LXI),WORK(NK))
   10 CONTINUE
  return
  100 CONTINUE
  call XERMSG ('SLATEC', 'DBSPPP', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'DBSPPP', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'DBSPPP', 'LDC DOES NOT SATISFY LDC >= K', &
     2, 1)
  return
end
