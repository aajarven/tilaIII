subroutine BSPDR (T, A, N, K, NDERIV, AD)
!
!! BSPDR uses the B-representation to construct a divided difference
!  table preparatory to a (right) derivative calculation.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3
!***TYPE      SINGLE PRECISION (BSPDR-S, DBSPDR-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, DIFFERENTIATION OF SPLINES,
!             INTERPOLATION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract
!         BSPDR is the BSPLDR routine of the reference.
!
!         BSPDR uses the B-representation (T,A,N,K) to construct a
!         divided difference table ADIF preparatory to a (right)
!         derivative calculation in BSPEV.  The lower triangular matrix
!         ADIF is stored in vector AD by columns.  The arrays are
!         related by
!
!           ADIF(I,J) = AD(I-J+1 + (2*N-J+2)*(J-1)/2)
!
!         I = J,N , J = 1,NDERIV .
!
!     Description of Arguments
!         Input
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the spline, K  >=  1
!          NDERIV  - number of derivatives, 1  <=  NDERIV  <=  K.
!                    NDERIV=1 gives the zero-th derivative = function
!                    value
!
!         Output
!          AD      - table of differences in a vector of length
!                    (2*N-NDERIV+1)*NDERIV/2 for input to BSPEV
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BSPDR
!
  INTEGER I, ID, II, IPKMID, JJ, JM, K, KMID, N, NDERIV
  REAL A, AD, DIFF, FKMID, T
!     DIMENSION T(N+K), AD((2*N-NDERIV+1)*NDERIV/2)
  DIMENSION T(*), A(*), AD(*)
!***FIRST EXECUTABLE STATEMENT  BSPDR
  if ( K < 1) go to 100
  if ( N < K) go to 105
  if ( NDERIV < 1 .OR. NDERIV > K) go to 110
  DO 10 I=1,N
    AD(I) = A(I)
   10 CONTINUE
  if (NDERIV == 1) RETURN
  KMID = K
  JJ = N
  JM = 0
  DO 30 ID=2,NDERIV
    KMID = KMID - 1
    FKMID = KMID
    II = 1
    DO 20 I=ID,N
      IPKMID = I + KMID
      DIFF = T(IPKMID) - T(I)
      if (DIFF /= 0.0E0) AD(II+JJ) = (AD(II+JM+1)-AD(II+JM))/ &
       DIFF*FKMID
      II = II + 1
   20   CONTINUE
    JM = JJ
    JJ = JJ + N - ID + 1
   30 CONTINUE
  return
!
!
  100 CONTINUE
  call XERMSG ('SLATEC', 'BSPDR', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'BSPDR', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'BSPDR', &
     'NDERIV DOES NOT SATISFY 1 <= NDERIV <= K', 2, 1)
  return
end
