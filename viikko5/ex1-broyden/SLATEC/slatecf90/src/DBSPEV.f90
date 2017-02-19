subroutine DBSPEV (T, AD, N, K, NDERIV, X, INEV, SVALUE, WORK)
!
!! DBSPEV calculates the value of the spline and its derivatives from...
!  the B-representation.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BSPEV-S, DBSPEV-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!         DBSPEV is the BSPLEV routine of the reference.
!
!         DBSPEV calculates the value of the spline and its derivatives
!         at X from the B-representation (T,A,N,K) and returns them in
!         SVALUE(I),I=1,NDERIV, T(K)  <=  X  <=  T(N+1).  AD(I) can be
!         the B-spline coefficients A(I), I=1,N) if NDERIV=1.  Otherwise
!         AD must be computed before hand by a call to DBSPDR (T,A,N,K,
!         NDERIV,AD).  If X=T(I),I=K,N), right limiting values are
!         obtained.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         DBSPEV calls DINTRV, DBSPVN
!
!     Description of Arguments
!
!         Input      T,AD,X, are double precision
!          T       - knot vector of length N+K
!          AD      - vector of length (2*N-NDERIV+1)*NDERIV/2 containing
!                    the difference table from DBSPDR.
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K  >=  1
!          NDERIV  - number of derivatives, 1  <=  NDERIV  <=  K.
!                    NDERIV=1 gives the zero-th derivative =
!                    function value
!          X       - argument, T(K)  <=  X  <=  T(N+1)
!          INEV    - an initialization parameter which must be set
!                    to 1 the first time DBSPEV is called.
!
!         Output     SVALUE,WORK are double precision
!          INEV    - INEV contains information for efficient process-
!                    ing after the initial call and INEV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INEV parameters.
!          SVALUE  - vector of length NDERIV containing the spline
!                    value in SVALUE(1) and the NDERIV-1 derivatives
!                    in the remaining components.
!          WORK    - work vector of length 3*K
!
!     Error Conditions
!         Improper input is a fatal error.
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  DBSPVN, DINTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSPEV
!
  INTEGER I,ID,INEV,IWORK,JJ,K,KP1,KP1MN,L,LEFT,LL,MFLAG, &
   N, NDERIV
  DOUBLE PRECISION AD, SVALUE, SUM, T, WORK, X
!     DIMENSION T(N+K)
  DIMENSION T(*), AD(*), SVALUE(*), WORK(*)
!***FIRST EXECUTABLE STATEMENT  DBSPEV
  if ( K < 1) go to 100
  if ( N < K) go to 105
  if ( NDERIV < 1 .OR. NDERIV > K) go to 115
  ID = NDERIV
  call DINTRV(T, N+1, X, INEV, I, MFLAG)
  if (X < T(K)) go to 110
  if (MFLAG == 0) go to 30
  if (X > T(I)) go to 110
   20 if (I == K) go to 120
  I = I - 1
  if (X == T(I)) go to 20
!
! *I* HAS BEEN FOUND IN (K,N) SO THAT T(I)  <=  X  <  T(I+1)
!     (OR  <=  T(I+1), if T(I)  <  T(I+1) = T(N+1) ).
   30 KP1MN = K + 1 - ID
  KP1 = K + 1
  call DBSPVN(T, KP1MN, K, 1, X, I, WORK(1),WORK(KP1),IWORK)
  JJ = (N+N-ID+2)*(ID-1)/2
!     ADIF(LEFTPL,ID) = AD(LEFTPL-ID+1 + (2*N-ID+2)*(ID-1)/2)
!     LEFTPL = LEFT + L
   40 LEFT = I - KP1MN
  SUM = 0.0D0
  LL = LEFT + JJ + 2 - ID
  DO 50 L=1,KP1MN
    SUM = SUM + WORK(L)*AD(LL)
    LL = LL + 1
   50 CONTINUE
  SVALUE(ID) = SUM
  ID = ID - 1
  if (ID == 0) go to 60
  JJ = JJ-(N-ID+1)
  KP1MN = KP1MN + 1
  call DBSPVN(T, KP1MN, K, 2, X, I, WORK(1), WORK(KP1),IWORK)
  go to 40
!
   60 RETURN
!
!
  100 CONTINUE
  call XERMSG ('SLATEC', 'DBSPEV', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'DBSPEV', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'DBSPEV', &
     'X IS NOT IN T(K) <= X <= T(N+1)', 2, 1)
  return
  115 CONTINUE
  call XERMSG ('SLATEC', 'DBSPEV', &
     'NDERIV DOES NOT SATISFY 1 <= NDERIV <= K', 2, 1)
  return
  120 CONTINUE
  call XERMSG ('SLATEC', 'DBSPEV', &
     'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
  return
end
