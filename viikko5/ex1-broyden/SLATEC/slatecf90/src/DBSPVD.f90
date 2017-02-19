subroutine DBSPVD (T, K, NDERIV, X, ILEFT, LDVNIK, VNIKX, WORK)
!
!! DBSPVD calculates the value and all derivatives of order less than ...
!  NDERIV of all basis functions which do not vanish at X.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BSPVD-S, DBSPVD-D)
!***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!
!         DBSPVD is the BSPLVD routine of the reference.
!
!         DBSPVD calculates the value and all derivatives of order
!         less than NDERIV of all basis functions which do not
!         (possibly) vanish at X.  ILEFT is input such that
!         T(ILEFT)  <=  X  <  T(ILEFT+1).  A call to INTRV(T,N+1,X,
!         ILO,ILEFT,MFLAG) will produce the proper ILEFT.  The output of
!         DBSPVD is a matrix VNIKX(I,J) of dimension at least (K,NDERIV)
!         whose columns contain the K nonzero basis functions and
!         their NDERIV-1 right derivatives at X, I=1,K, J=1,NDERIV.
!         These basis functions have indices ILEFT-K+I, I=1,K,
!         K  <=  ILEFT  <=  N.  The nonzero part of the I-th basis
!         function lies in (T(I),T(I+K)), I=1,N).
!
!         If X=T(ILEFT+1) then VNIKX contains left limiting values
!         (left derivatives) at T(ILEFT+1).  In particular, ILEFT = N
!         produces left limiting values at the right end point
!         X=T(N+1).  To obtain left limiting values at T(I), I=K+1,N+1,
!         set X= next lower distinct knot, call INTRV to get ILEFT,
!         set X=T(I), and then call DBSPVD.
!
!     Description of Arguments
!         Input      T,X are double precision
!          T       - knot vector of length N+K, where
!                    N = number of B-spline basis functions
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K  >=  1
!          NDERIV  - number of derivatives = NDERIV-1,
!                    1  <=  NDERIV  <=  K
!          X       - argument of basis functions,
!                    T(K)  <=  X  <=  T(N+1)
!          ILEFT   - largest integer such that
!                    T(ILEFT)  <=  X  <   T(ILEFT+1)
!          LDVNIK  - leading dimension of matrix VNIKX
!
!         Output     VNIKX,WORK are double precision
!          VNIKX   - matrix of dimension at least (K,NDERIV) contain-
!                    ing the nonzero basis functions at X and their
!                    derivatives columnwise.
!          WORK    - a work vector of length (K+1)*(K+2)/2
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  DBSPVN, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSPVD
!
  INTEGER I,IDERIV,ILEFT,IPKMD,J,JJ,JLOW,JM,JP1MID,K,KMD, KP1, L, &
   LDUMMY, M, MHIGH, NDERIV
  DOUBLE PRECISION FACTOR, FKMD, T, V, VNIKX, WORK, X
!     DIMENSION T(ILEFT+K), WORK((K+1)*(K+2)/2)
!     A(I,J) = WORK(I+J*(J+1)/2),  I=1,J+1  J=1,K-1
!     A(I,K) = W0RK(I+K*(K-1)/2)  I=1.K
!     WORK(1) AND WORK((K+1)*(K+2)/2) ARE NOT USED.
  DIMENSION T(*), VNIKX(LDVNIK,*), WORK(*)
!***FIRST EXECUTABLE STATEMENT  DBSPVD
  if ( K < 1) go to 200
  if ( NDERIV < 1 .OR. NDERIV > K) go to 205
  if ( LDVNIK < K) go to 210
  IDERIV = NDERIV
  KP1 = K + 1
  JJ = KP1 - IDERIV
  call DBSPVN(T, JJ, K, 1, X, ILEFT, VNIKX, WORK, IWORK)
  if (IDERIV == 1) go to 100
  MHIGH = IDERIV
  DO 20 M=2,MHIGH
    JP1MID = 1
    DO 10 J=IDERIV,K
      VNIKX(J,IDERIV) = VNIKX(JP1MID,1)
      JP1MID = JP1MID + 1
   10   CONTINUE
    IDERIV = IDERIV - 1
    JJ = KP1 - IDERIV
    call DBSPVN(T, JJ, K, 2, X, ILEFT, VNIKX, WORK, IWORK)
   20 CONTINUE
!
  JM = KP1*(KP1+1)/2
  DO 30 L = 1,JM
    WORK(L) = 0.0D0
   30 CONTINUE
!     A(I,I) = WORK(I*(I+3)/2) = 1.0       I = 1,K
  L = 2
  J = 0
  DO 40 I = 1,K
    J = J + L
    WORK(J) = 1.0D0
    L = L + 1
   40 CONTINUE
  KMD = K
  DO 90 M=2,MHIGH
    KMD = KMD - 1
    FKMD = KMD
    I = ILEFT
    J = K
    JJ = J*(J+1)/2
    JM = JJ - J
    DO 60 LDUMMY=1,KMD
      IPKMD = I + KMD
      FACTOR = FKMD/(T(IPKMD)-T(I))
      DO 50 L=1,J
        WORK(L+JJ) = (WORK(L+JJ)-WORK(L+JM))*FACTOR
   50     CONTINUE
      I = I - 1
      J = J - 1
      JJ = JM
      JM = JM - J
   60   CONTINUE
!
    DO 80 I=1,K
      V = 0.0D0
      JLOW = MAX(I,M)
      JJ = JLOW*(JLOW+1)/2
      DO 70 J=JLOW,K
        V = WORK(I+JJ)*VNIKX(J,M) + V
        JJ = JJ + J + 1
   70     CONTINUE
      VNIKX(I,M) = V
   80   CONTINUE
   90 CONTINUE
  100 RETURN
!
!
  200 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVD', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  205 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVD', &
     'NDERIV DOES NOT SATISFY 1 <= NDERIV <= K', 2, 1)
  return
  210 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVD', &
     'LDVNIK DOES NOT SATISFY LDVNIK >= K', 2, 1)
  return
end
