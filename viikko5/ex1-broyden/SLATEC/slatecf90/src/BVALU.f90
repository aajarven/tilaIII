function BVALU (T, A, N, K, IDERIV, X, INBV, WORK)
!
!! BVALU evaluates the B-representation of a B-spline at X for the
!  function value or any of its derivatives.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      SINGLE PRECISION (BVALU-S, DBVALU-D)
!***KEYWORDS  DIFFERENTIATION OF B-SPLINE, EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract
!         BVALU is the BVALUE function of the reference.
!
!         BVALU evaluates the B-representation (T,A,N,K) of a B-spline
!         at X for the function value on IDERIV = 0 or any of its
!         derivatives on IDERIV = 1,2,...,K-1.  Right limiting values
!         (right derivatives) are returned except at the right end
!         point X=T(N+1) where left limiting values are computed.  The
!         spline is defined on T(K)  <=  X  <=  T(N+1).  BVALU returns
!         a fatal error message when X is outside of this interval.
!
!         To compute left derivatives or left limiting values at a
!         knot T(I), replace N by I-1 and set X=T(I), I=K+1,N+1.
!
!         BVALU calls INTRV
!
!     Description of Arguments
!         Input
!          T       - knot vector of length N+K
!          A       - B-spline coefficient vector of length N
!          N       - number of B-spline coefficients
!                    N = sum of knot multiplicities-K
!          K       - order of the B-spline, K  >=  1
!          IDERIV  - order of the derivative, 0  <=  IDERIV  <=  K-1
!                    IDERIV=0 returns the B-spline value
!          X       - argument, T(K)  <=  X  <=  T(N+1)
!          INBV    - an initialization parameter which must be set
!                    to 1 the first time BVALU is called.
!
!         Output
!          INBV    - INBV contains information for efficient process-
!                    ing after the initial call and INBV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INBV parameters.
!          WORK    - work vector of length 3*K.
!          BVALU   - value of the IDERIV-th derivative at X
!
!     Error Conditions
!         An improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  INTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BVALU
!
  INTEGER I,IDERIV,IDERP1,IHI,IHMKMJ,ILO,IMK,IMKPJ, INBV, IPJ, &
   IP1, IP1MJ, J, JJ, J1, J2, K, KMIDER, KMJ, KM1, KPK, MFLAG, N
  REAL A, FKMJ, T, WORK, X
!     DIMENSION T(N+K), WORK(3*K)
  DIMENSION T(*), A(*), WORK(*)
!***FIRST EXECUTABLE STATEMENT  BVALU
  BVALU = 0.0E0
  if ( K < 1) go to 102
  if ( N < K) go to 101
  if ( IDERIV < 0 .OR. IDERIV >= K) go to 110
  KMIDER = K - IDERIV
!
! *** FIND *I* IN (K,N) SUCH THAT T(I)  <=  X  <  T(I+1)
!     (OR,  <=  T(I+1) if T(I)  <  T(I+1) = T(N+1)).
  KM1 = K - 1
  call INTRV(T, N+1, X, INBV, I, MFLAG)
  if (X < T(K)) go to 120
  if (MFLAG == 0) go to 20
  if (X > T(I)) go to 130
   10 if (I == K) go to 140
  I = I - 1
  if (X == T(I)) go to 10
!
! *** DIFFERENCE THE COEFFICIENTS *IDERIV* TIMES
!     WORK(I) = AJ(I), WORK(K+I) = DP(I), WORK(K+K+I) = DM(I), I=1.K
!
   20 IMK = I - K
  DO 30 J=1,K
    IMKPJ = IMK + J
    WORK(J) = A(IMKPJ)
   30 CONTINUE
  if (IDERIV == 0) go to 60
  DO 50 J=1,IDERIV
    KMJ = K - J
    FKMJ = KMJ
    DO 40 JJ=1,KMJ
      IHI = I + JJ
      IHMKMJ = IHI - KMJ
      WORK(JJ) = (WORK(JJ+1)-WORK(JJ))/(T(IHI)-T(IHMKMJ))*FKMJ
   40   CONTINUE
   50 CONTINUE
!
! *** COMPUTE VALUE AT *X* IN (T(I),(T(I+1)) OF IDERIV-TH DERIVATIVE,
!     GIVEN ITS RELEVANT B-SPLINE COEFF. IN AJ(1),...,AJ(K-IDERIV).
   60 if (IDERIV == KM1) go to 100
  IP1 = I + 1
  KPK = K + K
  J1 = K + 1
  J2 = KPK + 1
  DO 70 J=1,KMIDER
    IPJ = I + J
    WORK(J1) = T(IPJ) - X
    IP1MJ = IP1 - J
    WORK(J2) = X - T(IP1MJ)
    J1 = J1 + 1
    J2 = J2 + 1
   70 CONTINUE
  IDERP1 = IDERIV + 1
  DO 90 J=IDERP1,KM1
    KMJ = K - J
    ILO = KMJ
    DO 80 JJ=1,KMJ
      WORK(JJ) = (WORK(JJ+1)*WORK(KPK+ILO)+WORK(JJ) &
                *WORK(K+JJ))/(WORK(KPK+ILO)+WORK(K+JJ))
      ILO = ILO - 1
   80   CONTINUE
   90 CONTINUE
  100 BVALU = WORK(1)
  return
!
!
  101 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
  102 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', &
     'IDERIV DOES NOT SATISFY 0 <= IDERIV < K', 2, 1)
  return
  120 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', &
     'X IS N0T GREATER THAN OR EQUAL TO T(K)', 2, 1)
  return
  130 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', &
     'X IS NOT LESS THAN OR EQUAL TO T(N+1)', 2, 1)
  return
  140 CONTINUE
  call XERMSG ('SLATEC', 'BVALU', &
     'A LEFT LIMITING VALUE CANNOT BE OBTAINED AT T(K)', 2, 1)
  return
end
