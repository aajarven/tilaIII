subroutine DBKISR (X, N, SUM, IERR)
!
!! DBKISR is subsidiary to DBSKIN.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BKISR-S, DBKISR-D)
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     DBKISR computes repeated integrals of the K0 Bessel function
!     by the series for N=0,1, and 2.
!
!***SEE ALSO  DBSKIN
!***ROUTINES CALLED  D1MACH, DPSIXN
!***REVISION HISTORY  (YYMMDD)
!   820601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DBKISR
  INTEGER I, IERR, K, KK, KKN, K1, N, NP
  DOUBLE PRECISION AK, ATOL, BK, C, FK, FN, HX, HXS, POL, PR, SUM, &
   TKP, TOL, TRM, X, XLN
  DOUBLE PRECISION DPSIXN, D1MACH
  DIMENSION C(2)
  SAVE C
!
  DATA C(1), C(2) /1.57079632679489662D+00,1.0D0/
!***FIRST EXECUTABLE STATEMENT  DBKISR
  IERR=0
  TOL = MAX(D1MACH(4),1.0D-18)
  if (X < TOL) go to 50
  PR = 1.0D0
  POL = 0.0D0
  if (N == 0) go to 20
  DO 10 I=1,N
    POL = -POL*X + C(I)
    PR = PR*X/I
   10 CONTINUE
   20 CONTINUE
  HX = X*0.5D0
  HXS = HX*HX
  XLN = LOG(HX)
  NP = N + 1
  TKP = 3.0D0
  FK = 2.0D0
  FN = N
  BK = 4.0D0
  AK = 2.0D0/((FN+1.0D0)*(FN+2.0D0))
  SUM = AK*(DPSIXN(N+3)-DPSIXN(3)+DPSIXN(2)-XLN)
  ATOL = SUM*TOL*0.75D0
  DO 30 K=2,20
    AK = AK*(HXS/BK)*((TKP+1.0D0)/(TKP+FN+1.0D0))*(TKP/(TKP+FN))
    K1 = K + 1
    KK = K1 + K
    KKN = KK + N
    TRM = (DPSIXN(K1)+DPSIXN(KKN)-DPSIXN(KK)-XLN)*AK
    SUM = SUM + TRM
    if (ABS(TRM) <= ATOL) go to 40
    TKP = TKP + 2.0D0
    BK = BK + TKP
    FK = FK + 1.0D0
   30 CONTINUE
  go to 80
   40 CONTINUE
  SUM = (SUM*HXS+DPSIXN(NP)-XLN)*PR
  if (N == 1) SUM = -SUM
  SUM = POL + SUM
  return
!-----------------------------------------------------------------------
!     SMALL X CASE, X < WORD TOLERANCE
!-----------------------------------------------------------------------
   50 CONTINUE
  if (N > 0) go to 60
  HX = X*0.5D0
  SUM = DPSIXN(1) - LOG(HX)
  return
   60 CONTINUE
  SUM = C(N)
  return
   80 CONTINUE
  IERR=2
  return
end
