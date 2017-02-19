subroutine ZKSCL (ZRR, ZRI, FNU, N, YR, YI, NZ, RZR, RZI, ASCLE, &
     TOL, ELIM)
!
!! ZKSCL is subsidiary to ZBESK.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (CKSCL-A, ZKSCL-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     SET K FUNCTIONS TO ZERO ON UNDERFLOW, CONTINUE RECURRENCE
!     ON SCALED FUNCTIONS UNTIL TWO MEMBERS COME ON SCALE, THEN
!     return WITH MIN(NZ+2,N) VALUES SCALED BY 1/TOL.
!
!***SEE ALSO  ZBESK
!***ROUTINES CALLED  ZABS, ZLOG, ZUCHK
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!   930122  Added ZLOG to EXTERNAL statement.  (RWC)
!***END PROLOGUE  ZKSCL
!     COMPLEX CK,CS,CY,CZERO,RZ,S1,S2,Y,ZR,ZD,CELM
  DOUBLE PRECISION ACS, AS, ASCLE, CKI, CKR, CSI, CSR, CYI, &
   CYR, ELIM, FN, FNU, RZI, RZR, STR, S1I, S1R, S2I, &
   S2R, TOL, YI, YR, ZEROI, ZEROR, ZRI, ZRR, ZABS, &
   ZDR, ZDI, CELMR, ELM, HELIM, ALAS
  INTEGER I, IC, IDUM, KK, N, NN, NW, NZ
  DIMENSION YR(N), YI(N), CYR(2), CYI(2)
  EXTERNAL ZABS, ZLOG
  DATA ZEROR,ZEROI / 0.0D0 , 0.0D0 /
!***FIRST EXECUTABLE STATEMENT  ZKSCL
  NZ = 0
  IC = 0
  NN = MIN(2,N)
  DO 10 I=1,NN
    S1R = YR(I)
    S1I = YI(I)
    CYR(I) = S1R
    CYI(I) = S1I
    AS = ZABS(S1R,S1I)
    ACS = -ZRR + LOG(AS)
    NZ = NZ + 1
    YR(I) = ZEROR
    YI(I) = ZEROI
    if (ACS < (-ELIM)) go to 10
    call ZLOG(S1R, S1I, CSR, CSI, IDUM)
    CSR = CSR - ZRR
    CSI = CSI - ZRI
    STR = EXP(CSR)/TOL
    CSR = STR*COS(CSI)
    CSI = STR*SIN(CSI)
    call ZUCHK(CSR, CSI, NW, ASCLE, TOL)
    if (NW /= 0) go to 10
    YR(I) = CSR
    YI(I) = CSI
    IC = I
    NZ = NZ - 1
   10 CONTINUE
  if (N == 1) RETURN
  if (IC > 1) go to 20
  YR(1) = ZEROR
  YI(1) = ZEROI
  NZ = 2
   20 CONTINUE
  if (N == 2) RETURN
  if (NZ == 0) RETURN
  FN = FNU + 1.0D0
  CKR = FN*RZR
  CKI = FN*RZI
  S1R = CYR(1)
  S1I = CYI(1)
  S2R = CYR(2)
  S2I = CYI(2)
  HELIM = 0.5D0*ELIM
  ELM = EXP(-ELIM)
  CELMR = ELM
  ZDR = ZRR
  ZDI = ZRI
!
!     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
!     S2 GETS LARGER THAN EXP(ELIM/2)
!
  DO 30 I=3,N
    KK = I
    CSR = S2R
    CSI = S2I
    S2R = CKR*CSR - CKI*CSI + S1R
    S2I = CKI*CSR + CKR*CSI + S1I
    S1R = CSR
    S1I = CSI
    CKR = CKR + RZR
    CKI = CKI + RZI
    AS = ZABS(S2R,S2I)
    ALAS = LOG(AS)
    ACS = -ZDR + ALAS
    NZ = NZ + 1
    YR(I) = ZEROR
    YI(I) = ZEROI
    if (ACS < (-ELIM)) go to 25
    call ZLOG(S2R, S2I, CSR, CSI, IDUM)
    CSR = CSR - ZDR
    CSI = CSI - ZDI
    STR = EXP(CSR)/TOL
    CSR = STR*COS(CSI)
    CSI = STR*SIN(CSI)
    call ZUCHK(CSR, CSI, NW, ASCLE, TOL)
    if (NW /= 0) go to 25
    YR(I) = CSR
    YI(I) = CSI
    NZ = NZ - 1
    if (IC == KK-1) go to 40
    IC = KK
    go to 30
   25   CONTINUE
    if ( ALAS < HELIM) go to 30
    ZDR = ZDR - ELIM
    S1R = S1R*CELMR
    S1I = S1I*CELMR
    S2R = S2R*CELMR
    S2I = S2I*CELMR
   30 CONTINUE
  NZ = N
  if ( IC == N) NZ=N-1
  go to 45
   40 CONTINUE
  NZ = KK - 2
   45 CONTINUE
  DO 50 I=1,NZ
    YR(I) = ZEROR
    YI(I) = ZEROI
   50 CONTINUE
  return
end
