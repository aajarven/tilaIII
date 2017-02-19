subroutine CKSCL (ZR, FNU, N, Y, NZ, RZ, ASCLE, TOL, ELIM)
!
!! CKSCL is subsidiary to CBKNU, CUNK1 and CUNK2.
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
!***SEE ALSO  CBKNU, CUNK1, CUNK2
!***ROUTINES CALLED  CUCHK
!***REVISION HISTORY  (YYMMDD)
!   ??????  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CKSCL
  COMPLEX CK, CS, CY, CZERO, RZ, S1, S2, Y, ZR, ZD, CELM
  REAL AA, ASCLE, ACS, AS, CSI, CSR, ELIM, FN, FNU, TOL, XX, ZRI, &
   ELM, ALAS, HELIM
  INTEGER I, IC, K, KK, N, NN, NW, NZ
  DIMENSION Y(N), CY(2)
  DATA CZERO / (0.0E0,0.0E0) /
!***FIRST EXECUTABLE STATEMENT  CUCHK
  NZ = 0
  IC = 0
  XX = REAL(ZR)
  NN = MIN(2,N)
  DO 10 I=1,NN
    S1 = Y(I)
    CY(I) = S1
    AS = ABS(S1)
    ACS = -XX + ALOG(AS)
    NZ = NZ + 1
    Y(I) = CZERO
    if (ACS < (-ELIM)) go to 10
    CS = -ZR + CLOG(S1)
    CSR = REAL(CS)
    CSI = AIMAG(CS)
    AA = EXP(CSR)/TOL
    CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
    call CUCHK(CS, NW, ASCLE, TOL)
    if (NW /= 0) go to 10
    Y(I) = CS
    NZ = NZ - 1
    IC = I
   10 CONTINUE
  if (N == 1) RETURN
  if (IC > 1) go to 20
  Y(1) = CZERO
  NZ = 2
   20 CONTINUE
  if (N == 2) RETURN
  if (NZ == 0) RETURN
  FN = FNU + 1.0E0
  CK = CMPLX(FN,0.0E0)*RZ
  S1 = CY(1)
  S2 = CY(2)
  HELIM = 0.5E0*ELIM
  ELM = EXP(-ELIM)
  CELM = CMPLX(ELM,0.0E0)
  ZRI =AIMAG(ZR)
  ZD = ZR
!
!     FIND TWO CONSECUTIVE Y VALUES ON SCALE. SCALE RECURRENCE IF
!     S2 GETS LARGER THAN EXP(ELIM/2)
!
  DO 30 I=3,N
    KK = I
    CS = S2
    S2 = CK*S2 + S1
    S1 = CS
    CK = CK + RZ
    AS = ABS(S2)
    ALAS = ALOG(AS)
    ACS = -XX + ALAS
    NZ = NZ + 1
    Y(I) = CZERO
    if (ACS < (-ELIM)) go to 25
    CS = -ZD + CLOG(S2)
    CSR = REAL(CS)
    CSI = AIMAG(CS)
    AA = EXP(CSR)/TOL
    CS = CMPLX(AA,0.0E0)*CMPLX(COS(CSI),SIN(CSI))
    call CUCHK(CS, NW, ASCLE, TOL)
    if (NW /= 0) go to 25
    Y(I) = CS
    NZ = NZ - 1
    if (IC == (KK-1)) go to 40
    IC = KK
    go to 30
   25   CONTINUE
    if ( ALAS < HELIM) go to 30
    XX = XX-ELIM
    S1 = S1*CELM
    S2 = S2*CELM
    ZD = CMPLX(XX,ZRI)
   30 CONTINUE
  NZ = N
  if ( IC == N) NZ=N-1
  go to 45
   40 CONTINUE
  NZ = KK - 2
   45 CONTINUE
  DO 50 K=1,NZ
    Y(K) = CZERO
   50 CONTINUE
  return
end
