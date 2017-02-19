subroutine BESK (X, FNU, KODE, N, Y, NZ)
!
!! BESK implements forward recursion on the three term recursion ...
!  relation for a sequence of non-negative order Bessel
!  functions K/SUB(FNU+I-1)/(X), or scaled Bessel functions
!  EXP(X)*K/SUB(FNU+I-1)/(X), I=1,...,N for real, positive
!  X and non-negative orders FNU.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10B3
!***TYPE      SINGLE PRECISION (BESK-S, DBESK-D)
!***KEYWORDS  K BESSEL FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         BESK implements forward recursion on the three term
!         recursion relation for a sequence of non-negative order Bessel
!         functions K/sub(FNU+I-1)/(X), or scaled Bessel functions
!         EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N for real X  >  0.0E0 and
!         non-negative orders FNU.  If FNU  <  NULIM, orders FNU and
!         FNU+1 are obtained from BESKNU to start the recursion.  If
!         FNU  >=  NULIM, the uniform asymptotic expansion is used for
!         orders FNU and FNU+1 to start the recursion.  NULIM is 35 or
!         70 depending on whether N=1 or N  >=  2.  Under and overflow
!         tests are made on the leading term of the asymptotic expansion
!         before any extensive computation is done.
!
!     Description of Arguments
!
!         Input
!           X      - X  >  0.0E0
!           FNU    - order of the initial K function, FNU  >=  0.0E0
!           KODE   - a parameter to indicate the scaling option
!                    KODE=1 returns Y(I)=       K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!                    KODE=2 returns Y(I)=EXP(X)*K/sub(FNU+I-1)/(X),
!                                        I=1,...,N
!           N      - number of members in the sequence, N  >=  1
!
!         Output
!           y      - a vector whose first n components contain values
!                    for the sequence
!                    Y(I)=       K/sub(FNU+I-1)/(X), I=1,...,N  or
!                    Y(I)=EXP(X)*K/sub(FNU+I-1)/(X), I=1,...,N
!                    depending on KODE
!           NZ     - number of components of Y set to zero due to
!                    underflow with KODE=1,
!                    NZ=0   , normal return, computation completed
!                    NZ  /=  0, first NZ components of Y set to zero
!                             due to underflow, Y(I)=0.0E0, I=1,...,NZ
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!         Underflow with KODE=1 -  a non-fatal error (NZ  /=  0)
!
!***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!***ROUTINES CALLED  ASYIK, BESK0, BESK0E, BESK1, BESK1E, BESKNU,
!                    I1MACH, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BESK
!
  INTEGER I, J, K, KODE, MZ, N, NB, ND, NN, NUD, NULIM, NZ
  INTEGER I1MACH
  REAL CN, DNU, ELIM, ETX, FLGIK,FN, FNN, FNU,GLN,GNU,RTZ,S,S1,S2, &
   T, TM, TRX, W, X, XLIM, Y, ZN
  REAL BESK0, BESK1, BESK1E, BESK0E, R1MACH
  DIMENSION W(2), NULIM(2), Y(*)
  SAVE NULIM
  DATA NULIM(1),NULIM(2) / 35 , 70 /
!***FIRST EXECUTABLE STATEMENT  BESK
  NN = -I1MACH(12)
  ELIM = 2.303E0*(NN*R1MACH(5)-3.0E0)
  XLIM = R1MACH(1)*1.0E+3
  if (KODE < 1 .OR. KODE > 2) go to 280
  if (FNU < 0.0E0) go to 290
  if (X <= 0.0E0) go to 300
  if (X < XLIM) go to 320
  if (N < 1) go to 310
  ETX = KODE - 1
!
!     ND IS A DUMMY VARIABLE FOR N
!     GNU IS A DUMMY VARIABLE FOR FNU
!     NZ = NUMBER OF UNDERFLOWS ON KODE=1
!
  ND = N
  NZ = 0
  NUD = INT(FNU)
  DNU = FNU - NUD
  GNU = FNU
  NN = MIN(2,ND)
  FN = FNU + N - 1
  FNN = FN
  if (FN < 2.0E0) go to 150
!
!     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE LAST ORDER, FNU+N-1 >= NULIM
!
  ZN = X/FN
  if (ZN == 0.0E0) go to 320
  RTZ = SQRT(1.0E0+ZN*ZN)
  GLN = LOG((1.0E0+RTZ)/ZN)
  T = RTZ*(1.0E0-ETX) + ETX/(ZN+RTZ)
  CN = -FN*(T-GLN)
  if (CN > ELIM) go to 320
  if (NUD < NULIM(NN)) go to 30
  if (NN == 1) go to 20
   10 CONTINUE
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE FIRST ORDER, FNU >= NULIM
!
  FN = GNU
  ZN = X/FN
  RTZ = SQRT(1.0E0+ZN*ZN)
  GLN = LOG((1.0E0+RTZ)/ZN)
  T = RTZ*(1.0E0-ETX) + ETX/(ZN+RTZ)
  CN = -FN*(T-GLN)
   20 CONTINUE
  if (CN < -ELIM) go to 230
!
!     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1 >= NULIM
!
  FLGIK = -1.0E0
  call ASYIK(X,GNU,KODE,FLGIK,RTZ,CN,NN,Y)
  if (NN == 1) go to 240
  TRX = 2.0E0/X
  TM = (GNU+GNU+2.0E0)/X
  go to 130
!
   30 CONTINUE
  if (KODE == 2) go to 40
!
!     UNDERFLOW TEST (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION IN X)
!     FOR ORDER DNU
!
  if (X > ELIM) go to 230
   40 CONTINUE
  if (DNU /= 0.0E0) go to 80
  if (KODE == 2) go to 50
  S1 = BESK0(X)
  go to 60
   50 S1 = BESK0E(X)
   60 CONTINUE
  if (NUD == 0 .AND. ND == 1) go to 120
  if (KODE == 2) go to 70
  S2 = BESK1(X)
  go to 90
   70 S2 = BESK1E(X)
  go to 90
   80 CONTINUE
  NB = 2
  if (NUD == 0 .AND. ND == 1) NB = 1
  call BESKNU(X, DNU, KODE, NB, W, NZ)
  S1 = W(1)
  if (NB == 1) go to 120
  S2 = W(2)
   90 CONTINUE
  TRX = 2.0E0/X
  TM = (DNU+DNU+2.0E0)/X
!     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
  if (ND == 1) NUD = NUD - 1
  if (NUD > 0) go to 100
  if (ND > 1) go to 120
  S1 = S2
  go to 120
  100 CONTINUE
  DO 110 I=1,NUD
    S = S2
    S2 = TM*S2 + S1
    S1 = S
    TM = TM + TRX
  110 CONTINUE
  if (ND == 1) S1 = S2
  120 CONTINUE
  Y(1) = S1
  if (ND == 1) go to 240
  Y(2) = S2
  130 CONTINUE
  if (ND == 2) go to 240
!     FORWARD RECUR FROM FNU+2 TO FNU+N-1
  DO 140 I=3,ND
    Y(I) = TM*Y(I-1) + Y(I-2)
    TM = TM + TRX
  140 CONTINUE
  go to 240
!
  150 CONTINUE
!     UNDERFLOW TEST FOR KODE=1
  if (KODE == 2) go to 160
  if (X > ELIM) go to 230
  160 CONTINUE
!     OVERFLOW TEST
  if (FN <= 1.0E0) go to 170
  if (-FN*(LOG(X)-0.693E0) > ELIM) go to 320
  170 CONTINUE
  if (DNU == 0.0E0) go to 180
  call BESKNU(X, FNU, KODE, ND, Y, MZ)
  go to 240
  180 CONTINUE
  J = NUD
  if (J == 1) go to 210
  J = J + 1
  if (KODE == 2) go to 190
  Y(J) = BESK0(X)
  go to 200
  190 Y(J) = BESK0E(X)
  200 if (ND == 1) go to 240
  J = J + 1
  210 if (KODE == 2) go to 220
  Y(J) = BESK1(X)
  go to 240
  220 Y(J) = BESK1E(X)
  go to 240
!
!     UPDATE PARAMETERS ON UNDERFLOW
!
  230 CONTINUE
  NUD = NUD + 1
  ND = ND - 1
  if (ND == 0) go to 240
  NN = MIN(2,ND)
  GNU = GNU + 1.0E0
  if (FNN < 2.0E0) go to 230
  if (NUD < NULIM(NN)) go to 230
  go to 10
  240 CONTINUE
  NZ = N - ND
  if (NZ == 0) RETURN
  if (ND == 0) go to 260
  DO 250 I=1,ND
    J = N - I + 1
    K = ND - I + 1
    Y(J) = Y(K)
  250 CONTINUE
  260 CONTINUE
  DO 270 I=1,NZ
    Y(I) = 0.0E0
  270 CONTINUE
  return
!
!
!
  280 CONTINUE
  call XERMSG ('SLATEC', 'BESK', 'SCALING OPTION, KODE, NOT 1 OR 2' &
     , 2, 1)
  return
  290 CONTINUE
  call XERMSG ('SLATEC', 'BESK', 'ORDER, FNU, LESS THAN ZERO', 2, &
     1)
  return
  300 CONTINUE
  call XERMSG ('SLATEC', 'BESK', 'X LESS THAN OR EQUAL TO ZERO', 2, &
     1)
  return
  310 CONTINUE
  call XERMSG ('SLATEC', 'BESK', 'N LESS THAN ONE', 2, 1)
  return
  320 CONTINUE
  call XERMSG ('SLATEC', 'BESK', &
     'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
  return
end
