subroutine DBESY (X, FNU, N, Y)
!
!! DBESY implements forward recursion on the three term recursion ...
!            relation for a sequence of non-negative order Bessel ...
!            functions Y/SUB(FNU+I-1)/(X), I=1,...,N for real, positive ...
!            X and non-negative orders FNU.
!
!***LIBRARY   SLATEC
!***CATEGORY  C10A3
!***TYPE      DOUBLE PRECISION (BESY-S, DBESY-D)
!***KEYWORDS  SPECIAL FUNCTIONS, Y BESSEL FUNCTION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  **** a double precision routine ****
!         DBESY implements forward recursion on the three term
!         recursion relation for a sequence of non-negative order Bessel
!         functions Y/sub(FNU+I-1)/(X), I=1,N for real X  >  0.0D0 and
!         non-negative orders FNU.  If FNU  <  NULIM, orders FNU and
!         FNU+1 are obtained from DBSYNU which computes by a power
!         series for X  <=  2, the K Bessel function of an imaginary
!         argument for 2  <  X  <=  20 and the asymptotic expansion for
!         X  >  20.
!
!         If FNU  >=  NULIM, the uniform asymptotic expansion is coded
!         in DASYJY for orders FNU and FNU+1 to start the recursion.
!         NULIM is 70 or 100 depending on whether N=1 or N  >=  2.  An
!         overflow test is made on the leading term of the asymptotic
!         expansion before any extensive computation is done.
!
!         The maximum number of significant digits obtainable
!         is the smaller of 14 and the number of digits carried in
!         double precision arithmetic.
!
!     Description of Arguments
!
!         Input
!           X      - X  >  0.0D0
!           FNU    - order of the initial Y function, FNU  >=  0.0D0
!           N      - number of members in the sequence, N  >=  1
!
!         Output
!           Y      - a vector whose first N components contain values
!                    for the sequence Y(I)=Y/sub(FNU+I-1)/(X), I=1,N.
!
!     Error Conditions
!         Improper input arguments - a fatal error
!         Overflow - a fatal error
!
!***REFERENCES  F. W. J. Olver, Tables of Bessel Functions of Moderate
!                 or Large Orders, NPL Mathematical Tables 6, Her
!                 Majesty's Stationery Office, London, 1962.
!               N. M. Temme, On the numerical evaluation of the modified
!                 Bessel function of the third kind, Journal of
!                 Computational Physics 19, (1975), pp. 324-337.
!               N. M. Temme, On the numerical evaluation of the ordinary
!                 Bessel function of the second kind, Journal of
!                 Computational Physics 21, (1976), pp. 343-350.
!***ROUTINES CALLED  D1MACH, DASYJY, DBESY0, DBESY1, DBSYNU, DYAIRY,
!                    I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBESY
!
  EXTERNAL DYAIRY
  INTEGER I, IFLW, J, N, NB, ND, NN, NUD, NULIM
  INTEGER I1MACH
  DOUBLE PRECISION AZN,CN,DNU,ELIM,FLGJY,FN,FNU,RAN,S,S1,S2,TM,TRX, &
             W,WK,W2N,X,XLIM,XXN,Y
  DOUBLE PRECISION DBESY0, DBESY1, D1MACH
  DIMENSION W(2), NULIM(2), Y(*), WK(7)
  SAVE NULIM
  DATA NULIM(1),NULIM(2) / 70 , 100 /
!***FIRST EXECUTABLE STATEMENT  DBESY
  NN = -I1MACH(15)
  ELIM = 2.303D0*(NN*D1MACH(5)-3.0D0)
  XLIM = D1MACH(1)*1.0D+3
  if (FNU < 0.0D0) go to 140
  if (X <= 0.0D0) go to 150
  if (X < XLIM) go to 170
  if (N < 1) go to 160
!
!     ND IS A DUMMY VARIABLE FOR N
!
  ND = N
  NUD = INT(FNU)
  DNU = FNU - NUD
  NN = MIN(2,ND)
  FN = FNU + N - 1
  if (FN < 2.0D0) go to 100
!
!     OVERFLOW TEST  (LEADING EXPONENTIAL OF ASYMPTOTIC EXPANSION)
!     FOR THE LAST ORDER, FNU+N-1 >= NULIM
!
  XXN = X/FN
  W2N = 1.0D0-XXN*XXN
  if ( W2N <= 0.0D0) go to 10
  RAN = SQRT(W2N)
  AZN = LOG((1.0D0+RAN)/XXN) - RAN
  CN = FN*AZN
  if ( CN > ELIM) go to 170
   10 CONTINUE
  if (NUD < NULIM(NN)) go to 20
!
!     ASYMPTOTIC EXPANSION FOR ORDERS FNU AND FNU+1 >= NULIM
!
  FLGJY = -1.0D0
  call DASYJY(DYAIRY,X,FNU,FLGJY,NN,Y,WK,IFLW)
  if ( IFLW /= 0) go to 170
  if (NN == 1) RETURN
  TRX = 2.0D0/X
  TM = (FNU+FNU+2.0D0)/X
  go to 80
!
   20 CONTINUE
  if (DNU /= 0.0D0) go to 30
  S1 = DBESY0(X)
  if (NUD == 0 .AND. ND == 1) go to 70
  S2 = DBESY1(X)
  go to 40
   30 CONTINUE
  NB = 2
  if (NUD == 0 .AND. ND == 1) NB = 1
  call DBSYNU(X, DNU, NB, W)
  S1 = W(1)
  if (NB == 1) go to 70
  S2 = W(2)
   40 CONTINUE
  TRX = 2.0D0/X
  TM = (DNU+DNU+2.0D0)/X
!     FORWARD RECUR FROM DNU TO FNU+1 TO GET Y(1) AND Y(2)
  if (ND == 1) NUD = NUD - 1
  if (NUD > 0) go to 50
  if (ND > 1) go to 70
  S1 = S2
  go to 70
   50 CONTINUE
  DO 60 I=1,NUD
    S = S2
    S2 = TM*S2 - S1
    S1 = S
    TM = TM + TRX
   60 CONTINUE
  if (ND == 1) S1 = S2
   70 CONTINUE
  Y(1) = S1
  if (ND == 1) RETURN
  Y(2) = S2
   80 CONTINUE
  if (ND == 2) RETURN
!     FORWARD RECUR FROM FNU+2 TO FNU+N-1
  DO 90 I=3,ND
    Y(I) = TM*Y(I-1) - Y(I-2)
    TM = TM + TRX
   90 CONTINUE
  return
!
  100 CONTINUE
!     OVERFLOW TEST
  if (FN <= 1.0D0) go to 110
  if (-FN*(LOG(X)-0.693D0) > ELIM) go to 170
  110 CONTINUE
  if (DNU == 0.0D0) go to 120
  call DBSYNU(X, FNU, ND, Y)
  return
  120 CONTINUE
  J = NUD
  if (J == 1) go to 130
  J = J + 1
  Y(J) = DBESY0(X)
  if (ND == 1) RETURN
  J = J + 1
  130 CONTINUE
  Y(J) = DBESY1(X)
  if (ND == 1) RETURN
  TRX = 2.0D0/X
  TM = TRX
  go to 80
!
!
!
  140 CONTINUE
  call XERMSG ('SLATEC', 'DBESY', 'ORDER, FNU, LESS THAN ZERO', 2, &
     1)
  return
  150 CONTINUE
  call XERMSG ('SLATEC', 'DBESY', 'X LESS THAN OR EQUAL TO ZERO', &
     2, 1)
  return
  160 CONTINUE
  call XERMSG ('SLATEC', 'DBESY', 'N LESS THAN ONE', 2, 1)
  return
  170 CONTINUE
  call XERMSG ('SLATEC', 'DBESY', &
     'OVERFLOW, FNU OR N TOO LARGE OR X TOO SMALL', 6, 1)
  return
end
