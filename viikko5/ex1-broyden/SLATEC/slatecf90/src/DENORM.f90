  DOUBLE PRECISION FUNCTION DENORM (N, X)
!
!! DENORM is subsidiary to DNSQ and DNSQE.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (ENORM-S, DENORM-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an N-vector X, this function calculates the
!     Euclidean norm of X.
!
!     The Euclidean norm is computed by accumulating the sum of
!     squares in three different sums. The sums of squares for the
!     small and large components are scaled so that no overflows
!     occur. Non-destructive underflows are permitted. Underflows
!     and overflows do not occur in the computation of the unscaled
!     sum of squares for the intermediate components.
!     The definitions of small, intermediate and large components
!     depend on two constants, RDWARF and RGIANT. The main
!     restrictions on these constants are that RDWARF**2 not
!     underflow and RGIANT**2 not overflow. The constants
!     given here are suitable for every known computer.
!
!     The function statement is
!
!       DOUBLE PRECISION FUNCTION DENORM(N,X)
!
!     where
!
!       N is a positive integer input variable.
!
!       X is an input array of length N.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DENORM
  INTEGER I, N
  DOUBLE PRECISION AGIANT, FLOATN, ONE, RDWARF, RGIANT, S1, S2, S3, &
       X(*), X1MAX, X3MAX, XABS, ZERO
  SAVE ONE, ZERO, RDWARF, RGIANT
  DATA ONE,ZERO,RDWARF,RGIANT /1.0D0,0.0D0,3.834D-20,1.304D19/
!***FIRST EXECUTABLE STATEMENT  DENORM
  S1 = ZERO
  S2 = ZERO
  S3 = ZERO
  X1MAX = ZERO
  X3MAX = ZERO
  FLOATN = N
  AGIANT = RGIANT/FLOATN
  DO 90 I = 1, N
     XABS = ABS(X(I))
     if (XABS  >  RDWARF .AND. XABS  <  AGIANT) go to 70
        if (XABS  <=  RDWARF) go to 30
!
!              SUM FOR LARGE COMPONENTS.
!
           if (XABS  <=  X1MAX) go to 10
              S1 = ONE + S1*(X1MAX/XABS)**2
              X1MAX = XABS
              go to 20
   10          CONTINUE
              S1 = S1 + (XABS/X1MAX)**2
   20          CONTINUE
           go to 60
   30       CONTINUE
!
!              SUM FOR SMALL COMPONENTS.
!
           if (XABS  <=  X3MAX) go to 40
              S3 = ONE + S3*(X3MAX/XABS)**2
              X3MAX = XABS
              go to 50
   40          CONTINUE
              if (XABS  /=  ZERO) S3 = S3 + (XABS/X3MAX)**2
   50          CONTINUE
   60       CONTINUE
        go to 80
   70    CONTINUE
!
!           SUM FOR INTERMEDIATE COMPONENTS.
!
        S2 = S2 + XABS**2
   80    CONTINUE
   90    CONTINUE
!
!     CALCULATION OF NORM.
!
  if (S1  ==  ZERO) go to 100
     DENORM = X1MAX*SQRT(S1+(S2/X1MAX)/X1MAX)
     go to 130
  100 CONTINUE
     if (S2  ==  ZERO) go to 110
        if (S2  >=  X3MAX) &
           DENORM = SQRT(S2*(ONE+(X3MAX/S2)*(X3MAX*S3)))
        if (S2  <  X3MAX) &
           DENORM = SQRT(X3MAX*((S2/X3MAX)+(X3MAX*S3)))
        go to 120
  110    CONTINUE
        DENORM = X3MAX*SQRT(S3)
  120    CONTINUE
  130 CONTINUE
  return
!
!     LAST CARD OF FUNCTION DENORM.
!
end
