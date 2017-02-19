DOUBLE PRECISION FUNCTION DBESK0 (X)
!
!! DBESK0 computes the modified (hyperbolic) Bessel function of the ...
!            third kind of order zero.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B1
!***TYPE      DOUBLE PRECISION (BESK0-S, DBESK0-D)
!***KEYWORDS  FNLIB, HYPERBOLIC BESSEL FUNCTION,
!             MODIFIED BESSEL FUNCTION, ORDER ZERO, SPECIAL FUNCTIONS,
!             THIRD KIND
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESK0(X) calculates the double precision modified (hyperbolic)
! Bessel function of the third kind of order zero for double
! precision argument X.  The argument must be greater than zero
! but not so large that the result underflows.
!
! Series for BK0        on the interval  0.          to  4.00000E+00
!                                        with weighted error   3.08E-33
!                                         log weighted error  32.51
!                               significant figures required  32.05
!                                    decimal places required  33.11
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DBESI0, DBSK0E, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESK0
  DOUBLE PRECISION X, BK0CS(16), XMAX, XMAXT, XSML, Y, &
    D1MACH, DCSEVL, DBESI0, DBSK0E
  LOGICAL FIRST
  SAVE BK0CS, NTK0, XSML, XMAX, FIRST
  DATA BK0CS(  1) / -.353273932339027687201140060063153D-1    /
  DATA BK0CS(  2) / +.344289899924628486886344927529213D+0    /
  DATA BK0CS(  3) / +.359799365153615016265721303687231D-1    /
  DATA BK0CS(  4) / +.126461541144692592338479508673447D-2    /
  DATA BK0CS(  5) / +.228621210311945178608269830297585D-4    /
  DATA BK0CS(  6) / +.253479107902614945730790013428354D-6    /
  DATA BK0CS(  7) / +.190451637722020885897214059381366D-8    /
  DATA BK0CS(  8) / +.103496952576336245851008317853089D-10   /
  DATA BK0CS(  9) / +.425981614279108257652445327170133D-13   /
  DATA BK0CS( 10) / +.137446543588075089694238325440000D-15   /
  DATA BK0CS( 11) / +.357089652850837359099688597333333D-18   /
  DATA BK0CS( 12) / +.763164366011643737667498666666666D-21   /
  DATA BK0CS( 13) / +.136542498844078185908053333333333D-23   /
  DATA BK0CS( 14) / +.207527526690666808319999999999999D-26   /
  DATA BK0CS( 15) / +.271281421807298560000000000000000D-29   /
  DATA BK0CS( 16) / +.308259388791466666666666666666666D-32   /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESK0
  if (FIRST) THEN
     NTK0 = INITDS (BK0CS, 16, 0.1*REAL(D1MACH(3)))
     XSML = SQRT(4.0D0*D1MACH(3))
     XMAXT = -LOG(D1MACH(1))
     XMAX = XMAXT - 0.5D0*XMAXT*LOG(XMAXT)/(XMAXT+0.5D0)
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBESK0', &
     'X IS ZERO OR NEGATIVE', 2, 2)
  if (X > 2.0D0) go to 20
!
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBESK0 = -LOG(0.5D0*X)*DBESI0(X) - 0.25D0 + DCSEVL (.5D0*Y-1.D0, &
    BK0CS, NTK0)
  return
!
 20   DBESK0 = 0.D0
  if (X  >  XMAX) call XERMSG ('SLATEC', 'DBESK0', &
     'X SO BIG K0 UNDERFLOWS', 1, 1)
  if (X > XMAX) RETURN
!
  DBESK0 = EXP(-X) * DBSK0E(X)
!
  return
end
