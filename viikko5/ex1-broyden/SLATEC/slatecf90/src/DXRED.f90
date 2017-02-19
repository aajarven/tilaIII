subroutine DXRED (X, IX, IERROR)
!
!! DXRED provides double-precision floating-point arithmetic ...
!            with an extended exponent range.
!
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      DOUBLE PRECISION (XRED-S, DXRED-D)
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     DOUBLE PRECISION X
!     INTEGER IX
!
!                  IF
!                  RADIX**(-2L)  <=  (ABS(X),IX)  <=  RADIX**(2L)
!                  THEN DXRED TRANSFORMS (X,IX) SO THAT IX=0.
!                  if (X,IX) IS OUTSIDE THE ABOVE RANGE,
!                  THEN DXRED TAKES NO ACTION.
!                  THIS SUBROUTINE IS USEFUL if THE
!                  RESULTS OF EXTENDED-RANGE CALCULATIONS
!                  ARE TO BE USED IN SUBSEQUENT ORDINARY
!                  DOUBLE-PRECISION CALCULATIONS.
!
!***SEE ALSO  DXSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DXBLK2
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXRED
  DOUBLE PRECISION X
  INTEGER IX
  DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R, XA
  INTEGER L, L2, KMAX
  COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
  SAVE /DXBLK2/
!
!***FIRST EXECUTABLE STATEMENT  DXRED
  IERROR=0
  if (X == 0.0D0) go to 90
  XA = ABS(X)
  if (IX == 0) go to 70
  IXA = ABS(IX)
  IXA1 = IXA/L2
  IXA2 = MOD(IXA,L2)
  if (IX > 0) go to 40
   10 CONTINUE
  if (XA > 1.0D0) go to 20
  XA = XA*RAD2L
  IXA1 = IXA1 + 1
  go to 10
   20 XA = XA/RADIX**IXA2
  if (IXA1 == 0) go to 70
  DO 30 I=1,IXA1
    if (XA < 1.0D0) go to 100
    XA = XA/RAD2L
   30 CONTINUE
  go to 70
!
   40 CONTINUE
  if (XA < 1.0D0) go to 50
  XA = XA/RAD2L
  IXA1 = IXA1 + 1
  go to 40
   50 XA = XA*RADIX**IXA2
  if (IXA1 == 0) go to 70
  DO 60 I=1,IXA1
    if (XA > 1.0D0) go to 100
    XA = XA*RAD2L
   60 CONTINUE
   70 if (XA > RAD2L) go to 100
  if (XA > 1.0D0) go to 80
  if (RAD2L*XA < 1.0D0) go to 100
   80 X = SIGN(XA,X)
   90 IX = 0
  100 RETURN
end
