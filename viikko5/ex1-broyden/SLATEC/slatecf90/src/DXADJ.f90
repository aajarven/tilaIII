subroutine DXADJ (X, IX, IERROR)
!
!! DXADJ provides double-precision floating-point arithmetic ...
!            with an extended exponent range.
!
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      DOUBLE PRECISION (XADJ-S, DXADJ-D)
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     DOUBLE PRECISION X
!     INTEGER IX
!
!                  TRANSFORMS (X,IX) SO THAT
!                  RADIX**(-L)  <=  ABS(X)  <  RADIX**L.
!                  ON MOST COMPUTERS THIS TRANSFORMATION DOES
!                  NOT CHANGE THE MANTISSA OF X PROVIDED RADIX IS
!                  THE NUMBER BASE OF DOUBLE-PRECISION ARITHMETIC.
!
!***SEE ALSO  DXSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***COMMON BLOCKS    DXBLK2
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!           CALLs to XERROR changed to CALLs to XERMSG.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXADJ
  DOUBLE PRECISION X
  INTEGER IX
  DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
  INTEGER L, L2, KMAX
  COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
  SAVE /DXBLK2/
!
!   THE CONDITION IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! IS
!     2*L  <=  KMAX
!
! THIS CONDITION MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE DXSET.
!
!***FIRST EXECUTABLE STATEMENT  DXADJ
  IERROR=0
  if (X == 0.0D0) go to 50
  if (ABS(X) >= 1.0D0) go to 20
  if (RADIXL*ABS(X) >= 1.0D0) go to 60
  X = X*RAD2L
  if (IX < 0) go to 10
  IX = IX - L2
  go to 70
   10 if (IX < -KMAX+L2) go to 40
  IX = IX - L2
  go to 70
   20 if (ABS(X) < RADIXL) go to 60
  X = X/RAD2L
  if (IX > 0) go to 30
  IX = IX + L2
  go to 70
   30 if (IX > KMAX-L2) go to 40
  IX = IX + L2
  go to 70
   40 call XERMSG ('SLATEC', 'DXADJ', 'overflow in auxiliary index', &
               207, 1)
  IERROR=207
  return
   50 IX = 0
   60 if (ABS(IX) > KMAX) go to 40
   70 RETURN
end
