FUNCTION PCHST (ARG1, ARG2)
!
!! PCHST is the PCHIP Sign-Testing Routine
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHST-S, DPCHST-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         PCHST:  PCHIP Sign-Testing Routine.
!
!     Returns:
!        -1. if ARG1 and ARG2 are of opposite sign.
!         0. if either argument is zero.
!        +1. if ARG1 and ARG2 are of the same sign.
!
!     The object is to do this without multiplying ARG1*ARG2, to avoid
!     possible over/underflow problems.
!
!  Fortran intrinsics used:  SIGN.
!
!***SEE ALSO  PCHCE, PCHCI, PCHCS, PCHIM
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHST
!
!**End
!
!  DECLARE ARGUMENTS.
!
  REAL PCHST
  REAL  ARG1, ARG2
!
!  DECLARE LOCAL VARIABLES.
!
  REAL  ONE, ZERO
  SAVE ZERO, ONE
  DATA  ZERO /0./,  ONE /1./
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  PCHST
  PCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
  if ((ARG1 == ZERO) .OR. (ARG2 == ZERO))  PCHST = ZERO
!
  return
!------------- LAST LINE OF PCHST FOLLOWS ------------------------------
end
