  DOUBLE PRECISION FUNCTION DPCHST (ARG1, ARG2)
!
!! DPCHST is the DPCHIP Sign-Testing Routine.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (PCHST-S, DPCHST-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         DPCHST:  DPCHIP Sign-Testing Routine.
!
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
!***SEE ALSO  DPCHCE, DPCHCI, DPCHCS, DPCHIM
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   811103  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHST
!
!**End
!
!  DECLARE ARGUMENTS.
!
  DOUBLE PRECISION  ARG1, ARG2
!
!  DECLARE LOCAL VARIABLES.
!
  DOUBLE PRECISION  ONE, ZERO
  SAVE ZERO, ONE
  DATA  ZERO /0.D0/,  ONE/1.D0/
!
!  PERFORM THE TEST.
!
!***FIRST EXECUTABLE STATEMENT  DPCHST
  DPCHST = SIGN(ONE,ARG1) * SIGN(ONE,ARG2)
  if ((ARG1 == ZERO) .OR. (ARG2 == ZERO))  DPCHST = ZERO
!
  return
end
