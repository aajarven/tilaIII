FUNCTION DCHFIE (X1, X2, F1, F2, D1, D2, A, B)
!
!! DCHFIE evaluates the integral of a single cubic for DPCHIA.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (CHFIE-S, DCHFIE-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          DCHFIE:  Cubic Hermite Function Integral Evaluator.
!
!     Called by  DPCHIA  to evaluate the integral of a single cubic (in
!     Hermite form) over an arbitrary interval (A,B).
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
!        DOUBLE PRECISION  VALUE, DCHFIE
!
!        VALUE = DCHFIE (X1, X2, F1, F2, D1, D2, A, B)
!
!   Parameters:
!
!     VALUE -- (output) value of the requested integral.
!
!     X1,X2 -- (input) endpoints if interval of definition of cubic.
!
!     F1,F2 -- (input) function values at the ends of the interval.
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!
!     A,B -- (input) endpoints of interval of integration.
!
!***SEE ALSO  DPCHIA
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820730  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected subroutine name from DCHIV to DCHFIV.
!   870813  Minor cosmetic changes.
!   890411  1. Added SAVE statements (Vers. 3.2).
!           2. Added SIX to DOUBLE PRECISION declaration.
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR section in prologue.  (WRB)
!   930503  Corrected to set VALUE=0 when IERR.ne.0.  (FNF)
!   930504  Eliminated IERR and changed name DCHFIV to DCHFIE.  (FNF)
!***END PROLOGUE  DCHFIE
!
!  Programming notes:
!  1. There is no error return from this routine because zero is
!     indeed the mathematically correct answer when X1 == X2 .
!**End
!
!  DECLARE ARGUMENTS.
!
  DOUBLE PRECISION DCHFIE
  DOUBLE PRECISION  X1, X2, F1, F2, D1, D2, A, B
!
!  DECLARE LOCAL VARIABLES.
!
  DOUBLE PRECISION  DTERM, FOUR, FTERM, H, HALF, PHIA1, PHIA2, &
        PHIB1, PHIB2, PSIA1, PSIA2, PSIB1, PSIB2, SIX, TA1, TA2, &
        TB1, TB2, THREE, TWO, UA1, UA2, UB1, UB2
  SAVE HALF, TWO, THREE, FOUR, SIX
!
!  INITIALIZE.
!
  DATA  HALF/.5D0/, TWO/2.D0/, THREE/3.D0/, FOUR/4.D0/, SIX/6.D0/
!
!  VALIDITY CHECK INPUT.
!
!***FIRST EXECUTABLE STATEMENT  DCHFIE
  if (X1  ==  X2)  THEN
     DCHFIE = 0
  ELSE
     H = X2 - X1
     TA1 = (A - X1) / H
     TA2 = (X2 - A) / H
     TB1 = (B - X1) / H
     TB2 = (X2 - B) / H
!
     UA1 = TA1**3
     PHIA1 = UA1 * (TWO - TA1)
     PSIA1 = UA1 * (THREE*TA1 - FOUR)
     UA2 = TA2**3
     PHIA2 =  UA2 * (TWO - TA2)
     PSIA2 = -UA2 * (THREE*TA2 - FOUR)
!
     UB1 = TB1**3
     PHIB1 = UB1 * (TWO - TB1)
     PSIB1 = UB1 * (THREE*TB1 - FOUR)
     UB2 = TB2**3
     PHIB2 =  UB2 * (TWO - TB2)
     PSIB2 = -UB2 * (THREE*TB2 - FOUR)
!
     FTERM =   F1*(PHIA2 - PHIB2) + F2*(PHIB1 - PHIA1)
     DTERM = ( D1*(PSIA2 - PSIB2) + D2*(PSIB1 - PSIA1) )*(H/SIX)
!
     DCHFIE = (HALF*H) * (FTERM + DTERM)
  end if
!
  return
end
