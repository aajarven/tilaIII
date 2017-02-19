FUNCTION QWGTC (X, C, P2, P3, P4, KP)
!
!! QWGTC defines the WEIGHT function for QAWC.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (QWGTC-S, DQWGTC-D)
!***KEYWORDS  CAUCHY PRINCIPAL VALUE, WEIGHT FUNCTION
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***SEE ALSO  QK15W
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   830518  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  QWGTC
!
  REAL QWGTC
  REAL C,P2,P3,P4,X
  INTEGER KP
!***FIRST EXECUTABLE STATEMENT  QWGTC
  QWGTC = 0.1E+01/(X-C)
  return
end
