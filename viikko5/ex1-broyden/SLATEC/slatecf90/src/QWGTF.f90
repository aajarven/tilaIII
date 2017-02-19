FUNCTION QWGTF (X, OMEGA, P2, P3, P4, INTEGR)
!
!! QWGTF defines the weight function for QAWF.
!
!***SUBSIDIARY
!***PURPOSE  This function subprogram is used together with the
!            routine QAWF and defines the WEIGHT function.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (QWGTF-S, DQWGTF-D)
!***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
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
!***END PROLOGUE  QWGTF
!
  REAL QWGTF
  REAL OMEGA,OMX,P2,P3,P4,X
  INTEGER INTEGR
!***FIRST EXECUTABLE STATEMENT  QWGTF
  OMX = OMEGA*X
  go to(10,20),INTEGR
   10 QWGTF = COS(OMX)
  go to 30
   20 QWGTF = SIN(OMX)
   30 RETURN
end
