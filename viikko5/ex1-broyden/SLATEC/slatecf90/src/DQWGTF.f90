  DOUBLE PRECISION FUNCTION DQWGTF (X, OMEGA, P2, P3, P4, INTEGR)
!
!! DQWGTF defines the weight function for DQAWF.
!
!***SUBSIDIARY
!***PURPOSE  This function subprogram is used together with the
!            routine DQAWF and defines the WEIGHT function.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QWGTF-S, DQWGTF-D)
!***KEYWORDS  COS OR SIN IN WEIGHT FUNCTION
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***SEE ALSO  DQK15W
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQWGTF
!
  DOUBLE PRECISION OMEGA,OMX,P2,P3,P4,X
  INTEGER INTEGR
!***FIRST EXECUTABLE STATEMENT  DQWGTF
  OMX = OMEGA*X
  go to(10,20),INTEGR
   10 DQWGTF = COS(OMX)
  go to 30
   20 DQWGTF = SIN(OMX)
   30 RETURN
end
