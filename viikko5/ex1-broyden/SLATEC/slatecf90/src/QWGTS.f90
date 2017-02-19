FUNCTION QWGTS (X, A, B, ALFA, BETA, INTEGR)
!
!! QWGTS defines the weight function for QAWS.
!
!***SUBSIDIARY
!***PURPOSE  This function subprogram is used together with the
!            routine QAWS and defines the WEIGHT function.
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (QWGTS-S, DQWGTS-D)
!***KEYWORDS  ALGEBRAICO-LOGARITHMIC, END POINT SINGULARITIES,
!             WEIGHT FUNCTION
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
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  QWGTS
!
  REAL QWGTS
  REAL A,ALFA,B,BETA,BMX,X,XMA
  INTEGER INTEGR
!***FIRST EXECUTABLE STATEMENT  QWGTS
  XMA = X-A
  BMX = B-X
  QWGTS = XMA**ALFA*BMX**BETA
  go to (40,10,20,30),INTEGR
   10 QWGTS = QWGTS*LOG(XMA)
  go to 40
   20 QWGTS = QWGTS*LOG(BMX)
  go to 40
   30 QWGTS = QWGTS*LOG(XMA)*LOG(BMX)
   40 RETURN
end
