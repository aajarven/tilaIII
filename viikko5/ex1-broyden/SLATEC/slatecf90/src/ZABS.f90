  DOUBLE PRECISION FUNCTION ZABS (ZR, ZI)
!
!! ZABS is subsidiary to ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZAIRY and ZBIRY.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (ZABS-A)
!***AUTHOR  Amos, D. E., (SNL)
!***DESCRIPTION
!
!     ZABS COMPUTES THE ABSOLUTE VALUE OR MAGNITUDE OF A DOUBLE
!     PRECISION COMPLEX VARIABLE CMPLX(ZR,ZI)
!
!***SEE ALSO  ZAIRY, ZBESH, ZBESI, ZBESJ, ZBESK, ZBESY, ZBIRY
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830501  DATE WRITTEN
!   910415  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  ZABS
  DOUBLE PRECISION ZR, ZI, U, V, Q, S
!***FIRST EXECUTABLE STATEMENT  ZABS
  U = ABS(ZR)
  V = ABS(ZI)
  S = U + V
!-----------------------------------------------------------------------
!     S*1.0D0 MAKES AN UNNORMALIZED UNDERFLOW ON CDC MACHINES INTO A
!     TRUE FLOATING ZERO
!-----------------------------------------------------------------------
  S = S*1.0D+0
  if (S == 0.0D+0) go to 20
  if (U > V) go to 10
  Q = U/V
  ZABS = V*SQRT(1.D+0+Q*Q)
  return
   10 Q = V/U
  ZABS = U*SQRT(1.D+0+Q*Q)
  return
   20 ZABS = 0.0D+0
  return
end
