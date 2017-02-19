function R9PAK (Y, N)
!
!! R9PAK packs a base 2 exponent into a floating point number.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  A6B
!***TYPE      SINGLE PRECISION (R9PAK-S, D9PAK-D)
!***KEYWORDS  FNLIB, PACK
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Pack a base 2 exponent into floating point number Y.  This
! routine is almost the inverse of R9UPAK.  It is not exactly
! the inverse, because ABS(X) need not be between 0.5 and
! 1.0.  If both R9PAK and 2.0**N were known to be in range, we
! could compute
!       R9PAK = Y * 2.0**N .
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  I1MACH, R1MACH, R9UPAK, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901009  Routine used I1MACH(7) where it should use I1MACH(10),
!           Corrected (RWC)
!***END PROLOGUE  R9PAK
  LOGICAL FIRST
  SAVE NMIN, NMAX, A1N210, FIRST
  DATA A1N210 / 3.321928094887362E0/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9PAK
  if (FIRST) THEN
     A1N2B = 1.0
     if (I1MACH(10) /= 2) A1N2B = R1MACH(5)*A1N210
     NMIN = A1N2B*I1MACH(12)
     NMAX = A1N2B*I1MACH(13)
  end if
  FIRST = .FALSE.
!
  call R9UPAK(Y,R9PAK,NY)
!
  NSUM = N + NY
  if (NSUM < NMIN) go to 40
  if (NSUM  >  NMAX) call XERMSG ('SLATEC', 'R9PAK', &
     'PACKED NUMBER OVERFLOWS', 2, 2)
!
  if (NSUM == 0) RETURN
  if (NSUM > 0) go to 30
!
 20   R9PAK = 0.5*R9PAK
  NSUM = NSUM + 1
  if ( NSUM /= 0) go to 20
  return
!
30    R9PAK = 2.0*R9PAK
  NSUM = NSUM - 1
  if ( NSUM /= 0) go to 30
  return
!
40    call XERMSG ('SLATEC', 'R9PAK', 'PACKED NUMBER UNDERFLOWS', 1, 1)
  R9PAK = 0.0
  return
!
end
