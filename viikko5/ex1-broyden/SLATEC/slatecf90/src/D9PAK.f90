DOUBLE PRECISION FUNCTION D9PAK (Y, N)
!
!! D9PAK packs a base 2 exponent into a floating point number.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  A6B
!***TYPE      DOUBLE PRECISION (R9PAK-S, D9PAK-D)
!***KEYWORDS  FNLIB, PACK
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Pack a base 2 exponent into floating point number X.  This routine is
! almost the inverse of D9UPAK.  It is not exactly the inverse, because
! ABS(X) need not be between 0.5 and 1.0.  If both D9PAK and 2.d0**N
! were known to be in range we could compute
!               D9PAK = X *2.0d0**N
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9UPAK, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891009  Corrected error when XERROR called.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901009  Routine used I1MACH(7) where it should use I1MACH(10),
!           Corrected (RWC)
!***END PROLOGUE  D9PAK
  DOUBLE PRECISION Y, A1N2B,A1N210,D1MACH
  LOGICAL FIRST
  SAVE NMIN, NMAX, A1N210, FIRST
  DATA A1N210 / 3.321928094887362347870319429489D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9PAK
  if (FIRST) THEN
     A1N2B = 1.0D0
     if ( I1MACH(10) /= 2) A1N2B=D1MACH(5)*A1N210
     NMIN = A1N2B*I1MACH(15)
     NMAX = A1N2B*I1MACH(16)
  end if
  FIRST = .FALSE.
!
  call D9UPAK(Y,D9PAK,NY)
!
  NSUM=N+NY
  if ( NSUM < NMIN)go to 40
  if (NSUM  >  NMAX) call XERMSG ('SLATEC', 'D9PAK', &
     'PACKED NUMBER OVERFLOWS', 1, 2)
!
  if (NSUM == 0) RETURN
  if ( NSUM > 0) go to 30
!
 20   D9PAK = 0.5D0*D9PAK
  NSUM=NSUM+1
  if ( NSUM /= 0) go to 20
  return
!
 30   D9PAK = 2.0D0*D9PAK
  NSUM=NSUM - 1
  if (NSUM /= 0) go to 30
  return
!
 40   call XERMSG ('SLATEC', 'D9PAK', 'PACKED NUMBER UNDERFLOWS', 1, 1)
  D9PAK = 0.0D0
  return
!
end
