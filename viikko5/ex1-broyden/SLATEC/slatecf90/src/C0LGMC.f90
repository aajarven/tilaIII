FUNCTION C0LGMC (Z)
!
!! C0LGMC evaluates (Z+0.5)*LOG((Z+1.)/Z) - 1.0 with relative accuracy.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (C0LGMC-C)
!***KEYWORDS  FNLIB, GAMMA FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  (Z+0.5)*LOG((Z+1.0)/Z) - 1.0  with relative error accuracy
! Let Q = 1.0/Z so that
!     (Z+0.5)*LOG(1+1/Z) - 1 = (Z+0.5)*(LOG(1+Q) - Q + Q*Q/2) - Q*Q/4
!        = (Z+0.5)*Q**3*C9LN2R(Q) - Q**2/4,
! where  C9LN2R  is (LOG(1+Q) - Q + 0.5*Q**2) / Q**3.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  C9LN2R, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  C0LGMC
  COMPLEX C0LGMC
  COMPLEX Z, Q, C9LN2R
  SAVE RBIG
  DATA RBIG / 0.0 /
!***FIRST EXECUTABLE STATEMENT  C0LGMC
  if (RBIG == 0.0) RBIG = 1.0/R1MACH(3)

  CABSZ = ABS(Z)

  if ( CABSZ > RBIG ) then
    C0LGMC = -(Z+0.5)*LOG(Z) - Z
    RETURN
  end if

  Q = 1.0/Z
  if (CABSZ <= 1.23) C0LGMC = (Z+0.5)*LOG(1.0+Q) - 1.0
  if (CABSZ > 1.23) C0LGMC = ((1.+.5*Q)*C9LN2R(Q) - .25) * Q**2

  return
end
