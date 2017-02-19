function ERF (X)
!
!! ERF computes the error function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C8A, L5A1E
!***TYPE      SINGLE PRECISION (ERF-S, DERF-D)
!***KEYWORDS  ERF, ERROR FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ERF(X) calculates the single precision error function for
! single precision argument X.
!
! Series for ERF        on the interval  0.          to  1.00000D+00
!                                        with weighted error   7.10E-18
!                                         log weighted error  17.15
!                               significant figures required  16.31
!                                    decimal places required  17.71
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, ERFC, INITS, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900727  Added EXTERNAL statement.  (WRB)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  ERF
  DIMENSION ERFCS(13)
  LOGICAL FIRST
  EXTERNAL ERFC
  SAVE ERFCS, SQRTPI, NTERF, XBIG, SQEPS, FIRST
  DATA ERFCS( 1) /   -.049046121234691808E0 /
  DATA ERFCS( 2) /   -.14226120510371364E0 /
  DATA ERFCS( 3) /    .010035582187599796E0 /
  DATA ERFCS( 4) /   -.000576876469976748E0 /
  DATA ERFCS( 5) /    .000027419931252196E0 /
  DATA ERFCS( 6) /   -.000001104317550734E0 /
  DATA ERFCS( 7) /    .000000038488755420E0 /
  DATA ERFCS( 8) /   -.000000001180858253E0 /
  DATA ERFCS( 9) /    .000000000032334215E0 /
  DATA ERFCS(10) /   -.000000000000799101E0 /
  DATA ERFCS(11) /    .000000000000017990E0 /
  DATA ERFCS(12) /   -.000000000000000371E0 /
  DATA ERFCS(13) /    .000000000000000007E0 /
  DATA SQRTPI /1.7724538509055160E0/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ERF
  if (FIRST) THEN
     NTERF = INITS (ERFCS, 13, 0.1*R1MACH(3))
     XBIG = SQRT(-LOG(SQRTPI*R1MACH(3)))
     SQEPS = SQRT(2.0*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 1.) go to 20
!
! ERF(X) = 1. - ERFC(X) FOR -1.  <=  X  <=  1.
!
  if (Y <= SQEPS) ERF = 2.0*X/SQRTPI
  if (Y > SQEPS) ERF = X*(1.0 + CSEVL(2.*X**2-1., ERFCS, NTERF))
  return
!
! ERF(X) = 1. - ERFC(X) FOR  ABS(X)  >  1.
!
 20   if (Y <= XBIG) ERF = SIGN (1.0-ERFC(Y), X)
  if (Y > XBIG) ERF = SIGN (1.0, X)
!
  return
end
