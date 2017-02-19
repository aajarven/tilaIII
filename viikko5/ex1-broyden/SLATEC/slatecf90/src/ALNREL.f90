function ALNREL (X)
!
!! ALNREL evaluates ln(1+X) accurate in the sense of relative error.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      SINGLE PRECISION (ALNREL-S, DLNREL-D, CLNREL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ALNREL(X) evaluates ln(1+X) accurately in the sense of relative
! error when X is very small.  This routine must be used to
! maintain relative error accuracy whenever X is small and
! accurately known.
!
! Series for ALNR       on the interval -3.75000D-01 to  3.75000D-01
!                                        with weighted error   1.93E-17
!                                         log weighted error  16.72
!                               significant figures required  16.44
!                                    decimal places required  17.40
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  ALNREL
  DIMENSION ALNRCS(23)
  LOGICAL FIRST
  SAVE ALNRCS, NLNREL, XMIN, FIRST
  DATA ALNRCS( 1) /   1.0378693562743770E0 /
  DATA ALNRCS( 2) /   -.13364301504908918E0 /
  DATA ALNRCS( 3) /    .019408249135520563E0 /
  DATA ALNRCS( 4) /   -.003010755112753577E0 /
  DATA ALNRCS( 5) /    .000486946147971548E0 /
  DATA ALNRCS( 6) /   -.000081054881893175E0 /
  DATA ALNRCS( 7) /    .000013778847799559E0 /
  DATA ALNRCS( 8) /   -.000002380221089435E0 /
  DATA ALNRCS( 9) /    .000000416404162138E0 /
  DATA ALNRCS(10) /   -.000000073595828378E0 /
  DATA ALNRCS(11) /    .000000013117611876E0 /
  DATA ALNRCS(12) /   -.000000002354670931E0 /
  DATA ALNRCS(13) /    .000000000425227732E0 /
  DATA ALNRCS(14) /   -.000000000077190894E0 /
  DATA ALNRCS(15) /    .000000000014075746E0 /
  DATA ALNRCS(16) /   -.000000000002576907E0 /
  DATA ALNRCS(17) /    .000000000000473424E0 /
  DATA ALNRCS(18) /   -.000000000000087249E0 /
  DATA ALNRCS(19) /    .000000000000016124E0 /
  DATA ALNRCS(20) /   -.000000000000002987E0 /
  DATA ALNRCS(21) /    .000000000000000554E0 /
  DATA ALNRCS(22) /   -.000000000000000103E0 /
  DATA ALNRCS(23) /    .000000000000000019E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ALNREL
  if (FIRST) THEN
     NLNREL = INITS (ALNRCS, 23, 0.1*R1MACH(3))
     XMIN = -1.0 + SQRT(R1MACH(4))
  end if
  FIRST = .FALSE.
!
  if (X  <=  (-1.0)) call XERMSG ('SLATEC', 'ALNREL', 'X IS LE -1', &
     2, 2)
  if (X  <  XMIN) call XERMSG ('SLATEC', 'ALNREL', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 1, 1)
!
  if (ABS(X) <= 0.375) ALNREL = X*(1. - &
    X*CSEVL (X/.375, ALNRCS, NLNREL))
  if (ABS(X) > 0.375) ALNREL = LOG (1.0+X)
!
  return
end
