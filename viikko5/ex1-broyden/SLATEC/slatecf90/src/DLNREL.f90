  DOUBLE PRECISION FUNCTION DLNREL (X)
!
!! DLNREL evaluates ln(1+X) accurate in the sense of relative error.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      DOUBLE PRECISION (ALNREL-S, DLNREL-D, CLNREL-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DLNREL(X) calculates the double precision natural logarithm of
! (1.0+X) for double precision argument X.  This routine should
! be used when X is small and accurate to calculate the logarithm
! accurately (in the relative error sense) in the neighborhood
! of 1.0.
!
! Series for ALNR       on the interval -3.75000E-01 to  3.75000E-01
!                                        with weighted error   6.35E-32
!                                         log weighted error  31.20
!                               significant figures required  30.93
!                                    decimal places required  32.01
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DLNREL
  DOUBLE PRECISION ALNRCS(43), X, XMIN,  DCSEVL, D1MACH
  LOGICAL FIRST
  SAVE ALNRCS, NLNREL, XMIN, FIRST
  DATA ALNRCS(  1) / +.10378693562743769800686267719098D+1     /
  DATA ALNRCS(  2) / -.13364301504908918098766041553133D+0     /
  DATA ALNRCS(  3) / +.19408249135520563357926199374750D-1     /
  DATA ALNRCS(  4) / -.30107551127535777690376537776592D-2     /
  DATA ALNRCS(  5) / +.48694614797154850090456366509137D-3     /
  DATA ALNRCS(  6) / -.81054881893175356066809943008622D-4     /
  DATA ALNRCS(  7) / +.13778847799559524782938251496059D-4     /
  DATA ALNRCS(  8) / -.23802210894358970251369992914935D-5     /
  DATA ALNRCS(  9) / +.41640416213865183476391859901989D-6     /
  DATA ALNRCS( 10) / -.73595828378075994984266837031998D-7     /
  DATA ALNRCS( 11) / +.13117611876241674949152294345011D-7     /
  DATA ALNRCS( 12) / -.23546709317742425136696092330175D-8     /
  DATA ALNRCS( 13) / +.42522773276034997775638052962567D-9     /
  DATA ALNRCS( 14) / -.77190894134840796826108107493300D-10    /
  DATA ALNRCS( 15) / +.14075746481359069909215356472191D-10    /
  DATA ALNRCS( 16) / -.25769072058024680627537078627584D-11    /
  DATA ALNRCS( 17) / +.47342406666294421849154395005938D-12    /
  DATA ALNRCS( 18) / -.87249012674742641745301263292675D-13    /
  DATA ALNRCS( 19) / +.16124614902740551465739833119115D-13    /
  DATA ALNRCS( 20) / -.29875652015665773006710792416815D-14    /
  DATA ALNRCS( 21) / +.55480701209082887983041321697279D-15    /
  DATA ALNRCS( 22) / -.10324619158271569595141333961932D-15    /
  DATA ALNRCS( 23) / +.19250239203049851177878503244868D-16    /
  DATA ALNRCS( 24) / -.35955073465265150011189707844266D-17    /
  DATA ALNRCS( 25) / +.67264542537876857892194574226773D-18    /
  DATA ALNRCS( 26) / -.12602624168735219252082425637546D-18    /
  DATA ALNRCS( 27) / +.23644884408606210044916158955519D-19    /
  DATA ALNRCS( 28) / -.44419377050807936898878389179733D-20    /
  DATA ALNRCS( 29) / +.83546594464034259016241293994666D-21    /
  DATA ALNRCS( 30) / -.15731559416479562574899253521066D-21    /
  DATA ALNRCS( 31) / +.29653128740247422686154369706666D-22    /
  DATA ALNRCS( 32) / -.55949583481815947292156013226666D-23    /
  DATA ALNRCS( 33) / +.10566354268835681048187284138666D-23    /
  DATA ALNRCS( 34) / -.19972483680670204548314999466666D-24    /
  DATA ALNRCS( 35) / +.37782977818839361421049855999999D-25    /
  DATA ALNRCS( 36) / -.71531586889081740345038165333333D-26    /
  DATA ALNRCS( 37) / +.13552488463674213646502024533333D-26    /
  DATA ALNRCS( 38) / -.25694673048487567430079829333333D-27    /
  DATA ALNRCS( 39) / +.48747756066216949076459519999999D-28    /
  DATA ALNRCS( 40) / -.92542112530849715321132373333333D-29    /
  DATA ALNRCS( 41) / +.17578597841760239233269760000000D-29    /
  DATA ALNRCS( 42) / -.33410026677731010351377066666666D-30    /
  DATA ALNRCS( 43) / +.63533936180236187354180266666666D-31    /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DLNREL
  if (FIRST) THEN
     NLNREL = INITDS (ALNRCS, 43, 0.1*REAL(D1MACH(3)))
     XMIN = -1.0D0 + SQRT(D1MACH(4))
  end if
  FIRST = .FALSE.
!
  if (X  <=  (-1.D0)) call XERMSG ('SLATEC', 'DLNREL', 'X IS LE -1' &
     , 2, 2)
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DLNREL', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR -1', 1, 1)
!
  if (ABS(X) <= 0.375D0) DLNREL = X*(1.D0 - &
    X*DCSEVL (X/.375D0, ALNRCS, NLNREL))
!
  if (ABS(X) > 0.375D0) DLNREL = LOG (1.0D0+X)
!
  return
end