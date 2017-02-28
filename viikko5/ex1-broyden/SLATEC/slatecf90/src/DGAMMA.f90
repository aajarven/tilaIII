  DOUBLE PRECISION FUNCTION DGAMMA (X)
!
!! DGAMMA computes the complete Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      DOUBLE PRECISION (GAMMA-S, DGAMMA-D, CGAMMA-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DGAMMA(X) calculates the double precision complete Gamma function
! for double precision argument X.
!
! Series for GAM        on the interval  0.          to  1.00000E+00
!                                        with weighted error   5.79E-32
!                                         log weighted error  31.24
!                               significant figures required  30.00
!                                    decimal places required  32.05
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9LGMC, DCSEVL, DGAMLM, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable name.  (RWC, WRB)
!***END PROLOGUE  DGAMMA
  DOUBLE PRECISION X, GAMCS(42), DXREL, PI, SINPIY, SQ2PIL, XMAX, &
    XMIN, Y, D9LGMC, DCSEVL, D1MACH
  LOGICAL FIRST
!
  SAVE GAMCS, PI, SQ2PIL, NGAM, XMIN, XMAX, DXREL, FIRST
  DATA GAMCS(  1) / +.8571195590989331421920062399942D-2      /
  DATA GAMCS(  2) / +.4415381324841006757191315771652D-2      /
  DATA GAMCS(  3) / +.5685043681599363378632664588789D-1      /
  DATA GAMCS(  4) / -.4219835396418560501012500186624D-2      /
  DATA GAMCS(  5) / +.1326808181212460220584006796352D-2      /
  DATA GAMCS(  6) / -.1893024529798880432523947023886D-3      /
  DATA GAMCS(  7) / +.3606925327441245256578082217225D-4      /
  DATA GAMCS(  8) / -.6056761904460864218485548290365D-5      /
  DATA GAMCS(  9) / +.1055829546302283344731823509093D-5      /
  DATA GAMCS( 10) / -.1811967365542384048291855891166D-6      /
  DATA GAMCS( 11) / +.3117724964715322277790254593169D-7      /
  DATA GAMCS( 12) / -.5354219639019687140874081024347D-8      /
  DATA GAMCS( 13) / +.9193275519859588946887786825940D-9      /
  DATA GAMCS( 14) / -.1577941280288339761767423273953D-9      /
  DATA GAMCS( 15) / +.2707980622934954543266540433089D-10     /
  DATA GAMCS( 16) / -.4646818653825730144081661058933D-11     /
  DATA GAMCS( 17) / +.7973350192007419656460767175359D-12     /
  DATA GAMCS( 18) / -.1368078209830916025799499172309D-12     /
  DATA GAMCS( 19) / +.2347319486563800657233471771688D-13     /
  DATA GAMCS( 20) / -.4027432614949066932766570534699D-14     /
  DATA GAMCS( 21) / +.6910051747372100912138336975257D-15     /
  DATA GAMCS( 22) / -.1185584500221992907052387126192D-15     /
  DATA GAMCS( 23) / +.2034148542496373955201026051932D-16     /
  DATA GAMCS( 24) / -.3490054341717405849274012949108D-17     /
  DATA GAMCS( 25) / +.5987993856485305567135051066026D-18     /
  DATA GAMCS( 26) / -.1027378057872228074490069778431D-18     /
  DATA GAMCS( 27) / +.1762702816060529824942759660748D-19     /
  DATA GAMCS( 28) / -.3024320653735306260958772112042D-20     /
  DATA GAMCS( 29) / +.5188914660218397839717833550506D-21     /
  DATA GAMCS( 30) / -.8902770842456576692449251601066D-22     /
  DATA GAMCS( 31) / +.1527474068493342602274596891306D-22     /
  DATA GAMCS( 32) / -.2620731256187362900257328332799D-23     /
  DATA GAMCS( 33) / +.4496464047830538670331046570666D-24     /
  DATA GAMCS( 34) / -.7714712731336877911703901525333D-25     /
  DATA GAMCS( 35) / +.1323635453126044036486572714666D-25     /
  DATA GAMCS( 36) / -.2270999412942928816702313813333D-26     /
  DATA GAMCS( 37) / +.3896418998003991449320816639999D-27     /
  DATA GAMCS( 38) / -.6685198115125953327792127999999D-28     /
  DATA GAMCS( 39) / +.1146998663140024384347613866666D-28     /
  DATA GAMCS( 40) / -.1967938586345134677295103999999D-29     /
  DATA GAMCS( 41) / +.3376448816585338090334890666666D-30     /
  DATA GAMCS( 42) / -.5793070335782135784625493333333D-31     /
  DATA PI / 3.14159265358979323846264338327950D0 /
  DATA SQ2PIL / 0.91893853320467274178032973640562D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DGAMMA
  if (FIRST) THEN
     NGAM = INITDS (GAMCS, 42, 0.1*REAL(D1MACH(3)) )
!
     call DGAMLM (XMIN, XMAX)
     DXREL = SQRT(D1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 10.D0) go to 50
!
! COMPUTE GAMMA(X) FOR -XBND  <=  X  <=  XBND.  REDUCE INTERVAL AND FIND
! GAMMA(1+Y) FOR 0.0  <=  Y  <  1.0 FIRST OF ALL.
!
  N = X
  if (X < 0.D0) N = N - 1
  Y = X - N
  N = N - 1
  DGAMMA = 0.9375D0 + DCSEVL (2.D0*Y-1.D0, GAMCS, NGAM)
  if (N == 0) RETURN
!
  if (N > 0) go to 30
!
! COMPUTE GAMMA(X) FOR X  <  1.0
!
  N = -N
  if (X  ==  0.D0) call XERMSG ('SLATEC', 'DGAMMA', 'X IS 0', 4, 2)
  if (X  <  0.0 .AND. X+N-2  ==  0.D0) call XERMSG ('SLATEC', &
     'DGAMMA', 'X IS A NEGATIVE INTEGER', 4, 2)
  if (X  <  (-0.5D0) .AND. ABS((X-AINT(X-0.5D0))/X)  <  DXREL) &
     call XERMSG ('SLATEC', 'DGAMMA', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
     1, 1)
!
  DO 20 I=1,N
    DGAMMA = DGAMMA/(X+I-1 )
 20   CONTINUE
  return
!
! GAMMA(X) FOR X  >=  2.0 AND X  <=  10.0
!
 30   DO 40 I=1,N
    DGAMMA = (Y+I) * DGAMMA
 40   CONTINUE
  return
!
! GAMMA(X) FOR ABS(X)  >  10.0.  RECALL Y = ABS(X).
!
 50   if (X  >  XMAX) call XERMSG ('SLATEC', 'DGAMMA', &
     'X SO BIG GAMMA OVERFLOWS', 3, 2)
!
  DGAMMA = 0.D0
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DGAMMA', &
     'X SO SMALL GAMMA UNDERFLOWS', 2, 1)
  if (X < XMIN) RETURN
!
  DGAMMA = EXP ((Y-0.5D0)*LOG(Y) - Y + SQ2PIL + D9LGMC(Y) )
  if (X > 0.D0) RETURN
!
  if (ABS((X-AINT(X-0.5D0))/X)  <  DXREL) call XERMSG ('SLATEC', &
     'DGAMMA', &
     'ANSWER LT HALF PRECISION, X TOO NEAR NEGATIVE INTEGER', 1, 1)
!
  SINPIY = SIN (PI*Y)
  if (SINPIY  ==  0.D0) call XERMSG ('SLATEC', 'DGAMMA', &
     'X IS A NEGATIVE INTEGER', 4, 2)
!
  DGAMMA = -PI/(Y*SINPIY*DGAMMA)
!
  return
end