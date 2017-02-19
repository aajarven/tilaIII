DOUBLE PRECISION FUNCTION DASINH (X)
!
!! DASINH computes the arc hyperbolic sine.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4C
!***TYPE      DOUBLE PRECISION (ASINH-S, DASINH-D, CASINH-C)
!***KEYWORDS  ARC HYPERBOLIC SINE, ASINH, ELEMENTARY FUNCTIONS, FNLIB,
!             INVERSE HYPERBOLIC SINE
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DASINH(X) calculates the double precision arc hyperbolic
! sine for double precision argument X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DASINH
  DOUBLE PRECISION X, ASNHCS(39), ALN2, SQEPS, XMAX, Y, &
    DCSEVL, D1MACH
  LOGICAL FIRST
  SAVE ASNHCS, ALN2, NTERMS, XMAX, SQEPS, FIRST
  DATA ASNHCS(  1) / -.12820039911738186343372127359268D+0     /
  DATA ASNHCS(  2) / -.58811761189951767565211757138362D-1     /
  DATA ASNHCS(  3) / +.47274654322124815640725249756029D-2     /
  DATA ASNHCS(  4) / -.49383631626536172101360174790273D-3     /
  DATA ASNHCS(  5) / +.58506207058557412287494835259321D-4     /
  DATA ASNHCS(  6) / -.74669983289313681354755069217188D-5     /
  DATA ASNHCS(  7) / +.10011693583558199265966192015812D-5     /
  DATA ASNHCS(  8) / -.13903543858708333608616472258886D-6     /
  DATA ASNHCS(  9) / +.19823169483172793547317360237148D-7     /
  DATA ASNHCS( 10) / -.28847468417848843612747272800317D-8     /
  DATA ASNHCS( 11) / +.42672965467159937953457514995907D-9     /
  DATA ASNHCS( 12) / -.63976084654366357868752632309681D-10    /
  DATA ASNHCS( 13) / +.96991686089064704147878293131179D-11    /
  DATA ASNHCS( 14) / -.14844276972043770830246658365696D-11    /
  DATA ASNHCS( 15) / +.22903737939027447988040184378983D-12    /
  DATA ASNHCS( 16) / -.35588395132732645159978942651310D-13    /
  DATA ASNHCS( 17) / +.55639694080056789953374539088554D-14    /
  DATA ASNHCS( 18) / -.87462509599624678045666593520162D-15    /
  DATA ASNHCS( 19) / +.13815248844526692155868802298129D-15    /
  DATA ASNHCS( 20) / -.21916688282900363984955142264149D-16    /
  DATA ASNHCS( 21) / +.34904658524827565638313923706880D-17    /
  DATA ASNHCS( 22) / -.55785788400895742439630157032106D-18    /
  DATA ASNHCS( 23) / +.89445146617134012551050882798933D-19    /
  DATA ASNHCS( 24) / -.14383426346571317305551845239466D-19    /
  DATA ASNHCS( 25) / +.23191811872169963036326144682666D-20    /
  DATA ASNHCS( 26) / -.37487007953314343674570604543999D-21    /
  DATA ASNHCS( 27) / +.60732109822064279404549242880000D-22    /
  DATA ASNHCS( 28) / -.98599402764633583177370173440000D-23    /
  DATA ASNHCS( 29) / +.16039217452788496315232638293333D-23    /
  DATA ASNHCS( 30) / -.26138847350287686596716134399999D-24    /
  DATA ASNHCS( 31) / +.42670849606857390833358165333333D-25    /
  DATA ASNHCS( 32) / -.69770217039185243299730773333333D-26    /
  DATA ASNHCS( 33) / +.11425088336806858659812693333333D-26    /
  DATA ASNHCS( 34) / -.18735292078860968933021013333333D-27    /
  DATA ASNHCS( 35) / +.30763584414464922794065920000000D-28    /
  DATA ASNHCS( 36) / -.50577364031639824787046399999999D-29    /
  DATA ASNHCS( 37) / +.83250754712689142224213333333333D-30    /
  DATA ASNHCS( 38) / -.13718457282501044163925333333333D-30    /
  DATA ASNHCS( 39) / +.22629868426552784104106666666666D-31    /
  DATA ALN2 / 0.69314718055994530941723212145818D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DASINH
  if (FIRST) THEN
     NTERMS = INITDS (ASNHCS, 39, 0.1*REAL(D1MACH(3)) )
     SQEPS = SQRT(D1MACH(3))
     XMAX = 1.0D0/SQEPS
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y > 1.0D0) go to 20
!
  DASINH = X
  if (Y > SQEPS) DASINH = X*(1.0D0 + DCSEVL (2.D0*X*X-1.D0, &
    ASNHCS, NTERMS) )
  return
 20   if (Y < XMAX) DASINH = LOG (Y+SQRT(Y*Y+1.D0))
  if (Y >= XMAX) DASINH = ALN2 + LOG(Y)
  DASINH = SIGN (DASINH, X)
  return
!
end
