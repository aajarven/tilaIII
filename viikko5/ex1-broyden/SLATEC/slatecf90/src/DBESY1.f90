  DOUBLE PRECISION FUNCTION DBESY1 (X)
!
!! DBESY1 computes the Bessel function of the second kind of order one.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10A1
!***TYPE      DOUBLE PRECISION (BESY1-S, DBESY1-D)
!***KEYWORDS  BESSEL FUNCTION, FNLIB, ORDER ONE, SECOND KIND,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! DBESY1(X) calculates the double precision Bessel function of the
! second kind of order for double precision argument X.
!
! Series for BY1        on the interval  0.          to  1.60000E+01
!                                        with weighted error   8.65E-33
!                                         log weighted error  32.06
!                               significant figures required  32.17
!                                    decimal places required  32.71
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, D9B1MP, DBESJ1, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DBESY1
  DOUBLE PRECISION X, BY1CS(20), AMPL, THETA, TWODPI, XMIN, XSML, &
    Y, D1MACH, DCSEVL, DBESJ1
  LOGICAL FIRST
  SAVE BY1CS, TWODPI, NTY1, XMIN, XSML, FIRST
  DATA BY1CS(  1) / +.320804710061190862932352018628015D-1    /
  DATA BY1CS(  2) / +.126270789743350044953431725999727D+1    /
  DATA BY1CS(  3) / +.649996189992317500097490637314144D-2    /
  DATA BY1CS(  4) / -.893616452886050411653144160009712D-1    /
  DATA BY1CS(  5) / +.132508812217570954512375510370043D-1    /
  DATA BY1CS(  6) / -.897905911964835237753039508298105D-3    /
  DATA BY1CS(  7) / +.364736148795830678242287368165349D-4    /
  DATA BY1CS(  8) / -.100137438166600055549075523845295D-5    /
  DATA BY1CS(  9) / +.199453965739017397031159372421243D-7    /
  DATA BY1CS( 10) / -.302306560180338167284799332520743D-9    /
  DATA BY1CS( 11) / +.360987815694781196116252914242474D-11   /
  DATA BY1CS( 12) / -.348748829728758242414552947409066D-13   /
  DATA BY1CS( 13) / +.278387897155917665813507698517333D-15   /
  DATA BY1CS( 14) / -.186787096861948768766825352533333D-17   /
  DATA BY1CS( 15) / +.106853153391168259757070336000000D-19   /
  DATA BY1CS( 16) / -.527472195668448228943872000000000D-22   /
  DATA BY1CS( 17) / +.227019940315566414370133333333333D-24   /
  DATA BY1CS( 18) / -.859539035394523108693333333333333D-27   /
  DATA BY1CS( 19) / +.288540437983379456000000000000000D-29   /
  DATA BY1CS( 20) / -.864754113893717333333333333333333D-32   /
  DATA TWODPI / 0.636619772367581343075535053490057D0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DBESY1
  if (FIRST) THEN
     NTY1 = INITDS (BY1CS, 20, 0.1*REAL(D1MACH(3)))
!
     XMIN = 1.571D0 * EXP (MAX(LOG(D1MACH(1)), -LOG(D1MACH(2))) + &
       0.01D0)
     XSML = SQRT(4.0D0*D1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X  <=  0.D0) call XERMSG ('SLATEC', 'DBESY1', &
     'X IS ZERO OR NEGATIVE', 1, 2)
  if (X > 4.0D0) go to 20
!
  if (X  <  XMIN) call XERMSG ('SLATEC', 'DBESY1', &
     'X SO SMALL Y1 OVERFLOWS', 3, 2)
  Y = 0.D0
  if (X > XSML) Y = X*X
  DBESY1 = TWODPI * LOG(0.5D0*X)*DBESJ1(X) + (0.5D0 + &
    DCSEVL (.125D0*Y-1.D0, BY1CS, NTY1))/X
  return
!
 20   call D9B1MP (X, AMPL, THETA)
  DBESY1 = AMPL * SIN(THETA)
  return
!
end