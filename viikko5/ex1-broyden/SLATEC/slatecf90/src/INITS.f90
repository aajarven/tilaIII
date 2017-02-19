function INITS (OS, NOS, ETA)
!
!! INITS determines the number of terms needed in an orthogonal ...
!            polynomial series so that it meets a specified accuracy.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C3A2
!***TYPE      SINGLE PRECISION (INITS-S, INITDS-D)
!***KEYWORDS  CHEBYSHEV, FNLIB, INITIALIZE, ORTHOGONAL POLYNOMIAL,
!             ORTHOGONAL SERIES, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!  Initialize the orthogonal series, represented by the array OS, so
!  that INITS is the number of terms needed to insure the error is no
!  larger than ETA.  Ordinarily, ETA will be chosen to be one-tenth
!  machine precision.
!
!             Input Arguments --
!   OS     single precision array of NOS coefficients in an orthogonal
!          series.
!   NOS    number of coefficients in OS.
!   ETA    single precision scalar containing requested accuracy of
!          series.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891115  Modified error message.  (WRB)
!   891115  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  INITS
  REAL OS(*)
!***FIRST EXECUTABLE STATEMENT  INITS
  if (NOS  <  1) call XERMSG ('SLATEC', 'INITS', &
     'Number of coefficients is less than 1', 2, 1)
!
  ERR = 0.
  DO 10 II = 1,NOS
    I = NOS + 1 - II
    ERR = ERR + ABS(OS(I))
    if (ERR > ETA) go to 20
   10 CONTINUE
!
   20 if (I  ==  NOS) call XERMSG ('SLATEC', 'INITS', &
     'Chebyshev series too short for specified accuracy', 1, 1)
  INITS = I
!
  return
end
