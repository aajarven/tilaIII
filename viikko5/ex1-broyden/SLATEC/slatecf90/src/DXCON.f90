subroutine DXCON (X, IX, IERROR)
!
!! DXCON provides double-precision floating-point arithmetic ...
!            with an extended exponent range.
!
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      DOUBLE PRECISION (XCON-S, DXCON-D)
!***KEYWORDS  EXTENDED-RANGE DOUBLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     DOUBLE PRECISION X
!     INTEGER IX
!
!                  CONVERTS (X,IX) = X*RADIX**IX
!                  TO DECIMAL FORM IN PREPARATION FOR
!                  PRINTING, SO THAT (X,IX) = X*10**IX
!                  WHERE 1/10  <=  ABS(X)  <  1
!                  IS RETURNED, EXCEPT THAT IF
!                  (ABS(X),IX) IS BETWEEN RADIX**(-2L)
!                  AND RADIX**(2L) THEN THE REDUCED
!                  FORM WITH IX = 0 IS RETURNED.
!
!***SEE ALSO  DXSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  DXADJ, DXC210, DXRED
!***COMMON BLOCKS    DXBLK2
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXCON
  DOUBLE PRECISION X
  INTEGER IX
!
!   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! ARE
!    (1) 4  <=  L  <=  2**NBITS - 1 - KMAX
!
!    (2) KMAX  <=  ((2**NBITS)-2)/LOG10R - L
!
! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE DXSET.
!
  DOUBLE PRECISION RADIX, RADIXL, RAD2L, DLG10R
  INTEGER L, L2, KMAX
  COMMON /DXBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
  SAVE /DXBLK2/, ISPACE
!
  DOUBLE PRECISION A, B, Z
!
  DATA ISPACE /1/
!   THE PARAMETER ISPACE IS THE INCREMENT USED IN FORM-
! ING THE AUXILIARY INDEX OF THE DECIMAL EXTENDED-RANGE
! FORM. THE RETURNED VALUE OF IX WILL BE AN INTEGER MULT-
! IPLE OF ISPACE. ISPACE MUST SATISFY 1  <=  ISPACE  <=
! L/2. if A VALUE GREATER THAN 1 IS TAKEN, THE RETURNED
! VALUE OF X WILL SATISFY 10**(-ISPACE)  <=  ABS(X)  <=  1
! WHEN (ABS(X),IX)  <  RADIX**(-2L), AND 1/10  <=  ABS(X)
!  <  10**(ISPACE-1) WHEN (ABS(X),IX)  >  RADIX**(2L).
!
!***FIRST EXECUTABLE STATEMENT  DXCON
  IERROR=0
  call DXRED(X, IX,IERROR)
  if (IERROR /= 0) RETURN
  if (IX == 0) go to 150
  call DXADJ(X, IX,IERROR)
  if (IERROR /= 0) RETURN
!
! CASE 1 IS WHEN (X,IX) IS LESS THAN RADIX**(-2L) IN MAGNITUDE,
! CASE 2 IS WHEN (X,IX) IS GREATER THAN RADIX**(2L) IN MAGNITUDE.
  ITEMP = 1
  ICASE = (3+SIGN(ITEMP,IX))/2
  go to (10, 20), ICASE
   10 if (ABS(X) < 1.0D0) go to 30
  X = X/RADIXL
  IX = IX + L
  go to 30
   20 if (ABS(X) >= 1.0D0) go to 30
  X = X*RADIXL
  IX = IX - L
   30 CONTINUE
!
! AT THIS POINT, RADIX**(-L)  <=  ABS(X)  <  1.0D0     IN CASE 1,
!                      1.0D0  <=  ABS(X)  <  RADIX**L  IN CASE 2.
  I = LOG10(ABS(X))/DLG10R
  A = RADIX**I
  go to (40, 60), ICASE
   40 if (A <= RADIX*ABS(X)) go to 50
  I = I - 1
  A = A/RADIX
  go to 40
   50 if (ABS(X) < A) go to 80
  I = I + 1
  A = A*RADIX
  go to 50
   60 if (A <= ABS(X)) go to 70
  I = I - 1
  A = A/RADIX
  go to 60
   70 if (ABS(X) < RADIX*A) go to 80
  I = I + 1
  A = A*RADIX
  go to 70
   80 CONTINUE
!
! AT THIS POINT I IS SUCH THAT
! RADIX**(I-1)  <=  ABS(X)  <  RADIX**I      IN CASE 1,
!     RADIX**I  <=  ABS(X)  <  RADIX**(I+1)  IN CASE 2.
  ITEMP = ISPACE/DLG10R
  A = RADIX**ITEMP
  B = 10.0D0**ISPACE
   90 if (A <= B) go to 100
  ITEMP = ITEMP - 1
  A = A/RADIX
  go to 90
  100 if (B < A*RADIX) go to 110
  ITEMP = ITEMP + 1
  A = A*RADIX
  go to 100
  110 CONTINUE
!
! AT THIS POINT ITEMP IS SUCH THAT
! RADIX**ITEMP  <=  10**ISPACE  <  RADIX**(ITEMP+1).
  if (ITEMP > 0) go to 120
! ITEMP = 0 IF, AND ONLY IF, ISPACE = 1 AND RADIX = 16.0D0
  X = X*RADIX**(-I)
  IX = IX + I
  call DXC210(IX, Z, J,IERROR)
  if (IERROR /= 0) RETURN
  X = X*Z
  IX = J
  go to (130, 140), ICASE
  120 CONTINUE
  I1 = I/ITEMP
  X = X*RADIX**(-I1*ITEMP)
  IX = IX + I1*ITEMP
!
! AT THIS POINT,
! RADIX**(-ITEMP)  <=  ABS(X)  <  1.0D0        IN CASE 1,
!           1.0D0  <=  ABS(X)  <  RADIX**ITEMP IN CASE 2.
  call DXC210(IX, Z, J,IERROR)
  if (IERROR /= 0) RETURN
  J1 = J/ISPACE
  J2 = J - J1*ISPACE
  X = X*Z*10.0D0**J2
  IX = J1*ISPACE
!
! AT THIS POINT,
!  10.0D0**(-2*ISPACE)  <=  ABS(X)  <  1.0D0                IN CASE 1,
!           10.0D0**-1  <=  ABS(X)  <  10.0D0**(2*ISPACE-1) IN CASE 2.
  go to (130, 140), ICASE
  130 if (B*ABS(X) >= 1.0D0) go to 150
  X = X*B
  IX = IX - ISPACE
  go to 130
  140 if (10.0D0*ABS(X) < B) go to 150
  X = X/B
  IX = IX + ISPACE
  go to 140
  150 RETURN
end
