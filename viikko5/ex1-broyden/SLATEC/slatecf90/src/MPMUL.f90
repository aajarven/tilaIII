subroutine MPMUL (X, Y, Z)
!
!! MPMUL multiples two 'mp' numbers.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DQDOTA and DQDOTI
!***LIBRARY   SLATEC
!***TYPE      ALL (MPMUL-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Multiplies X and Y, returning result in Z, for 'mp' X, Y and Z.
!  The simple o(t**2) algorithm is used, with four guard digits and
!  R*-rounding. Advantage is taken of zero digits in X, but not in Y.
!  Asymptotically faster algorithms are known (see Knuth, VOL. 2),
!  but are difficult to implement in FORTRAN in an efficient and
!  machine-independent manner. In comments to other 'mp' routines,
!  M(t) is the time to perform t-digit 'mp' multiplication. Thus
!  M(t) = o(t**2) with the present version of MPMUL, but
!  M(t) = o(t.log(t).log(log(t))) is theoretically possible.
!
!  The arguments X(*), Y(*), and Z(*), and the variable R in COMMON are
!  all INTEGER arrays of size 30.  See the comments in the routine
!  MPBLAS for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPCHK, MPERR, MPMLP, MPNZR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPMUL
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Y(*), Z(*), RS, RE, XI, C, RI
!***FIRST EXECUTABLE STATEMENT  MPMUL
  call MPCHK (1, 4)
  I2 = T + 4
  I2P = I2 + 1
! FORM SIGN OF PRODUCT
  RS = X(1)*Y(1)
  if (RS /= 0) go to 10
! SET RESULT TO ZERO
  Z(1) = 0
  return
! FORM EXPONENT OF PRODUCT
   10 RE = X(2) + Y(2)
! CLEAR ACCUMULATOR
  DO 20 I = 1, I2
   20 R(I) = 0
! PERFORM MULTIPLICATION
  C = 8
  DO 40 I = 1, T
  XI = X(I+2)
! FOR SPEED, PUT THE NUMBER WITH MANY ZEROS FIRST
  if (XI == 0) go to 40
  call MPMLP (R(I+1), Y(3), XI, MIN (T, I2 - I))
  C = C - 1
  if (C > 0) go to 40
! CHECK FOR LEGAL BASE B DIGIT
  if ((XI < 0).OR.(XI >= B)) go to 90
! PROPAGATE CARRIES AT END AND EVERY EIGHTH TIME,
! FASTER THAN DOING IT EVERY TIME.
  DO 30 J = 1, I2
  J1 = I2P - J
  RI = R(J1) + C
  if (RI < 0) go to 70
  C = RI/B
   30 R(J1) = RI - B*C
  if (C /= 0) go to 90
  C = 8
   40 CONTINUE
  if (C == 8) go to 60
  if ((XI < 0).OR.(XI >= B)) go to 90
  C = 0
  DO 50 J = 1, I2
  J1 = I2P - J
  RI = R(J1) + C
  if (RI < 0) go to 70
  C = RI/B
   50 R(J1) = RI - B*C
  if (C /= 0) go to 90
! NORMALIZE AND ROUND RESULT
   60 call MPNZR (RS, RE, Z, 0)
  return
   70 WRITE (LUN, 80)
   80 FORMAT (' *** INTEGER OVERFLOW IN MPMUL, B TOO LARGE ***')
  go to 110
   90 WRITE (LUN, 100)
  100 FORMAT (' *** ILLEGAL BASE B DIGIT IN call TO MPMUL,', &
          ' POSSIBLE OVERWRITING PROBLEM ***')
  110 call MPERR
  Z(1) = 0
  return
end
