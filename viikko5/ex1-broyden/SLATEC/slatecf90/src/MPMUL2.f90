subroutine MPMUL2 (X, IY, Z, TRUNC)
!
!! MPMUL2 multiplies an 'mp' number by a single precision integer.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DQDOTA and DQDOTI
!***LIBRARY   SLATEC
!***TYPE      ALL (MPMUL2-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Multiplies 'mp' X by single-precision integer IY giving 'mp' Z.
!  Multiplication by 1 may be used to normalize a number even if some
!  digits are greater than B-1. Result is rounded if TRUNC == 0,
!  otherwise truncated.
!
!  The arguments X(*) and Z(*), and the variable R in COMMON are all
!  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
!  for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPOVFL, MPSTR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPMUL2
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Z(*), TRUNC, RE, RS
  INTEGER C, C1, C2, RI, T1, T3, T4
!***FIRST EXECUTABLE STATEMENT  MPMUL2
  RS = X(1)
  if (RS == 0) go to 10
  J = IY
  if (J) 20, 10, 50
! RESULT ZERO
   10 Z(1) = 0
  return
   20 J = -J
  RS = -RS
! CHECK FOR MULTIPLICATION BY B
  if (J /= B) go to 50
  if (X(2) < M) go to 40
  call MPCHK (1, 4)
  WRITE (LUN, 30)
   30 FORMAT (' *** OVERFLOW OCCURRED IN MPMUL2 ***')
  call MPOVFL (Z)
  return
   40 call MPSTR (X, Z)
  Z(1) = RS
  Z(2) = X(2) + 1
  return
! SET EXPONENT TO EXPONENT(X) + 4
   50 RE = X(2) + 4
! FORM PRODUCT IN ACCUMULATOR
  C = 0
  T1 = T + 1
  T3 = T + 3
  T4 = T + 4
! if J*B NOT REPRESENTABLE AS AN INTEGER WE HAVE TO SIMULATE
! DOUBLE-PRECISION MULTIPLICATION.
  if (J >= MAX(8*B, 32767/B)) go to 110
  DO 60 IJ = 1, T
  I = T1 - IJ
  RI = J*X(I+2) + C
  C = RI/B
   60 R(I+4) = RI - B*C
! CHECK FOR INTEGER OVERFLOW
  if (RI < 0) go to 130
! HAVE TO TREAT FIRST FOUR WORDS OF R SEPARATELY
  DO 70 IJ = 1, 4
  I = 5 - IJ
  RI = C
  C = RI/B
   70 R(I) = RI - B*C
  if (C == 0) go to 100
! HAVE TO SHIFT RIGHT HERE AS CARRY OFF END
   80 DO 90 IJ = 1, T3
  I = T4 - IJ
   90 R(I+1) = R(I)
  RI = C
  C = RI/B
  R(1) = RI - B*C
  RE = RE + 1
  if (C) 130, 100, 80
! NORMALIZE AND ROUND OR TRUNCATE RESULT
  100 call MPNZR (RS, RE, Z, TRUNC)
  return
! HERE J IS TOO LARGE FOR SINGLE-PRECISION MULTIPLICATION
  110 J1 = J/B
  J2 = J - J1*B
! FORM PRODUCT
  DO 120 IJ = 1, T4
  C1 = C/B
  C2 = C - B*C1
  I = T1 - IJ
  IX = 0
  if (I > 0) IX = X(I+2)
  RI = J2*IX + C2
  IS = RI/B
  C = J1*IX + C1 + IS
  120 R(I+4) = RI - B*IS
  if (C) 130, 100, 80
! CAN ONLY GET HERE if INTEGER OVERFLOW OCCURRED
  130 call MPCHK (1, 4)
  WRITE (LUN, 140)
  140 FORMAT (' *** INTEGER OVERFLOW IN MPMUL2, B TOO LARGE ***')
  call MPERR
  go to 10
end
