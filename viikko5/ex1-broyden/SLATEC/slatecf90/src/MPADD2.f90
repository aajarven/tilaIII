subroutine MPADD2 (X, Y, Z, Y1, TRUNC)
!
!! MPADD2 is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPADD2-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Called by MPADD, MPSUB etc.
!  X, Y and Z are MP numbers, Y1 and TRUNC are integers.
!  To force call by reference rather than value/result, Y1 is
!  declared as an array, but only Y1(1) is ever used.
!  Sets Z = X + Y1(1)*ABS(Y), where Y1(1) = +- Y(1).
!  If TRUNC  ==  0, R*-rounding is used;  otherwise, truncation.
!  R*-rounding is defined in the Kuki and Cody reference.
!
!  The arguments X(*), Y(*), and Z(*) are all INTEGER arrays of size
!  30.  See the comments in the routine MPBLAS for the reason for this
!  choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***REFERENCES  H. Kuki and W. J. Cody, A statistical study of floating
!                 point number systems, Communications of the ACM 16, 4
!                 (April 1973), pp. 223-230.
!               R. P. Brent, On the precision attainable with various
!                 floating-point number systems, IEEE Transactions on
!                 Computers C-22, 6 (June 1973), pp. 601-607.
!               R. P. Brent, A Fortran multiple-precision arithmetic
!                 package, ACM Transactions on Mathematical Software 4,
!                 1 (March 1978), pp. 57-70.
!               R. P. Brent, MP, a Fortran multiple-precision arithmetic
!                 package, Algorithm 524, ACM Transactions on Mathema-
!                 tical Software 4, 1 (March 1978), pp. 71-81.
!***ROUTINES CALLED  MPADD3, MPCHK, MPERR, MPNZR, MPSTR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920528  Added a REFERENCES section revised.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPADD2
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Y(*), Z(*), Y1(*), TRUNC
  INTEGER S, ED, RS, RE
!***FIRST EXECUTABLE STATEMENT  MPADD2
  if (X(1) /= 0) go to 20
   10 call MPSTR(Y, Z)
  Z(1) = Y1(1)
  return
   20 if (Y1(1) /= 0) go to 40
   30 call MPSTR (X, Z)
  return
! COMPARE SIGNS
   40 S = X(1)*Y1(1)
  if (ABS(S) <= 1) go to 60
  call MPCHK (1, 4)
  WRITE (LUN, 50)
   50 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN call TO MPADD2,', &
          ' POSSIBLE OVERWRITING PROBLEM ***')
  call MPERR
  Z(1) = 0
  return
! COMPARE EXPONENTS
   60 ED = X(2) - Y(2)
  MED = ABS(ED)
  if (ED) 90, 70, 120
! EXPONENTS EQUAL SO COMPARE SIGNS, THEN FRACTIONS if NEC.
   70 if (S > 0) go to 100
  DO 80 J = 1, T
  if (X(J+2) - Y(J+2)) 100, 80, 130
   80 CONTINUE
! RESULT IS ZERO
  Z(1) = 0
  return
! HERE EXPONENT(Y)  >=  EXPONENT(X)
   90 if (MED > T) go to 10
  100 RS = Y1(1)
  RE = Y(2)
  call MPADD3 (X, Y, S, MED, RE)
! NORMALIZE, ROUND OR TRUNCATE, AND RETURN
  110 call MPNZR (RS, RE, Z, TRUNC)
  return
! ABS(X)  >  ABS(Y)
  120 if (MED > T) go to 30
  130 RS = X(1)
  RE = X(2)
  call MPADD3 (Y, X, S, MED, RE)
  go to 110
end
