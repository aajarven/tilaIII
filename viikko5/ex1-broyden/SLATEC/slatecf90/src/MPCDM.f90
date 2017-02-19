subroutine MPCDM (DX, Z)
!
!! MPCDM is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPCDM-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! Converts double-precision number DX to multiple-precision Z.
! Some numbers will not convert exactly on machines with base
! other than two, four or sixteen. This routine is not called
! by any other routine in 'mp', so may be omitted if double-
! precision is not available.
!
! The argument Z(*) and the variable R in COMMON are both INTEGER
! arrays of size 30.  See the comments in the routine MPBLAS for the
! for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPCHK, MPDIVI, MPMULI, MPNZR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPCDM
  DOUBLE PRECISION DB, DJ, DX
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, Z(*), RS, RE, TP
!***FIRST EXECUTABLE STATEMENT  MPCDM
  call MPCHK (1, 4)
  I2 = T + 4
! CHECK SIGN
  if (DX) 20, 10, 30
! if DX = 0D0 RETURN 0
   10 Z(1) = 0
  return
! DX  <  0D0
   20 RS = -1
  DJ = -DX
  go to 40
! DX  >  0D0
   30 RS = 1
  DJ = DX
   40 IE = 0
   50 if (DJ < 1D0) go to 60
! INCREASE IE AND DIVIDE DJ BY 16.
  IE = IE + 1
  DJ = 0.0625D0*DJ
  go to 50
   60 if (DJ >= 0.0625D0) go to 70
  IE = IE - 1
  DJ = 16D0*DJ
  go to 60
! NOW DJ IS DY DIVIDED BY SUITABLE POWER OF 16
! SET EXPONENT TO 0
   70 RE = 0
  DB = DBLE(B)
! CONVERSION LOOP (ASSUME DOUBLE-PRECISION OPS. EXACT)
  DO 80 I = 1, I2
  DJ = DB*DJ
  R(I) = INT(DJ)
   80 DJ = DJ - DBLE(R(I))
! NORMALIZE RESULT
  call MPNZR (RS, RE, Z, 0)
  IB = MAX(7*B*B, 32767)/16
  TP = 1
! NOW MULTIPLY BY 16**IE
  if (IE) 90, 130, 110
   90 K = -IE
  DO 100 I = 1, K
  TP = 16*TP
  if ((TP <= IB).AND.(TP /= B).AND.(I < K)) go to 100
  call MPDIVI (Z, TP, Z)
  TP = 1
  100 CONTINUE
  return
  110 DO 120 I = 1, IE
  TP = 16*TP
  if ((TP <= IB).AND.(TP /= B).AND.(I < IE)) go to 120
  call MPMULI (Z, TP, Z)
  TP = 1
  120 CONTINUE
  130 RETURN
end
