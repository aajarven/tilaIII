subroutine MPCMD (X, DZ)
!
!! MPCMD is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPCMD-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Converts multiple-precision X to double-precision DZ. Assumes
!  X is in allowable range for double-precision numbers. There is
!  some loss of accuracy if the exponent is large.
!
!  The argument X(*) is INTEGER array of size 30.  See the comments in
!  the routine MPBLAS for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPCHK, MPERR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPCMD
  DOUBLE PRECISION DB, DZ, DZ2
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), TM
!***FIRST EXECUTABLE STATEMENT  MPCMD
  call MPCHK (1, 4)
  DZ = 0D0
  if (X(1) == 0) RETURN
  DB = DBLE(B)
  DO 10 I = 1, T
  DZ = DB*DZ + DBLE(X(I+2))
  TM = I
! CHECK if FULL DOUBLE-PRECISION ACCURACY ATTAINED
  DZ2 = DZ + 1D0
! TEST BELOW NOT ALWAYS EQUIVALENT TO - if (DZ2 <= DZ) go to 20,
! FOR EXAMPLE ON CYBER 76.
  if ((DZ2-DZ) <= 0D0) go to 20
   10 CONTINUE
! NOW ALLOW FOR EXPONENT
   20 DZ = DZ*(DB**(X(2)-TM))
! CHECK REASONABLENESS OF RESULT.
  if (DZ <= 0D0) go to 30
! LHS SHOULD BE  <=  0.5 BUT ALLOW FOR SOME ERROR IN LOG
  if (ABS(DBLE(X(2))-(LOG(DZ)/ &
      LOG(DBLE(B))+0.5D0)) > 0.6D0) go to 30
  if (X(1) < 0) DZ = -DZ
  return
! FOLLOWING MESSAGE INDICATES THAT X IS TOO LARGE OR SMALL -
! TRY USING MPCMDE INSTEAD.
   30 WRITE (LUN, 40)
   40 FORMAT (' *** FLOATING-POINT OVER/UNDER-FLOW IN MPCMD ***')
  call MPERR
  return
end
