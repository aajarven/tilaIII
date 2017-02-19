subroutine MPCHK (I, J)
!
!! MPCHK is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPCHK-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Checks legality of B, T, M, MXR and LUN which should be set
!  in COMMON. The condition on MXR (the dimension of the EP arrays)
!  is that  MXR  >=  (I*T + J)
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  I1MACH, MPERR
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   891009  Removed unreferenced statement label.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPCHK
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R
!***FIRST EXECUTABLE STATEMENT  MPCHK
  LUN = I1MACH(4)
! NOW CHECK LEGALITY OF B, T AND M
  if (B > 1) go to 40
  WRITE (LUN, 30) B
   30 FORMAT (' *** B =', I10, ' ILLEGAL IN call TO MPCHK,'/ &
   ' PERHAPS NOT SET BEFORE call TO AN MP ROUTINE ***')
  call MPERR
   40 if (T > 1) go to 60
  WRITE (LUN, 50) T
   50 FORMAT (' *** T =', I10, ' ILLEGAL IN call TO MPCHK,'/ &
   ' PERHAPS NOT SET BEFORE call TO AN MP ROUTINE ***')
  call MPERR
   60 if (M > T) go to 80
  WRITE (LUN, 70)
   70 FORMAT (' *** M  <=  T IN call TO MPCHK,'/ &
   ' PERHAPS NOT SET BEFORE call TO AN MP ROUTINE ***')
  call MPERR
! 8*B*B-1 SHOULD BE REPRESENTABLE, if NOT WILL OVERFLOW
! AND MAY BECOME NEGATIVE, SO CHECK FOR THIS
   80 IB = 4*B*B - 1
  if ((IB > 0).AND.((2*IB+1) > 0)) go to 100
  WRITE (LUN, 90)
   90 FORMAT (' *** B TOO LARGE IN call TO MPCHK ***')
  call MPERR
! CHECK THAT SPACE IN COMMON IS SUFFICIENT
  100 MX = I*T + J
  if (MXR >= MX) RETURN
! HERE COMMON IS TOO SMALL, SO GIVE ERROR MESSAGE.
  WRITE (LUN, 110) I, J, MX, MXR, T
  110 FORMAT (' *** MXR TOO SMALL OR NOT SET TO DIM(R) BEFORE CALL', &
   ' TO AN MP ROUTINE *** ' / &
   ' *** MXR SHOULD BE AT LEAST', I3, '*T +', I4, ' =', I6, '  ***' &
   / ' *** ACTUALLY MXR =', I10, ', AND T =', I10, '  ***')
  call MPERR
  return
end
