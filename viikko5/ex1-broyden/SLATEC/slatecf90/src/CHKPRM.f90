subroutine CHKPRM (INTL, IORDER, A, B, M, MBDCND, C, D, N, NBDCND, &
     COFX, COFY, IDMN, IERROR)
!
!! CHKPRM is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CHKPRM-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This program checks the input parameters for errors.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CHKPRM
!
  EXTERNAL        COFX       ,COFY
!***FIRST EXECUTABLE STATEMENT  CHKPRM
  IERROR = 1
  if (A >= B .OR. C >= D) RETURN
!
!     CHECK BOUNDARY SWITCHES
!
  IERROR = 2
  if (MBDCND < 0 .OR. MBDCND > 4) RETURN
  IERROR = 3
  if (NBDCND < 0 .OR. NBDCND > 4) RETURN
!
!     CHECK FIRST DIMENSION IN CALLING ROUTINE
!
  IERROR = 5
  if (IDMN  <  7) RETURN
!
!     CHECK M
!
  IERROR = 6
  if (M > (IDMN-1) .OR. M < 6) RETURN
!
!     CHECK N
!
  IERROR = 7
  if (N  <  5) RETURN
!
!     CHECK IORDER
!
  IERROR = 8
  if (IORDER /= 2 .AND. IORDER /= 4) RETURN
!
!     CHECK INTL
!
  IERROR = 9
  if (INTL /= 0 .AND. INTL /= 1) RETURN
!
!     CHECK THAT EQUATION IS ELLIPTIC
!
  DLX = (B-A)/M
  DLY = (D-C)/N
  DO  30 I=2,M
     XI = A+(I-1)*DLX
     call COFX (XI,AI,BI,CI)
     DO  20 J=2,N
        YJ = C+(J-1)*DLY
        call COFY (YJ,DJ,EJ,FJ)
        if (AI*DJ  >  0.0) go to  10
        IERROR = 10
        return
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
!
!     NO ERROR FOUND
!
  IERROR = 0
  return
end
