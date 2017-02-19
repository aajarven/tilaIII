subroutine CHKSN4 (MBDCND, NBDCND, ALPHA, BETA, COFX, SINGLR)
!
!! CHKSN4 is subsidiary to SEPX4.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CHKSN4-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine checks if the PDE SEPX4
!     must solve is a singular operator.
!
!***SEE ALSO  SEPX4
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPL4
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CHKSN4
!
  COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  LOGICAL         SINGLR
  EXTERNAL COFX
!***FIRST EXECUTABLE STATEMENT  CHKSN4
  SINGLR = .FALSE.
!
!     CHECK if THE BOUNDARY CONDITIONS ARE
!     ENTIRELY PERIODIC AND/OR MIXED
!
  if ((MBDCND /= 0 .AND. MBDCND /= 3) .OR. &
      (NBDCND /= 0 .AND. NBDCND /= 3)) RETURN
!
!     CHECK THAT MIXED CONDITIONS ARE PURE NEUMAN
!
  if (MBDCND  /=  3) go to  10
  if (ALPHA /= 0.0 .OR. BETA /= 0.0) RETURN
   10 CONTINUE
!
!     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
!     ARE ZERO
!
  DO  30 I=IS,MS
     XI = AIT+(I-1)*DLX
     call COFX (XI,AI,BI,CI)
     if (CI  /=  0.0) RETURN
   30 CONTINUE
!
!     THE OPERATOR MUST BE SINGULAR if THIS POINT IS REACHED
!
  SINGLR = .TRUE.
  return
end
