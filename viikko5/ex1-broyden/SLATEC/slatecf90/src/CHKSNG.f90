subroutine CHKSNG (MBDCND, NBDCND, ALPHA, BETA, GAMA, XNU, COFX, &
     COFY, SINGLR)
!
!! CHKSNG is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CHKSNG-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine checks if the PDE SEPELI
!     must solve is a singular operator.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CHKSNG
!
  COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  LOGICAL         SINGLR
!***FIRST EXECUTABLE STATEMENT  CHKSNG
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
   10 if (NBDCND  /=  3) go to  20
  if (GAMA /= 0.0 .OR. XNU /= 0.0) RETURN
   20 CONTINUE
!
!     CHECK THAT NON-DERIVATIVE COEFFICIENT FUNCTIONS
!     ARE ZERO
!
  DO  30 I=IS,MS
     XI = AIT+(I-1)*DLX
     call COFX (XI,AI,BI,CI)
     if (CI  /=  0.0) RETURN
   30 CONTINUE
  DO  40 J=JS,NS
     YJ = CIT+(J-1)*DLY
     call COFY (YJ,DJ,EJ,FJ)
     if (FJ  /=  0.0) RETURN
   40 CONTINUE
!
!     THE OPERATOR MUST BE SINGULAR if THIS POINT IS REACHED
!
  SINGLR = .TRUE.
  return
end
