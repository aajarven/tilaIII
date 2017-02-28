subroutine ORTHO4 (USOL, IDMN, ZN, ZM, PERTRB)
!
!! ORTHO4 is subsidiary to SEPX4.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (ORTHO4-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine orthogonalizes the array USOL with respect to
!     the constant array in a weighted least squares norm.
!
!***SEE ALSO  SEPX4
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPL4
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  ORTHO4
!
  COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  DIMENSION       USOL(IDMN,*)           ,ZN(*)      ,ZM(*)
!***FIRST EXECUTABLE STATEMENT  ORTHO4
  ISTR = IS
  IFNL = MS
  JSTR = JS
  JFNL = NS
!
!     COMPUTE WEIGHTED INNER PRODUCTS
!
  UTE = 0.0
  ETE = 0.0
  DO  20 I=IS,MS
     II = I-IS+1
     DO  10 J=JS,NS
        JJ = J-JS+1
        ETE = ETE+ZM(II)*ZN(JJ)
        UTE = UTE+USOL(I,J)*ZM(II)*ZN(JJ)
   10    CONTINUE
   20 CONTINUE
!
!     SET PERTURBATION PARAMETER
!
  PERTRB = UTE/ETE
!
!     SUBTRACT OFF CONSTANT PERTRB
!
  DO  40 I=ISTR,IFNL
     DO  30 J=JSTR,JFNL
        USOL(I,J) = USOL(I,J)-PERTRB
   30    CONTINUE
   40 CONTINUE
  return
end