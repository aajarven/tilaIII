subroutine MINSOL (USOL, IDMN, ZN, ZM, PERTB)
!
!! MINSOL is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (MINSOL-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine orthogonalizes the array USOL with respect to
!     the constant array in a weighted least squares norm.
!
!     Entry at MINSOL occurs when the final solution is
!     to be minimized with respect to the weighted
!     least squares norm.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  MINSOL
!
  COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  DIMENSION       USOL(IDMN,*)           ,ZN(*)      ,ZM(*)
!***FIRST EXECUTABLE STATEMENT  MINSOL
  ISTR = 1
  IFNL = K
  JSTR = 1
  JFNL = L
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
