subroutine DEFER (COFX, COFY, IDMN, USOL, GRHS)
!
!! DEFER is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DEFER-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine first approximates the truncation error given by
!     TRUN1(X,Y)=DLX**2*TX+DLY**2*TY where
!     TX=AFUN(X)*UXXXX/12.0+BFUN(X)*UXXX/6.0 on the interior and
!     at the boundaries if periodic (here UXXX,UXXXX are the third
!     and fourth partial derivatives of U with respect to X).
!     TX is of the form AFUN(X)/3.0*(UXXXX/4.0+UXXX/DLX)
!     at X=A or X=B if the boundary condition there is mixed.
!     TX=0.0 along specified boundaries.  TY has symmetric form
!     in Y with X,AFUN(X),BFUN(X) replaced by Y,DFUN(Y),EFUN(Y).
!     The second order solution in USOL is used to approximate
!     (via second order finite differencing) the truncation error
!     and the result is added to the right hand side in GRHS
!     and then transferred to USOL to be used as a new right
!     hand side when calling BLKTRI for a fourth order solution.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  DX, DY
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  DEFER
!
  COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  DIMENSION       GRHS(IDMN,*)           ,USOL(IDMN,*)
  EXTERNAL        COFX       ,COFY
!***FIRST EXECUTABLE STATEMENT  DEFER
  DO  40 J=JS,NS
     YJ = CIT+(J-1)*DLY
     call COFY (YJ,DJ,EJ,FJ)
     DO  30 I=IS,MS
        XI = AIT+(I-1)*DLX
        call COFX (XI,AI,BI,CI)
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT (XI,YJ)
!
        call DX (USOL,IDMN,I,J,UXXX,UXXXX)
        call DY (USOL,IDMN,I,J,UYYY,UYYYY)
        TX = AI*UXXXX/12.0+BI*UXXX/6.0
        TY = DJ*UYYYY/12.0+EJ*UYYY/6.0
!
!     RESET FORM OF TRUNCATION if AT BOUNDARY WHICH IS NON-PERIODIC
!
        if (KSWX == 1 .OR. (I > 1 .AND. I < K)) go to  10
        TX = AI/3.0*(UXXXX/4.0+UXXX/DLX)
   10       if (KSWY == 1 .OR. (J > 1 .AND. J < L)) go to  20
        TY = DJ/3.0*(UYYYY/4.0+UYYY/DLY)
   20       GRHS(I,J) = GRHS(I,J)+DLX**2*TX+DLY**2*TY
   30    CONTINUE
   40 CONTINUE
!
!     RESET THE RIGHT HAND SIDE IN USOL
!
  DO  60 I=IS,MS
     DO  50 J=JS,NS
        USOL(I,J) = GRHS(I,J)
   50    CONTINUE
   60 CONTINUE
  return
end
