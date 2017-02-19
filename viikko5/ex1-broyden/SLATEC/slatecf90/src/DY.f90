subroutine DY (U, IDMN, I, J, UYYY, UYYYY)
!
!! DY is subsidiary to SEPELI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DY-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This program computes second order finite difference
!     approximations to the third and fourth Y
!     partial derivatives of U at the (I,J) mesh point.
!
!***SEE ALSO  SEPELI
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPLPCM
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  DY
!
  COMMON /SPLPCM/ KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  DIMENSION       U(IDMN,*)
!***FIRST EXECUTABLE STATEMENT  DY
  if (J > 2 .AND. J < (L-1)) go to  50
  if (J  ==  1) go to  10
  if (J  ==  2) go to  30
  if (J  ==  L-1) go to  60
  if (J  ==  L) go to  80
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C
!
   10 if (KSWY  ==  1) go to  20
  UYYY = (-5.0*U(I,1)+18.0*U(I,2)-24.0*U(I,3)+14.0*U(I,4)- &
                                                   3.0*U(I,5))/TDLY3
  UYYYY = (3.0*U(I,1)-14.0*U(I,2)+26.0*U(I,3)-24.0*U(I,4)+ &
                                        11.0*U(I,5)-2.0*U(I,6))/DLY4
  return
!
!     PERIODIC AT X=A
!
   20 UYYY = (-U(I,L-2)+2.0*U(I,L-1)-2.0*U(I,2)+U(I,3))/TDLY3
  UYYYY = (U(I,L-2)-4.0*U(I,L-1)+6.0*U(I,1)-4.0*U(I,2)+U(I,3))/DLY4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=C+DLY
!
   30 if (KSWY  ==  1) go to  40
  UYYY = (-3.0*U(I,1)+10.0*U(I,2)-12.0*U(I,3)+6.0*U(I,4)-U(I,5))/ &
         TDLY3
  UYYYY = (2.0*U(I,1)-9.0*U(I,2)+16.0*U(I,3)-14.0*U(I,4)+6.0*U(I,5)- &
                                                        U(I,6))/DLY4
  return
!
!     PERIODIC AT Y=C+DLY
!
   40 UYYY = (-U(I,L-1)+2.0*U(I,1)-2.0*U(I,3)+U(I,4))/TDLY3
  UYYYY = (U(I,L-1)-4.0*U(I,1)+6.0*U(I,2)-4.0*U(I,3)+U(I,4))/DLY4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
!
   50 CONTINUE
  UYYY = (-U(I,J-2)+2.0*U(I,J-1)-2.0*U(I,J+1)+U(I,J+2))/TDLY3
  UYYYY = (U(I,J-2)-4.0*U(I,J-1)+6.0*U(I,J)-4.0*U(I,J+1)+U(I,J+2))/ &
          DLY4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D-DLY
!
   60 if (KSWY  ==  1) go to  70
  UYYY = (U(I,L-4)-6.0*U(I,L-3)+12.0*U(I,L-2)-10.0*U(I,L-1)+ &
                                                   3.0*U(I,L))/TDLY3
  UYYYY = (-U(I,L-5)+6.0*U(I,L-4)-14.0*U(I,L-3)+16.0*U(I,L-2)- &
                                       9.0*U(I,L-1)+2.0*U(I,L))/DLY4
  return
!
!     PERIODIC AT Y=D-DLY
!
   70 CONTINUE
  UYYY = (-U(I,L-3)+2.0*U(I,L-2)-2.0*U(I,1)+U(I,2))/TDLY3
  UYYYY = (U(I,L-3)-4.0*U(I,L-2)+6.0*U(I,L-1)-4.0*U(I,1)+U(I,2))/ &
          DLY4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT Y=D
!
   80 UYYY = -(3.0*U(I,L-4)-14.0*U(I,L-3)+24.0*U(I,L-2)-18.0*U(I,L-1)+ &
                                                   5.0*U(I,L))/TDLY3
  UYYYY = (-2.0*U(I,L-5)+11.0*U(I,L-4)-24.0*U(I,L-3)+26.0*U(I,L-2)- &
                                      14.0*U(I,L-1)+3.0*U(I,L))/DLY4
  return
end
