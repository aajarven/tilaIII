subroutine DX4 (U, IDMN, I, J, UXXX, UXXXX)
!
!! DX4 is subsidiary to SEPX4.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DX4-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This program computes second order finite difference
!     approximations to the third and fourth X
!     partial derivatives of U at the (I,J) mesh point.
!
!***SEE ALSO  SEPX4
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    SPL4
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  DX4
!
  COMMON /SPL4/   KSWX       ,KSWY       ,K          ,L          , &
                  AIT        ,BIT        ,CIT        ,DIT        , &
                  MIT        ,NIT        ,IS         ,MS         , &
                  JS         ,NS         ,DLX        ,DLY        , &
                  TDLX3      ,TDLY3      ,DLX4       ,DLY4
  DIMENSION       U(IDMN,*)
!***FIRST EXECUTABLE STATEMENT  DX4
  if (I > 2 .AND. I < (K-1)) go to  50
  if (I  ==  1) go to  10
  if (I  ==  2) go to  30
  if (I  ==  K-1) go to  60
  if (I  ==  K) go to  80
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A
!
   10 if (KSWX  ==  1) go to  20
  UXXX = (-5.0*U(1,J)+18.0*U(2,J)-24.0*U(3,J)+14.0*U(4,J)- &
                                                 3.0*U(5,J))/(TDLX3)
  UXXXX = (3.0*U(1,J)-14.0*U(2,J)+26.0*U(3,J)-24.0*U(4,J)+ &
                                        11.0*U(5,J)-2.0*U(6,J))/DLX4
  return
!
!     PERIODIC AT X=A
!
   20 UXXX = (-U(K-2,J)+2.0*U(K-1,J)-2.0*U(2,J)+U(3,J))/(TDLX3)
  UXXXX = (U(K-2,J)-4.0*U(K-1,J)+6.0*U(1,J)-4.0*U(2,J)+U(3,J))/DLX4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=A+DLX
!
   30 if (KSWX  ==  1) go to  40
  UXXX = (-3.0*U(1,J)+10.0*U(2,J)-12.0*U(3,J)+6.0*U(4,J)-U(5,J))/ &
         TDLX3
  UXXXX = (2.0*U(1,J)-9.0*U(2,J)+16.0*U(3,J)-14.0*U(4,J)+6.0*U(5,J)- &
                                                        U(6,J))/DLX4
  return
!
!     PERIODIC AT X=A+DLX
!
   40 UXXX = (-U(K-1,J)+2.0*U(1,J)-2.0*U(3,J)+U(4,J))/(TDLX3)
  UXXXX = (U(K-1,J)-4.0*U(1,J)+6.0*U(2,J)-4.0*U(3,J)+U(4,J))/DLX4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS ON THE INTERIOR
!
   50 CONTINUE
  UXXX = (-U(I-2,J)+2.0*U(I-1,J)-2.0*U(I+1,J)+U(I+2,J))/TDLX3
  UXXXX = (U(I-2,J)-4.0*U(I-1,J)+6.0*U(I,J)-4.0*U(I+1,J)+U(I+2,J))/ &
          DLX4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B-DLX
!
   60 if (KSWX  ==  1) go to  70
  UXXX = (U(K-4,J)-6.0*U(K-3,J)+12.0*U(K-2,J)-10.0*U(K-1,J)+ &
                                                   3.0*U(K,J))/TDLX3
  UXXXX = (-U(K-5,J)+6.0*U(K-4,J)-14.0*U(K-3,J)+16.0*U(K-2,J)- &
                                       9.0*U(K-1,J)+2.0*U(K,J))/DLX4
  return
!
!     PERIODIC AT X=B-DLX
!
   70 UXXX = (-U(K-3,J)+2.0*U(K-2,J)-2.0*U(1,J)+U(2,J))/TDLX3
  UXXXX = (U(K-3,J)-4.0*U(K-2,J)+6.0*U(K-1,J)-4.0*U(1,J)+U(2,J))/ &
          DLX4
  return
!
!     COMPUTE PARTIAL DERIVATIVE APPROXIMATIONS AT X=B
!
   80 UXXX = -(3.0*U(K-4,J)-14.0*U(K-3,J)+24.0*U(K-2,J)-18.0*U(K-1,J)+ &
                                                   5.0*U(K,J))/TDLX3
  UXXXX = (-2.0*U(K-5,J)+11.0*U(K-4,J)-24.0*U(K-3,J)+26.0*U(K-2,J)- &
                                      14.0*U(K-1,J)+3.0*U(K,J))/DLX4
  return
end
