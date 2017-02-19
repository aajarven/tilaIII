subroutine BSPLVD (T, K, X, ILEFT, VNIKX, NDERIV)
!
!! BSPLVD is subsidiary to FC.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BSPLVD-S, DFSPVD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! Calculates value and deriv.s of all B-splines which do not vanish at X
!
!  Fill VNIKX(J,IDERIV), J=IDERIV, ... ,K  with nonzero values of
!  B-splines of order K+1-IDERIV , IDERIV=NDERIV, ... ,1, by repeated
!  calls to BSPLVN
!
!***SEE ALSO  FC
!***ROUTINES CALLED  BSPLVN
!***REVISION HISTORY  (YYMMDD)
!   780801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  BSPLVD
  DIMENSION T(*),VNIKX(K,*)
  DIMENSION A(20,20)
!***FIRST EXECUTABLE STATEMENT  BSPLVD
  call BSPLVN(T,K+1-NDERIV,1,X,ILEFT,VNIKX(NDERIV,NDERIV))
  if (NDERIV  <=  1)               go to 99
  IDERIV = NDERIV
  DO 15 I=2,NDERIV
     IDERVM = IDERIV-1
     DO 11 J=IDERIV,K
   11       VNIKX(J-1,IDERVM) = VNIKX(J,IDERIV)
     IDERIV = IDERVM
     call BSPLVN(T,0,2,X,ILEFT,VNIKX(IDERIV,IDERIV))
   15    CONTINUE
!
  DO 20 I=1,K
     DO 19 J=1,K
   19       A(I,J) = 0.
   20    A(I,I) = 1.
  KMD = K
  DO 40 M=2,NDERIV
     KMD = KMD-1
     FKMD = KMD
     I = ILEFT
     J = K
   21       JM1 = J-1
        IPKMD = I + KMD
        DIFF = T(IPKMD) - T(I)
        if (JM1  ==  0)            go to 26
        if (DIFF  ==  0.)          go to 25
        DO 24 L=1,J
   24          A(L,J) = (A(L,J) - A(L,J-1))/DIFF*FKMD
   25       J = JM1
        I = I - 1
                                   go to 21
   26    if (DIFF  ==  0.)             go to 30
     A(1,1) = A(1,1)/DIFF*FKMD
!
   30    DO 40 I=1,K
        V = 0.
        JLOW = MAX(I,M)
        DO 35 J=JLOW,K
   35          V = A(I,J)*VNIKX(J,M) + V
   40       VNIKX(I,M) = V
   99                                  return
end
