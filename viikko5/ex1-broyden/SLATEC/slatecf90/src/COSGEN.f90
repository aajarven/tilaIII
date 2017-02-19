subroutine COSGEN (N, IJUMP, FNUM, FDEN, A)
!
!! COSGEN is subsidiary to GENBUN.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (COSGEN-S, CMPCSG-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine computes required cosine values in ascending
!     order.  When IJUMP  >  1 the routine computes values
!
!        2*COS(J*PI/L) , J=1,2,...,L and J  /=  0(MOD N/IJUMP+1)
!
!     where L = IJUMP*(N/IJUMP+1).
!
!
!     when IJUMP = 1 it computes
!
!            2*COS((J-FNUM)*PI/(N+FDEN)) ,  J=1, 2, ... ,N
!
!     where
!        FNUM = 0.5, FDEN = 0.0, for regular reduction values.
!        FNUM = 0.0, FDEN = 1.0, for B-R and C-R when ISTAG = 1
!        FNUM = 0.0, FDEN = 0.5, for B-R and C-R when ISTAG = 2
!        FNUM = 0.5, FDEN = 0.5, for B-R and C-R when ISTAG = 2
!                                in POISN2 only.
!
!***SEE ALSO  GENBUN
!***ROUTINES CALLED  PIMACH
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  COSGEN
  DIMENSION       A(*)
!
!
!***FIRST EXECUTABLE STATEMENT  COSGEN
  PI = PIMACH(DUM)
  if (N  ==  0) go to 105
  if (IJUMP  ==  1) go to 103
  K3 = N/IJUMP+1
  K4 = K3-1
  PIBYN = PI/(N+IJUMP)
  DO 102 K=1,IJUMP
     K1 = (K-1)*K3
     K5 = (K-1)*K4
     DO 101 I=1,K4
        X = K1+I
        K2 = K5+I
        A(K2) = -2.*COS(X*PIBYN)
  101    CONTINUE
  102 CONTINUE
  go to 105
  103 CONTINUE
  NP1 = N+1
  Y = PI/(N+FDEN)
  DO 104 I=1,N
     X = NP1-I-FNUM
     A(I) = 2.*COS(X*Y)
  104 CONTINUE
  105 CONTINUE
  return
end
