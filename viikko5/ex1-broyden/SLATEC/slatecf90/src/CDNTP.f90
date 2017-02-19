subroutine CDNTP (H, K, N, NQ, T, TOUT, YH, Y)
!
!! CDNTP interpolates the K-th derivative of Y at TOUT, using the data ...
!  in the YH array.  If K has a value greater than NQ, the NQ-th derivative ...
!  is calculated.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDNTP-S, DDNTP-D, CDNTP-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDNTP
  INTEGER I, J, JJ, K, KK, KUSED, N, NQ
  COMPLEX Y(*), YH(N,*)
  REAL FACTOR, H, R, T, TOUT
!***FIRST EXECUTABLE STATEMENT  CDNTP
  if (K  ==  0) THEN
    DO 10 I = 1,N
 10       Y(I) = YH(I,NQ+1)
    R = ((TOUT - T)/H)
    DO 20 JJ = 1,NQ
      J = NQ + 1 - JJ
      DO 20 I = 1,N
 20         Y(I) = YH(I,J) + R*Y(I)
  ELSE
    KUSED = MIN(K, NQ)
    FACTOR = 1.E0
    DO 40 KK = 1,KUSED
 40       FACTOR = FACTOR*(NQ+1-KK)
    DO 50 I = 1,N
 50       Y(I) = FACTOR*YH(I,NQ+1)
    R = ((TOUT - T)/H)
    DO 80 JJ = KUSED+1,NQ
      J = KUSED + 1 + NQ - JJ
      FACTOR = 1.E0
      DO 60 KK = 1,KUSED
 60         FACTOR = FACTOR*(J-KK)
      DO 70 I = 1,N
 70         Y(I) = FACTOR*YH(I,J) + R*Y(I)
 80       CONTINUE
    DO 100 I = 1,N
 100      Y(I) = Y(I)*H**(-KUSED)
  end if
  return
end
