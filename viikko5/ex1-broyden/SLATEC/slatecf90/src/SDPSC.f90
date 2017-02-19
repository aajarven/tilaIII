subroutine SDPSC (KSGN, N, NQ, YH)
!
!! SDPSC computes the predicted YH values by effectively multiplying ...
!  the YH array by the Pascal triangle
!            matrix when KSGN is +1, and performs the inverse function
!            when KSGN is -1.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      SINGLE PRECISION (SDPSC-S, DDPSC-D, CDPSC-C)
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
!***END PROLOGUE  SDPSC
  INTEGER I, J, J1, J2, KSGN, N, NQ
  REAL YH(N,*)
!***FIRST EXECUTABLE STATEMENT  SDPSC
  if (KSGN  >  0) THEN
    DO 10 J1 = 1,NQ
      DO 10 J2 = J1,NQ
        J = NQ - J2 + J1
        DO 10 I = 1,N
 10           YH(I,J) = YH(I,J) + YH(I,J+1)
  ELSE
    DO 30 J1 = 1,NQ
      DO 30 J2 = J1,NQ
        J = NQ - J2 + J1
        DO 30 I = 1,N
 30           YH(I,J) = YH(I,J) - YH(I,J+1)
  end if
  return
end
