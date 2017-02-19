subroutine SDSCL (HMAX, N, NQ, RMAX, H, RC, RH, YH)
!
!! SDSCL rescales the YH array whenever the step size is changed.
!!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      SINGLE PRECISION (SDSCL-S, DDSCL-D, CDSCL-C)
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
!***END PROLOGUE  SDSCL
  INTEGER I, J, N, NQ
  REAL H, HMAX, RC, RH, RMAX, R1, YH(N,*)
!***FIRST EXECUTABLE STATEMENT  SDSCL
  if (H  <  1.E0) THEN
    RH = MIN(ABS(H)*RH, ABS(H)*RMAX, HMAX)/ABS(H)
  ELSE
    RH = MIN(RH, RMAX, HMAX/ABS(H))
  end if
  R1 = 1.E0
  DO 10 J = 1,NQ
    R1 = R1*RH
    DO 10 I = 1,N
 10       YH(I,J+1) = YH(I,J+1)*R1
  H = H*RH
  RC = RC*RH
  return
end
