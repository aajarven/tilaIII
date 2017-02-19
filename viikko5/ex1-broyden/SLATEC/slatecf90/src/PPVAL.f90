function PPVAL (LDC, C, XI, LXI, K, IDERIV, X, INPPV)
!
!! PPVAL calculates the value of the IDERIV-th derivative of the B-spline ...
!  from the PP-representation.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      SINGLE PRECISION (PPVAL-S, DPPVAL-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract
!         PPVAL is the PPVALU function of the reference.
!
!         PPVAL calculates (at X) the value of the IDERIV-th
!         derivative of the B-spline from the PP-representation
!         (C,XI,LXI,K).  The Taylor expansion about XI(J) for X in
!         the interval XI(J)  <=  X  <  XI(J+1) is evaluated, J=1,LXI.
!         Right limiting values at X=XI(J) are obtained.  PPVAL will
!         extrapolate beyond XI(1) and XI(LXI+1).
!
!         To obtain left limiting values (left derivatives) at XI(J),
!         replace LXI by J-1 and set X=XI(J),J=2,LXI+1.
!
!     Description of Arguments
!         Input
!          LDC     - leading dimension of C matrix, LDC  >=  K
!          C       - matrix of dimension at least (K,LXI) containing
!                    right derivatives at break points XI(*).
!          XI      - break point vector of length LXI+1
!          LXI     - number of polynomial pieces
!          K       - order of B-spline, K  >=  1
!          IDERIV  - order of the derivative, 0  <=  IDERIV  <=  K-1
!                    IDERIV=0 gives the B-spline value
!          X       - argument, XI(1)  <=  X  <=  XI(LXI+1)
!          INPPV   - an initialization parameter which must be set
!                    to 1 the first time PPVAL is called.
!
!         Output
!          INPPV   - INPPV contains information for efficient process-
!                    ing after the initial call and INPPV must not
!                    be changed by the user.  Distinct splines require
!                    distinct INPPV parameters.
!          PPVAL   - value of the IDERIV-th derivative at X
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  INTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  PPVAL
!
  INTEGER I, IDERIV, INPPV, J, K, LDC, LXI, NDUMMY
  REAL C, DX, FLTK, X, XI
  DIMENSION XI(*), C(LDC,*)
!***FIRST EXECUTABLE STATEMENT  PPVAL
  PPVAL = 0.0E0
  if ( K < 1) go to 90
  if ( LDC < K) go to 80
  if ( LXI < 1) go to 85
  if ( IDERIV < 0 .OR. IDERIV >= K) go to 95
  I = K - IDERIV
  FLTK = I
  call INTRV(XI, LXI, X, INPPV, I, NDUMMY)
  DX = X - XI(I)
  J = K
   10 PPVAL = (PPVAL/FLTK)*DX + C(J,I)
  J = J - 1
  FLTK = FLTK - 1.0E0
  if (FLTK > 0.0E0) go to 10
  return
!
!
   80 CONTINUE
  call XERMSG ('SLATEC', 'PPVAL', 'LDC DOES NOT SATISFY LDC >= K', &
     2, 1)
  return
   85 CONTINUE
  call XERMSG ('SLATEC', 'PPVAL', 'LXI DOES NOT SATISFY LXI >= 1', &
     2, 1)
  return
   90 CONTINUE
  call XERMSG ('SLATEC', 'PPVAL', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
   95 CONTINUE
  call XERMSG ('SLATEC', 'PPVAL', &
     'IDERIV DOES NOT SATISFY 0 <= IDERIV < K', 2, 1)
  return
end
