FUNCTION PCHDF (K, X, S, IERR)
!
!! PCHDF computes divided differences for PCHCE and PCHSP.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      SINGLE PRECISION (PCHDF-S, DPCHDF-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!          PCHDF:   PCHIP Finite Difference Formula
!
!     Uses a divided difference formulation to compute a K-point approx-
!     imation to the derivative at X(K) based on the data in X and S.
!
!     Called by  PCHCE  and  PCHSP  to compute 3- and 4-point boundary
!     derivative approximations.
!
! ----------------------------------------------------------------------
!
!     On input:
!        K      is the order of the desired derivative approximation.
!               K must be at least 3 (error return if not).
!        X      contains the K values of the independent variable.
!               X need not be ordered, but the values **MUST** be
!               distinct.  (Not checked here.)
!        S      contains the associated slope values:
!                  S(I) = (F(I+1)-F(I))/(X(I+1)-X(I)), I=1(1)K-1.
!               (Note that S need only be of length K-1.)
!
!     On return:
!        S      will be destroyed.
!        IERR   will be set to -1 if K < 2 .
!        PCHDF  will be set to the desired derivative approximation if
!               IERR=0 or to zero if IERR=-1.
!
! ----------------------------------------------------------------------
!
!***SEE ALSO  PCHCE, PCHSP
!***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!                 Verlag, New York, 1978, pp. 10-16.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820503  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870813  Minor cosmetic changes.
!   890411  Added SAVE statements (Vers. 3.2).
!   890411  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   920429  Revised format and order of references.  (WRB,FNF)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  PCHDF
!
!**End
!
!  DECLARE ARGUMENTS.
!
  REAL PCHDF
  INTEGER  K, IERR
  REAL  X(K), S(K)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  I, J
  REAL  VALUE, ZERO
  SAVE ZERO
  DATA  ZERO /0./
!
!  CHECK FOR LEGAL VALUE OF K.
!
!***FIRST EXECUTABLE STATEMENT  PCHDF
  if (K  <  3)  go to 5001
!
!  COMPUTE COEFFICIENTS OF INTERPOLATING POLYNOMIAL.
!
  DO 10  J = 2, K-1
     DO 9  I = 1, K-J
        S(I) = (S(I+1)-S(I))/(X(I+J)-X(I))
    9    CONTINUE
   10 CONTINUE
!
!  EVALUATE DERIVATIVE AT X(K).
!
  VALUE = S(1)
  DO 20  I = 2, K-1
     VALUE = S(I) + VALUE*(X(K)-X(I))
   20 CONTINUE
!
!  NORMAL RETURN.
!
  IERR = 0
  PCHDF = VALUE
  return
!
!  ERROR RETURN.
!
 5001 CONTINUE
!     K < 3 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'PCHDF', 'K LESS THAN THREE', IERR, 1)
  PCHDF = ZERO
  return
!------------- LAST LINE OF PCHDF FOLLOWS ------------------------------
end
