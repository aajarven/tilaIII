subroutine PCHKT (N, X, KNOTYP, T)
!
!! PCHKT computes B-spline knot sequence for PCHBS.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E3
!***TYPE      SINGLE PRECISION (PCHKT-S, DPCHKT-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!     Set a knot sequence for the B-spline representation of a PCH
!     function with breakpoints X.  All knots will be at least double.
!     Endknots are set as:
!        (1) quadruple knots at endpoints if KNOTYP=0;
!        (2) extrapolate the length of end interval if KNOTYP=1;
!        (3) periodic if KNOTYP=2.
!
!  Input arguments:  N, X, KNOTYP.
!  Output arguments:  T.
!
!  Restrictions/assumptions:
!     1. N >= 2 .  (not checked)
!     2. X(i) < X(i+1), i=1,...,N .  (not checked)
!     3. 0 <= KNOTYP <= 2 .  (Acts like KNOTYP=0 for any other value.)
!
!***SEE ALSO  PCHBS
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   870701  DATE WRITTEN
!   900405  Converted Fortran to upper case.
!   900410  Converted prologue to SLATEC 4.0 format.
!   900410  Minor cosmetic changes.
!   930514  Changed NKNOTS from an output to an input variable.  (FNF)
!   930604  Removed unused variable NKNOTS from argument list.  (FNF)
!***END PROLOGUE  PCHKT
!
!*Internal Notes:
!
!  Since this is subsidiary to PCHBS, which validates its input before
!  calling, it is unnecessary for such validation to be done here.
!
!**End
!
!  Declare arguments.
!
  INTEGER  N, KNOTYP
  REAL  X(*), T(*)
!
!  Declare local variables.
!
  INTEGER  J, K, NDIM
  REAL  HBEG, HEND
!***FIRST EXECUTABLE STATEMENT  PCHKT
!
!  Initialize.
!
  NDIM = 2*N
!
!  Set interior knots.
!
  J = 1
  DO 20  K = 1, N
     J = J + 2
     T(J) = X(K)
     T(J+1) = T(J)
   20 CONTINUE
!     Assertion:  At this point T(3),...,T(NDIM+2) have been set and
!                 J=NDIM+1.
!
!  Set end knots according to KNOTYP.
!
  HBEG = X(2) - X(1)
  HEND = X(N) - X(N-1)
  if (KNOTYP == 1 )  THEN
!          Extrapolate.
     T(2) = X(1) - HBEG
     T(NDIM+3) = X(N) + HEND
  ELSE if ( KNOTYP == 2 )  THEN
!          Periodic.
     T(2) = X(1) - HEND
     T(NDIM+3) = X(N) + HBEG
  ELSE
!          Quadruple end knots.
     T(2) = X(1)
     T(NDIM+3) = X(N)
  end if
  T(1) = T(2)
  T(NDIM+4) = T(NDIM+3)
!
!  Terminate.
!
  return
!------------- LAST LINE OF PCHKT FOLLOWS ------------------------------
end
