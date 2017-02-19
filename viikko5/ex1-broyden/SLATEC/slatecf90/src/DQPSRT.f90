subroutine DQPSRT (LIMIT, LAST, MAXERR, ERMAX, ELIST, IORD, NRMAX)
!
!! DQPSRT sorts the local error estimates for a quadrature routine.
!
!***SUBSIDIARY
!***PURPOSE  This routine maintains the descending ordering in the
!            list of the local error estimated resulting from the
!            interval subdivision process. At each call two error
!            estimates are inserted using the sequential search
!            method, top-down for the largest error estimate and
!            bottom-up for the smallest error estimate.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QPSRT-S, DQPSRT-D)
!***KEYWORDS  SEQUENTIAL SORTING
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!           Ordering routine
!           Standard fortran subroutine
!           Double precision version
!
!           PARAMETERS (MEANING AT OUTPUT)
!              LIMIT  - Integer
!                       Maximum number of error estimates the list
!                       can contain
!
!              LAST   - Integer
!                       Number of error estimates currently in the list
!
!              MAXERR - Integer
!                       MAXERR points to the NRMAX-th largest error
!                       estimate currently in the list
!
!              ERMAX  - Double precision
!                       NRMAX-th largest error estimate
!                       ERMAX = ELIST(MAXERR)
!
!              ELIST  - Double precision
!                       Vector of dimension LAST containing
!                       the error estimates
!
!              IORD   - Integer
!                       Vector of dimension LAST, the first K elements
!                       of which contain pointers to the error
!                       estimates, such that
!                       ELIST(IORD(1)),...,  ELIST(IORD(K))
!                       form a decreasing sequence, with
!                       K = LAST if LAST <= (LIMIT/2+2), and
!                       K = LIMIT+1-LAST otherwise
!
!              NRMAX  - Integer
!                       MAXERR = IORD(NRMAX)
!
!***SEE ALSO  DQAGE, DQAGIE, DQAGPE, DQAWSE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQPSRT
!
  DOUBLE PRECISION ELIST,ERMAX,ERRMAX,ERRMIN
  INTEGER I,IBEG,IDO,IORD,ISUCC,J,JBND,JUPBN,K,LAST,LIMIT,MAXERR, &
    NRMAX
  DIMENSION ELIST(*),IORD(*)
!
!           CHECK WHETHER THE LIST CONTAINS MORE THAN
!           TWO ERROR ESTIMATES.
!
!***FIRST EXECUTABLE STATEMENT  DQPSRT
  if ( LAST > 2) go to 10
  IORD(1) = 1
  IORD(2) = 2
  go to 90
!
!           THIS PART OF THE ROUTINE IS ONLY EXECUTED IF, DUE TO A
!           DIFFICULT INTEGRAND, SUBDIVISION INCREASED THE ERROR
!           ESTIMATE. IN THE NORMAL CASE THE INSERT PROCEDURE SHOULD
!           START AFTER THE NRMAX-TH LARGEST ERROR ESTIMATE.
!
   10 ERRMAX = ELIST(MAXERR)
  if ( NRMAX == 1) go to 30
  IDO = NRMAX-1
  DO 20 I = 1,IDO
    ISUCC = IORD(NRMAX-1)
! ***JUMP OUT OF DO-LOOP
    if ( ERRMAX <= ELIST(ISUCC)) go to 30
    IORD(NRMAX) = ISUCC
    NRMAX = NRMAX-1
   20    CONTINUE
!
!           COMPUTE THE NUMBER OF ELEMENTS IN THE LIST TO BE MAINTAINED
!           IN DESCENDING ORDER. THIS NUMBER DEPENDS ON THE NUMBER OF
!           SUBDIVISIONS STILL ALLOWED.
!
   30 JUPBN = LAST
  if ( LAST > (LIMIT/2+2)) JUPBN = LIMIT+3-LAST
  ERRMIN = ELIST(LAST)
!
!           INSERT ERRMAX BY TRAVERSING THE LIST TOP-DOWN,
!           STARTING COMPARISON FROM THE ELEMENT ELIST(IORD(NRMAX+1)).
!
  JBND = JUPBN-1
  IBEG = NRMAX+1
  if ( IBEG > JBND) go to 50
  DO 40 I=IBEG,JBND
    ISUCC = IORD(I)
! ***JUMP OUT OF DO-LOOP
    if ( ERRMAX >= ELIST(ISUCC)) go to 60
    IORD(I-1) = ISUCC
   40 CONTINUE
   50 IORD(JBND) = MAXERR
  IORD(JUPBN) = LAST
  go to 90
!
!           INSERT ERRMIN BY TRAVERSING THE LIST BOTTOM-UP.
!
   60 IORD(I-1) = MAXERR
  K = JBND
  DO 70 J=I,JBND
    ISUCC = IORD(K)
! ***JUMP OUT OF DO-LOOP
    if ( ERRMIN < ELIST(ISUCC)) go to 80
    IORD(K+1) = ISUCC
    K = K-1
   70 CONTINUE
  IORD(I) = LAST
  go to 90
   80 IORD(K+1) = LAST
!
!           SET MAXERR AND ERMAX.
!
   90 MAXERR = IORD(NRMAX)
  ERMAX = ELIST(MAXERR)
  return
end
