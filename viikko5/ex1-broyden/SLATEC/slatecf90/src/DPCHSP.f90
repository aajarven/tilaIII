subroutine DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!! DPCHSP sets derivatives needed to determine the Hermite representation ...
!  of the cubic spline interpolant to given data, with specified boundary
!  conditions.
!
!***LIBRARY   SLATEC (PCHIP)
!***CATEGORY  E1A
!***TYPE      DOUBLE PRECISION (PCHSP-S, DPCHSP-D)
!***KEYWORDS  CUBIC HERMITE INTERPOLATION, PCHIP,
!             PIECEWISE CUBIC INTERPOLATION, SPLINE INTERPOLATION
!***AUTHOR  Fritsch, F. N., (LLNL)
!             Lawrence Livermore National Laboratory
!             P.O. Box 808  (L-316)
!             Livermore, CA  94550
!             FTS 532-4275, (510) 422-4275
!***DESCRIPTION
!
!          DPCHSP:   Piecewise Cubic Hermite Spline
!
!     Computes the Hermite representation of the cubic spline inter-
!     polant to the data given in X and F satisfying the boundary
!     conditions specified by IC and VC.
!
!     To facilitate two-dimensional applications, includes an increment
!     between successive values of the F- and D arrays.
!
!     The resulting piecewise cubic Hermite function may be evaluated
!     by DPCHFE or DPCHFD.
!
!     NOTE:  This is a modified version of C. de Boor's cubic spline
!            routine CUBSPL.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        PARAMETER  (INCFD = ...)
!        INTEGER  IC(2), N, NWK, IERR
!        DOUBLE PRECISION  VC(2), X(N), F(INCFD,N), D(INCFD,N), WK(NWK)
!
!        call  DPCHSP (IC, VC, N, X, F, D, INCFD, WK, NWK, IERR)
!
!   Parameters:
!
!     IC -- (input) integer array of length 2 specifying desired
!           boundary conditions:
!           IC(1) = IBEG, desired condition at beginning of data.
!           IC(2) = IEND, desired condition at end of data.
!
!           IBEG = 0  to set D(1) so that the third derivative is con-
!              tinuous at X(2).  This is the "not a knot" condition
!              provided by de Boor's cubic spline routine CUBSPL.
!              < This is the default boundary condition. >
!           IBEG = 1  if first derivative at X(1) is given in VC(1).
!           IBEG = 2  if second derivative at X(1) is given in VC(1).
!           IBEG = 3  to use the 3-point difference formula for D(1).
!                     (Reverts to the default b.c. if N < 3 .)
!           IBEG = 4  to use the 4-point difference formula for D(1).
!                     (Reverts to the default b.c. if N < 4 .)
!          NOTES:
!           1. An error return is taken if IBEG is out of range.
!           2. For the "natural" boundary condition, use IBEG=2 and
!              VC(1)=0.
!
!           IEND may take on the same values as IBEG, but applied to
!           derivative at X(N).  In case IEND = 1 or 2, the value is
!           given in VC(2).
!
!          NOTES:
!           1. An error return is taken if IEND is out of range.
!           2. For the "natural" boundary condition, use IEND=2 and
!              VC(2)=0.
!
!     VC -- (input) real*8 array of length 2 specifying desired boundary
!           values, as indicated above.
!           VC(1) need be set only if IC(1) = 1 or 2 .
!           VC(2) need be set only if IC(2) = 1 or 2 .
!
!     N -- (input) number of data points.  (Error return if N < 2 .)
!
!     X -- (input) real*8 array of independent variable values.  The
!           elements of X must be strictly increasing:
!                X(I-1)  <  X(I),  I = 2(1)N.
!           (Error return if not.)
!
!     F -- (input) real*8 array of dependent variable values to be
!           interpolated.  F(1+(I-1)*INCFD) is value corresponding to
!           X(I).
!
!     D -- (output) real*8 array of derivative values at the data
!           points.  These values will determine the cubic spline
!           interpolant with the requested boundary conditions.
!           The value corresponding to X(I) is stored in
!                D(1+(I-1)*INCFD),  I=1(1)N.
!           No other entries in D are changed.
!
!     INCFD -- (input) increment between successive values in F and D.
!           This argument is provided primarily for 2-D applications.
!           (Error return if  INCFD < 1 .)
!
!     WK -- (scratch) real*8 array of working storage.
!
!     NWK -- (input) length of work array.
!           (Error return if NWK < 2*N .)
!
!     IERR -- (output) error flag.
!           Normal return:
!              IERR = 0  (no errors).
!           "Recoverable" errors:
!              IERR = -1  if N < 2 .
!              IERR = -2  if INCFD < 1 .
!              IERR = -3  if the X-array is not strictly increasing.
!              IERR = -4  if IBEG < 0 or IBEG > 4 .
!              IERR = -5  if IEND < 0 of IEND > 4 .
!              IERR = -6  if both of the above are true.
!              IERR = -7  if NWK is too small.
!               NOTE:  The above errors are checked in the order listed,
!                   and following arguments have **NOT** been validated.
!             (The D array has not been changed in any of these cases.)
!              IERR = -8  in case of trouble solving the linear system
!                         for the interior derivative values.
!             (The D array may have been changed in this case.)
!             (             Do **NOT** use it!                )
!
!***REFERENCES  Carl de Boor, A Practical Guide to Splines, Springer-
!                 Verlag, New York, 1978, pp. 53-59.
!***ROUTINES CALLED  DPCHDF, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820503  DATE WRITTEN
!   820804  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   890206  Corrected XERROR calls.
!   890411  Added SAVE statements (Vers. 3.2).
!   890703  Corrected category record.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920429  Revised format and order of references.  (WRB,FNF)
!***END PROLOGUE  DPCHSP
!  Programming notes:
!
!     To produce a single precision version, simply:
!        a. Change DPCHSP to PCHSP wherever it occurs,
!        b. Change the double precision declarations to real, and
!        c. Change the constants ZERO, HALF, ... to single precision.
!
!  DECLARE ARGUMENTS.
!
  INTEGER  IC(2), N, INCFD, NWK, IERR
  DOUBLE PRECISION  VC(2), X(*), F(INCFD,*), D(INCFD,*), WK(2,*)
!
!  DECLARE LOCAL VARIABLES.
!
  INTEGER  IBEG, IEND, INDEX, J, NM1
  DOUBLE PRECISION  G, HALF, ONE, STEMP(3), THREE, TWO, XTEMP(4), &
    ZERO
  SAVE ZERO, HALF, ONE, TWO, THREE
  DOUBLE PRECISION  DPCHDF
!
  DATA  ZERO /0.D0/, HALF/.5D0/, ONE/1.D0/, TWO/2.D0/, THREE/3.D0/
!
!  VALIDITY-CHECK ARGUMENTS.
!
!***FIRST EXECUTABLE STATEMENT  DPCHSP
  if ( N < 2 )  go to 5001
  if ( INCFD < 1 )  go to 5002
  DO 1  J = 2, N
     if ( X(J) <= X(J-1) )  go to 5003
    1 CONTINUE
!
  IBEG = IC(1)
  IEND = IC(2)
  IERR = 0
  if ( (IBEG < 0).OR.(IBEG > 4) )  IERR = IERR - 1
  if ( (IEND < 0).OR.(IEND > 4) )  IERR = IERR - 2
  if ( IERR < 0 )  go to 5004
!
!  FUNCTION DEFINITION IS OK -- GO ON.
!
  if ( NWK  <  2*N )  go to 5007
!
!  COMPUTE FIRST DIFFERENCES OF X SEQUENCE AND STORE IN WK(1,.). ALSO,
!  COMPUTE FIRST DIVIDED DIFFERENCE OF DATA AND STORE IN WK(2,.).
  DO 5  J=2,N
     WK(1,J) = X(J) - X(J-1)
     WK(2,J) = (F(1,J) - F(1,J-1))/WK(1,J)
    5 CONTINUE
!
!  SET TO DEFAULT BOUNDARY CONDITIONS if N IS TOO SMALL.
!
  if ( IBEG > N )  IBEG = 0
  if ( IEND > N )  IEND = 0
!
!  SET UP FOR BOUNDARY CONDITIONS.
!
  if ( (IBEG == 1).OR.(IBEG == 2) )  THEN
     D(1,1) = VC(1)
  ELSE if (IBEG  >  2)  THEN
!        PICK UP FIRST IBEG POINTS, IN REVERSE ORDER.
     DO 10  J = 1, IBEG
        INDEX = IBEG-J+1
!           INDEX RUNS FROM IBEG DOWN TO 1.
        XTEMP(J) = X(INDEX)
        if (J  <  IBEG)  STEMP(J) = WK(2,INDEX)
   10    CONTINUE
!                 --------------------------------
     D(1,1) = DPCHDF (IBEG, XTEMP, STEMP, IERR)
!                 --------------------------------
     if (IERR  /=  0)  go to 5009
     IBEG = 1
  end if
!
  if ( (IEND == 1).OR.(IEND == 2) )  THEN
     D(1,N) = VC(2)
  ELSE if (IEND  >  2)  THEN
!        PICK UP LAST IEND POINTS.
     DO 15  J = 1, IEND
        INDEX = N-IEND+J
!           INDEX RUNS FROM N+1-IEND UP TO N.
        XTEMP(J) = X(INDEX)
        if (J  <  IEND)  STEMP(J) = WK(2,INDEX+1)
   15    CONTINUE
!                 --------------------------------
     D(1,N) = DPCHDF (IEND, XTEMP, STEMP, IERR)
!                 --------------------------------
     if (IERR  /=  0)  go to 5009
     IEND = 1
  end if
!
! --------------------( BEGIN CODING FROM CUBSPL )--------------------
!
!  **** A TRIDIAGONAL LINEAR SYSTEM FOR THE UNKNOWN SLOPES S(J) OF
!  F  AT X(J), J=1,...,N, IS GENERATED AND THEN SOLVED BY GAUSS ELIM-
!  INATION, WITH S(J) ENDING UP IN D(1,J), ALL J.
!     WK(1,.) AND WK(2,.) ARE USED FOR TEMPORARY STORAGE.
!
!  CONSTRUCT FIRST EQUATION FROM FIRST BOUNDARY CONDITION, OF THE FORM
!             WK(2,1)*S(1) + WK(1,1)*S(2) = D(1,1)
!
  if (IBEG  ==  0)  THEN
     if (N  ==  2)  THEN
!           NO CONDITION AT LEFT END AND N = 2.
        WK(2,1) = ONE
        WK(1,1) = ONE
        D(1,1) = TWO*WK(2,2)
     ELSE
!           NOT-A-KNOT CONDITION AT LEFT END AND N  >  2.
        WK(2,1) = WK(1,3)
        WK(1,1) = WK(1,2) + WK(1,3)
        D(1,1) =((WK(1,2) + TWO*WK(1,1))*WK(2,2)*WK(1,3) &
                          + WK(1,2)**2*WK(2,3)) / WK(1,1)
     ENDIF
  ELSE if (IBEG  ==  1)  THEN
!        SLOPE PRESCRIBED AT LEFT END.
     WK(2,1) = ONE
     WK(1,1) = ZERO
  ELSE
!        SECOND DERIVATIVE PRESCRIBED AT LEFT END.
     WK(2,1) = TWO
     WK(1,1) = ONE
     D(1,1) = THREE*WK(2,2) - HALF*WK(1,2)*D(1,1)
  end if
!
!  if THERE ARE INTERIOR KNOTS, GENERATE THE CORRESPONDING EQUATIONS AND
!  CARRY OUT THE FORWARD PASS OF GAUSS ELIMINATION, AFTER WHICH THE J-TH
!  EQUATION READS    WK(2,J)*S(J) + WK(1,J)*S(J+1) = D(1,J).
!
  NM1 = N-1
  if (NM1  >  1)  THEN
     DO 20 J=2,NM1
        if (WK(2,J-1)  ==  ZERO)  go to 5008
        G = -WK(1,J+1)/WK(2,J-1)
        D(1,J) = G*D(1,J-1) &
                    + THREE*(WK(1,J)*WK(2,J+1) + WK(1,J+1)*WK(2,J))
        WK(2,J) = G*WK(1,J-1) + TWO*(WK(1,J) + WK(1,J+1))
   20    CONTINUE
  end if
!
!  CONSTRUCT LAST EQUATION FROM SECOND BOUNDARY CONDITION, OF THE FORM
!           (-G*WK(2,N-1))*S(N-1) + WK(2,N)*S(N) = D(1,N)
!
!     if SLOPE IS PRESCRIBED AT RIGHT END, ONE CAN GO DIRECTLY TO BACK-
!     SUBSTITUTION, SINCE ARRAYS HAPPEN TO BE SET UP JUST RIGHT FOR IT
!     AT THIS POINT.
  if (IEND  ==  1)  go to 30
!
  if (IEND  ==  0)  THEN
     if (N == 2 .AND. IBEG == 0)  THEN
!           NOT-A-KNOT AT RIGHT ENDPOINT AND AT LEFT ENDPOINT AND N = 2.
        D(1,2) = WK(2,2)
        go to 30
     ELSE if ((N == 2) .OR. (N == 3 .AND. IBEG == 0))  THEN
!           EITHER (N=3 AND NOT-A-KNOT ALSO AT LEFT) OR (N=2 AND *NOT*
!           NOT-A-KNOT AT LEFT END POINT).
        D(1,N) = TWO*WK(2,N)
        WK(2,N) = ONE
        if (WK(2,N-1)  ==  ZERO)  go to 5008
        G = -ONE/WK(2,N-1)
     ELSE
!           NOT-A-KNOT AND N  >=  3, AND EITHER N > 3 OR  ALSO NOT-A-
!           KNOT AT LEFT END POINT.
        G = WK(1,N-1) + WK(1,N)
!           DO NOT NEED TO CHECK FOLLOWING DENOMINATORS (X-DIFFERENCES).
        D(1,N) = ((WK(1,N)+TWO*G)*WK(2,N)*WK(1,N-1) &
                    + WK(1,N)**2*(F(1,N-1)-F(1,N-2))/WK(1,N-1))/G
        if (WK(2,N-1)  ==  ZERO)  go to 5008
        G = -G/WK(2,N-1)
        WK(2,N) = WK(1,N-1)
     ENDIF
  ELSE
!        SECOND DERIVATIVE PRESCRIBED AT RIGHT ENDPOINT.
     D(1,N) = THREE*WK(2,N) + HALF*WK(1,N)*D(1,N)
     WK(2,N) = TWO
     if (WK(2,N-1)  ==  ZERO)  go to 5008
     G = -ONE/WK(2,N-1)
  end if
!
!  COMPLETE FORWARD PASS OF GAUSS ELIMINATION.
!
  WK(2,N) = G*WK(1,N-1) + WK(2,N)
  if (WK(2,N)  ==  ZERO)   go to 5008
  D(1,N) = (G*D(1,N-1) + D(1,N))/WK(2,N)
!
!  CARRY OUT BACK SUBSTITUTION
!
   30 CONTINUE
  DO 40 J=NM1,1,-1
     if (WK(2,J)  ==  ZERO)  go to 5008
     D(1,J) = (D(1,J) - WK(1,J)*D(1,J+1))/WK(2,J)
   40 CONTINUE
! --------------------(  END  CODING FROM CUBSPL )--------------------
!
!  NORMAL RETURN.
!
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     N < 2 RETURN.
  IERR = -1
  call XERMSG ('SLATEC', 'DPCHSP', &
     'NUMBER OF DATA POINTS LESS THAN TWO', IERR, 1)
  return
!
 5002 CONTINUE
!     INCFD < 1 RETURN.
  IERR = -2
  call XERMSG ('SLATEC', 'DPCHSP', 'INCREMENT LESS THAN ONE', IERR, &
     1)
  return
!
 5003 CONTINUE
!     X-ARRAY NOT STRICTLY INCREASING.
  IERR = -3
  call XERMSG ('SLATEC', 'DPCHSP', &
     'X-ARRAY NOT STRICTLY INCREASING', IERR, 1)
  return
!
 5004 CONTINUE
!     IC OUT OF RANGE RETURN.
  IERR = IERR - 3
  call XERMSG ('SLATEC', 'DPCHSP', 'IC OUT OF RANGE', IERR, 1)
  return
!
 5007 CONTINUE
!     NWK TOO SMALL RETURN.
  IERR = -7
  call XERMSG ('SLATEC', 'DPCHSP', 'WORK ARRAY TOO SMALL', IERR, 1)
  return
!
 5008 CONTINUE
!     SINGULAR SYSTEM.
!   *** THEORETICALLY, THIS CAN ONLY OCCUR if SUCCESSIVE X-VALUES   ***
!   *** ARE EQUAL, WHICH SHOULD ALREADY HAVE BEEN CAUGHT (IERR=-3). ***
  IERR = -8
  call XERMSG ('SLATEC', 'DPCHSP', 'SINGULAR LINEAR SYSTEM', IERR, &
     1)
  return
!
 5009 CONTINUE
!     ERROR RETURN FROM DPCHDF.
!   *** THIS CASE SHOULD NEVER OCCUR ***
  IERR = -9
  call XERMSG ('SLATEC', 'DPCHSP', 'ERROR RETURN FROM DPCHDF', &
     IERR, 1)
  return
end
