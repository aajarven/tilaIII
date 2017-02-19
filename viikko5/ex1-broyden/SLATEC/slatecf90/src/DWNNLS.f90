subroutine DWNNLS (W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, &
     IWORK, WORK)
!
!! DWNNLS solves a linearly constrained least squares problem with ...
!            equality constraints and nonnegativity constraints on
!            selected variables.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1A2A
!***TYPE      DOUBLE PRECISION (WNNLS-S, DWNNLS-D)
!***KEYWORDS  CONSTRAINED LEAST SQUARES, CURVE FITTING, DATA FITTING,
!             EQUALITY CONSTRAINTS, INEQUALITY CONSTRAINTS,
!             NONNEGATIVITY CONSTRAINTS, QUADRATIC PROGRAMMING
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     Abstract
!
!     This subprogram solves a linearly constrained least squares
!     problem.  Suppose there are given matrices E and A of
!     respective dimensions ME by N and MA by N, and vectors F
!     and B of respective lengths ME and MA.  This subroutine
!     solves the problem
!
!               EX = F, (equations to be exactly satisfied)
!
!               AX = B, (equations to be approximately satisfied,
!                        in the least squares sense)
!
!               subject to components L+1,...,N nonnegative
!
!     Any values ME >= 0, MA >= 0 and 0 <=  L  <= N are permitted.
!
!     The problem is reposed as problem DWNNLS
!
!               (WT*E)X = (WT*F)
!               (   A)    (   B), (least squares)
!               subject to components L+1,...,N nonnegative.
!
!     The subprogram chooses the heavy weight (or penalty parameter) WT.
!
!     The parameters for DWNNLS are
!
!     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     W(*,*),MDW,  The array W(*,*) is double subscripted with first
!     ME,MA,N,L    dimensioning parameter equal to MDW.  For this
!                  discussion let us call M = ME + MA.  Then MDW
!                  must satisfy MDW >= M.  The condition MDW < M
!                  is an error.
!
!                  The array W(*,*) contains the matrices and vectors
!
!                       (E  F)
!                       (A  B)
!
!                  in rows and columns 1,...,M and 1,...,N+1
!                  respectively.  Columns 1,...,L correspond to
!                  unconstrained variables X(1),...,X(L).  The
!                  remaining variables are constrained to be
!                  nonnegative. The condition L < 0 or L > N is
!                  an error.
!
!     PRGOPT(*)    This double precision array is the option vector.
!                  If the user is satisfied with the nominal
!                  subprogram features set
!
!                  PRGOPT(1)=1 (or PRGOPT(1)=1.0)
!
!                  Otherwise PRGOPT(*) is a linked list consisting of
!                  groups of data of the following form
!
!                  LINK
!                  KEY
!                  DATA SET
!
!                  The parameters LINK and KEY are each one word.
!                  The DATA SET can be comprised of several words.
!                  The number of items depends on the value of KEY.
!                  The value of LINK points to the first
!                  entry of the next group of data within
!                  PRGOPT(*).  The exception is when there are
!                  no more options to change.  In that
!                  case LINK=1 and the values KEY and DATA SET
!                  are not referenced. The general layout of
!                  PRGOPT(*) is as follows.
!
!               ...PRGOPT(1)=LINK1 (link to first entry of next group)
!               .  PRGOPT(2)=KEY1 (key to the option change)
!               .  PRGOPT(3)=DATA VALUE (data value for this change)
!               .       .
!               .       .
!               .       .
!               ...PRGOPT(LINK1)=LINK2 (link to the first entry of
!               .                       next group)
!               .  PRGOPT(LINK1+1)=KEY2 (key to the option change)
!               .  PRGOPT(LINK1+2)=DATA VALUE
!               ...     .
!               .       .
!               .       .
!               ...PRGOPT(LINK)=1 (no more options to change)
!
!                  Values of LINK that are nonpositive are errors.
!                  A value of LINK > NLINK=100000 is also an error.
!                  This helps prevent using invalid but positive
!                  values of LINK that will probably extend
!                  beyond the program limits of PRGOPT(*).
!                  Unrecognized values of KEY are ignored.  The
!                  order of the options is arbitrary and any number
!                  of options can be changed with the following
!                  restriction.  To prevent cycling in the
!                  processing of the option array a count of the
!                  number of options changed is maintained.
!                  Whenever this count exceeds NOPT=1000 an error
!                  message is printed and the subprogram returns.
!
!                  OPTIONS..
!
!                  KEY=6
!                         Scale the nonzero columns of the
!                  entire data matrix
!                  (E)
!                  (A)
!                  to have length one. The DATA SET for
!                  this option is a single value.  It must
!                  be nonzero if unit length column scaling is
!                  desired.
!
!                  KEY=7
!                         Scale columns of the entire data matrix
!                  (E)
!                  (A)
!                  with a user-provided diagonal matrix.
!                  The DATA SET for this option consists
!                  of the N diagonal scaling factors, one for
!                  each matrix column.
!
!                  KEY=8
!                         Change the rank determination tolerance from
!                  the nominal value of SQRT(SRELPR).  This quantity
!                  can be no smaller than SRELPR, The arithmetic-
!                  storage precision.  The quantity used
!                  here is internally restricted to be at
!                  least SRELPR.  The DATA SET for this option
!                  is the new tolerance.
!
!                  KEY=9
!                         Change the blow-up parameter from the
!                  nominal value of SQRT(SRELPR).  The reciprocal of
!                  this parameter is used in rejecting solution
!                  components as too large when a variable is
!                  first brought into the active set.  Too large
!                  means that the proposed component times the
!                  reciprocal of the parameter is not less than
!                  the ratio of the norms of the right-side
!                  vector and the data matrix.
!                  This parameter can be no smaller than SRELPR,
!                  the arithmetic-storage precision.
!
!                  For example, suppose we want to provide
!                  a diagonal matrix to scale the problem
!                  matrix and change the tolerance used for
!                  determining linear dependence of dropped col
!                  vectors.  For these options the dimensions of
!                  PRGOPT(*) must be at least N+6.  The FORTRAN
!                  statements defining these options would
!                  be as follows.
!
!                  PRGOPT(1)=N+3 (link to entry N+3 in PRGOPT(*))
!                  PRGOPT(2)=7 (user-provided scaling key)
!
!                  call DCOPY(N,D,1,PRGOPT(3),1) (copy the N
!                  scaling factors from a user array called D(*)
!                  into PRGOPT(3)-PRGOPT(N+2))
!
!                  PRGOPT(N+3)=N+6 (link to entry N+6 of PRGOPT(*))
!                  PRGOPT(N+4)=8 (linear dependence tolerance key)
!                  PRGOPT(N+5)=... (new value of the tolerance)
!
!                  PRGOPT(N+6)=1 (no more options to change)
!
!
!     IWORK(1),    The amounts of working storage actually allocated
!     IWORK(2)     for the working arrays WORK(*) and IWORK(*),
!                  respectively.  These quantities are compared with
!                  the actual amounts of storage needed for DWNNLS( ).
!                  Insufficient storage allocated for either WORK(*)
!                  or IWORK(*) is considered an error.  This feature
!                  was included in DWNNLS( ) because miscalculating
!                  the storage formulas for WORK(*) and IWORK(*)
!                  might very well lead to subtle and hard-to-find
!                  execution errors.
!
!                  The length of WORK(*) must be at least
!
!                  LW = ME+MA+5*N
!                  This test will not be made if IWORK(1) <= 0.
!
!                  The length of IWORK(*) must be at least
!
!                  LIW = ME+MA+N
!                  This test will not be made if IWORK(2) <= 0.
!
!     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     X(*)         An array dimensioned at least N, which will
!                  contain the N components of the solution vector
!                  on output.
!
!     RNORM        The residual norm of the solution.  The value of
!                  RNORM contains the residual vector length of the
!                  equality constraints and least squares equations.
!
!     MODE         The value of MODE indicates the success or failure
!                  of the subprogram.
!
!                  MODE = 0  Subprogram completed successfully.
!
!                       = 1  Max. number of iterations (equal to
!                            3*(N-L)) exceeded. Nearly all problems
!                            should complete in fewer than this
!                            number of iterations. An approximate
!                            solution and its corresponding residual
!                            vector length are in X(*) and RNORM.
!
!                       = 2  Usage error occurred.  The offending
!                            condition is noted with the error
!                            processing subprogram, XERMSG( ).
!
!     User-designated
!     Working arrays..
!
!     WORK(*)      A double precision working array of length at least
!                  M + 5*N.
!
!     IWORK(*)     An integer-valued working array of length at least
!                  M+N.
!
!***REFERENCES  K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Report SAND77-0552, Sandia
!                 Laboratories, June 1978.
!               K. H. Haskell and R. J. Hanson, Selected algorithms for
!                 the linearly constrained least squares problem - a
!                 users guide, Report SAND78-1290, Sandia Laboratories,
!                 August 1979.
!               K. H. Haskell and R. J. Hanson, An algorithm for
!                 linear least squares problems with equality and
!                 nonnegativity constraints, Mathematical Programming
!                 21 (1981), pp. 98-118.
!               R. J. Hanson and K. H. Haskell, Two algorithms for the
!                 linearly constrained least squares problem, ACM
!                 Transactions on Mathematical Software, September 1982.
!               C. L. Lawson and R. J. Hanson, Solving Least Squares
!                 Problems, Prentice-Hall, Inc., 1974.
!***ROUTINES CALLED  DWNLSM, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890618  Completely restructured and revised.  (WRB & RWC)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls, change Prologue
!           comments to agree with WNNLS.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DWNNLS
  INTEGER IWORK(*), L, L1, L2, L3, L4, L5, LIW, LW, MA, MDW, ME, &
       MODE, N
  DOUBLE PRECISION  PRGOPT(*), RNORM, W(MDW,*), WORK(*), X(*)
  CHARACTER*8 XERN1
!***FIRST EXECUTABLE STATEMENT  DWNNLS
  MODE = 0
  if (MA+ME <= 0 .OR. N <= 0) RETURN
!
  if (IWORK(1) > 0) THEN
     LW = ME + MA + 5*N
     if (IWORK(1) < LW) THEN
        WRITE (XERN1, '(I8)') LW
        call XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' // &
           'ALLOCATED FOR WORK(*), NEED LW = ' // XERN1, 2, 1)
        MODE = 2
        return
     ENDIF
  end if
!
  if (IWORK(2) > 0) THEN
     LIW = ME + MA + N
     if (IWORK(2) < LIW) THEN
        WRITE (XERN1, '(I8)') LIW
        call XERMSG ('SLATEC', 'DWNNLS', 'INSUFFICIENT STORAGE ' // &
           'ALLOCATED FOR IWORK(*), NEED LIW = ' // XERN1, 2, 1)
        MODE = 2
        return
     ENDIF
  end if
!
  if (MDW < ME+MA) THEN
     call XERMSG ('SLATEC', 'DWNNLS', &
        'THE VALUE MDW < ME+MA IS AN ERROR', 1, 1)
     MODE = 2
     return
  end if
!
  if (L < 0 .OR. L > N) THEN
     call XERMSG ('SLATEC', 'DWNNLS', &
        'L >= 0 .AND. L <= N IS REQUIRED', 2, 1)
     MODE = 2
     return
  end if
!
!     THE PURPOSE OF THIS SUBROUTINE IS TO BREAK UP THE ARRAYS
!     WORK(*) AND IWORK(*) INTO SEPARATE WORK ARRAYS
!     REQUIRED BY THE MAIN SUBROUTINE DWNLSM( ).
!
  L1 = N + 1
  L2 = L1 + N
  L3 = L2 + ME + MA
  L4 = L3 + N
  L5 = L4 + N
!
  call DWNLSM(W, MDW, ME, MA, N, L, PRGOPT, X, RNORM, MODE, IWORK, &
              IWORK(L1), WORK(1), WORK(L1), WORK(L2), WORK(L3), &
              WORK(L4), WORK(L5))
  return
end
