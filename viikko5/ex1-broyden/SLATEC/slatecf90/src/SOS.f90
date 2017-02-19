subroutine SOS (FNC, NEQ, X, RTOLX, ATOLX, TOLF, IFLAG, RW, LRW, &
     IW, LIW)
!
!! SOS solves a square system of nonlinear equations.
!
!***LIBRARY   SLATEC
!***CATEGORY  F2A
!***TYPE      SINGLE PRECISION (SOS-S, DSOS-D)
!***KEYWORDS  BROWN'S METHOD, NEWTON'S METHOD, NONLINEAR EQUATIONS,
!             ROOTS, SOLUTIONS
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     SOS solves a system of NEQ simultaneous nonlinear equations in
!     NEQ unknowns.  That is, it solves the problem   F(X)=0
!     where X is a vector with components  X(1),...,X(NEQ)  and  F
!     is a vector of nonlinear functions.  Each equation is of the form
!
!               F (X(1),...,X(NEQ))=0     for K=1,...,NEQ.
!                K
!
!     The algorithm is based on an iterative method which is a
!     variation of Newton's method using Gaussian elimination
!     in a manner similar to the Gauss-Seidel process.  Convergence
!     is roughly quadratic.  All partial derivatives required by
!     the algorithm are approximated by first difference quotients.
!     The convergence behavior of this code is affected by the
!     ordering of the equations, and it is advantageous to place linear
!     and mildly nonlinear equations first in the ordering.
!
!     Actually, SOS is merely an interfacing routine for
!     calling subroutine SOSEQS which embodies the solution
!     algorithm.  The purpose of this is to add greater
!     flexibility and ease of use for the prospective user.
!
!     SOSEQS calls the accompanying routine SOSSOL, which solves special
!     triangular linear systems by back-substitution.
!
!     The user must supply a function subprogram which evaluates the
!     K-th equation only (K specified by SOSEQS) for each call
!     to the subprogram.
!
!     SOS represents an implementation of the mathematical algorithm
!     described in the references below.  It is a modification of the
!     code SOSNLE written by H. A. Watts in 1973.
!
! **********************************************************************
!   -Input-
!
!     FNC -Name of the function program which evaluates the equations.
!          This name must be in an EXTERNAL statement in the calling
!          program.  The user must supply FNC in the form FNC(X,K),
!          where X is the solution vector (which must be dimensioned
!          in FNC) and FNC returns the value of the K-th function.
!
!     NEQ -Number of equations to be solved.
!
!     X   -Solution vector.  Initial guesses must be supplied.
!
!     RTOLX -Relative error tolerance used in the convergence criteria.
!          Each solution component X(I) is checked by an accuracy test
!          of the form   ABS(X(I)-XOLD(I))  <=  RTOLX*ABS(X(I))+ATOLX,
!          where XOLD(I) represents the previous iteration value.
!          RTOLX must be non-negative.
!
!     ATOLX -Absolute error tolerance used in the convergence criteria.
!          ATOLX must be non-negative.  If the user suspects some
!          solution component may be zero, he should set ATOLX to an
!          appropriate (depends on the scale of the remaining variables)
!          positive value for better efficiency.
!
!     TOLF -Residual error tolerance used in the convergence criteria.
!          Convergence will be indicated if all residuals (values of the
!          functions or equations) are not bigger than TOLF in
!          magnitude.  Note that extreme care must be given in assigning
!          an appropriate value for TOLF because this convergence test
!          is dependent on the scaling of the equations.  An
!          inappropriate value can cause premature termination of the
!          iteration process.
!
!     IFLAG -Optional input indicator.  You must set  IFLAG=-1  if you
!          want to use any of the optional input items listed below.
!          Otherwise set it to zero.
!
!     RW  -A REAL work array which is split apart by SOS and used
!          internally by SOSEQS.
!
!     LRW -Dimension of the RW array.  LRW must be at least
!                    1 + 6*NEQ + NEQ*(NEQ+1)/2
!
!     IW  -An INTEGER work array which is split apart by SOS and used
!          internally by SOSEQS.
!
!     LIW -Dimension of the IW array. LIW must be at least  3 + NEQ.
!
!   -Optional Input-
!
!     IW(1) -Internal printing parameter.  You must set  IW(1)=-1  if
!          you want the intermediate solution iterates to be printed.
!
!     IW(2) -Iteration limit.  The maximum number of allowable
!          iterations can be specified, if desired.  To override the
!          default value of 50, set IW(2) to the number wanted.
!
!     Remember, if you tell the code that you are using one of the
!               options (by setting IFLAG=-1), you must supply values
!               for both IW(1) and IW(2).
!
! **********************************************************************
!   -Output-
!
!     X   -Solution vector.
!
!     IFLAG -Status indicator
!
!                         *** Convergence to a Solution ***
!
!          1 Means satisfactory convergence to a solution was achieved.
!            Each solution component X(I) satisfies the error tolerance
!            test   ABS(X(I)-XOLD(I))  <=  RTOLX*ABS(X(I))+ATOLX.
!
!          2 Means procedure converged to a solution such that all
!            residuals are at most TOLF in magnitude,
!            ABS(FNC(X,I))  <=  TOLF.
!
!          3 Means that conditions for both IFLAG=1 and IFLAG=2 hold.
!
!          4 Means possible numerical convergence.  Behavior indicates
!            limiting precision calculations as a result of user asking
!            for too much accuracy or else convergence is very slow.
!            Residual norms and solution increment norms have
!            remained roughly constant over several consecutive
!            iterations.
!
!                         *** Task Interrupted ***
!
!          5 Means the allowable number of iterations has been met
!            without obtaining a solution to the specified accuracy.
!            Very slow convergence may be indicated.  Examine the
!            approximate solution returned and see if the error
!            tolerances seem appropriate.
!
!          6 Means the allowable number of iterations has been met and
!            the iterative process does not appear to be converging.
!            A local minimum may have been encountered or there may be
!            limiting precision difficulties.
!
!          7 Means that the iterative scheme appears to be diverging.
!            Residual norms and solution increment norms have
!            increased over several consecutive iterations.
!
!                         *** Task Cannot Be Continued ***
!
!          8 Means that a Jacobian-related matrix was singular.
!
!          9 Means improper input parameters.
!
!          *** IFLAG should be examined after each call to   ***
!          *** SOS with the appropriate action being taken.  ***
!
!
!     RW(1) -Contains a norm of the residual.
!
!     IW(3) -Contains the number of iterations used by the process.
!
! **********************************************************************
!***REFERENCES  K. M. Brown, Solution of simultaneous nonlinear
!                 equations, Algorithm 316, Communications of the
!                 A.C.M. 10, (1967), pp. 728-729.
!               K. M. Brown, A quadratically convergent Newton-like
!                 method based upon Gaussian elimination, SIAM Journal
!                 on Numerical Analysis 6, (1969), pp. 560-569.
!***ROUTINES CALLED  SOSEQS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Convert XERRWV calls to XERMSG calls, changed Prologue
!           comments to agree with DSOS.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SOS
  DIMENSION X(*), RW(*), IW(*)
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3, XERN4
  EXTERNAL FNC
!***FIRST EXECUTABLE STATEMENT  SOS
  INPFLG = IFLAG
!
!     CHECK FOR VALID INPUT
!
  if (NEQ  <=  0) THEN
     WRITE (XERN1, '(I8)') NEQ
     call XERMSG ('SLATEC', 'SOS', 'THE NUMBER OF EQUATIONS ' // &
        'MUST BE A POSITIVE INTEGER.  YOU HAVE CALLED THE ' // &
        'CODE WITH NEQ = ' // XERN1, 1, 1)
     IFLAG = 9
  end if
!
  if (RTOLX  <  0.0D0 .OR. ATOLX  <  0.0D0) THEN
     WRITE (XERN3, '(1PE15.6)') ATOLX
     WRITE (XERN4, '(1PE15.6)') RTOLX
     call XERMSG ('SLATEC', 'SOS', 'THE ERROR TOLERANCES FOR ' // &
        'THE SOLUTION ITERATES CANNOT BE NEGATIVE. YOU HAVE ' // &
        'CALLED THE CODE WITH  RTOLX = ' // XERN3 // &
        ' AND ATOLX = ' // XERN4,2, 1)
        IFLAG = 9
  end if
!
  if (TOLF  <  0.0D0) THEN
     WRITE (XERN3, '(1PE15.6)') TOLF
     call XERMSG ('SLATEC', 'SOS', 'THE RESIDUAL ERROR ' // &
        'TOLERANCE MUST BE NON-NEGATIVE.  YOU HAVE CALLED THE ' // &
        'CODE WITH TOLF = ' // XERN3, 3, 1)
        IFLAG = 9
  end if
!
  IPRINT = 0
  MXIT = 50
  if (INPFLG  ==  (-1)) THEN
     if (IW(1)  ==  (-1)) IPRINT = -1
     MXIT = IW(2)
     if (MXIT  <=  0) THEN
        WRITE (XERN1, '(I8)') MXIT
        call XERMSG ('SLATEC', 'SOS', 'YOU HAVE TOLD THE CODE ' // &
           'TO USE OPTIONAL IN PUT ITEMS BY SETTING  IFLAG=-1. ' // &
           'HOWEVER YOU HAVE CALLED THE CODE WITH THE MAXIMUM ' // &
           'ALLOWABLE NUMBER OF ITERATIONS SET TO  IW(2) = ' // &
           XERN1, 4, 1)
        IFLAG = 9
     ENDIF
  end if
!
  NC = (NEQ*(NEQ+1))/2
  if (LRW  <  1 + 6*NEQ + NC) THEN
     WRITE (XERN1, '(I8)') LRW
     call XERMSG ('SLATEC', 'SOS', 'DIMENSION OF THE RW ARRAY ' // &
        'MUST BE AT LEAST 1 + 6*NEQ + NEQ*(NEQ+1)/2 .  YOU HAVE ' // &
        'CALLED THE CODE WITH LRW = ' // XERN1, 5, 1)
     IFLAG = 9
  end if
!
  if (LIW  <  3 + NEQ) THEN
     WRITE (XERN1, '(I8)') LIW
     call XERMSG ('SLATEC', 'SOS', 'DIMENSION OF THE IW ARRAY ' // &
        'MUST BE AT LEAST   3 + NEQ.  YOU HAVE CALLED THE CODE ' // &
        'WITH  LIW = ' // XERN1, 6, 1)
     IFLAG = 9
  end if
!
  if (IFLAG  /=  9) THEN
     NCJS = 6
     NSRRC = 4
     NSRI = 5
!
     K1 = NC + 2
     K2 = K1 + NEQ
     K3 = K2 + NEQ
     K4 = K3 + NEQ
     K5 = K4 + NEQ
     K6 = K5 + NEQ
!
     call SOSEQS(FNC, NEQ, X, RTOLX, ATOLX, TOLF, IFLAG, MXIT, NCJS, &
                 NSRRC, NSRI, IPRINT, RW(1), RW(2), NC, RW(K1), &
                 RW(K2), RW(K3), RW(K4), RW(K5), RW(K6), IW(4))
!
     IW(3) = MXIT
  end if
  return
end
