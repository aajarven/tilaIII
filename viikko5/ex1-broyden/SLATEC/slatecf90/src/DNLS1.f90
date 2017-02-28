subroutine DNLS1 (FCN, IOPT, M, N, X, FVEC, FJAC, LDFJAC, FTOL, &
     XTOL, GTOL, MAXFEV, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, &
     NFEV, NJEV, IPVT, QTF, WA1, WA2, WA3, WA4)
!
!! DNLS1 minimizes the sum of the squares of M nonlinear functions ...
!  in N variables by a modification of the Levenberg-Marquardt algorithm.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1B1A1, K1B1A2
!***TYPE      DOUBLE PRECISION (SNLS1-S, DNLS1-D)
!***KEYWORDS  LEVENBERG-MARQUARDT, NONLINEAR DATA FITTING,
!             NONLINEAR LEAST SQUARES
!***AUTHOR  Hiebert, K. L., (SNLA)
!***DESCRIPTION
!
! 1. Purpose.
!
!       The purpose of DNLS1 is to minimize the sum of the squares of M
!       nonlinear functions in N variables by a modification of the
!       Levenberg-Marquardt algorithm.  The user must provide a subrou-
!       tine which calculates the functions.  The user has the option
!       of how the Jacobian will be supplied.  The user can supply the
!       full Jacobian, or the rows of the Jacobian (to avoid storing
!       the full Jacobian), or let the code approximate the Jacobian by
!       forward-differencing.   This code is the combination of the
!       MINPACK codes (Argonne) LMDER, LMDIF, and LMSTR.
!
!
! 2. Subroutine and Type Statements.
!
!       SUBROUTINE DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
!      *                 GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO
!      *                 ,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
!       INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
!       INTEGER IPVT(N)
!       DOUBLE PRECISION FTOL,XTOL,GTOL,EPSFCN,FACTOR
!       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),DIAG(N),QTF(N),
!      *     WA1(N),WA2(N),WA3(N),WA4(M)
!
!
! 3. Parameters.
!
!       Parameters designated as input parameters must be specified on
!       entry to DNLS1 and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from DNLS1.
!
!      FCN is the name of the user-supplied subroutine which calculate
!         the functions.  If the user wants to supply the Jacobian
!         (IOPT=2 or 3), then FCN must be written to calculate the
!         Jacobian, as well as the functions.  See the explanation
!         of the IOPT argument below.
!         If the user wants the iterates printed (NPRINT positive), then
!         FCN must do the printing.  See the explanation of NPRINT
!         below.  FCN must be declared in an EXTERNAL statement in the
!         calling program and should be written as follows.
!
!
!         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!         INTEGER IFLAG,LDFJAC,M,N
!         DOUBLE PRECISION X(N),FVEC(M)
!         ----------
!         FJAC and LDFJAC may be ignored       , if IOPT=1.
!         DOUBLE PRECISION FJAC(LDFJAC,N)      , if IOPT=2.
!         DOUBLE PRECISION FJAC(N)             , if IOPT=3.
!         ----------
!           If IFLAG=0, the values in X and FVEC are available
!           for printing.  See the explanation of NPRINT below.
!           IFLAG will never be zero unless NPRINT is positive.
!           The values of X and FVEC must not be changed.
!         return
!         ----------
!           If IFLAG=1, calculate the functions at X and return
!           this vector in FVEC.
!         return
!         ----------
!           If IFLAG=2, calculate the full Jacobian at X and return
!           this matrix in FJAC.  Note that IFLAG will never be 2 unless
!           IOPT=2.  FVEC contains the function values at X and must
!           not be altered.  FJAC(I,J) must be set to the derivative
!           of FVEC(I) with respect to X(J).
!         return
!         ----------
!           If IFLAG=3, calculate the LDFJAC-th row of the Jacobian
!           and return this vector in FJAC.  Note that IFLAG will
!           never be 3 unless IOPT=3.  FVEC contains the function
!           values at X and must not be altered.  FJAC(J) must be
!           set to the derivative of FVEC(LDFJAC) with respect to X(J).
!         return
!         ----------
!         END
!
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of DNLS1.  In this case, set
!         IFLAG to a negative integer.
!
!
!       IOPT is an input variable which specifies how the Jacobian will
!         be calculated.  If IOPT=2 or 3, then the user must supply the
!         Jacobian, as well as the function values, through the
!         subroutine FCN.  If IOPT=2, the user supplies the full
!         Jacobian with one call to FCN.  If IOPT=3, the user supplies
!         one row of the Jacobian with each call.  (In this manner,
!         storage can be saved because the full Jacobian is not stored.)
!         If IOPT=1, the code will approximate the Jacobian by forward
!         differencing.
!
!       M is a positive integer input variable set to the number of
!         functions.
!
!       N is a positive integer input variable set to the number of
!         variables.  N must not exceed M.
!
!       X is an array of length N.  On input, X must contain an initial
!         estimate of the solution vector.  On output, X contains the
!         final estimate of the solution vector.
!
!       FVEC is an output array of length M which contains the functions
!         evaluated at the output X.
!
!       FJAC is an output array.  For IOPT=1 and 2, FJAC is an M by N
!         array.  For IOPT=3, FJAC is an N by N array.  The upper N by N
!         submatrix of FJAC contains an upper triangular matrix R with
!         diagonal elements of nonincreasing magnitude such that
!
!                T     T           T
!               P *(JAC *JAC)*P = R *R,
!
!         where P is a permutation matrix and JAC is the final calcu-
!         lated Jacobian.  Column J of P is column IPVT(J) (see below)
!         of the identity matrix.  The lower part of FJAC contains
!         information generated during the computation of R.
!
!       LDFJAC is a positive integer input variable which specifies
!         the leading dimension of the array FJAC.  For IOPT=1 and 2,
!         LDFJAC must not be less than M.  For IOPT=3, LDFJAC must not
!         be less than N.
!
!       FTOL is a non-negative input variable.  Termination occurs when
!         both the actual and predicted relative reductions in the sum
!         of squares are at most FTOL.  Therefore, FTOL measures the
!         relative error desired in the sum of squares.  Section 4 con-
!         tains more details about FTOL.
!
!       XTOL is a non-negative input variable.  Termination occurs when
!         the relative error between two consecutive iterates is at most
!         XTOL.  Therefore, XTOL measures the relative error desired in
!         the approximate solution.  Section 4 contains more details
!         about XTOL.
!
!       GTOL is a non-negative input variable.  Termination occurs when
!         the cosine of the angle between FVEC and any column of the
!         Jacobian is at most GTOL in absolute value.  Therefore, GTOL
!         measures the orthogonality desired between the function vector
!         and the columns of the Jacobian.  Section 4 contains more
!         details about GTOL.
!
!       MAXFEV is a positive integer input variable.  Termination occurs
!         when the number of calls to FCN to evaluate the functions
!         has reached MAXFEV.
!
!       EPSFCN is an input variable used in determining a suitable step
!         for the forward-difference approximation.  This approximation
!         assumes that the relative errors in the functions are of the
!         order of EPSFCN.  If EPSFCN is less than the machine preci-
!         sion, it is assumed that the relative errors in the functions
!         are of the order of the machine precision.  If IOPT=2 or 3,
!         then EPSFCN can be ignored (treat it as a dummy argument).
!
!       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
!         internally set.  If MODE = 2, DIAG must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       MODE is an integer input variable.  If MODE = 1, the variables
!         will be scaled internally.  If MODE = 2, the scaling is speci-
!         fied by the input DIAG.  Other values of MODE are equivalent
!         to MODE = 1.
!
!       FACTOR is a positive input variable used in determining the ini-
!         tial step bound.  This bound is set to the product of FACTOR
!         and the Euclidean norm of DIAG*X if nonzero, or else to FACTOR
!         itself.  In most cases FACTOR should lie in the interval
!         (.1,100.).  100. is a generally recommended value.
!
!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive.  In this case, FCN is
!         called with IFLAG = 0 at the beginning of the first iteration
!         and every NPRINT iterations thereafter and immediately prior
!         to return, with X and FVEC available for printing. Appropriate
!         print statements must be added to FCN (see example) and
!         FVEC should not be altered.  If NPRINT is not positive, no
!         special calls to FCN with IFLAG = 0 are made.
!
!       INFO is an integer output variable.  If the user has terminated
!        execution, INFO is set to the (negative) value of IFLAG.  See
!        description of FCN and JAC. Otherwise, INFO is set as follows
!
!         INFO = 0  improper input parameters.
!
!         INFO = 1  both actual and predicted relative reductions in the
!                   sum of squares are at most FTOL.
!
!         INFO = 2  relative error between two consecutive iterates is
!                   at most XTOL.
!
!         INFO = 3  conditions for INFO = 1 and INFO = 2 both hold.
!
!         INFO = 4  the cosine of the angle between FVEC and any column
!                   of the Jacobian is at most GTOL in absolute value.
!
!         INFO = 5  number of calls to FCN for function evaluation
!                   has reached MAXFEV.
!
!         INFO = 6  FTOL is too small.  No further reduction in the sum
!                   of squares is possible.
!
!         INFO = 7  XTOL is too small.  No further improvement in the
!                   approximate solution X is possible.
!
!         INFO = 8  GTOL is too small.  FVEC is orthogonal to the
!                   columns of the Jacobian to machine precision.
!
!         Sections 4 and 5 contain more details about INFO.
!
!       NFEV is an integer output variable set to the number of calls to
!         FCN for function evaluation.
!
!       NJEV is an integer output variable set to the number of
!         evaluations of the full Jacobian.  If IOPT=2, only one call to
!         FCN is required for each evaluation of the full Jacobian.
!         If IOPT=3, the M calls to FCN are required.
!         If IOPT=1, then NJEV is set to zero.
!
!       IPVT is an integer output array of length N.  IPVT defines a
!         permutation matrix P such that JAC*P = Q*R, where JAC is the
!         final calculated Jacobian, Q is orthogonal (not stored), and R
!         is upper triangular with diagonal elements of nonincreasing
!         magnitude.  Column J of P is column IPVT(J) of the identity
!         matrix.
!
!       QTF is an output array of length N which contains the first N
!         elements of the vector (Q transpose)*FVEC.
!
!       WA1, WA2, and WA3 are work arrays of length N.
!
!       WA4 is a work array of length M.
!
!
! 4. Successful Completion.
!
!       The accuracy of DNLS1 is controlled by the convergence parame-
!       ters FTOL, XTOL, and GTOL.  These parameters are used in tests
!       which make three types of comparisons between the approximation
!       X and a solution XSOL.  DNLS1 terminates when any of the tests
!       is satisfied.  If any of the convergence parameters is less than
!       the machine precision (as defined by the function R1MACH(4)),
!       then DNLS1 only attempts to satisfy the test defined by the
!       machine precision.  Further progress is not usually possible.
!
!       The tests assume that the functions are reasonably well behaved,
!       and, if the Jacobian is supplied by the user, that the functions
!       and the Jacobian are coded consistently.  If these conditions
!       are not satisfied, then DNLS1 may incorrectly indicate conver-
!       gence.  If the Jacobian is coded correctly or IOPT=1,
!       then the validity of the answer can be checked, for example, by
!       rerunning DNLS1 with tighter tolerances.
!
!       First Convergence Test.  If ENORM(Z) denotes the Euclidean norm
!         of a vector Z, then this test attempts to guarantee that
!
!               ENORM(FVEC)  <=  (1+FTOL)*ENORM(FVECS),
!
!         where FVECS denotes the functions evaluated at XSOL.  If this
!         condition is satisfied with FTOL = 10**(-K), then the final
!         residual norm ENORM(FVEC) has K significant decimal digits and
!         INFO is set to 1 (or to 3 if the second test is also satis-
!         fied).  Unless high precision solutions are required, the
!         recommended value for FTOL is the square root of the machine
!         precision.
!
!       Second Convergence Test.  If D is the diagonal matrix whose
!         entries are defined by the array DIAG, then this test attempts
!         to guarantee that
!
!               ENORM(D*(X-XSOL))  <=  XTOL*ENORM(D*XSOL).
!
!         If this condition is satisfied with XTOL = 10**(-K), then the
!         larger components of D*X have K significant decimal digits and
!         INFO is set to 2 (or to 3 if the first test is also satis-
!         fied).  There is a danger that the smaller components of D*X
!         may have large relative errors, but if MODE = 1, then the
!         accuracy of the components of X is usually related to their
!         sensitivity.  Unless high precision solutions are required,
!         the recommended value for XTOL is the square root of the
!         machine precision.
!
!       Third Convergence Test.  This test is satisfied when the cosine
!         of the angle between FVEC and any column of the Jacobian at X
!         is at most GTOL in absolute value.  There is no clear rela-
!         tionship between this test and the accuracy of DNLS1, and
!         furthermore, the test is equally well satisfied at other crit-
!         ical points, namely maximizers and saddle points.  Therefore,
!         termination caused by this test (INFO = 4) should be examined
!         carefully.  The recommended value for GTOL is zero.
!
!
! 5. Unsuccessful Completion.
!
!       Unsuccessful termination of DNLS1 can be due to improper input
!       parameters, arithmetic interrupts, or an excessive number of
!       function evaluations.
!
!       Improper Input Parameters.  INFO is set to 0 if IOPT  <  1
!         or IOPT  >  3, or N  <=  0, or M  <  N, or for IOPT=1 or 2
!         LDFJAC  <  M, or for IOPT=3 LDFJAC  <  N, or FTOL  <  0.E0,
!         or XTOL  <  0.E0, or GTOL  <  0.E0, or MAXFEV  <=  0, or
!         FACTOR  <=  0.E0.
!
!       Arithmetic Interrupts.  If these interrupts occur in the FCN
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of X by DNLS1.  In this
!         case, it may be possible to remedy the situation by rerunning
!         DNLS1 with a smaller value of FACTOR.
!
!       Excessive Number of Function Evaluations.  A reasonable value
!         for MAXFEV is 100*(N+1) for IOPT=2 or 3 and 200*(N+1) for
!         IOPT=1.  If the number of calls to FCN reaches MAXFEV, then
!         this indicates that the routine is converging very slowly
!         as measured by the progress of FVEC, and INFO is set to 5.
!         In this case, it may be helpful to restart DNLS1 with MODE
!         set to 1.
!
!
! 6. Characteristics of the Algorithm.
!
!       DNLS1 is a modification of the Levenberg-Marquardt algorithm.
!       Two of its main characteristics involve the proper use of
!       implicitly scaled variables (if MODE = 1) and an optimal choice
!       for the correction.  The use of implicitly scaled variables
!       achieves scale invariance of DNLS1 and limits the size of the
!       correction in any direction where the functions are changing
!       rapidly.  The optimal choice of the correction guarantees (under
!       reasonable conditions) global convergence from starting points
!       far from the solution and a fast rate of convergence for
!       problems with small residuals.
!
!       Timing.  The time required by DNLS1 to solve a given problem
!         depends on M and N, the behavior of the functions, the accu-
!         racy requested, and the starting point.  The number of arith-
!         metic operations needed by DNLS1 is about N**3 to process each
!         evaluation of the functions (call to FCN) and to process each
!         evaluation of the Jacobian it takes M*N**2 for IOPT=2 (one
!         call to FCN), M*N**2 for IOPT=1 (N calls to FCN) and
!         1.5*M*N**2 for IOPT=3 (M calls to FCN).  Unless FCN
!         can be evaluated quickly, the timing of DNLS1 will be
!         strongly influenced by the time spent in FCN.
!
!       Storage.  DNLS1 requires (M*N + 2*M + 6*N) for IOPT=1 or 2 and
!         (N**2 + 2*M + 6*N) for IOPT=3 single precision storage
!         locations and N integer storage locations, in addition to
!         the storage required by the program.  There are no internally
!         declared storage arrays.
!
! *Long Description:
!
! 7. Example.
!
!       The problem is to determine the values of X(1), X(2), and X(3)
!       which provide the best fit (in the least squares sense) of
!
!             X(1) + U(I)/(V(I)*X(2) + W(I)*X(3)),  I = 1, 15
!
!       to the data
!
!             Y = (0.14,0.18,0.22,0.25,0.29,0.32,0.35,0.39,
!                  0.37,0.58,0.73,0.96,1.34,2.10,4.39),
!
!       where U(I) = I, V(I) = 16 - I, and W(I) = MIN(U(I),V(I)).  The
!       I-th component of FVEC is thus defined by
!
!             Y(I) - (X(1) + U(I)/(V(I)*X(2) + W(I)*X(3))).
!
!       **********
!
!       PROGRAM TEST
! C
! C     Driver for DNLS1 example.
! C
!       INTEGER J,IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV,
!      *        NWRITE
!       INTEGER IPVT(3)
!       DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,FNORM,EPSFCN
!       DOUBLE PRECISION X(3),FVEC(15),FJAC(15,3),DIAG(3),QTF(3),
!      *     WA1(3),WA2(3),WA3(3),WA4(15)
!       DOUBLE PRECISION DENORM,D1MACH
!       EXTERNAL FCN
!       DATA NWRITE /6/
! C
!       IOPT = 1
!       M = 15
!       N = 3
! C
! C     The following starting values provide a rough fit.
! C
!       X(1) = 1.E0
!       X(2) = 1.E0
!       X(3) = 1.E0
! C
!       LDFJAC = 15
! C
! C     Set FTOL and XTOL to the square root of the machine precision
! C     and GTOL to zero.  Unless high precision solutions are
! C     required, these are the recommended settings.
! C
!       FTOL = SQRT(R1MACH(4))
!       XTOL = SQRT(R1MACH(4))
!       GTOL = 0.E0
! C
!       MAXFEV = 400
!       EPSFCN = 0.0
!       MODE = 1
!       FACTOR = 1.E2
!       NPRINT = 0
! C
!       call DNLS1(FCN,IOPT,M,N,X,FVEC,FJAC,LDFJAC,FTOL,XTOL,
!      *           GTOL,MAXFEV,EPSFCN,DIAG,MODE,FACTOR,NPRINT,
!      *           INFO,NFEV,NJEV,IPVT,QTF,WA1,WA2,WA3,WA4)
!       FNORM = ENORM(M,FVEC)
!       WRITE (NWRITE,1000) FNORM,NFEV,NJEV,INFO,(X(J),J=1,N)
!       STOP
!  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
!      *        5X,' NUMBER OF JACOBIAN EVALUATIONS',I10 //
!      *        5X,' EXIT PARAMETER',16X,I10 //
!      *        5X,' FINAL APPROXIMATE SOLUTION' // 5X,3E15.7)
!       END
!       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,DUM,IDUM)
! C     This is the form of the FCN routine if IOPT=1,
! C     that is, if the user does not calculate the Jacobian.
!       INTEGER I,M,N,IFLAG
!       DOUBLE PRECISION X(N),FVEC(M),Y(15)
!       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
! C
!       if (IFLAG  /=  0) go to 5
! C
! C     Insert print statements here when NPRINT is positive.
! C
!       return
!     5 CONTINUE
!       DO 10 I = 1, M
!          TMP1 = I
!          TMP2 = 16 - I
!          TMP3 = TMP1
!          if (I  >  8) TMP3 = TMP2
!          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!    10    CONTINUE
!       return
!       END
!
!
!       Results obtained with different compilers or machines
!       may be slightly different.
!
!       FINAL L2 NORM OF THE RESIDUALS  0.9063596E-01
!
!       NUMBER OF FUNCTION EVALUATIONS        25
!
!       NUMBER OF JACOBIAN EVALUATIONS         0
!
!       EXIT PARAMETER                         1
!
!       FINAL APPROXIMATE SOLUTION
!
!        0.8241058E-01  0.1133037E+01  0.2343695E+01
!
!
!       For IOPT=2, FCN would be modified as follows to also
!       calculate the full Jacobian when IFLAG=2.
!
!       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
! C
! C     This is the form of the FCN routine if IOPT=2,
! C     that is, if the user calculates the full Jacobian.
! C
!       INTEGER I,LDFJAC,M,N,IFLAG
!       DOUBLE PRECISION X(N),FVEC(M),FJAC(LDFJAC,N),Y(15)
!       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
! C
!       if (IFLAG  /=  0) go to 5
! C
! C     Insert print statements here when NPRINT is positive.
! C
!       return
!     5 CONTINUE
!       if ( IFLAG /= 1) go to 20
!       DO 10 I = 1, M
!          TMP1 = I
!          TMP2 = 16 - I
!          TMP3 = TMP1
!          if (I  >  8) TMP3 = TMP2
!          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!    10    CONTINUE
!       return
! C
! C     Below, calculate the full Jacobian.
! C
!    20    CONTINUE
! C
!       DO 30 I = 1, M
!          TMP1 = I
!          TMP2 = 16 - I
!          TMP3 = TMP1
!          if (I  >  8) TMP3 = TMP2
!          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!          FJAC(I,1) = -1.E0
!          FJAC(I,2) = TMP1*TMP2/TMP4
!          FJAC(I,3) = TMP1*TMP3/TMP4
!    30    CONTINUE
!       return
!       END
!
!
!       For IOPT = 3, FJAC would be dimensioned as FJAC(3,3),
!         LDFJAC would be set to 3, and FCN would be written as
!         follows to calculate a row of the Jacobian when IFLAG=3.
!
!       SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
! C     This is the form of the FCN routine if IOPT=3,
! C     that is, if the user calculates the Jacobian row by row.
!       INTEGER I,M,N,IFLAG
!       DOUBLE PRECISION X(N),FVEC(M),FJAC(N),Y(15)
!       DOUBLE PRECISION TMP1,TMP2,TMP3,TMP4
!       DATA Y(1),Y(2),Y(3),Y(4),Y(5),Y(6),Y(7),Y(8),
!      *     Y(9),Y(10),Y(11),Y(12),Y(13),Y(14),Y(15)
!      *     /1.4E-1,1.8E-1,2.2E-1,2.5E-1,2.9E-1,3.2E-1,3.5E-1,3.9E-1,
!      *      3.7E-1,5.8E-1,7.3E-1,9.6E-1,1.34E0,2.1E0,4.39E0/
! C
!       if (IFLAG  /=  0) go to 5
! C
! C     Insert print statements here when NPRINT is positive.
! C
!       return
!     5 CONTINUE
!       if (  IFLAG /= 1) go to 20
!       DO 10 I = 1, M
!          TMP1 = I
!          TMP2 = 16 - I
!          TMP3 = TMP1
!          if (I  >  8) TMP3 = TMP2
!          FVEC(I) = Y(I) - (X(1) + TMP1/(X(2)*TMP2 + X(3)*TMP3))
!    10    CONTINUE
!       return
! C
! C     Below, calculate the LDFJAC-th row of the Jacobian.
! C
!    20 CONTINUE
!
!       I = LDFJAC
!          TMP1 = I
!          TMP2 = 16 - I
!          TMP3 = TMP1
!          if (I  >  8) TMP3 = TMP2
!          TMP4 = (X(2)*TMP2 + X(3)*TMP3)**2
!          FJAC(1) = -1.E0
!          FJAC(2) = TMP1*TMP2/TMP4
!          FJAC(3) = TMP1*TMP3/TMP4
!       return
!       END
!
!***REFERENCES  Jorge J. More, The Levenberg-Marquardt algorithm:
!                 implementation and theory.  In Numerical Analysis
!                 Proceedings (Dundee, June 28 - July 1, 1977, G. A.
!                 Watson, Editor), Lecture Notes in Mathematics 630,
!                 Springer-Verlag, 1978.
!***ROUTINES CALLED  D1MACH, DCKDER, DENORM, DFDJC3, DMPAR, DQRFAC,
!                    DWUPDT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   920205  Corrected XERN1 declaration.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNLS1
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER IOPT,M,N,LDFJAC,MAXFEV,MODE,NPRINT,INFO,NFEV,NJEV
  INTEGER IJUNK,NROW,IPVT(*)
  DOUBLE PRECISION FTOL,XTOL,GTOL,FACTOR,EPSFCN
  DOUBLE PRECISION X(*),FVEC(*),FJAC(LDFJAC,*),DIAG(*),QTF(*), &
       WA1(*),WA2(*),WA3(*),WA4(*)
  LOGICAL SING
  EXTERNAL FCN
  INTEGER I,IFLAG,ITER,J,L,MODECH
  DOUBLE PRECISION ACTRED,DELTA,DIRDER,EPSMCH,FNORM,FNORM1,GNORM, &
       ONE,PAR,PNORM,PRERED,P1,P5,P25,P75,P0001,RATIO,SUM,TEMP, &
       TEMP1,TEMP2,XNORM,ZERO
  DOUBLE PRECISION D1MACH,DENORM,ERR,CHKLIM
  CHARACTER*8 XERN1
  CHARACTER*16 XERN3
  SAVE CHKLIM, ONE, P1, P5, P25, P75, P0001, ZERO
!
  DATA CHKLIM/.1D0/
  DATA ONE,P1,P5,P25,P75,P0001,ZERO &
       /1.0D0,1.0D-1,5.0D-1,2.5D-1,7.5D-1,1.0D-4,0.0D0/
!***FIRST EXECUTABLE STATEMENT  DNLS1
  EPSMCH = D1MACH(4)
!
  INFO = 0
  IFLAG = 0
  NFEV = 0
  NJEV = 0
!
!     CHECK THE INPUT PARAMETERS FOR ERRORS.
!
  if (IOPT  <  1 .OR. IOPT  >  3 .OR. N  <=  0 .OR. &
      M  <  N .OR. LDFJAC  <  N .OR. FTOL  <  ZERO &
      .OR. XTOL  <  ZERO .OR. GTOL  <  ZERO &
      .OR. MAXFEV  <=  0 .OR. FACTOR  <=  ZERO) go to 300
  if (IOPT  <  3 .AND. LDFJAC  <  M) go to 300
  if (MODE  /=  2) go to 20
  DO 10 J = 1, N
     if (DIAG(J)  <=  ZERO) go to 300
   10    CONTINUE
   20 CONTINUE
!
!     EVALUATE THE FUNCTION AT THE STARTING POINT
!     AND CALCULATE ITS NORM.
!
  IFLAG = 1
  IJUNK = 1
  call FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
  NFEV = 1
  if (IFLAG  <  0) go to 300
  FNORM = DENORM(M,FVEC)
!
!     INITIALIZE LEVENBERG-MARQUARDT PARAMETER AND ITERATION COUNTER.
!
  PAR = ZERO
  ITER = 1
!
!     BEGINNING OF THE OUTER LOOP.
!
   30 CONTINUE
!
!        if REQUESTED, call FCN TO ENABLE PRINTING OF ITERATES.
!
     if (NPRINT  <=  0) go to 40
     IFLAG = 0
     if (MOD(ITER-1,NPRINT)  ==  0) &
        call FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
     if (IFLAG  <  0) go to 300
   40    CONTINUE
!
!        CALCULATE THE JACOBIAN MATRIX.
!
  if (IOPT  ==  3) go to 475
!
!     STORE THE FULL JACOBIAN USING M*N STORAGE
!
  if (IOPT  ==  1) go to 410
!
!     THE USER SUPPLIES THE JACOBIAN
!
     IFLAG = 2
     call FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
     NJEV = NJEV + 1
!
!             ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN
!
     if (ITER  <=  1) THEN
        if (IFLAG  <  0) go to 300
!
!           GET THE INCREMENTED X-VALUES INTO WA1(*).
!
        MODECH = 1
        call DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
!
!           EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT IN WA4(*).
!
        IFLAG = 1
        call FCN(IFLAG,M,N,WA1,WA4,FJAC,LDFJAC)
        NFEV = NFEV + 1
        if ( IFLAG  <  0) go to 300
        DO 350 I = 1, M
           MODECH = 2
           call DCKDER(1,N,X,FVEC(I),FJAC(I,1),LDFJAC,WA1, &
                WA4(I),MODECH,ERR)
           if (ERR  <  CHKLIM) THEN
              WRITE (XERN1, '(I8)') I
              WRITE (XERN3, '(1PE15.6)') ERR
              call XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF ' // &
                 'FUNCTION ' // XERN1 // ' MAY BE WRONG, ERR = ' // &
                 XERN3 // ' TOO CLOSE TO 0.', 7, 0)
           ENDIF
  350       CONTINUE
     ENDIF
!
     go to 420
!
!     THE CODE APPROXIMATES THE JACOBIAN
!
410      IFLAG = 1
     call DFDJC3(FCN,M,N,X,FVEC,FJAC,LDFJAC,IFLAG,EPSFCN,WA4)
     NFEV = NFEV + N
  420    if (IFLAG  <  0) go to 300
!
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
     call DQRFAC(M,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
!
!        FORM (Q TRANSPOSE)*FVEC AND STORE THE FIRST N COMPONENTS IN
!        QTF.
!
     DO 430 I = 1, M
        WA4(I) = FVEC(I)
  430         CONTINUE
     DO 470 J = 1, N
        if (FJAC(J,J)  ==  ZERO) go to 460
        SUM = ZERO
        DO 440 I = J, M
           SUM = SUM + FJAC(I,J)*WA4(I)
  440          CONTINUE
        TEMP = -SUM/FJAC(J,J)
        DO 450 I = J, M
           WA4(I) = WA4(I) + FJAC(I,J)*TEMP
  450          CONTINUE
  460       CONTINUE
        FJAC(J,J) = WA1(J)
        QTF(J) = WA4(J)
  470       CONTINUE
     go to 560
!
!        ACCUMULATE THE JACOBIAN BY ROWS IN ORDER TO SAVE STORAGE.
!        COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX
!        CALCULATED ONE ROW AT A TIME, WHILE SIMULTANEOUSLY
!        FORMING (Q TRANSPOSE)*FVEC AND STORING THE FIRST
!        N COMPONENTS IN QTF.
!
  475    DO 490 J = 1, N
        QTF(J) = ZERO
        DO 480 I = 1, N
           FJAC(I,J) = ZERO
  480          CONTINUE
  490        CONTINUE
     DO 500 I = 1, M
        NROW = I
        IFLAG = 3
        call FCN(IFLAG,M,N,X,FVEC,WA3,NROW)
        if (IFLAG  <  0) go to 300
!
!            ON THE FIRST ITERATION, CHECK THE USER SUPPLIED JACOBIAN.
!
        if ( ITER  >  1) go to 498
!
!            GET THE INCREMENTED X-VALUES INTO WA1(*).
!
        MODECH = 1
        call DCKDER(M,N,X,FVEC,FJAC,LDFJAC,WA1,WA4,MODECH,ERR)
!
!            EVALUATE AT INCREMENTED VALUES, if NOT ALREADY EVALUATED.
!
        if ( I  /=  1) go to 495
!
!            EVALUATE FUNCTION AT INCREMENTED VALUE AND PUT INTO WA4(*).
!
        IFLAG = 1
        call FCN(IFLAG,M,N,WA1,WA4,FJAC,NROW)
        NFEV = NFEV + 1
        if ( IFLAG  <  0) go to 300
495         CONTINUE
        MODECH = 2
        call DCKDER(1,N,X,FVEC(I),WA3,1,WA1,WA4(I),MODECH,ERR)
        if (ERR  <  CHKLIM) THEN
           WRITE (XERN1, '(I8)') I
           WRITE (XERN3, '(1PE15.6)') ERR
           call XERMSG ('SLATEC', 'DNLS1', 'DERIVATIVE OF FUNCTION ' &
              // XERN1 // ' MAY BE WRONG, ERR = ' // XERN3 // &
              ' TOO CLOSE TO 0.', 7, 0)
        ENDIF
498         CONTINUE
!
        TEMP = FVEC(I)
        call DWUPDT(N,FJAC,LDFJAC,WA3,QTF,TEMP,WA1,WA2)
  500       CONTINUE
     NJEV = NJEV + 1
!
!        if THE JACOBIAN IS RANK DEFICIENT, call DQRFAC TO
!        REORDER ITS COLUMNS AND UPDATE THE COMPONENTS OF QTF.
!
     SING = .FALSE.
     DO 510 J = 1, N
        if (FJAC(J,J)  ==  ZERO) SING = .TRUE.
        IPVT(J) = J
        WA2(J) = DENORM(J,FJAC(1,J))
  510       CONTINUE
     if (.NOT.SING) go to 560
     call DQRFAC(N,N,FJAC,LDFJAC,.TRUE.,IPVT,N,WA1,WA2,WA3)
     DO 550 J = 1, N
        if (FJAC(J,J)  ==  ZERO) go to 540
        SUM = ZERO
        DO 520 I = J, N
           SUM = SUM + FJAC(I,J)*QTF(I)
  520         CONTINUE
        TEMP = -SUM/FJAC(J,J)
        DO 530 I = J, N
           QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  530          CONTINUE
  540       CONTINUE
        FJAC(J,J) = WA1(J)
  550       CONTINUE
  560    CONTINUE
!
!        ON THE FIRST ITERATION AND if MODE IS 1, SCALE ACCORDING
!        TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
     if (ITER  /=  1) go to 80
     if (MODE  ==  2) go to 60
     DO 50 J = 1, N
        DIAG(J) = WA2(J)
        if (WA2(J)  ==  ZERO) DIAG(J) = ONE
   50       CONTINUE
   60    CONTINUE
!
!        ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED X
!        AND INITIALIZE THE STEP BOUND DELTA.
!
     DO 70 J = 1, N
        WA3(J) = DIAG(J)*X(J)
   70       CONTINUE
     XNORM = DENORM(N,WA3)
     DELTA = FACTOR*XNORM
     if (DELTA  ==  ZERO) DELTA = FACTOR
   80    CONTINUE
!
!        COMPUTE THE NORM OF THE SCALED GRADIENT.
!
     GNORM = ZERO
     if (FNORM  ==  ZERO) go to 170
     DO 160 J = 1, N
        L = IPVT(J)
        if (WA2(L)  ==  ZERO) go to 150
        SUM = ZERO
        DO 140 I = 1, J
           SUM = SUM + FJAC(I,J)*(QTF(I)/FNORM)
  140          CONTINUE
        GNORM = MAX(GNORM,ABS(SUM/WA2(L)))
  150       CONTINUE
  160       CONTINUE
  170    CONTINUE
!
!        TEST FOR CONVERGENCE OF THE GRADIENT NORM.
!
     if (GNORM  <=  GTOL) INFO = 4
     if (INFO  /=  0) go to 300
!
!        RESCALE if NECESSARY.
!
     if (MODE  ==  2) go to 190
     DO 180 J = 1, N
        DIAG(J) = MAX(DIAG(J),WA2(J))
  180       CONTINUE
  190    CONTINUE
!
!        BEGINNING OF THE INNER LOOP.
!
  200    CONTINUE
!
!           DETERMINE THE LEVENBERG-MARQUARDT PARAMETER.
!
        call DMPAR(N,FJAC,LDFJAC,IPVT,DIAG,QTF,DELTA,PAR,WA1,WA2, &
                   WA3,WA4)
!
!           STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
        DO 210 J = 1, N
           WA1(J) = -WA1(J)
           WA2(J) = X(J) + WA1(J)
           WA3(J) = DIAG(J)*WA1(J)
  210          CONTINUE
        PNORM = DENORM(N,WA3)
!
!           ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
        if (ITER  ==  1) DELTA = MIN(DELTA,PNORM)
!
!           EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!
        IFLAG = 1
        call FCN(IFLAG,M,N,WA2,WA4,FJAC,IJUNK)
        NFEV = NFEV + 1
        if (IFLAG  <  0) go to 300
        FNORM1 = DENORM(M,WA4)
!
!           COMPUTE THE SCALED ACTUAL REDUCTION.
!
        ACTRED = -ONE
        if (P1*FNORM1  <  FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!
!           COMPUTE THE SCALED PREDICTED REDUCTION AND
!           THE SCALED DIRECTIONAL DERIVATIVE.
!
        DO 230 J = 1, N
           WA3(J) = ZERO
           L = IPVT(J)
           TEMP = WA1(L)
           DO 220 I = 1, J
              WA3(I) = WA3(I) + FJAC(I,J)*TEMP
  220             CONTINUE
  230          CONTINUE
        TEMP1 = DENORM(N,WA3)/FNORM
        TEMP2 = (SQRT(PAR)*PNORM)/FNORM
        PRERED = TEMP1**2 + TEMP2**2/P5
        DIRDER = -(TEMP1**2 + TEMP2**2)
!
!           COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!           REDUCTION.
!
        RATIO = ZERO
        if (PRERED  /=  ZERO) RATIO = ACTRED/PRERED
!
!           UPDATE THE STEP BOUND.
!
        if (RATIO  >  P25) go to 240
           if (ACTRED  >=  ZERO) TEMP = P5
           if (ACTRED  <  ZERO) &
              TEMP = P5*DIRDER/(DIRDER + P5*ACTRED)
           if (P1*FNORM1  >=  FNORM .OR. TEMP  <  P1) TEMP = P1
           DELTA = TEMP*MIN(DELTA,PNORM/P1)
           PAR = PAR/TEMP
           go to 260
  240       CONTINUE
           if (PAR  /=  ZERO .AND. RATIO  <  P75) go to 250
           DELTA = PNORM/P5
           PAR = P5*PAR
  250          CONTINUE
  260       CONTINUE
!
!           TEST FOR SUCCESSFUL ITERATION.
!
        if (RATIO  <  P0001) go to 290
!
!           SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!
        DO 270 J = 1, N
           X(J) = WA2(J)
           WA2(J) = DIAG(J)*X(J)
  270          CONTINUE
        DO 280 I = 1, M
           FVEC(I) = WA4(I)
  280          CONTINUE
        XNORM = DENORM(N,WA2)
        FNORM = FNORM1
        ITER = ITER + 1
  290       CONTINUE
!
!           TESTS FOR CONVERGENCE.
!
        if (ABS(ACTRED)  <=  FTOL .AND. PRERED  <=  FTOL &
            .AND. P5*RATIO  <=  ONE) INFO = 1
        if (DELTA  <=  XTOL*XNORM) INFO = 2
        if (ABS(ACTRED)  <=  FTOL .AND. PRERED  <=  FTOL &
            .AND. P5*RATIO  <=  ONE .AND. INFO  ==  2) INFO = 3
        if (INFO  /=  0) go to 300
!
!           TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!
        if (NFEV  >=  MAXFEV) INFO = 5
        if (ABS(ACTRED)  <=  EPSMCH .AND. PRERED  <=  EPSMCH &
            .AND. P5*RATIO  <=  ONE) INFO = 6
        if (DELTA  <=  EPSMCH*XNORM) INFO = 7
        if (GNORM  <=  EPSMCH) INFO = 8
        if (INFO  /=  0) go to 300
!
!           END OF THE INNER LOOP. REPEAT if ITERATION UNSUCCESSFUL.
!
        if (RATIO  <  P0001) go to 200
!
!        END OF THE OUTER LOOP.
!
     go to 30
  300 CONTINUE
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
  if (IFLAG  <  0) INFO = IFLAG
  IFLAG = 0
  if (NPRINT  >  0) call FCN(IFLAG,M,N,X,FVEC,FJAC,IJUNK)
  if (INFO  <  0) call XERMSG ('SLATEC', 'DNLS1', &
     'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
  if (INFO  ==  0) call XERMSG ('SLATEC', 'DNLS1', &
     'INVALID INPUT PARAMETER.', 2, 1)
  if (INFO  ==  4) call XERMSG ('SLATEC', 'DNLS1', &
     'THIRD CONVERGENCE CONDITION, CHECK RESULTS BEFORE ACCEPTING.', &
     1, 1)
  if (INFO  ==  5) call XERMSG ('SLATEC', 'DNLS1', &
     'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
  if (INFO  >=  6) call XERMSG ('SLATEC', 'DNLS1', &
     'TOLERANCES TOO SMALL, NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
  return
end