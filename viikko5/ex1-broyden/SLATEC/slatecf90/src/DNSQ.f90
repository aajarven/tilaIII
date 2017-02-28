subroutine DNSQ (FCN, JAC, IOPT, N, X, FVEC, FJAC, LDFJAC, XTOL, &
     MAXFEV, ML, MU, EPSFCN, DIAG, MODE, FACTOR, NPRINT, INFO, NFEV, &
     NJEV, R, LR, QTF, WA1, WA2, WA3, WA4)
!
!! DNSQ finds a zero of a system of a N nonlinear functions in N variables ...
!  by a modification of the Powell hybrid method.
!
!***LIBRARY   SLATEC
!***CATEGORY  F2A
!***TYPE      DOUBLE PRECISION (SNSQ-S, DNSQ-D)
!***KEYWORDS  NONLINEAR SQUARE SYSTEM, POWELL HYBRID METHOD, ZEROS
!***AUTHOR  Hiebert, K. L. (SNLA)
!***DESCRIPTION
!
! 1. Purpose.
!
!       The purpose of DNSQ is to find a zero of a system of N nonlinear
!       functions in N variables by a modification of the Powell
!       hybrid method.  The user must provide a subroutine which
!       calculates the functions.  The user has the option of either to
!       provide a subroutine which calculates the Jacobian or to let the
!       code calculate it by a forward-difference approximation.
!       This code is the combination of the MINPACK codes (Argonne)
!       HYBRD and HYBRDJ.
!
! 2. Subroutine and Type Statements.
!
!       SUBROUTINE DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,
!      *                 ML,MU,EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,
!      *                 NJEV,R,LR,QTF,WA1,WA2,WA3,WA4)
!       INTEGER IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,NJEV,LR
!       DOUBLE PRECISION XTOL,EPSFCN,FACTOR
!       DOUBLE PRECISION
!       X(N),FVEC(N),DIAG(N),FJAC(LDFJAC,N),R(LR),QTF(N),
!      *     WA1(N),WA2(N),WA3(N),WA4(N)
!       EXTERNAL FCN,JAC
!
! 3. Parameters.
!
!       Parameters designated as input parameters must be specified on
!       entry to DNSQ and are not changed on exit, while parameters
!       designated as output parameters need not be specified on entry
!       and are set to appropriate values on exit from DNSQ.
!
!       FCN is the name of the user-supplied subroutine which calculates
!         the functions.  FCN must be declared in an EXTERNAL statement
!         in the user calling program, and should be written as follows.
!
!         SUBROUTINE FCN(N,X,FVEC,IFLAG)
!         INTEGER N,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N)
!         ----------
!         CALCULATE THE FUNCTIONS AT X AND
!         return THIS VECTOR IN FVEC.
!         ----------
!         return
!         END
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of DNSQ.  In this case set
!         IFLAG to a negative integer.
!
!       JAC is the name of the user-supplied subroutine which calculates
!         the Jacobian.  If IOPT=1, then JAC must be declared in an
!         EXTERNAL statement in the user calling program, and should be
!         written as follows.
!
!         SUBROUTINE JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
!         INTEGER N,LDFJAC,IFLAG
!         DOUBLE PRECISION X(N),FVEC(N),FJAC(LDFJAC,N)
!         ----------
!         Calculate the Jacobian at X and return this
!         matrix in FJAC.  FVEC contains the function
!         values at X and should not be altered.
!         ----------
!         return
!         END
!
!         The value of IFLAG should not be changed by JAC unless the
!         user wants to terminate execution of DNSQ.  In this case set
!         IFLAG to a negative integer.
!
!         If IOPT=2, JAC can be ignored (treat it as a dummy argument).
!
!       IOPT is an input variable which specifies how the Jacobian will
!         be calculated.  If IOPT=1, then the user must supply the
!         Jacobian through the subroutine JAC.  If IOPT=2, then the
!         code will approximate the Jacobian by forward-differencing.
!
!       N is a positive integer input variable set to the number of
!         functions and variables.
!
!       X is an array of length N.  On input X must contain an initial
!         estimate of the solution vector.  On output X contains the
!         final estimate of the solution vector.
!
!       FVEC is an output array of length N which contains the functions
!         evaluated at the output X.
!
!       FJAC is an output N by N array which contains the orthogonal
!         matrix Q produced by the QR factorization of the final
!         approximate Jacobian.
!
!       LDFJAC is a positive integer input variable not less than N
!         which specifies the leading dimension of the array FJAC.
!
!       XTOL is a nonnegative input variable.  Termination occurs when
!         the relative error between two consecutive iterates is at most
!         XTOL.  Therefore, XTOL measures the relative error desired in
!         the approximate solution.  Section 4 contains more details
!         about XTOL.
!
!       MAXFEV is a positive integer input variable.  Termination occurs
!         when the number of calls to FCN is at least MAXFEV by the end
!         of an iteration.
!
!       ML is a nonnegative integer input variable which specifies the
!         number of subdiagonals within the band of the Jacobian matrix.
!         If the Jacobian is not banded or IOPT=1, set ML to at
!         least N - 1.
!
!       MU is a nonnegative integer input variable which specifies the
!         number of superdiagonals within the band of the Jacobian
!         matrix.  If the Jacobian is not banded or IOPT=1, set MU to at
!         least N - 1.
!
!       EPSFCN is an input variable used in determining a suitable step
!         for the forward-difference approximation.  This approximation
!         assumes that the relative errors in the functions are of the
!         order of EPSFCN.  If EPSFCN is less than the machine
!         precision, it is assumed that the relative errors in the
!         functions are of the order of the machine precision.  If
!         IOPT=1, then EPSFCN can be ignored (treat it as a dummy
!         argument).
!
!       DIAG is an array of length N.  If MODE = 1 (see below), DIAG is
!         internally set.  If MODE = 2, DIAG must contain positive
!         entries that serve as implicit (multiplicative) scale factors
!         for the variables.
!
!       MODE is an integer input variable.  If MODE = 1, the variables
!         will be scaled internally.  If MODE = 2, the scaling is
!         specified by the input DIAG.  Other values of MODE are
!         equivalent to MODE = 1.
!
!       FACTOR is a positive input variable used in determining the
!         initial step bound.  This bound is set to the product of
!         FACTOR and the Euclidean norm of DIAG*X if nonzero, or else to
!         FACTOR itself.  In most cases FACTOR should lie in the
!         interval (.1,100.).  100. is a generally recommended value.
!
!       NPRINT is an integer input variable that enables controlled
!         printing of iterates if it is positive.  In this case, FCN is
!         called with IFLAG = 0 at the beginning of the first iteration
!         and every NPRINT iterations thereafter and immediately prior
!         to return, with X and FVEC available for printing. appropriate
!         print statements must be added to FCN(see example).  If NPRINT
!         is not positive, no special calls of FCN with IFLAG = 0 are
!         made.
!
!       INFO is an integer output variable.  If the user has terminated
!         execution, INFO is set to the (negative) value of IFLAG.  See
!         description of FCN and JAC. Otherwise, INFO is set as follows.
!
!         INFO = 0  Improper input parameters.
!
!         INFO = 1  Relative error between two consecutive iterates is
!                   at most XTOL.
!
!         INFO = 2  Number of calls to FCN has reached or exceeded
!                   MAXFEV.
!
!         INFO = 3  XTOL is too small.  No further improvement in the
!                   approximate solution X is possible.
!
!         INFO = 4  Iteration is not making good progress, as measured
!                   by the improvement from the last five Jacobian
!                   evaluations.
!
!         INFO = 5  Iteration is not making good progress, as measured
!                   by the improvement from the last ten iterations.
!
!         Sections 4 and 5 contain more details about INFO.
!
!       NFEV is an integer output variable set to the number of calls to
!         FCN.
!
!       NJEV is an integer output variable set to the number of calls to
!         JAC. (If IOPT=2, then NJEV is set to zero.)
!
!       R is an output array of length LR which contains the upper
!         triangular matrix produced by the QR factorization of the
!         final approximate Jacobian, stored rowwise.
!
!       LR is a positive integer input variable not less than
!         (N*(N+1))/2.
!
!       QTF is an output array of length N which contains the vector
!         (Q transpose)*FVEC.
!
!       WA1, WA2, WA3, and WA4 are work arrays of length N.
!
!
! 4. Successful completion.
!
!       The accuracy of DNSQ is controlled by the convergence parameter
!       XTOL.  This parameter is used in a test which makes a comparison
!       between the approximation X and a solution XSOL.  DNSQ
!       terminates when the test is satisfied.  If the convergence
!       parameter is less than the machine precision (as defined by the
!       function D1MACH(4)), then DNSQ only attempts to satisfy the test
!       defined by the machine precision.  Further progress is not
!       usually possible.
!
!       The test assumes that the functions are reasonably well behaved,
!       and, if the Jacobian is supplied by the user, that the functions
!       and the Jacobian are coded consistently.  If these conditions
!       are not satisfied, then DNSQ may incorrectly indicate
!       convergence.  The coding of the Jacobian can be checked by the
!       subroutine DCKDER. If the Jacobian is coded correctly or IOPT=2,
!       then the validity of the answer can be checked, for example, by
!       rerunning DNSQ with a tighter tolerance.
!
!       Convergence Test.  If DENORM(Z) denotes the Euclidean norm of a
!         vector Z and D is the diagonal matrix whose entries are
!         defined by the array DIAG, then this test attempts to
!         guarantee that
!
!               DENORM(D*(X-XSOL))  <=  XTOL*DENORM(D*XSOL).
!
!         If this condition is satisfied with XTOL = 10**(-K), then the
!         larger components of D*X have K significant decimal digits and
!         INFO is set to 1.  There is a danger that the smaller
!         components of D*X may have large relative errors, but the fast
!         rate of convergence of DNSQ usually avoids this possibility.
!         Unless high precision solutions are required, the recommended
!         value for XTOL is the square root of the machine precision.
!
!
! 5. Unsuccessful Completion.
!
!       Unsuccessful termination of DNSQ can be due to improper input
!       parameters, arithmetic interrupts, an excessive number of
!       function evaluations, or lack of good progress.
!
!       Improper Input Parameters.  INFO is set to 0 if IOPT .LT .1,
!         or IOPT  >  2, or N  <=  0, or LDFJAC  <  N, or
!         XTOL  <  0.E0, or MAXFEV  <=  0, or ML  <  0, or MU  <  0,
!         or FACTOR  <=  0.E0, or LR  <  (N*(N+1))/2.
!
!       Arithmetic Interrupts.  If these interrupts occur in the FCN
!         subroutine during an early stage of the computation, they may
!         be caused by an unacceptable choice of X by DNSQ.  In this
!         case, it may be possible to remedy the situation by rerunning
!         DNSQ with a smaller value of FACTOR.
!
!       Excessive Number of Function Evaluations.  A reasonable value
!         for MAXFEV is 100*(N+1) for IOPT=1 and 200*(N+1) for IOPT=2.
!         If the number of calls to FCN reaches MAXFEV, then this
!         indicates that the routine is converging very slowly as
!         measured by the progress of FVEC, and INFO is set to 2. This
!         situation should be unusual because, as indicated below, lack
!         of good progress is usually diagnosed earlier by DNSQ,
!         causing termination with info = 4 or INFO = 5.
!
!       Lack of Good Progress.  DNSQ searches for a zero of the system
!         by minimizing the sum of the squares of the functions.  In so
!         doing, it can become trapped in a region where the minimum
!         does not correspond to a zero of the system and, in this
!         situation, the iteration eventually fails to make good
!         progress.  In particular, this will happen if the system does
!         not have a zero.  If the system has a zero, rerunning DNSQ
!         from a different starting point may be helpful.
!
!
! 6. Characteristics of The Algorithm.
!
!       DNSQ is a modification of the Powell Hybrid method.  Two of its
!       main characteristics involve the choice of the correction as a
!       convex combination of the Newton and scaled gradient directions,
!       and the updating of the Jacobian by the rank-1 method of
!       Broyden.  The choice of the correction guarantees (under
!       reasonable conditions) global convergence for starting points
!       far from the solution and a fast rate of convergence.  The
!       Jacobian is calculated at the starting point by either the
!       user-supplied subroutine or a forward-difference approximation,
!       but it is not recalculated until the rank-1 method fails to
!       produce satisfactory progress.
!
!       Timing.  The time required by DNSQ to solve a given problem
!         depends on N, the behavior of the functions, the accuracy
!         requested, and the starting point.  The number of arithmetic
!         operations needed by DNSQ is about 11.5*(N**2) to process
!         each evaluation of the functions (call to FCN) and 1.3*(N**3)
!         to process each evaluation of the Jacobian (call to JAC,
!         if IOPT = 1).  Unless FCN and JAC can be evaluated quickly,
!         the timing of DNSQ will be strongly influenced by the time
!         spent in FCN and JAC.
!
!       Storage.  DNSQ requires (3*N**2 + 17*N)/2 single precision
!         storage locations, in addition to the storage required by the
!         program.  There are no internally declared storage arrays.
!
! *Long Description:
!
! 7. Example.
!
!       The problem is to determine the values of X(1), X(2), ..., X(9),
!       which solve the system of tridiagonal equations
!
!       (3-2*X(1))*X(1)           -2*X(2)                   = -1
!               -X(I-1) + (3-2*X(I))*X(I)         -2*X(I+1) = -1, I=2-8
!                                   -X(8) + (3-2*X(9))*X(9) = -1
! C     **********
!
!       PROGRAM TEST
! C
! C     Driver for DNSQ example.
! C
!       INTEGER J,IOPT,N,MAXFEV,ML,MU,MODE,NPRINT,INFO,NFEV,LDFJAC,LR,
!      *        NWRITE
!       DOUBLE PRECISION XTOL,EPSFCN,FACTOR,FNORM
!       DOUBLE PRECISION X(9),FVEC(9),DIAG(9),FJAC(9,9),R(45),QTF(9),
!      *     WA1(9),WA2(9),WA3(9),WA4(9)
!       DOUBLE PRECISION DENORM,D1MACH
!       EXTERNAL FCN
!       DATA NWRITE /6/
! C
!       IOPT = 2
!       N = 9
! C
! C     THE FOLLOWING STARTING VALUES PROVIDE A ROUGH SOLUTION.
! C
!       DO 10 J = 1, 9
!          X(J) = -1.E0
!    10    CONTINUE
! C
!       LDFJAC = 9
!       LR = 45
! C
! C     SET XTOL TO THE SQUARE ROOT OF THE MACHINE PRECISION.
! C     UNLESS HIGH PRECISION SOLUTIONS ARE REQUIRED,
! C     THIS IS THE RECOMMENDED SETTING.
! C
!       XTOL = SQRT(D1MACH(4))
! C
!       MAXFEV = 2000
!       ML = 1
!       MU = 1
!       EPSFCN = 0.E0
!       MODE = 2
!       DO 20 J = 1, 9
!          DIAG(J) = 1.E0
!    20    CONTINUE
!       FACTOR = 1.E2
!       NPRINT = 0
! C
!       call DNSQ(FCN,JAC,IOPT,N,X,FVEC,FJAC,LDFJAC,XTOL,MAXFEV,ML,MU,
!      *           EPSFCN,DIAG,MODE,FACTOR,NPRINT,INFO,NFEV,NJEV,
!      *           R,LR,QTF,WA1,WA2,WA3,WA4)
!       FNORM = DENORM(N,FVEC)
!       WRITE (NWRITE,1000) FNORM,NFEV,INFO,(X(J),J=1,N)
!       STOP
!  1000 FORMAT (5X,' FINAL L2 NORM OF THE RESIDUALS',E15.7 //
!      *        5X,' NUMBER OF FUNCTION EVALUATIONS',I10 //
!      *        5X,' EXIT PARAMETER',16X,I10 //
!      *        5X,' FINAL APPROXIMATE SOLUTION' // (5X,3E15.7))
!       END
!       SUBROUTINE FCN(N,X,FVEC,IFLAG)
!       INTEGER N,IFLAG
!       DOUBLE PRECISION X(N),FVEC(N)
!       INTEGER K
!       DOUBLE PRECISION ONE,TEMP,TEMP1,TEMP2,THREE,TWO,ZERO
!       DATA ZERO,ONE,TWO,THREE /0.E0,1.E0,2.E0,3.E0/
! C
!       if (IFLAG  /=  0) go to 5
! C
! C     INSERT PRINT STATEMENTS HERE WHEN NPRINT IS POSITIVE.
! C
!       return
!     5 CONTINUE
!       DO 10 K = 1, N
!          TEMP = (THREE - TWO*X(K))*X(K)
!          TEMP1 = ZERO
!          if (K  /=  1) TEMP1 = X(K-1)
!          TEMP2 = ZERO
!          if (K  /=  N) TEMP2 = X(K+1)
!          FVEC(K) = TEMP - TEMP1 - TWO*TEMP2 + ONE
!    10    CONTINUE
!       return
!       END
!
!       Results obtained with different compilers or machines
!       may be slightly different.
!
!       Final L2 norm of the residuals  0.1192636E-07
!
!       Number of function evaluations        14
!
!       Exit parameter                         1
!
!       Final approximate solution
!
!       -0.5706545E+00 -0.6816283E+00 -0.7017325E+00
!       -0.7042129E+00 -0.7013690E+00 -0.6918656E+00
!       -0.6657920E+00 -0.5960342E+00 -0.4164121E+00
!
!***REFERENCES  M. J. D. Powell, A hybrid method for nonlinear equa-
!                 tions. In Numerical Methods for Nonlinear Algebraic
!                 Equations, P. Rabinowitz, Editor.  Gordon and Breach,
!                 1988.
!***ROUTINES CALLED  D1MACH, D1MPYQ, D1UPDT, DDOGLG, DENORM, DFDJC1,
!                    DQFORM, DQRFAC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DNSQ
  DOUBLE PRECISION D1MACH,DENORM
  INTEGER I, IFLAG, INFO, IOPT, ITER, IWA(1), J, JM1, L, LDFJAC, &
       LR, MAXFEV, ML, MODE, MU, N, NCFAIL, NCSUC, NFEV, NJEV, &
       NPRINT, NSLOW1, NSLOW2
  DOUBLE PRECISION ACTRED, DELTA, DIAG(*), EPSFCN, EPSMCH, FACTOR, &
       FJAC(LDFJAC,*), FNORM, FNORM1, FVEC(*), ONE, P0001, P001, &
       P1, P5, PNORM, PRERED, QTF(*), R(*), RATIO, SUM, TEMP, &
       WA1(*), WA2(*), WA3(*), WA4(*), X(*), XNORM, XTOL, ZERO
  EXTERNAL FCN
  LOGICAL JEVAL,SING
  SAVE ONE, P1, P5, P001, P0001, ZERO
  DATA ONE,P1,P5,P001,P0001,ZERO &
       /1.0D0,1.0D-1,5.0D-1,1.0D-3,1.0D-4,0.0D0/
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 320
!***FIRST EXECUTABLE STATEMENT  DNSQ
     EPSMCH = D1MACH(4)
!
     INFO = 0
     IFLAG = 0
     NFEV = 0
     NJEV = 0
!
!        CHECK THE INPUT PARAMETERS FOR ERRORS.
!
!     ...EXIT
     if (IOPT  <  1 .OR. IOPT  >  2 .OR. N  <=  0 &
         .OR. XTOL  <  ZERO .OR. MAXFEV  <=  0 .OR. ML  <  0 &
         .OR. MU  <  0 .OR. FACTOR  <=  ZERO .OR. LDFJAC  <  N &
         .OR. LR  <  (N*(N + 1))/2) go to 320
     if (MODE  /=  2) go to 20
        DO 10 J = 1, N
!     .........EXIT
           if (DIAG(J)  <=  ZERO) go to 320
   10       CONTINUE
   20    CONTINUE
!
!        EVALUATE THE FUNCTION AT THE STARTING POINT
!        AND CALCULATE ITS NORM.
!
     IFLAG = 1
     call FCN(N,X,FVEC,IFLAG)
     NFEV = 1
!     ...EXIT
     if (IFLAG  <  0) go to 320
     FNORM = DENORM(N,FVEC)
!
!        INITIALIZE ITERATION COUNTER AND MONITORS.
!
     ITER = 1
     NCSUC = 0
     NCFAIL = 0
     NSLOW1 = 0
     NSLOW2 = 0
!
!        BEGINNING OF THE OUTER LOOP.
!
   30    CONTINUE
!           BEGIN BLOCK PERMITTING ...EXITS TO 90
           JEVAL = .TRUE.
!
!              CALCULATE THE JACOBIAN MATRIX.
!
           if (IOPT  ==  2) go to 40
!
!                 USER SUPPLIES JACOBIAN
!
              call JAC(N,X,FVEC,FJAC,LDFJAC,IFLAG)
              NJEV = NJEV + 1
           go to 50
   40          CONTINUE
!
!                 CODE APPROXIMATES THE JACOBIAN
!
              IFLAG = 2
              call DFDJC1(FCN,N,X,FVEC,FJAC,LDFJAC,IFLAG,ML,MU, &
                          EPSFCN,WA1,WA2)
              NFEV = NFEV + MIN(ML+MU+1,N)
   50          CONTINUE
!
!     .........EXIT
           if (IFLAG  <  0) go to 320
!
!              COMPUTE THE QR FACTORIZATION OF THE JACOBIAN.
!
           call DQRFAC(N,N,FJAC,LDFJAC,.FALSE.,IWA,1,WA1,WA2,WA3)
!
!              ON THE FIRST ITERATION AND if MODE IS 1, SCALE ACCORDING
!              TO THE NORMS OF THE COLUMNS OF THE INITIAL JACOBIAN.
!
!           ...EXIT
           if (ITER  /=  1) go to 90
           if (MODE  ==  2) go to 70
              DO 60 J = 1, N
                 DIAG(J) = WA2(J)
                 if (WA2(J)  ==  ZERO) DIAG(J) = ONE
   60             CONTINUE
   70          CONTINUE
!
!              ON THE FIRST ITERATION, CALCULATE THE NORM OF THE SCALED
!              X AND INITIALIZE THE STEP BOUND DELTA.
!
           DO 80 J = 1, N
              WA3(J) = DIAG(J)*X(J)
   80          CONTINUE
           XNORM = DENORM(N,WA3)
           DELTA = FACTOR*XNORM
           if (DELTA  ==  ZERO) DELTA = FACTOR
   90       CONTINUE
!
!           FORM (Q TRANSPOSE)*FVEC AND STORE IN QTF.
!
        DO 100 I = 1, N
           QTF(I) = FVEC(I)
  100       CONTINUE
        DO 140 J = 1, N
           if (FJAC(J,J)  ==  ZERO) go to 130
              SUM = ZERO
              DO 110 I = J, N
                 SUM = SUM + FJAC(I,J)*QTF(I)
  110             CONTINUE
              TEMP = -SUM/FJAC(J,J)
              DO 120 I = J, N
                 QTF(I) = QTF(I) + FJAC(I,J)*TEMP
  120             CONTINUE
  130          CONTINUE
  140       CONTINUE
!
!           COPY THE TRIANGULAR FACTOR OF THE QR FACTORIZATION INTO R.
!
        SING = .FALSE.
        DO 170 J = 1, N
           L = J
           JM1 = J - 1
           if (JM1  <  1) go to 160
           DO 150 I = 1, JM1
              R(L) = FJAC(I,J)
              L = L + N - I
  150          CONTINUE
  160          CONTINUE
           R(L) = WA1(J)
           if (WA1(J)  ==  ZERO) SING = .TRUE.
  170       CONTINUE
!
!           ACCUMULATE THE ORTHOGONAL FACTOR IN FJAC.
!
        call DQFORM(N,N,FJAC,LDFJAC,WA1)
!
!           RESCALE if NECESSARY.
!
        if (MODE  ==  2) go to 190
           DO 180 J = 1, N
              DIAG(J) = MAX(DIAG(J),WA2(J))
  180          CONTINUE
  190       CONTINUE
!
!           BEGINNING OF THE INNER LOOP.
!
  200       CONTINUE
!
!              if REQUESTED, call FCN TO ENABLE PRINTING OF ITERATES.
!
           if (NPRINT  <=  0) go to 210
              IFLAG = 0
              if (MOD(ITER-1,NPRINT)  ==  0) &
                 call FCN(N,X,FVEC,IFLAG)
!     ............EXIT
              if (IFLAG  <  0) go to 320
  210          CONTINUE
!
!              DETERMINE THE DIRECTION P.
!
           call DDOGLG(N,R,LR,DIAG,QTF,DELTA,WA1,WA2,WA3)
!
!              STORE THE DIRECTION P AND X + P. CALCULATE THE NORM OF P.
!
           DO 220 J = 1, N
              WA1(J) = -WA1(J)
              WA2(J) = X(J) + WA1(J)
              WA3(J) = DIAG(J)*WA1(J)
  220          CONTINUE
           PNORM = DENORM(N,WA3)
!
!              ON THE FIRST ITERATION, ADJUST THE INITIAL STEP BOUND.
!
           if (ITER  ==  1) DELTA = MIN(DELTA,PNORM)
!
!              EVALUATE THE FUNCTION AT X + P AND CALCULATE ITS NORM.
!
           IFLAG = 1
           call FCN(N,WA2,WA4,IFLAG)
           NFEV = NFEV + 1
!     .........EXIT
           if (IFLAG  <  0) go to 320
           FNORM1 = DENORM(N,WA4)
!
!              COMPUTE THE SCALED ACTUAL REDUCTION.
!
           ACTRED = -ONE
           if (FNORM1  <  FNORM) ACTRED = ONE - (FNORM1/FNORM)**2
!
!              COMPUTE THE SCALED PREDICTED REDUCTION.
!
           L = 1
           DO 240 I = 1, N
              SUM = ZERO
              DO 230 J = I, N
                 SUM = SUM + R(L)*WA1(J)
                 L = L + 1
  230             CONTINUE
              WA3(I) = QTF(I) + SUM
  240          CONTINUE
           TEMP = DENORM(N,WA3)
           PRERED = ZERO
           if (TEMP  <  FNORM) PRERED = ONE - (TEMP/FNORM)**2
!
!              COMPUTE THE RATIO OF THE ACTUAL TO THE PREDICTED
!              REDUCTION.
!
           RATIO = ZERO
           if (PRERED  >  ZERO) RATIO = ACTRED/PRERED
!
!              UPDATE THE STEP BOUND.
!
           if (RATIO  >=  P1) go to 250
              NCSUC = 0
              NCFAIL = NCFAIL + 1
              DELTA = P5*DELTA
           go to 260
  250          CONTINUE
              NCFAIL = 0
              NCSUC = NCSUC + 1
              if (RATIO  >=  P5 .OR. NCSUC  >  1) &
                 DELTA = MAX(DELTA,PNORM/P5)
              if (ABS(RATIO-ONE)  <=  P1) DELTA = PNORM/P5
  260          CONTINUE
!
!              TEST FOR SUCCESSFUL ITERATION.
!
           if (RATIO  <  P0001) go to 280
!
!                 SUCCESSFUL ITERATION. UPDATE X, FVEC, AND THEIR NORMS.
!
              DO 270 J = 1, N
                 X(J) = WA2(J)
                 WA2(J) = DIAG(J)*X(J)
                 FVEC(J) = WA4(J)
  270             CONTINUE
              XNORM = DENORM(N,WA2)
              FNORM = FNORM1
              ITER = ITER + 1
  280          CONTINUE
!
!              DETERMINE THE PROGRESS OF THE ITERATION.
!
           NSLOW1 = NSLOW1 + 1
           if (ACTRED  >=  P001) NSLOW1 = 0
           if (JEVAL) NSLOW2 = NSLOW2 + 1
           if (ACTRED  >=  P1) NSLOW2 = 0
!
!              TEST FOR CONVERGENCE.
!
           if (DELTA  <=  XTOL*XNORM .OR. FNORM  ==  ZERO) INFO = 1
!     .........EXIT
           if (INFO  /=  0) go to 320
!
!              TESTS FOR TERMINATION AND STRINGENT TOLERANCES.
!
           if (NFEV  >=  MAXFEV) INFO = 2
           if (P1*MAX(P1*DELTA,PNORM)  <=  EPSMCH*XNORM) INFO = 3
           if (NSLOW2  ==  5) INFO = 4
           if (NSLOW1  ==  10) INFO = 5
!     .........EXIT
           if (INFO  /=  0) go to 320
!
!              CRITERION FOR RECALCULATING JACOBIAN
!
!           ...EXIT
           if (NCFAIL  ==  2) go to 310
!
!              CALCULATE THE RANK ONE MODIFICATION TO THE JACOBIAN
!              AND UPDATE QTF if NECESSARY.
!
           DO 300 J = 1, N
              SUM = ZERO
              DO 290 I = 1, N
                 SUM = SUM + FJAC(I,J)*WA4(I)
  290             CONTINUE
              WA2(J) = (SUM - WA3(J))/PNORM
              WA1(J) = DIAG(J)*((DIAG(J)*WA1(J))/PNORM)
              if (RATIO  >=  P0001) QTF(J) = SUM
  300          CONTINUE
!
!              COMPUTE THE QR FACTORIZATION OF THE UPDATED JACOBIAN.
!
           call D1UPDT(N,N,R,LR,WA1,WA2,WA3,SING)
           call D1MPYQ(N,N,FJAC,LDFJAC,WA2,WA3)
           call D1MPYQ(1,N,QTF,1,WA2,WA3)
!
!              END OF THE INNER LOOP.
!
           JEVAL = .FALSE.
        go to 200
  310       CONTINUE
!
!           END OF THE OUTER LOOP.
!
     go to 30
  320 CONTINUE
!
!     TERMINATION, EITHER NORMAL OR USER IMPOSED.
!
  if (IFLAG  <  0) INFO = IFLAG
  IFLAG = 0
  if (NPRINT  >  0) call FCN(N,X,FVEC,IFLAG)
  if (INFO  <  0) call XERMSG ('SLATEC', 'DNSQ', &
     'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
  if (INFO  ==  0) call XERMSG ('SLATEC', 'DNSQ', &
     'INVALID INPUT PARAMETER.', 2, 1)
  if (INFO  ==  2) call XERMSG ('SLATEC', 'DNSQ', &
     'TOO MANY FUNCTION EVALUATIONS.', 9, 1)
  if (INFO  ==  3) call XERMSG ('SLATEC', 'DNSQ', &
     'XTOL TOO SMALL. NO FURTHER IMPROVEMENT POSSIBLE.', 3, 1)
  if (INFO  >  4) call XERMSG ('SLATEC', 'DNSQ', &
     'ITERATION NOT MAKING GOOD PROGRESS.', 1, 1)
  return
!
!     LAST CARD OF SUBROUTINE DNSQ.
!
end