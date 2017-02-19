subroutine DSOSEQ (FNC, N, S, RTOLX, ATOLX, TOLF, IFLAG, MXIT, &
     NCJS, NSRRC, NSRI, IPRINT, FMAX, C, NC, B, P, TEMP, X, Y, FAC, &
     IS)
!
!! DSOSEQ is subsidiary to DSOS
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SOSEQS-S, DSOSEQ-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     DSOSEQ solves a system of N simultaneous nonlinear equations.
!     See the comments in the interfacing routine DSOS for a more
!     detailed description of some of the items in the calling list.
!
! **********************************************************************
!   -Input-
!
!     FNC- Function subprogram which evaluates the equations
!     N  -number of equations
!     S  -Solution vector of initial guesses
!     RTOLX-Relative error tolerance on solution components
!     ATOLX-Absolute error tolerance on solution components
!     TOLF-Residual error tolerance
!     MXIT-Maximum number of allowable iterations.
!     NCJS-Maximum number of consecutive iterative steps to perform
!          using the same triangular Jacobian matrix approximation.
!     NSRRC-Number of consecutive iterative steps for which the
!          limiting precision accuracy test must be satisfied
!          before the routine exits with IFLAG=4.
!     NSRI-Number of consecutive iterative steps for which the
!          diverging condition test must be satisfied before
!          the routine exits with IFLAG=7.
!     IPRINT-Internal printing parameter. You must set IPRINT=-1 if you
!          want the intermediate solution iterates and a residual norm
!          to be printed.
!     C   -Internal work array, dimensioned at least N*(N+1)/2.
!     NC  -Dimension of C array. NC   >=   N*(N+1)/2.
!     B   -Internal work array, dimensioned N.
!     P   -Internal work array, dimensioned N.
!     TEMP-Internal work array, dimensioned N.
!     X   -Internal work array, dimensioned N.
!     Y   -Internal work array, dimensioned N.
!     FAC -Internal work array, dimensioned N.
!     IS  -Internal work array, dimensioned N.
!
!   -Output-
!     S    -Solution vector
!     IFLAG-Status indicator flag
!     MXIT-The actual number of iterations performed
!     FMAX-Residual norm
!     C   -Upper unit triangular matrix which approximates the
!          forward triangularization of the full Jacobian matrix.
!          Stored in a vector with dimension at least N*(N+1)/2.
!     B   -Contains the residuals (function values) divided
!          by the corresponding components of the P vector
!     P   -Array used to store the partial derivatives. After
!          each iteration P(K) contains the maximal derivative
!          occurring in the K-th reduced equation.
!     TEMP-Array used to store the previous solution iterate.
!     X   -Solution vector. Contains the values achieved on the
!          last iteration loop upon exit from DSOS.
!     Y   -Array containing the solution increments.
!     FAC -Array containing factors used in computing numerical
!          derivatives.
!     IS  -Records the pivotal information (column interchanges)
!
! **********************************************************************
! *** Three machine dependent parameters appear in this subroutine.
!
! *** The smallest positive magnitude, zero, is defined by the function
! *** routine D1MACH(1).
!
! *** URO, the computer unit roundoff value, is defined by D1MACH(3) for
! *** machines that round or D1MACH(4) for machines that truncate.
! *** URO is the smallest positive number such that 1.+URO   >   1.
!
! *** The output tape unit number, LOUN, is defined by the function
! *** I1MACH(2).
! **********************************************************************
!
!***SEE ALSO  DSOS
!***ROUTINES CALLED  D1MACH, DSOSSL, I1MACH
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DSOSEQ
!
!
  INTEGER I1MACH
  DOUBLE PRECISION D1MACH
  INTEGER IC, ICR, IFLAG, IPRINT, IS(*), ISJ, ISV, IT, ITEM, ITRY, &
       J, JK, JS, K, KD, KJ, KK, KM1, KN, KSV, L, LOUN, LS, M, MIT, &
       MM, MXIT, N, NC, NCJS, NP1, NSRI, NSRRC
  DOUBLE PRECISION ATOLX, B(*), C(*), CSV, F, FAC(*), FACT, FDIF, &
       FMAX, FMIN, FMXS, FN1, FN2, FNC, FP, H, HX, P(*), PMAX, RE, &
       RTOLX, S(*), SRURO, TEMP(*), TEST, TOLF, URO, X(*), XNORM, &
       Y(*), YJ, YN1, YN2, YN3, YNORM, YNS, ZERO
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 430
!        BEGIN BLOCK PERMITTING ...EXITS TO 410
!           BEGIN BLOCK PERMITTING ...EXITS TO 390
!***FIRST EXECUTABLE STATEMENT  DSOSEQ
           URO = D1MACH(4)
           LOUN = I1MACH(2)
           ZERO = D1MACH(1)
           RE = MAX(RTOLX,URO)
           SRURO = SQRT(URO)
!
           IFLAG = 0
           NP1 = N + 1
           ICR = 0
           IC = 0
           ITRY = NCJS
           YN1 = 0.0D0
           YN2 = 0.0D0
           YN3 = 0.0D0
           YNS = 0.0D0
           MIT = 0
           FN1 = 0.0D0
           FN2 = 0.0D0
           FMXS = 0.0D0
!
!              INITIALIZE THE INTERCHANGE (PIVOTING) VECTOR AND
!              SAVE THE CURRENT SOLUTION APPROXIMATION FOR FUTURE USE.
!
           DO 10 K = 1, N
              IS(K) = K
              X(K) = S(K)
              TEMP(K) = X(K)
   10          CONTINUE
!
!
!              *********************************************************
!              **** BEGIN PRINCIPAL ITERATION LOOP  ****
!              *********************************************************
!
           DO 380 M = 1, MXIT
!                 BEGIN BLOCK PERMITTING ...EXITS TO 350
!                    BEGIN BLOCK PERMITTING ...EXITS TO 240
!
                    DO 20 K = 1, N
                       FAC(K) = SRURO
   20                   CONTINUE
!
   30                   CONTINUE
!                          BEGIN BLOCK PERMITTING ...EXITS TO 180
                          KN = 1
                          FMAX = 0.0D0
!
!
!                             ******** BEGIN SUBITERATION LOOP DEFINING
!                             THE LINEARIZATION OF EACH ********
!                             EQUATION WHICH RESULTS IN THE CONSTRUCTION
!                             OF AN UPPER ******** TRIANGULAR MATRIX
!                             APPROXIMATING THE FORWARD ********
!                             TRIANGULARIZATION OF THE FULL JACOBIAN
!                             MATRIX
!
                          DO 170 K = 1, N
!                                BEGIN BLOCK PERMITTING ...EXITS TO 160
                                KM1 = K - 1
!
!                                   BACK-SOLVE A TRIANGULAR LINEAR
!                                   SYSTEM OBTAINING IMPROVED SOLUTION
!                                   VALUES FOR K-1 OF THE VARIABLES FROM
!                                   THE FIRST K-1 EQUATIONS. THESE
!                                   VARIABLES ARE THEN ELIMINATED FROM
!                                   THE K-TH EQUATION.
!
                                if (KM1  ==  0) go to 50
                                   call DSOSSL(K,N,KM1,Y,C,B,KN)
                                   DO 40 J = 1, KM1
                                      JS = IS(J)
                                      X(JS) = TEMP(JS) + Y(J)
   40                                  CONTINUE
   50                               CONTINUE
!
!
!                                   EVALUATE THE K-TH EQUATION AND THE
!                                   INTERMEDIATE COMPUTATION FOR THE MAX
!                                   NORM OF THE RESIDUAL VECTOR.
!
                                F = FNC(X,K)
                                FMAX = MAX(FMAX,ABS(F))
!
!                                   if WE WISH TO PERFORM SEVERAL
!                                   ITERATIONS USING A FIXED
!                                   FACTORIZATION OF AN APPROXIMATE
!                                   JACOBIAN,WE NEED ONLY UPDATE THE
!                                   CONSTANT VECTOR.
!
!                                ...EXIT
                                if (ITRY  <  NCJS) go to 160
!
!
                                IT = 0
!
!                                   COMPUTE PARTIAL DERIVATIVES THAT ARE
!                                   REQUIRED IN THE LINEARIZATION OF THE
!                                   K-TH REDUCED EQUATION
!
                                DO 90 J = K, N
                                   ITEM = IS(J)
                                   HX = X(ITEM)
                                   H = FAC(ITEM)*HX
                                   if (ABS(H)  <=  ZERO) &
                                      H = FAC(ITEM)
                                   X(ITEM) = HX + H
                                   if (KM1  ==  0) go to 70
                                      Y(J) = H
                                      call DSOSSL(K,N,J,Y,C,B,KN)
                                      DO 60 L = 1, KM1
                                         LS = IS(L)
                                         X(LS) = TEMP(LS) + Y(L)
   60                                     CONTINUE
   70                                  CONTINUE
                                   FP = FNC(X,K)
                                   X(ITEM) = HX
                                   FDIF = FP - F
                                   if (ABS(FDIF)  >  URO*ABS(F)) &
                                      go to 80
                                      FDIF = 0.0D0
                                      IT = IT + 1
   80                                  CONTINUE
                                   P(J) = FDIF/H
   90                               CONTINUE
!
                                if (IT  <=  (N - K)) go to 110
!
!                                      ALL COMPUTED PARTIAL DERIVATIVES
!                                      OF THE K-TH EQUATION ARE
!                                      EFFECTIVELY ZERO.TRY LARGER
!                                      PERTURBATIONS OF THE INDEPENDENT
!                                      VARIABLES.
!
                                   DO 100 J = K, N
                                      ISJ = IS(J)
                                      FACT = 100.0D0*FAC(ISJ)
!           ..............................EXIT
                                      if (FACT  >  1.0D10) &
                                         go to 390
                                      FAC(ISJ) = FACT
  100                                  CONTINUE
!                          ............EXIT
                                   go to 180
  110                               CONTINUE
!
!                                ...EXIT
                                if (K  ==  N) go to 160
!
!                                   ACHIEVE A PIVOTING EFFECT BY
!                                   CHOOSING THE MAXIMAL DERIVATIVE
!                                   ELEMENT
!
                                PMAX = 0.0D0
                                DO 130 J = K, N
                                   TEST = ABS(P(J))
                                   if (TEST  <=  PMAX) go to 120
                                      PMAX = TEST
                                      ISV = J
  120                                  CONTINUE
  130                               CONTINUE
!           ........................EXIT
                                if (PMAX  ==  0.0D0) go to 390
!
!                                   SET UP THE COEFFICIENTS FOR THE K-TH
!                                   ROW OF THE TRIANGULAR LINEAR SYSTEM
!                                   AND SAVE THE PARTIAL DERIVATIVE OF
!                                   LARGEST MAGNITUDE
!
                                PMAX = P(ISV)
                                KK = KN
                                DO 140 J = K, N
                                   if (J  /=  ISV) &
                                      C(KK) = -P(J)/PMAX
                                   KK = KK + 1
  140                               CONTINUE
                                P(K) = PMAX
!
!
!                                ...EXIT
                                if (ISV  ==  K) go to 160
!
!                                   INTERCHANGE THE TWO COLUMNS OF C
!                                   DETERMINED BY THE PIVOTAL STRATEGY
!
                                KSV = IS(K)
                                IS(K) = IS(ISV)
                                IS(ISV) = KSV
!
                                KD = ISV - K
                                KJ = K
                                DO 150 J = 1, K
                                   CSV = C(KJ)
                                   JK = KJ + KD
                                   C(KJ) = C(JK)
                                   C(JK) = CSV
                                   KJ = KJ + N - J
  150                               CONTINUE
  160                            CONTINUE
!
                             KN = KN + NP1 - K
!
!                                STORE THE COMPONENTS FOR THE CONSTANT
!                                VECTOR
!
                             B(K) = -F/P(K)
!
  170                         CONTINUE
!                       ......EXIT
                          go to 190
  180                      CONTINUE
                    go to 30
  190                   CONTINUE
!
!                       ********
!                       ******** END OF LOOP CREATING THE TRIANGULAR
!                       LINEARIZATION MATRIX
!                       ********
!
!
!                        SOLVE THE RESULTING TRIANGULAR SYSTEM FOR A NEW
!                        SOLUTION APPROXIMATION AND OBTAIN THE SOLUTION
!                        INCREMENT NORM.
!
                    KN = KN - 1
                    Y(N) = B(N)
                    if (N  >  1) call DSOSSL(N,N,N,Y,C,B,KN)
                    XNORM = 0.0D0
                    YNORM = 0.0D0
                    DO 200 J = 1, N
                       YJ = Y(J)
                       YNORM = MAX(YNORM,ABS(YJ))
                       JS = IS(J)
                       X(JS) = TEMP(JS) + YJ
                       XNORM = MAX(XNORM,ABS(X(JS)))
  200                   CONTINUE
!
!
!                       PRINT INTERMEDIATE SOLUTION ITERATES AND
!                       RESIDUAL NORM if DESIRED
!
                    if (IPRINT  /=  (-1)) go to 220
                       MM = M - 1
                       WRITE (LOUN,210) FMAX,MM,(X(J), J = 1, N)
  210                      FORMAT ('0RESIDUAL NORM =', D9.2, / 1X, &
                               'SOLUTION ITERATE (', I3, ')', / &
                               (1X, 5D26.14))
  220                   CONTINUE
!
!                       TEST FOR CONVERGENCE TO A SOLUTION (RELATIVE
!                       AND/OR ABSOLUTE ERROR COMPARISON ON SUCCESSIVE
!                       APPROXIMATIONS OF EACH SOLUTION VARIABLE)
!
                    DO 230 J = 1, N
                       JS = IS(J)
!                    ......EXIT
                       if (ABS(Y(J))  >  RE*ABS(X(JS)) + ATOLX) &
                          go to 240
  230                   CONTINUE
                    if (FMAX  <=  FMXS) IFLAG = 1
  240                CONTINUE
!
!                    TEST FOR CONVERGENCE TO A SOLUTION BASED ON
!                    RESIDUALS
!
                 if (FMAX  <=  TOLF) IFLAG = IFLAG + 2
!        ............EXIT
                 if (IFLAG  >  0) go to 410
!
!
                 if (M  >  1) go to 250
                    FMIN = FMAX
                 go to 330
  250                CONTINUE
!                       BEGIN BLOCK PERMITTING ...EXITS TO 320
!
!                          SAVE SOLUTION HAVING MINIMUM RESIDUAL NORM.
!
                       if (FMAX  >=  FMIN) go to 270
                          MIT = M + 1
                          YN1 = YNORM
                          YN2 = YNS
                          FN1 = FMXS
                          FMIN = FMAX
                          DO 260 J = 1, N
                             S(J) = X(J)
  260                         CONTINUE
                          IC = 0
  270                      CONTINUE
!
!                          TEST FOR LIMITING PRECISION CONVERGENCE. VERY
!                          SLOWLY CONVERGENT PROBLEMS MAY ALSO BE
!                          DETECTED.
!
                       if (YNORM  >  SRURO*XNORM) go to 290
                       if (FMAX  <  0.2D0*FMXS &
                           .OR. FMAX  >  5.0D0*FMXS) go to 290
                       if (YNORM  <  0.2D0*YNS &
                           .OR. YNORM  >  5.0D0*YNS) go to 290
                          ICR = ICR + 1
                          if (ICR  >=  NSRRC) go to 280
                             IC = 0
!                       .........EXIT
                             go to 320
  280                         CONTINUE
                          IFLAG = 4
                          FMAX = FMIN
!     ........................EXIT
                          go to 430
  290                      CONTINUE
                       ICR = 0
!
!                          TEST FOR DIVERGENCE OF THE ITERATIVE SCHEME.
!
                       if (YNORM  >  2.0D0*YNS &
                           .OR. FMAX  >  2.0D0*FMXS) go to 300
                          IC = 0
                       go to 310
  300                      CONTINUE
                          IC = IC + 1
!                       ......EXIT
                          if (IC  <  NSRI) go to 320
                          IFLAG = 7
!        .....................EXIT
                          go to 410
  310                      CONTINUE
  320                   CONTINUE
  330                CONTINUE
!
!                    CHECK TO SEE if NEXT ITERATION CAN USE THE OLD
!                    JACOBIAN FACTORIZATION
!
                 ITRY = ITRY - 1
                 if (ITRY  ==  0) go to 340
                 if (20.0D0*YNORM  >  XNORM) go to 340
                 if (YNORM  >  2.0D0*YNS) go to 340
!                 ......EXIT
                    if (FMAX  <  2.0D0*FMXS) go to 350
  340                CONTINUE
                 ITRY = NCJS
  350             CONTINUE
!
!                 SAVE THE CURRENT SOLUTION APPROXIMATION AND THE
!                 RESIDUAL AND SOLUTION INCREMENT NORMS FOR USE IN THE
!                 NEXT ITERATION.
!
              DO 360 J = 1, N
                 TEMP(J) = X(J)
  360             CONTINUE
              if (M  /=  MIT) go to 370
                 FN2 = FMAX
                 YN3 = YNORM
  370             CONTINUE
              FMXS = FMAX
              YNS = YNORM
!
!
  380          CONTINUE
!
!              *********************************************************
!              **** END OF PRINCIPAL ITERATION LOOP ****
!              *********************************************************
!
!
!               TOO MANY ITERATIONS, CONVERGENCE WAS NOT ACHIEVED.
           M = MXIT
           IFLAG = 5
           if (YN1  >  10.0D0*YN2 .OR. YN3  >  10.0D0*YN1) &
              IFLAG = 6
           if (FN1  >  5.0D0*FMIN .OR. FN2  >  5.0D0*FMIN) &
              IFLAG = 6
           if (FMAX  >  5.0D0*FMIN) IFLAG = 6
!        ......EXIT
           go to 410
  390       CONTINUE
!
!
!           A JACOBIAN-RELATED MATRIX IS EFFECTIVELY SINGULAR.
        IFLAG = 8
        DO 400 J = 1, N
           S(J) = TEMP(J)
  400       CONTINUE
!     ......EXIT
        go to 430
  410    CONTINUE
!
!
     DO 420 J = 1, N
        S(J) = X(J)
  420    CONTINUE
  430 CONTINUE
!
!
  MXIT = M
  return
end
