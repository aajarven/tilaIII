subroutine SCOV (FCN, IOPT, M, N, X, FVEC, R, LDR, INFO, WA1, WA2, &
     WA3, WA4)
!
!! SCOV calculates the covariance matrix for a nonlinear data fitting problem.
!  It is intended to be used after a
!            successful return from either SNLS1 or SNLS1E.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1B1
!***TYPE      SINGLE PRECISION (SCOV-S, DCOV-D)
!***KEYWORDS  COVARIANCE MATRIX, NONLINEAR DATA FITTING,
!             NONLINEAR LEAST SQUARES
!***AUTHOR  Hiebert, K. L., (SNLA)
!***DESCRIPTION
!
!  1. Purpose.
!
!     SCOV calculates the covariance matrix for a nonlinear data
!     fitting problem.  It is intended to be used after a
!     successful return from either SNLS1 or SNLS1E. SCOV
!     and SNLS1 (and SNLS1E) have compatible parameters.  The
!     required external subroutine, FCN, is the same
!     for all three codes, SCOV, SNLS1, and SNLS1E.
!
!  2. Subroutine and Type Statements.
!
!     SUBROUTINE SCOV(FCN,IOPT,M,N,X,FVEC,R,LDR,INFO,
!                     WA1,WA2,WA3,WA4)
!     INTEGER IOPT,M,N,LDR,INFO
!     REAL X(N),FVEC(M),R(LDR,N),WA1(N),WA2(N),WA3(N),WA4(M)
!     EXTERNAL FCN
!
!  3. Parameters.
!
!       FCN is the name of the user-supplied subroutine which calculates
!         the functions.  If the user wants to supply the Jacobian
!         (IOPT=2 or 3), then FCN must be written to calculate the
!         Jacobian, as well as the functions.  See the explanation
!         of the IOPT argument below.  FCN must be declared in an
!         EXTERNAL statement in the calling program and should be
!         written as follows.
!
!         SUBROUTINE FCN(IFLAG,M,N,X,FVEC,FJAC,LDFJAC)
!         INTEGER IFLAG,LDFJAC,M,N
!         REAL X(N),FVEC(M)
!         ----------
!         FJAC and LDFJAC may be ignored     , if IOPT=1.
!         REAL FJAC(LDFJAC,N)                , if IOPT=2.
!         REAL FJAC(N)                       , if IOPT=3.
!         ----------
!           IFLAG will never be zero when FCN is called by SCOV.
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
!           never be 3 unless IOPT=3.  FJAC(J) must be set to
!           the derivative of FVEC(LDFJAC) with respect to X(J).
!         return
!         ----------
!         END
!
!
!         The value of IFLAG should not be changed by FCN unless the
!         user wants to terminate execution of SCOV.  In this case, set
!         IFLAG to a negative integer.
!
!
!    IOPT is an input variable which specifies how the Jacobian will
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
!       X is an array of length N.  On input X must contain the value
!         at which the covariance matrix is to be evaluated.  This is
!         usually the value for X returned from a successful run of
!         SNLS1 (or SNLS1E).  The value of X will not be changed.
!
!    FVEC is an output array of length M which contains the functions
!         evaluated at X.
!
!       R is an output array.  For IOPT=1 and 2, R is an M by N array.
!         For IOPT=3, R is an N by N array.  On output, if INFO=1,
!         the upper N by N submatrix of R contains the covariance
!         matrix evaluated at X.
!
!     LDR is a positive integer input variable which specifies
!         the leading dimension of the array R.  For IOPT=1 and 2,
!         LDR must not be less than M.  For IOPT=3, LDR must not
!         be less than N.
!
!    INFO is an integer output variable.  If the user has terminated
!         execution, INFO is set to the (negative) value of IFLAG.  See
!         description of FCN. Otherwise, INFO is set as follows.
!
!         INFO = 0 Improper input parameters (M <= 0 or N <= 0).
!
!         INFO = 1 Successful return.  The covariance matrix has been
!                  calculated and stored in the upper N by N
!                  submatrix of R.
!
!         INFO = 2 The Jacobian matrix is singular for the input value
!                  of X.  The covariance matrix cannot be calculated.
!                  The upper N by N submatrix of R contains the QR
!                  factorization of the Jacobian (probably not of
!                  interest to the user).
!
!     WA1 is a work array of length N.
!     WA2 is a work array of length N.
!     WA3 is a work array of length N.
!     WA4 is a work array of length M.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ENORM, FDJAC3, QRFAC, RWUPDT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810522  DATE WRITTEN
!   890505  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Fixed an error message.  (RWC)
!***END PROLOGUE  SCOV
!
!     REVISED 820707-1100
!     REVISED YYMMDD HHMM
!
  INTEGER I,IDUM,IFLAG,INFO,IOPT,J,K,KP1,LDR,M,N,NM1,NROW
  REAL X(*),R(LDR,*),FVEC(*),WA1(*),WA2(*),WA3(*),WA4(*)
  EXTERNAL FCN
  REAL ONE,SIGMA,TEMP,ZERO
  LOGICAL SING
  SAVE ZERO, ONE
  DATA ZERO/0.E0/,ONE/1.E0/
!***FIRST EXECUTABLE STATEMENT  SCOV
  SING=.FALSE.
  IFLAG=0
  if (M <= 0 .OR. N <= 0) go to 300
!
!     CALCULATE SIGMA = (SUM OF THE SQUARED RESIDUALS) / (M-N)
  IFLAG=1
  call FCN(IFLAG,M,N,X,FVEC,R,LDR)
  if (IFLAG < 0) go to 300
  TEMP=ENORM(M,FVEC)
  SIGMA=ONE
  if (M /= N) SIGMA=TEMP*TEMP/(M-N)
!
!     CALCULATE THE JACOBIAN
  if (IOPT == 3) go to 200
!
!     STORE THE FULL JACOBIAN USING M*N STORAGE
  if (IOPT == 1) go to 100
!
!     USER SUPPLIES THE JACOBIAN
  IFLAG=2
  call FCN(IFLAG,M,N,X,FVEC,R,LDR)
  go to 110
!
!     CODE APPROXIMATES THE JACOBIAN
100   call FDJAC3(FCN,M,N,X,FVEC,R,LDR,IFLAG,ZERO,WA4)
110   if (IFLAG < 0) go to 300
!
!     COMPUTE THE QR DECOMPOSITION
  call QRFAC(M,N,R,LDR,.FALSE.,IDUM,1,WA1,WA1,WA1)
  DO 120 I=1,N
120   R(I,I)=WA1(I)
  go to 225
!
!     COMPUTE THE QR FACTORIZATION OF THE JACOBIAN MATRIX CALCULATED ONE
!     ROW AT A TIME AND STORED IN THE UPPER TRIANGLE OF R.
!     ( (Q TRANSPOSE)*FVEC IS ALSO CALCULATED BUT NOT USED.)
200   CONTINUE
  DO 210 J=1,N
  WA2(J)=ZERO
  DO 205 I=1,N
  R(I,J)=ZERO
205   CONTINUE
210   CONTINUE
  IFLAG=3
  DO 220 I=1,M
  NROW = I
  call FCN(IFLAG,M,N,X,FVEC,WA1,NROW)
  if (IFLAG < 0) go to 300
  TEMP=FVEC(I)
  call RWUPDT(N,R,LDR,WA1,WA2,TEMP,WA3,WA4)
220   CONTINUE
!
!     CHECK if R IS SINGULAR.
225   CONTINUE
  DO 230 I=1,N
  if (R(I,I) == ZERO) SING=.TRUE.
230   CONTINUE
  if (SING) go to 300
!
!     R IS UPPER TRIANGULAR.  CALCULATE (R TRANSPOSE) INVERSE AND STORE
!     IN THE UPPER TRIANGLE OF R.
  if (N == 1) go to 275
  NM1=N-1
  DO 270 K=1,NM1
!
!     INITIALIZE THE RIGHT-HAND SIDE (WA1(*)) AS THE K-TH COLUMN OF THE
!     IDENTITY MATRIX.
  DO 240 I=1,N
  WA1(I)=ZERO
240   CONTINUE
  WA1(K)=ONE
!
  R(K,K)=WA1(K)/R(K,K)
  KP1=K+1
  DO 260 I=KP1,N
!
!     SUBTRACT R(K,I-1)*R(I-1,*) FROM THE RIGHT-HAND SIDE, WA1(*).
  DO 250 J=I,N
  WA1(J)=WA1(J)-R(K,I-1)*R(I-1,J)
250   CONTINUE
  R(K,I)=WA1(I)/R(I,I)
260   CONTINUE
270   CONTINUE
275   R(N,N)=ONE/R(N,N)
!
!     CALCULATE R-INVERSE * (R TRANSPOSE) INVERSE AND STORE IN THE UPPER
!     TRIANGLE OF R.
  DO 290 I=1,N
  DO 290 J=I,N
  TEMP=ZERO
  DO 280 K=J,N
  TEMP=TEMP+R(I,K)*R(J,K)
280   CONTINUE
  R(I,J)=TEMP*SIGMA
290   CONTINUE
  INFO=1
!
300   CONTINUE
  if (M <= 0 .OR. N <= 0) INFO=0
  if (IFLAG < 0) INFO=IFLAG
  if (SING) INFO=2
  if (INFO  <  0) call XERMSG ('SLATEC', 'SCOV', &
     'EXECUTION TERMINATED BECAUSE USER SET IFLAG NEGATIVE.', 1, 1)
  if (INFO  ==  0) call XERMSG ('SLATEC', 'SCOV', &
     'INVALID INPUT PARAMETER.', 2, 1)
  if (INFO  ==  2) call XERMSG ('SLATEC', 'SCOV', &
     'SINGULAR JACOBIAN MATRIX, COVARIANCE MATRIX CANNOT BE ' // &
     'CALCULATED.', 1, 1)
  return
end
