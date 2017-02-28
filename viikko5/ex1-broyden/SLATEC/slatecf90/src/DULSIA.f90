subroutine DULSIA (A, MDA, M, N, B, MDB, NB, RE, AE, KEY, MODE, &
     NP, KRANK, KSURE, RNORM, W, LW, IWORK, LIW, INFO)
!
!! DULSIA solves an underdetermined linear system of equations by ...
!            performing an LQ factorization of the matrix using
!            Householder transformations.  Emphasis is put on detecting
!            possible rank deficiency.
!
!***LIBRARY   SLATEC
!***CATEGORY  D9
!***TYPE      DOUBLE PRECISION (ULSIA-S, DULSIA-D)
!***KEYWORDS  LINEAR LEAST SQUARES, LQ FACTORIZATION,
!             UNDERDETERMINED LINEAR SYSTEM
!***AUTHOR  Manteuffel, T. A., (LANL)
!***DESCRIPTION
!
!     DULSIA computes the minimal length solution(s) to the problem AX=B
!     where A is an M by N matrix with M <= N and B is the M by NB
!     matrix of right hand sides.  User input bounds on the uncertainty
!     in the elements of A are used to detect numerical rank deficiency.
!     The algorithm employs a row and column pivot strategy to
!     minimize the growth of uncertainty and round-off errors.
!
!     DULSIA requires (MDA+1)*N + (MDB+1)*NB + 6*M dimensioned space
!
!   ******************************************************************
!   *                                                                *
!   *         WARNING - All input arrays are changed on exit.        *
!   *                                                                *
!   ******************************************************************
!
!     Input.. All TYPE REAL variables are DOUBLE PRECISION
!
!     A(,)          Linear coefficient matrix of AX=B, with MDA the
!      MDA,M,N      actual first dimension of A in the calling program.
!                   M is the row dimension (no. of EQUATIONS of the
!                   problem) and N the col dimension (no. of UNKNOWNS).
!                   Must have MDA >= M and M <= N.
!
!     B(,)          Right hand side(s), with MDB the actual first
!      MDB,NB       dimension of B in the calling program. NB is the
!                   number of M by 1 right hand sides.  Since the
!                   solution is returned in B, must have MDB >= N.  If
!                   NB = 0, B is never accessed.
!
!   ******************************************************************
!   *                                                                *
!   *         Note - Use of RE and AE are what make this             *
!   *                code significantly different from               *
!   *                other linear least squares solvers.             *
!   *                However, the inexperienced user is              *
!   *                advised to set RE=0.,AE=0.,KEY=0.               *
!   *                                                                *
!   ******************************************************************
!
!     RE(),AE(),KEY
!     RE()          RE() is a vector of length N such that RE(I) is
!                   the maximum relative uncertainty in row I of
!                   the matrix A. The values of RE() must be between
!                   0 and 1. A minimum of 10*machine precision will
!                   be enforced.
!
!     AE()          AE() is a vector of length N such that AE(I) is
!                   the maximum absolute uncertainty in row I of
!                   the matrix A. The values of AE() must be greater
!                   than or equal to 0.
!
!     KEY           For ease of use, RE and AE may be input as either
!                   vectors or scalars. If a scalar is input, the algo-
!                   rithm will use that value for each column of A.
!                   The parameter KEY indicates whether scalars or
!                   vectors are being input.
!                        KEY=0     RE scalar  AE scalar
!                        KEY=1     RE vector  AE scalar
!                        KEY=2     RE scalar  AE vector
!                        KEY=3     RE vector  AE vector
!
!
!     MODE          The integer MODE indicates how the routine
!                   is to react if rank deficiency is detected.
!                   If MODE = 0 return immediately, no solution
!                             1 compute truncated solution
!                             2 compute minimal length least squares sol
!                   The inexperienced user is advised to set MODE=0
!
!     NP            The first NP rows of A will not be interchanged
!                   with other rows even though the pivot strategy
!                   would suggest otherwise.
!                   The inexperienced user is advised to set NP=0.
!
!     WORK()        A real work array dimensioned 5*M.  However, if
!                   RE or AE have been specified as vectors, dimension
!                   WORK 4*M. If both RE and AE have been specified
!                   as vectors, dimension WORK 3*M.
!
!     LW            Actual dimension of WORK
!
!     IWORK()       Integer work array dimensioned at least N+M.
!
!     LIW           Actual dimension of IWORK.
!
!
!     INFO          Is a flag which provides for the efficient
!                   solution of subsequent problems involving the
!                   same A but different B.
!                   If INFO = 0 original call
!                      INFO = 1 subsequent calls
!                   On subsequent calls, the user must supply A, KRANK,
!                   LW, IWORK, LIW, and the first 2*M locations of WORK
!                   as output by the original call to DULSIA. MODE must
!                   be equal to the value of MODE in the original call.
!                   If MODE < 2, only the first N locations of WORK
!                   are accessed. AE, RE, KEY, and NP are not accessed.
!
!
!
!
!     Output..All TYPE REAL variables are DOUBLE PRECISION
!
!     A(,)          Contains the lower triangular part of the reduced
!                   matrix and the transformation information. It togeth
!                   with the first M elements of WORK (see below)
!                   completely specify the LQ factorization of A.
!
!     B(,)          Contains the N by NB solution matrix for X.
!
!     KRANK,KSURE   The numerical rank of A,  based upon the relative
!                   and absolute bounds on uncertainty, is bounded
!                   above by KRANK and below by KSURE. The algorithm
!                   returns a solution based on KRANK. KSURE provides
!                   an indication of the precision of the rank.
!
!     RNORM()       Contains the Euclidean length of the NB residual
!                   vectors  B(I)-AX(I), I=1,NB. If the matrix A is of
!                   full rank, then RNORM=0.0.
!
!     WORK()        The first M locations of WORK contain values
!                   necessary to reproduce the Householder
!                   transformation.
!
!     IWORK()       The first N locations contain the order in
!                   which the columns of A were used. The next
!                   M locations contain the order in which the
!                   rows of A were used.
!
!     INFO          Flag to indicate status of computation on completion
!                  -1   Parameter error(s)
!                   0 - Rank deficient, no solution
!                   1 - Rank deficient, truncated solution
!                   2 - Rank deficient, minimal length least squares sol
!                   3 - Numerical rank 0, zero solution
!                   4 - Rank  <  NP
!                   5 - Full rank
!
!***REFERENCES  T. Manteuffel, An interval analysis approach to rank
!                 determination in linear least squares problems,
!                 Report SAND80-0655, Sandia Laboratories, June 1980.
!***ROUTINES CALLED  D1MACH, DU11US, DU12US, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810801  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900510  Fixed an error message.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DULSIA
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DOUBLE PRECISION D1MACH
  DIMENSION A(MDA,*),B(MDB,*),RE(*),AE(*),RNORM(*),W(*)
  INTEGER IWORK(*)
!
!***FIRST EXECUTABLE STATEMENT  DULSIA
  if ( INFO < 0 .OR. INFO > 1) go to 514
  IT=INFO
  INFO=-1
  if ( NB == 0 .AND. IT == 1) go to 501
  if ( M < 1) go to 502
  if ( N < 1) go to 503
  if ( N < M) go to 504
  if ( MDA < M) go to 505
  if ( LIW < M+N) go to 506
  if ( MODE < 0 .OR. MODE > 3) go to 515
  if ( NB == 0) go to 4
  if ( NB < 0) go to 507
  if ( MDB < N) go to 508
  if ( IT == 0) go to 4
  go to 400
    4 if ( KEY < 0.OR.KEY > 3) go to 509
  if ( KEY == 0 .AND. LW < 5*M) go to 510
  if ( KEY == 1 .AND. LW < 4*M) go to 510
  if ( KEY == 2 .AND. LW < 4*M) go to 510
  if ( KEY == 3 .AND. LW < 3*M) go to 510
  if ( NP < 0 .OR. NP > M) go to 516
!
  EPS=10.*D1MACH(3)
  M1=1
  M2=M1+M
  M3=M2+M
  M4=M3+M
  M5=M4+M
!
  if ( KEY == 1) go to 100
  if ( KEY == 2) go to 200
  if ( KEY == 3) go to 300
!
  if ( RE(1) < 0.D00) go to 511
  if ( RE(1) > 1.0D0) go to 512
  if ( RE(1) < EPS) RE(1)=EPS
  if ( AE(1) < 0.0D0) go to 513
  DO 20 I=1,M
  W(M4-1+I)=RE(1)
  W(M5-1+I)=AE(1)
   20 CONTINUE
  call DU11US(A,MDA,M,N,W(M4),W(M5),MODE,NP,KRANK,KSURE, &
              W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
  go to 400
!
  100 CONTINUE
  if ( AE(1) < 0.0D0) go to 513
  DO 120 I=1,M
  if ( RE(I) < 0.0D0) go to 511
  if ( RE(I) > 1.0D0) go to 512
  if ( RE(I) < EPS) RE(I)=EPS
  W(M4-1+I)=AE(1)
  120 CONTINUE
  call DU11US(A,MDA,M,N,RE,W(M4),MODE,NP,KRANK,KSURE, &
              W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
  go to 400
!
  200 CONTINUE
  if ( RE(1) < 0.0D0) go to 511
  if ( RE(1) > 1.0D0) go to 512
  if ( RE(1) < EPS) RE(1)=EPS
  DO 220 I=1,M
  W(M4-1+I)=RE(1)
  if ( AE(I) < 0.0D0) go to 513
  220 CONTINUE
  call DU11US(A,MDA,M,N,W(M4),AE,MODE,NP,KRANK,KSURE, &
              W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
  go to 400
!
  300 CONTINUE
  DO 320 I=1,M
  if ( RE(I) < 0.0D0) go to 511
  if ( RE(I) > 1.0D0) go to 512
  if ( RE(I) < EPS) RE(I)=EPS
  if ( AE(I) < 0.0D0) go to 513
  320 CONTINUE
  call DU11US(A,MDA,M,N,RE,AE,MODE,NP,KRANK,KSURE, &
              W(M1),W(M2),W(M3),IWORK(M1),IWORK(M2))
!
!     DETERMINE INFO
!
  400 if ( KRANK /= M) go to 402
      INFO=5
      go to 410
  402 if ( KRANK /= 0) go to 404
      INFO=3
      go to 410
  404 if ( KRANK >= NP) go to 406
      INFO=4
      return
  406 INFO=MODE
  if ( MODE == 0) RETURN
  410 if ( NB == 0) RETURN
!
!
!     SOLUTION PHASE
!
  M1=1
  M2=M1+M
  M3=M2+M
  if ( INFO == 2) go to 420
  if ( LW < M2-1) go to 510
  call DU12US(A,MDA,M,N,B,MDB,NB,MODE,KRANK, &
              RNORM,W(M1),W(M1),IWORK(M1),IWORK(M2))
  return
!
  420 if ( LW < M3-1) go to 510
  call DU12US(A,MDA,M,N,B,MDB,NB,MODE,KRANK, &
              RNORM,W(M1),W(M2),IWORK(M1),IWORK(M2))
  return
!
!     ERROR MESSAGES
!
  501 call XERMSG ('SLATEC', 'DULSIA', &
     'SOLUTION ONLY (INFO=1) BUT NO RIGHT HAND SIDE (NB=0)', 1, 0)
  return
  502 call XERMSG ('SLATEC', 'DULSIA', 'M < 1', 2, 1)
  return
  503 call XERMSG ('SLATEC', 'DULSIA', 'N < 1', 2, 1)
  return
  504 call XERMSG ('SLATEC', 'DULSIA', 'N < M', 2, 1)
  return
  505 call XERMSG ('SLATEC', 'DULSIA', 'MDA < M', 2, 1)
  return
  506 call XERMSG ('SLATEC', 'DULSIA', 'LIW < M+N', 2, 1)
  return
  507 call XERMSG ('SLATEC', 'DULSIA', 'NB < 0', 2, 1)
  return
  508 call XERMSG ('SLATEC', 'DULSIA', 'MDB < N', 2, 1)
  return
  509 call XERMSG ('SLATEC', 'DULSIA', 'KEY OUT OF RANGE', 2, 1)
  return
  510 call XERMSG ('SLATEC', 'DULSIA', 'INSUFFICIENT WORK SPACE', 8, 1)
  INFO=-1
  return
  511 call XERMSG ('SLATEC', 'DULSIA', 'RE(I)  <  0', 2, 1)
  return
  512 call XERMSG ('SLATEC', 'DULSIA', 'RE(I)  >  1', 2, 1)
  return
  513 call XERMSG ('SLATEC', 'DULSIA', 'AE(I)  <  0', 2, 1)
  return
  514 call XERMSG ('SLATEC', 'DULSIA', 'INFO OUT OF RANGE', 2, 1)
  return
  515 call XERMSG ('SLATEC', 'DULSIA', 'MODE OUT OF RANGE', 2, 1)
  return
  516 call XERMSG ('SLATEC', 'DULSIA', 'NP OUT OF RANGE', 2, 1)
  return
end