subroutine DBNDSL (MODE, G, MDG, NB, IP, IR, X, N, RNORM)
!
!! DBNDSL solves the least squares problem for a banded matrix using ...
!  sequential accumulation of rows of the data matrix.
!  Exactly one right-hand side vector is permitted.
!
!***LIBRARY   SLATEC
!***CATEGORY  D9
!***TYPE      DOUBLE PRECISION (BNDSOL-S, DBNDSL-D)
!***KEYWORDS  BANDED MATRIX, CURVE FITTING, LEAST SQUARES
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     These subroutines solve the least squares problem Ax = b for
!     banded matrices A using sequential accumulation of rows of the
!     data matrix.  Exactly one right-hand side vector is permitted.
!
!     These subroutines are intended for the type of least squares
!     systems that arise in applications such as curve or surface
!     fitting of data.  The least squares equations are accumulated and
!     processed using only part of the data.  This requires a certain
!     user interaction during the solution of Ax = b.
!
!     Specifically, suppose the data matrix (A B) is row partitioned
!     into Q submatrices.  Let (E F) be the T-th one of these
!     submatrices where E = (0 C 0).  Here the dimension of E is MT by N
!     and the dimension of C is MT by NB.  The value of NB is the
!     bandwidth of A.  The dimensions of the leading block of zeros in E
!     are MT by JT-1.
!
!     The user of the subroutine DBNDAC provides MT,JT,C and F for
!     T=1,...,Q.  Not all of this data must be supplied at once.
!
!     Following the processing of the various blocks (E F), the matrix
!     (A B) has been transformed to the form (R D) where R is upper
!     triangular and banded with bandwidth NB.  The least squares
!     system Rx = d is then easily solved using back substitution by
!     executing the statement call DBNDSL(1,...). The sequence of
!     values for JT must be nondecreasing.  This may require some
!     preliminary interchanges of rows and columns of the matrix A.
!
!     The primary reason for these subroutines is that the total
!     processing can take place in a working array of dimension MU by
!     NB+1.  An acceptable value for MU is
!
!                       MU = MAX(MT + N + 1),
!
!     where N is the number of unknowns.
!
!     Here the maximum is taken over all values of MT for T=1,...,Q.
!     Notice that MT can be taken to be a small as one, showing that
!     MU can be as small as N+2.  The subprogram DBNDAC processes the
!     rows more efficiently if MU is large enough so that each new
!     block (C F) has a distinct value of JT.
!
!     The four principle parts of these algorithms are obtained by the
!     following call statements
!
!     call DBNDAC(...)  Introduce new blocks of data.
!
!     call DBNDSL(1,...)Compute solution vector and length of
!                       residual vector.
!
!     call DBNDSL(2,...)Given any row vector H solve YR = H for the
!                       row vector Y.
!
!     call DBNDSL(3,...)Given any column vector W solve RZ = W for
!                       the column vector Z.
!
!     The dots in the above call statements indicate additional
!     arguments that will be specified in the following paragraphs.
!
!     The user must dimension the array appearing in the call list..
!     G(MDG,NB+1)
!
!     Description of calling sequence for DBNDAC..
!
!     The entire set of parameters for DBNDAC are
!
!     Input.. All Type REAL variables are DOUBLE PRECISION
!
!     G(*,*)            The working array into which the user will
!                       place the MT by NB+1 block (C F) in rows IR
!                       through IR+MT-1, columns 1 through NB+1.
!                       See descriptions of IR and MT below.
!
!     MDG               The number of rows in the working array
!                       G(*,*).  The value of MDG should be  >=  MU.
!                       The value of MU is defined in the abstract
!                       of these subprograms.
!
!     NB                The bandwidth of the data matrix A.
!
!     IP                Set by the user to the value 1 before the
!                       first call to DBNDAC.  Its subsequent value
!                       is controlled by DBNDAC to set up for the
!                       next call to DBNDAC.
!
!     IR                Index of the row of G(*,*) where the user is
!                       the user to the value 1 before the first call
!                       to DBNDAC.  Its subsequent value is controlled
!                       by DBNDAC. A value of IR  >  MDG is considered
!                       an error.
!
!     MT,JT             Set by the user to indicate respectively the
!                       number of new rows of data in the block and
!                       the index of the first nonzero column in that
!                       set of rows (E F) = (0 C 0 F) being processed.
!     Output.. All Type REAL variables are DOUBLE PRECISION
!
!     G(*,*)            The working array which will contain the
!                       processed rows of that part of the data
!                       matrix which has been passed to DBNDAC.
!
!     IP,IR             The values of these arguments are advanced by
!                       DBNDAC to be ready for storing and processing
!                       a new block of data in G(*,*).
!
!     Description of calling sequence for DBNDSL..
!
!     The user must dimension the arrays appearing in the call list..
!
!     G(MDG,NB+1), X(N)
!
!     The entire set of parameters for DBNDSL are
!
!     Input..
!
!     MODE              Set by the user to one of the values 1, 2, or
!                       3.  These values respectively indicate that
!                       the solution of AX = B, YR = H or RZ = W is
!                       required.
!
!     G(*,*),MDG,       These arguments all have the same meaning and
!      NB,IP,IR         contents as following the last call to DBNDAC.
!
!     X(*)              With mode=2 or 3 this array contains,
!                       respectively, the right-side vectors H or W of
!                       the systems YR = H or RZ = W.
!
!     N                 The number of variables in the solution
!                       vector.  If any of the N diagonal terms are
!                       zero the subroutine DBNDSL prints an
!                       appropriate message.  This condition is
!                       considered an error.
!
!     Output..
!
!     X(*)              This array contains the solution vectors X,
!                       Y or Z of the systems AX = B, YR = H or
!                       RZ = W depending on the value of MODE=1,
!                       2 or 3.
!
!     RNORM             If MODE=1 RNORM is the Euclidean length of the
!                       residual vector AX-B.  When MODE=2 or 3 RNORM
!                       is set to zero.
!
!     Remarks..
!
!     To obtain the upper triangular matrix and transformed right-hand
!     side vector D so that the super diagonals of R form the columns
!     of G(*,*), execute the following Fortran statements.
!
!     NBP1=NB+1
!
!     DO 10 J=1, NBP1
!
!  10 G(IR,J) = 0.E0
!
!     MT=1
!
!     JT=N+1
!
!     call DBNDAC(G,MDG,NB,IP,IR,MT,JT)
!
!***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
!                 Problems, Prentice-Hall, Inc., 1974, Chapter 27.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBNDSL
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DIMENSION G(MDG,*),X(*)
!***FIRST EXECUTABLE STATEMENT  DBNDSL
  ZERO=0.D0
!
  RNORM=ZERO
  go to (10,90,50), MODE
!                                   ********************* MODE = 1
!                                   ALG. STEP 26
   10      DO 20 J=1,N
       X(J)=G(J,NB+1)
   20 CONTINUE
  RSQ=ZERO
  NP1=N+1
  IRM1=IR-1
  if (NP1 > IRM1) go to 40
       DO 30 J=NP1,IRM1
       RSQ=RSQ+G(J,NB+1)**2
   30 CONTINUE
  RNORM=SQRT(RSQ)
   40 CONTINUE
!                                   ********************* MODE = 3
!                                   ALG. STEP 27
   50      DO 80 II=1,N
       I=N+1-II
!                                   ALG. STEP 28
       S=ZERO
       L=MAX(0,I-IP)
!                                   ALG. STEP 29
       if (I == N) go to 70
!                                   ALG. STEP 30
       IE=MIN(N+1-I,NB)
            DO 60 J=2,IE
            JG=J+L
            IX=I-1+J
            S=S+G(I,JG)*X(IX)
   60 CONTINUE
!                                   ALG. STEP 31
   70      if (G(I,L+1)) 80,130,80
   80      X(I)=(X(I)-S)/G(I,L+1)
!                                   ALG. STEP 32
  return
!                                   ********************* MODE = 2
   90      DO 120 J=1,N
       S=ZERO
       if (J == 1) go to 110
       I1=MAX(1,J-NB+1)
       I2=J-1
            DO 100 I=I1,I2
            L=J-I+1+MAX(0,I-IP)
            S=S+X(I)*G(I,L)
  100 CONTINUE
  110      L=MAX(0,J-IP)
       if (G(J,L+1)) 120,130,120
  120      X(J)=(X(J)-S)/G(J,L+1)
  return
!
  130 CONTINUE
  NERR=1
  IOPT=2
  call XERMSG ('SLATEC', 'DBNDSL', &
     'A ZERO DIAGONAL TERM IS IN THE N BY N UPPER TRIANGULAR ' // &
     'MATRIX.', NERR, IOPT)
  return
end
