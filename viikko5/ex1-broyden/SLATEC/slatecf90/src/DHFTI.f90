subroutine DHFTI (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H, &
     G, IP)
!
!! DHFTI solves a least squares problem for banded matrices using ...
!            sequential accumulation of rows of the data matrix.
!            Exactly one right-hand side vector is permitted.
!
!***LIBRARY   SLATEC
!***CATEGORY  D9
!***TYPE      DOUBLE PRECISION (HFTI-S, DHFTI-D)
!***KEYWORDS  CURVE FITTING, LEAST SQUARES
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!***DESCRIPTION
!
!     DIMENSION A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N)
!
!     This subroutine solves a linear least squares problem or a set of
!     linear least squares problems having the same matrix but different
!     right-side vectors.  The problem data consists of an M by N matrix
!     A, an M by NB matrix B, and an absolute tolerance parameter TAU
!     whose usage is described below.  The NB column vectors of B
!     represent right-side vectors for NB distinct linear least squares
!     problems.
!
!     This set of problems can also be written as the matrix least
!     squares problem
!
!                       AX = B,
!
!     where X is the N by NB solution matrix.
!
!     Note that if B is the M by M identity matrix, then X will be the
!     pseudo-inverse of A.
!
!     This subroutine first transforms the augmented matrix (A B) to a
!     matrix (R C) using premultiplying Householder transformations with
!     column interchanges.  All subdiagonal elements in the matrix R are
!     zero and its diagonal elements satisfy
!
!                       ABS(R(I,I)) >= ABS(R(I+1,I+1)),
!
!                       I = 1,...,L-1, where
!
!                       L = MIN(M,N).
!
!     The subroutine will compute an integer, KRANK, equal to the number
!     of diagonal terms of R that exceed TAU in magnitude. Then a
!     solution of minimum Euclidean length is computed using the first
!     KRANK rows of (R C).
!
!     To be specific we suggest that the user consider an easily
!     computable matrix norm, such as, the maximum of all column sums of
!     magnitudes.
!
!     Now if the relative uncertainty of B is EPS, (norm of uncertainty/
!     norm of B), it is suggested that TAU be set approximately equal to
!     EPS*(norm of A).
!
!     The user must dimension all arrays appearing in the call list..
!     A(MDA,N),(B(MDB,NB) or B(M)),RNORM(NB),H(N),G(N),IP(N).  This
!     permits the solution of a range of problems in the same array
!     space.
!
!     The entire set of parameters for DHFTI are
!
!     INPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     A(*,*),MDA,M,N    The array A(*,*) initially contains the M by N
!                       matrix A of the least squares problem AX = B.
!                       The first dimensioning parameter of the array
!                       A(*,*) is MDA, which must satisfy MDA >= M
!                       Either M >= N or M < N is permitted.  There
!                       is no restriction on the rank of A.  The
!                       condition MDA < M is considered an error.
!
!     B(*),MDB,NB       If NB = 0 the subroutine will perform the
!                       orthogonal decomposition but will make no
!                       references to the array B(*).  If NB > 0
!                       the array B(*) must initially contain the M by
!                       NB matrix B of the least squares problem AX =
!                       B.  If NB >= 2 the array B(*) must be doubly
!                       subscripted with first dimensioning parameter
!                       MDB >= MAX(M,N).  If NB = 1 the array B(*) may
!                       be either doubly or singly subscripted.  In
!                       the latter case the value of MDB is arbitrary
!                       but it should be set to some valid integer
!                       value such as MDB = M.
!
!                       The condition of NB > 1.AND.MDB <  MAX(M,N)
!                       is considered an error.
!
!     TAU               Absolute tolerance parameter provided by user
!                       for pseudorank determination.
!
!     H(*),G(*),IP(*)   Arrays of working space used by DHFTI.
!
!     OUTPUT.. All TYPE REAL variables are DOUBLE PRECISION
!
!     A(*,*)            The contents of the array A(*,*) will be
!                       modified by the subroutine. These contents
!                       are not generally required by the user.
!
!     B(*)              On return the array B(*) will contain the N by
!                       NB solution matrix X.
!
!     KRANK             Set by the subroutine to indicate the
!                       pseudorank of A.
!
!     RNORM(*)          On return, RNORM(J) will contain the Euclidean
!                       norm of the residual vector for the problem
!                       defined by the J-th column vector of the array
!                       B(*,*) for J = 1,...,NB.
!
!     H(*),G(*)         On return these arrays respectively contain
!                       elements of the pre- and post-multiplying
!                       Householder transformations used to compute
!                       the minimum Euclidean length solution.
!
!     IP(*)             Array in which the subroutine records indices
!                       describing the permutation of column vectors.
!                       The contents of arrays H(*),G(*) and IP(*)
!                       are not generally required by the user.
!
!***REFERENCES  C. L. Lawson and R. J. Hanson, Solving Least Squares
!                 Problems, Prentice-Hall, Inc., 1974, Chapter 14.
!***ROUTINES CALLED  D1MACH, DH12, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901005  Replace usage of DDIFF with usage of D1MACH.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DHFTI
  INTEGER I, II, IOPT, IP(*), IP1, J, JB, JJ, K, KP1, KRANK, L, &
       LDIAG, LMAX, M, MDA, MDB, N, NB, NERR
  DOUBLE PRECISION A, B, D1MACH, DZERO, FACTOR, &
       G, H, HMAX, RELEPS, RNORM, SM, SM1, SZERO, TAU, TMP
  DIMENSION A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*)
  SAVE RELEPS
  DATA RELEPS /0.D0/
!     BEGIN BLOCK PERMITTING ...EXITS TO 360
!***FIRST EXECUTABLE STATEMENT  DHFTI
     if (RELEPS == 0.D0) RELEPS = D1MACH(4)
     SZERO = 0.0D0
     DZERO = 0.0D0
     FACTOR = 0.001D0
!
     K = 0
     LDIAG = MIN(M,N)
     if (LDIAG  <=  0) go to 350
!           BEGIN BLOCK PERMITTING ...EXITS TO 130
!              BEGIN BLOCK PERMITTING ...EXITS TO 120
              if (MDA  >=  M) go to 10
                 NERR = 1
                 IOPT = 2
                 call XERMSG ('SLATEC', 'DHFTI', &
                    'MDA < M, PROBABLE ERROR.', &
                    NERR, IOPT)
!     ...............EXIT
                 go to 360
   10             CONTINUE
!
              if (NB  <=  1 .OR. MAX(M,N)  <=  MDB) go to 20
                 NERR = 2
                 IOPT = 2
                 call XERMSG ('SLATEC', 'DHFTI', &
                    'MDB < MAX(M,N).AND.NB > 1. PROBABLE ERROR.', &
                    NERR, IOPT)
!     ...............EXIT
                 go to 360
   20             CONTINUE
!
              DO 100 J = 1, LDIAG
!                    BEGIN BLOCK PERMITTING ...EXITS TO 70
                    if (J  ==  1) go to 40
!
!                           UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
!                          ..
                       LMAX = J
                       DO 30 L = J, N
                          H(L) = H(L) - A(J-1,L)**2
                          if (H(L)  >  H(LMAX)) LMAX = L
   30                      CONTINUE
!                    ......EXIT
                       if (FACTOR*H(LMAX)  >  HMAX*RELEPS) go to 70
   40                   CONTINUE
!
!                        COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
!                       ..
                    LMAX = J
                    DO 60 L = J, N
                       H(L) = 0.0D0
                       DO 50 I = J, M
                          H(L) = H(L) + A(I,L)**2
   50                      CONTINUE
                       if (H(L)  >  H(LMAX)) LMAX = L
   60                   CONTINUE
                    HMAX = H(LMAX)
   70                CONTINUE
!                    ..
!                     LMAX HAS BEEN DETERMINED
!
!                     DO COLUMN INTERCHANGES if NEEDED.
!                    ..
                 IP(J) = LMAX
                 if (IP(J)  ==  J) go to 90
                    DO 80 I = 1, M
                       TMP = A(I,J)
                       A(I,J) = A(I,LMAX)
                       A(I,LMAX) = TMP
   80                   CONTINUE
                    H(LMAX) = H(J)
   90                CONTINUE
!
!                     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A
!                     AND B.
!                    ..
                 call DH12(1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA, &
                           N-J)
                 call DH12(2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
  100             CONTINUE
!
!                  DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE,
!                  TAU.
!                 ..
              DO 110 J = 1, LDIAG
!              ......EXIT
                 if (ABS(A(J,J))  <=  TAU) go to 120
  110             CONTINUE
              K = LDIAG
!           ......EXIT
              go to 130
  120          CONTINUE
           K = J - 1
  130       CONTINUE
        KP1 = K + 1
!
!           COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
!
        if (NB  <  1) go to 170
        DO 160 JB = 1, NB
           TMP = SZERO
           if (M  <  KP1) go to 150
           DO 140 I = KP1, M
              TMP = TMP + B(I,JB)**2
  140          CONTINUE
  150          CONTINUE
           RNORM(JB) = SQRT(TMP)
  160       CONTINUE
  170       CONTINUE
!           SPECIAL FOR PSEUDORANK = 0
        if (K  >  0) go to 210
           if (NB  <  1) go to 200
           DO 190 JB = 1, NB
              DO 180 I = 1, N
                 B(I,JB) = SZERO
  180             CONTINUE
  190          CONTINUE
  200          CONTINUE
        go to 340
  210       CONTINUE
!
!               if THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
!               DECOMPOSITION OF FIRST K ROWS.
!              ..
           if (K  ==  N) go to 230
              DO 220 II = 1, K
                 I = KP1 - II
                 call DH12(1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  220             CONTINUE
  230          CONTINUE
!
!
           if (NB  <  1) go to 330
           DO 320 JB = 1, NB
!
!                  SOLVE THE K BY K TRIANGULAR SYSTEM.
!                 ..
              DO 260 L = 1, K
                 SM = DZERO
                 I = KP1 - L
                 IP1 = I + 1
                 if (K  <  IP1) go to 250
                 DO 240 J = IP1, K
                    SM = SM + A(I,J)*B(J,JB)
  240                CONTINUE
  250                CONTINUE
                 SM1 = SM
                 B(I,JB) = (B(I,JB) - SM1)/A(I,I)
  260             CONTINUE
!
!                  COMPLETE COMPUTATION OF SOLUTION VECTOR.
!                 ..
              if (K  ==  N) go to 290
                 DO 270 J = KP1, N
                    B(J,JB) = SZERO
  270                CONTINUE
                 DO 280 I = 1, K
                    call DH12(2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1, &
                              MDB,1)
  280                CONTINUE
  290             CONTINUE
!
!                   RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
!                   COLUMN INTERCHANGES.
!                 ..
              DO 310 JJ = 1, LDIAG
                 J = LDIAG + 1 - JJ
                 if (IP(J)  ==  J) go to 300
                    L = IP(J)
                    TMP = B(L,JB)
                    B(L,JB) = B(J,JB)
                    B(J,JB) = TMP
  300                CONTINUE
  310             CONTINUE
  320          CONTINUE
  330          CONTINUE
  340       CONTINUE
  350    CONTINUE
!        ..
!         THE SOLUTION VECTORS, X, ARE NOW
!         IN THE FIRST  N  ROWS OF THE ARRAY B(,).
!
     KRANK = K
  360 CONTINUE
  return
end
