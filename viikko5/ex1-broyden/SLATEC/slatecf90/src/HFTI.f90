subroutine HFTI (A, MDA, M, N, B, MDB, NB, TAU, KRANK, RNORM, H, &
     G, IP)
!
!! HFTI solves a linear least squares problems by performing a QR ...
!  factorization of the matrix using Householder transformations.
!
!***LIBRARY   SLATEC
!***CATEGORY  D9
!***TYPE      SINGLE PRECISION (HFTI-S, DHFTI-D)
!***KEYWORDS  CURVE FITTING, LINEAR LEAST SQUARES, QR FACTORIZATION
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
!     The entire set of parameters for HFTI are
!
!     INPUT..
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
!     H(*),G(*),IP(*)   Arrays of working space used by HFTI.
!
!     OUTPUT..
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
!***ROUTINES CALLED  H12, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   901005  Replace usage of DIFF with usage of R1MACH.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HFTI
  DIMENSION A(MDA,*),B(MDB,*),H(*),G(*),RNORM(*)
  INTEGER IP(*)
  DOUBLE PRECISION SM,DZERO
  SAVE RELEPS
  DATA RELEPS /0.E0/
!***FIRST EXECUTABLE STATEMENT  HFTI
  if (RELEPS == 0) RELEPS = R1MACH(4)
  SZERO=0.
  DZERO=0.D0
  FACTOR=0.001
!
  K=0
  LDIAG=MIN(M,N)
  if (LDIAG <= 0) go to 270
  if (.NOT.MDA < M) go to 5
  NERR=1
  IOPT=2
  call XERMSG ('SLATEC', 'HFTI', 'MDA < M, PROBABLE ERROR.', &
     NERR, IOPT)
  return
    5 CONTINUE
!
  if (.NOT.(NB > 1.AND.MAX(M,N) > MDB)) go to 6
  NERR=2
  IOPT=2
  call XERMSG ('SLATEC', 'HFTI', &
     'MDB < MAX(M,N).AND.NB > 1. PROBABLE ERROR.', NERR, IOPT)
  return
    6 CONTINUE
!
      DO 80 J=1,LDIAG
      if (J == 1) go to 20
!
!     UPDATE SQUARED COLUMN LENGTHS AND FIND LMAX
!    ..
      LMAX=J
          DO 10 L=J,N
          H(L)=H(L)-A(J-1,L)**2
          if (H(L) > H(LMAX)) LMAX=L
   10         CONTINUE
      if (FACTOR*H(LMAX)  >  HMAX*RELEPS) go to 50
!
!     COMPUTE SQUARED COLUMN LENGTHS AND FIND LMAX
!    ..
   20     LMAX=J
          DO 40 L=J,N
          H(L)=0.
              DO 30 I=J,M
   30             H(L)=H(L)+A(I,L)**2
          if (H(L) > H(LMAX)) LMAX=L
   40         CONTINUE
      HMAX=H(LMAX)
!    ..
!     LMAX HAS BEEN DETERMINED
!
!     DO COLUMN INTERCHANGES if NEEDED.
!    ..
   50     CONTINUE
      IP(J)=LMAX
      if (IP(J) == J) go to 70
          DO 60 I=1,M
          TMP=A(I,J)
          A(I,J)=A(I,LMAX)
   60         A(I,LMAX)=TMP
      H(LMAX)=H(J)
!
!     COMPUTE THE J-TH TRANSFORMATION AND APPLY IT TO A AND B.
!    ..
   70     call H12 (1,J,J+1,M,A(1,J),1,H(J),A(1,J+1),1,MDA,N-J)
   80     call H12 (2,J,J+1,M,A(1,J),1,H(J),B,1,MDB,NB)
!
!     DETERMINE THE PSEUDORANK, K, USING THE TOLERANCE, TAU.
!    ..
      DO 90 J=1,LDIAG
      if (ABS(A(J,J)) <= TAU) go to 100
   90     CONTINUE
  K=LDIAG
  go to 110
  100 K=J-1
  110 KP1=K+1
!
!     COMPUTE THE NORMS OF THE RESIDUAL VECTORS.
!
  if (NB <= 0) go to 140
      DO 130 JB=1,NB
      TMP=SZERO
      if (KP1 > M) go to 130
          DO 120 I=KP1,M
  120         TMP=TMP+B(I,JB)**2
  130     RNORM(JB)=SQRT(TMP)
  140 CONTINUE
!                                           SPECIAL FOR PSEUDORANK = 0
  if (K > 0) go to 160
  if (NB <= 0) go to 270
      DO 150 JB=1,NB
          DO 150 I=1,N
  150         B(I,JB)=SZERO
  go to 270
!
!     if THE PSEUDORANK IS LESS THAN N COMPUTE HOUSEHOLDER
!     DECOMPOSITION OF FIRST K ROWS.
!    ..
  160 if (K == N) go to 180
      DO 170 II=1,K
      I=KP1-II
  170     call H12 (1,I,KP1,N,A(I,1),MDA,G(I),A,MDA,1,I-1)
  180 CONTINUE
!
!
  if (NB <= 0) go to 270
      DO 260 JB=1,NB
!
!     SOLVE THE K BY K TRIANGULAR SYSTEM.
!    ..
          DO 210 L=1,K
          SM=DZERO
          I=KP1-L
          if (I == K) go to 200
          IP1=I+1
              DO 190 J=IP1,K
  190             SM=SM+A(I,J)*DBLE(B(J,JB))
  200         SM1=SM
  210         B(I,JB)=(B(I,JB)-SM1)/A(I,I)
!
!     COMPLETE COMPUTATION OF SOLUTION VECTOR.
!    ..
      if (K == N) go to 240
          DO 220 J=KP1,N
  220         B(J,JB)=SZERO
          DO 230 I=1,K
  230         call H12 (2,I,KP1,N,A(I,1),MDA,G(I),B(1,JB),1,MDB,1)
!
!      RE-ORDER THE SOLUTION VECTOR TO COMPENSATE FOR THE
!      COLUMN INTERCHANGES.
!    ..
  240         DO 250 JJ=1,LDIAG
          J=LDIAG+1-JJ
          if (IP(J) == J) go to 250
          L=IP(J)
          TMP=B(L,JB)
          B(L,JB)=B(J,JB)
          B(J,JB)=TMP
  250         CONTINUE
  260     CONTINUE
!    ..
!     THE SOLUTION VECTORS, X, ARE NOW
!     IN THE FIRST  N  ROWS OF THE ARRAY B(,).
!
  270 KRANK=K
  return
end
