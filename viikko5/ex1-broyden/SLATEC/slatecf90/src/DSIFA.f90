subroutine DSIFA (A, LDA, N, KPVT, INFO)
!
!! DSIFA factors a real symmetric matrix by elimination with symmetric pivoting.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      DOUBLE PRECISION (SSIFA-S, DSIFA-D, CHIFA-C, CSIFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     DSIFA factors a double precision symmetric matrix by elimination
!     with symmetric pivoting.
!
!     To solve  A*X = B , follow DSIFA by DSISL.
!     To compute  INVERSE(A)*C , follow DSIFA by DSISL.
!     To compute  DETERMINANT(A) , follow DSIFA by DSIDI.
!     To compute  INERTIA(A) , follow DSIFA by DSIDI.
!     To compute  INVERSE(A) , follow DSIFA by DSIDI.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA,N)
!                the symmetric matrix to be factored.
!                Only the diagonal and upper triangle are used.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        A       a block diagonal matrix and the multipliers which
!                were used to obtain it.
!                The factorization can be written  A = U*D*TRANS(U)
!                where  U  is a product of permutation and unit
!                upper triangular matrices, TRANS(U) is the
!                transpose of  U , and  D  is block diagonal
!                with 1 by 1 and 2 by 2 blocks.
!
!        KPVT    INTEGER(N)
!                an integer vector of pivot indices.
!
!        INFO    INTEGER
!                = 0  normal value.
!                = K  if the K-th pivot block is singular.  This is
!                     not an error condition for this subroutine,
!                     but it does indicate that DSISL or DSIDI may
!                     divide by zero if called.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSWAP, IDAMAX
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Modified routine equivalence list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSIFA
  INTEGER LDA,N,KPVT(*),INFO
  DOUBLE PRECISION A(LDA,*)
!
  DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
  DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
  INTEGER IMAX,IMAXP1,J,JJ,JMAX,K,KM1,KM2,KSTEP,IDAMAX
  LOGICAL SWAP
!***FIRST EXECUTABLE STATEMENT  DSIFA
!
!     INITIALIZE
!
!     ALPHA IS USED IN CHOOSING PIVOT BLOCK SIZE.
!
  ALPHA = (1.0D0 + SQRT(17.0D0))/8.0D0
!
  INFO = 0
!
!     MAIN LOOP ON K, WHICH GOES FROM N TO 1.
!
  K = N
   10 CONTINUE
!
!        LEAVE THE LOOP if K=0 OR K=1.
!
     if (K  ==  0) go to 200
     if (K  >  1) go to 20
        KPVT(1) = 1
        if (A(1,1)  ==  0.0D0) INFO = 1
        go to 200
   20    CONTINUE
!
!        THIS SECTION OF CODE DETERMINES THE KIND OF
!        ELIMINATION TO BE PERFORMED.  WHEN IT IS COMPLETED,
!        KSTEP WILL BE SET TO THE SIZE OF THE PIVOT BLOCK, AND
!        SWAP WILL BE SET TO .TRUE. if AN INTERCHANGE IS
!        REQUIRED.
!
     KM1 = K - 1
     ABSAKK = ABS(A(K,K))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
     IMAX = IDAMAX(K-1,A(1,K),1)
     COLMAX = ABS(A(IMAX,K))
     if (ABSAKK  <  ALPHA*COLMAX) go to 30
        KSTEP = 1
        SWAP = .FALSE.
     go to 90
   30    CONTINUE
!
!           DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!           ROW IMAX.
!
        ROWMAX = 0.0D0
        IMAXP1 = IMAX + 1
        DO 40 J = IMAXP1, K
           ROWMAX = MAX(ROWMAX,ABS(A(IMAX,J)))
   40       CONTINUE
        if (IMAX  ==  1) go to 50
           JMAX = IDAMAX(IMAX-1,A(1,IMAX),1)
           ROWMAX = MAX(ROWMAX,ABS(A(JMAX,IMAX)))
   50       CONTINUE
        if (ABS(A(IMAX,IMAX))  <  ALPHA*ROWMAX) go to 60
           KSTEP = 1
           SWAP = .TRUE.
        go to 80
   60       CONTINUE
        if (ABSAKK  <  ALPHA*COLMAX*(COLMAX/ROWMAX)) go to 70
           KSTEP = 1
           SWAP = .FALSE.
        go to 80
   70       CONTINUE
           KSTEP = 2
           SWAP = IMAX  /=  KM1
   80       CONTINUE
   90    CONTINUE
     if (MAX(ABSAKK,COLMAX)  /=  0.0D0) go to 100
!
!  COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
        KPVT(K) = K
        INFO = K
     go to 190
  100    CONTINUE
     if (KSTEP  ==  2) go to 140
!
!  1 X 1 PIVOT BLOCK.
!
        if (.NOT.SWAP) go to 120
!
!  PERFORM AN INTERCHANGE.
!
           call DSWAP(IMAX,A(1,IMAX),1,A(1,K),1)
           DO 110 JJ = IMAX, K
              J = K + IMAX - JJ
              T = A(J,K)
              A(J,K) = A(IMAX,J)
              A(IMAX,J) = T
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
        DO 130 JJ = 1, KM1
           J = K - JJ
           MULK = -A(J,K)/A(K,K)
           T = MULK
           call DAXPY(J,T,A(1,K),1,A(1,J),1)
           A(J,K) = MULK
  130       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
        KPVT(K) = K
        if (SWAP) KPVT(K) = IMAX
     go to 190
  140    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
        if (.NOT.SWAP) go to 160
!
!              PERFORM AN INTERCHANGE.
!
           call DSWAP(IMAX,A(1,IMAX),1,A(1,K-1),1)
           DO 150 JJ = IMAX, KM1
              J = KM1 + IMAX - JJ
              T = A(J,K-1)
              A(J,K-1) = A(IMAX,J)
              A(IMAX,J) = T
  150          CONTINUE
           T = A(K-1,K)
           A(K-1,K) = A(IMAX,K)
           A(IMAX,K) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
        KM2 = K - 2
        if (KM2  ==  0) go to 180
           AK = A(K,K)/A(K-1,K)
           AKM1 = A(K-1,K-1)/A(K-1,K)
           DENOM = 1.0D0 - AK*AKM1
           DO 170 JJ = 1, KM2
              J = KM1 - JJ
              BK = A(J,K)/A(K-1,K)
              BKM1 = A(J,K-1)/A(K-1,K)
              MULK = (AKM1*BK - BKM1)/DENOM
              MULKM1 = (AK*BKM1 - BK)/DENOM
              T = MULK
              call DAXPY(J,T,A(1,K),1,A(1,J),1)
              T = MULKM1
              call DAXPY(J,T,A(1,K-1),1,A(1,J),1)
              A(J,K) = MULK
              A(J,K-1) = MULKM1
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
        KPVT(K) = 1 - K
        if (SWAP) KPVT(K) = -IMAX
        KPVT(K-1) = KPVT(K)
  190    CONTINUE
     K = K - KSTEP
  go to 10
  200 CONTINUE
  return
end
