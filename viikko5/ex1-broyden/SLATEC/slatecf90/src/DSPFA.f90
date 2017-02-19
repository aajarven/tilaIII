subroutine DSPFA (AP, N, KPVT, INFO)
!
!! DSPFA factors a real symmetric matrix stored in packed form by ...
!            elimination with symmetric pivoting.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      DOUBLE PRECISION (SSPFA-S, DSPFA-D, CHPFA-C, CSPFA-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX FACTORIZATION, PACKED,
!             SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     DSPFA factors a double precision symmetric matrix stored in
!     packed form by elimination with symmetric pivoting.
!
!     To solve  A*X = B , follow DSPFA by DSPSL.
!     To compute  INVERSE(A)*C , follow DSPFA by DSPSL.
!     To compute  DETERMINANT(A) , follow DSPFA by DSPDI.
!     To compute  INERTIA(A) , follow DSPFA by DSPDI.
!     To compute  INVERSE(A) , follow DSPFA by DSPDI.
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the packed form of a symmetric matrix  A .  The
!                columns of the upper triangle are stored sequentially
!                in a one-dimensional array of length  N*(N+1)/2 .
!                See comments below for details.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     Output
!
!        AP      a block diagonal matrix and the multipliers which
!                were used to obtain it stored in packed form.
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
!                     but it does indicate that DSPSL or DSPDI may
!                     divide by zero if called.
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a symmetric matrix.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K)  = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
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
!***END PROLOGUE  DSPFA
  INTEGER N,KPVT(*),INFO
  DOUBLE PRECISION AP(*)
!
  DOUBLE PRECISION AK,AKM1,BK,BKM1,DENOM,MULK,MULKM1,T
  DOUBLE PRECISION ABSAKK,ALPHA,COLMAX,ROWMAX
  INTEGER IDAMAX,IJ,IK,IKM1,IM,IMAX,IMAXP1,IMIM,IMJ,IMK
  INTEGER J,JJ,JK,JKM1,JMAX,JMIM,K,KK,KM1,KM1K,KM1KM1,KM2,KSTEP
  LOGICAL SWAP
!***FIRST EXECUTABLE STATEMENT  DSPFA
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
  IK = (N*(N - 1))/2
   10 CONTINUE
!
!        LEAVE THE LOOP if K=0 OR K=1.
!
     if (K  ==  0) go to 200
     if (K  >  1) go to 20
        KPVT(1) = 1
        if (AP(1)  ==  0.0D0) INFO = 1
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
     KK = IK + K
     ABSAKK = ABS(AP(KK))
!
!        DETERMINE THE LARGEST OFF-DIAGONAL ELEMENT IN
!        COLUMN K.
!
     IMAX = IDAMAX(K-1,AP(IK+1),1)
     IMK = IK + IMAX
     COLMAX = ABS(AP(IMK))
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
        IM = IMAX*(IMAX - 1)/2
        IMJ = IM + 2*IMAX
        DO 40 J = IMAXP1, K
           ROWMAX = MAX(ROWMAX,ABS(AP(IMJ)))
           IMJ = IMJ + J
   40       CONTINUE
        if (IMAX  ==  1) go to 50
           JMAX = IDAMAX(IMAX-1,AP(IM+1),1)
           JMIM = JMAX + IM
           ROWMAX = MAX(ROWMAX,ABS(AP(JMIM)))
   50       CONTINUE
        IMIM = IMAX + IM
        if (ABS(AP(IMIM))  <  ALPHA*ROWMAX) go to 60
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
!           COLUMN K IS ZERO.  SET INFO AND ITERATE THE LOOP.
!
        KPVT(K) = K
        INFO = K
     go to 190
  100    CONTINUE
     if (KSTEP  ==  2) go to 140
!
!           1 X 1 PIVOT BLOCK.
!
        if (.NOT.SWAP) go to 120
!
!              PERFORM AN INTERCHANGE.
!
           call DSWAP(IMAX,AP(IM+1),1,AP(IK+1),1)
           IMJ = IK + IMAX
           DO 110 JJ = IMAX, K
              J = K + IMAX - JJ
              JK = IK + J
              T = AP(JK)
              AP(JK) = AP(IMJ)
              AP(IMJ) = T
              IMJ = IMJ - (J - 1)
  110          CONTINUE
  120       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
        IJ = IK - (K - 1)
        DO 130 JJ = 1, KM1
           J = K - JJ
           JK = IK + J
           MULK = -AP(JK)/AP(KK)
           T = MULK
           call DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
           AP(JK) = MULK
           IJ = IJ - (J - 1)
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
        KM1K = IK + K - 1
        IKM1 = IK - (K - 1)
        if (.NOT.SWAP) go to 160
!
!              PERFORM AN INTERCHANGE.
!
           call DSWAP(IMAX,AP(IM+1),1,AP(IKM1+1),1)
           IMJ = IKM1 + IMAX
           DO 150 JJ = IMAX, KM1
              J = KM1 + IMAX - JJ
              JKM1 = IKM1 + J
              T = AP(JKM1)
              AP(JKM1) = AP(IMJ)
              AP(IMJ) = T
              IMJ = IMJ - (J - 1)
  150          CONTINUE
           T = AP(KM1K)
           AP(KM1K) = AP(IMK)
           AP(IMK) = T
  160       CONTINUE
!
!           PERFORM THE ELIMINATION.
!
        KM2 = K - 2
        if (KM2  ==  0) go to 180
           AK = AP(KK)/AP(KM1K)
           KM1KM1 = IKM1 + K - 1
           AKM1 = AP(KM1KM1)/AP(KM1K)
           DENOM = 1.0D0 - AK*AKM1
           IJ = IK - (K - 1) - (K - 2)
           DO 170 JJ = 1, KM2
              J = KM1 - JJ
              JK = IK + J
              BK = AP(JK)/AP(KM1K)
              JKM1 = IKM1 + J
              BKM1 = AP(JKM1)/AP(KM1K)
              MULK = (AKM1*BK - BKM1)/DENOM
              MULKM1 = (AK*BKM1 - BK)/DENOM
              T = MULK
              call DAXPY(J,T,AP(IK+1),1,AP(IJ+1),1)
              T = MULKM1
              call DAXPY(J,T,AP(IKM1+1),1,AP(IJ+1),1)
              AP(JK) = MULK
              AP(JKM1) = MULKM1
              IJ = IJ - (J - 1)
  170          CONTINUE
  180       CONTINUE
!
!           SET THE PIVOT ARRAY.
!
        KPVT(K) = 1 - K
        if (SWAP) KPVT(K) = -IMAX
        KPVT(K-1) = KPVT(K)
  190    CONTINUE
     IK = IK - (K - 1)
     if (KSTEP  ==  2) IK = IK - (K - 2)
     K = K - KSTEP
  go to 10
  200 CONTINUE
  return
end
