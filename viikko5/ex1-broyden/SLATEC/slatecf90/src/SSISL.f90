subroutine SSISL (A, LDA, N, KPVT, B)
!
!! SSISL solves a real symmetric system using the factors obtained from SSIFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      SINGLE PRECISION (SSISL-S, DSISL-D, CHISL-C, CSISL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, SOLVE, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     SSISL solves the real symmetric system
!     A * X = B
!     using the factors computed by SSIFA.
!
!     On Entry
!
!        A       REAL(LDA,N)
!                the output from SSIFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        KPVT    INTEGER(N)
!                the pivot vector from SSIFA.
!
!        B       REAL(N)
!                the right hand side vector.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero may occur if  SSICO  has set RCOND  ==  0.0
!        or  SSIFA  has set INFO  /=  0  .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call SSIFA(A,LDA,N,KPVT,INFO)
!           if (INFO  /=  0) go to ...
!           DO 10 J = 1, P
!              call SSISL(A,LDA,N,KPVT,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SDOT
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
!***END PROLOGUE  SSISL
  INTEGER LDA,N,KPVT(*)
  REAL A(LDA,*),B(*)
!
  REAL AK,AKM1,BK,BKM1,SDOT,DENOM,TEMP
  INTEGER K,KP
!
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
!     D INVERSE TO B.
!
!***FIRST EXECUTABLE STATEMENT  SSISL
  K = N
   10 if (K  ==  0) go to 80
     if (KPVT(K)  <  0) go to 40
!
!           1 X 1 PIVOT BLOCK.
!
        if (K  ==  1) go to 30
           KP = KPVT(K)
           if (KP  ==  K) go to 20
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
   20          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
           call SAXPY(K-1,B(K),A(1,K),1,B(1),1)
   30       CONTINUE
!
!           APPLY D INVERSE.
!
        B(K) = B(K)/A(K,K)
        K = K - 1
     go to 70
   40    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
        if (K  ==  2) go to 60
           KP = ABS(KPVT(K))
           if (KP  ==  K - 1) go to 50
!
!                 INTERCHANGE.
!
              TEMP = B(K-1)
              B(K-1) = B(KP)
              B(KP) = TEMP
   50          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
           call SAXPY(K-2,B(K),A(1,K),1,B(1),1)
           call SAXPY(K-2,B(K-1),A(1,K-1),1,B(1),1)
   60       CONTINUE
!
!           APPLY D INVERSE.
!
        AK = A(K,K)/A(K-1,K)
        AKM1 = A(K-1,K-1)/A(K-1,K)
        BK = B(K)/A(K-1,K)
        BKM1 = B(K-1)/A(K-1,K)
        DENOM = AK*AKM1 - 1.0E0
        B(K) = (AKM1*BK - BKM1)/DENOM
        B(K-1) = (AK*BKM1 - BK)/DENOM
        K = K - 2
   70    CONTINUE
  go to 10
   80 CONTINUE
!
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
!
  K = 1
   90 if (K  >  N) go to 160
     if (KPVT(K)  <  0) go to 120
!
!           1 X 1 PIVOT BLOCK.
!
        if (K  ==  1) go to 110
!
!              APPLY THE TRANSFORMATION.
!
           B(K) = B(K) + SDOT(K-1,A(1,K),1,B(1),1)
           KP = KPVT(K)
           if (KP  ==  K) go to 100
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
        K = K + 1
     go to 150
  120    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
        if (K  ==  1) go to 140
!
!              APPLY THE TRANSFORMATION.
!
           B(K) = B(K) + SDOT(K-1,A(1,K),1,B(1),1)
           B(K+1) = B(K+1) + SDOT(K-1,A(1,K+1),1,B(1),1)
           KP = ABS(KPVT(K))
           if (KP  ==  K) go to 130
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
        K = K + 2
  150    CONTINUE
  go to 90
  160 CONTINUE
  return
end
