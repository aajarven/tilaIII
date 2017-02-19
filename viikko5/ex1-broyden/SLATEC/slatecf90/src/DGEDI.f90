subroutine DGEDI (A, LDA, N, IPVT, DET, WORK, JOB)
!
!! DGEDI computes the determinant and inverse of a matrix using the ...
!            factors computed by DGECO or DGEFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D3A1, D2A1
!***TYPE      DOUBLE PRECISION (SGEDI-S, DGEDI-D, CGEDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DGEDI computes the determinant and inverse of a matrix
!     using the factors computed by DGECO or DGEFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA, N)
!                the output from DGECO or DGEFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        IPVT    INTEGER(N)
!                the pivot vector from DGECO or DGEFA.
!
!        WORK    DOUBLE PRECISION(N)
!                work vector.  Contents destroyed.
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        A       inverse of original matrix if requested.
!                Otherwise unchanged.
!
!        DET     DOUBLE PRECISION(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  ABS(DET(1))  <  10.0
!                or  DET(1)  ==  0.0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if DGECO has set RCOND  >  0.0 or DGEFA has set
!        INFO  ==  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL, DSWAP
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DGEDI
  INTEGER LDA,N,IPVT(*),JOB
  DOUBLE PRECISION A(LDA,*),DET(2),WORK(*)
!
  DOUBLE PRECISION T
  DOUBLE PRECISION TEN
  INTEGER I,J,K,KB,KP1,L,NM1
!***FIRST EXECUTABLE STATEMENT  DGEDI
!
!     COMPUTE DETERMINANT
!
  if (JOB/10  ==  0) go to 70
     DET(1) = 1.0D0
     DET(2) = 0.0D0
     TEN = 10.0D0
     DO 50 I = 1, N
        if (IPVT(I)  /=  I) DET(1) = -DET(1)
        DET(1) = A(I,I)*DET(1)
        if (DET(1)  ==  0.0D0) go to 60
   10       if (ABS(DET(1))  >=  1.0D0) go to 20
           DET(1) = TEN*DET(1)
           DET(2) = DET(2) - 1.0D0
        go to 10
   20       CONTINUE
   30       if (ABS(DET(1))  <  TEN) go to 40
           DET(1) = DET(1)/TEN
           DET(2) = DET(2) + 1.0D0
        go to 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(U)
!
  if (MOD(JOB,10)  ==  0) go to 150
     DO 100 K = 1, N
        A(K,K) = 1.0D0/A(K,K)
        T = -A(K,K)
        call DSCAL(K-1,T,A(1,K),1)
        KP1 = K + 1
        if (N  <  KP1) go to 90
        DO 80 J = KP1, N
           T = A(K,J)
           A(K,J) = 0.0D0
           call DAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM INVERSE(U)*INVERSE(L)
!
     NM1 = N - 1
     if (NM1  <  1) go to 140
     DO 130 KB = 1, NM1
        K = N - KB
        KP1 = K + 1
        DO 110 I = KP1, N
           WORK(I) = A(I,K)
           A(I,K) = 0.0D0
  110       CONTINUE
        DO 120 J = KP1, N
           T = WORK(J)
           call DAXPY(N,T,A(1,J),1,A(1,K),1)
  120       CONTINUE
        L = IPVT(K)
        if (L  /=  K) call DSWAP(N,A(1,K),1,A(1,L),1)
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
  return
end
