subroutine CPODI (A, LDA, N, DET, JOB)
!
!! CPODI computes the determinant and inverse of a certain complex ...
!            Hermitian positive definite matrix using the factors ...
!            computed by CPOCO, CPOFA, or CQRDC.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B, D3D1B
!***TYPE      COMPLEX (SPODI-S, DPODI-D, CPODI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPODI computes the determinant and inverse of a certain
!     complex Hermitian positive definite matrix (see below)
!     using the factors computed by CPOCO, CPOFA or CQRDC.
!
!     On Entry
!
!        A       COMPLEX(LDA, N)
!                the output  A  from CPOCO or CPOFA
!                or the output  X  from CQRDC.
!
!        LDA     INTEGER
!                the leading dimension of the array  A .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        JOB     INTEGER
!                = 11   both determinant and inverse.
!                = 01   inverse only.
!                = 10   determinant only.
!
!     On Return
!
!        A       If CPOCO or CPOFA was used to factor  A  then
!                CPODI produces the upper half of INVERSE(A) .
!                If CQRDC was used to decompose  X  then
!                CPODI produces the upper half of INVERSE(CTRANS(X)*X)
!                where CTRANS(X) is the conjugate transpose.
!                Elements of  A  below the diagonal are unchanged.
!                If the units digit of JOB is zero,  A  is unchanged.
!
!        DET     REAL(2)
!                determinant of  A  or of  CTRANS(X)*X  if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  DET(1)  <  10.0
!                or  DET(1)  ==  0.0 .
!
!     Error Condition
!
!        a division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if CPOCO or CPOFA has set INFO  ==  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CPODI
  INTEGER LDA,N,JOB
  COMPLEX A(LDA,*)
  REAL DET(2)
!
  COMPLEX T
  REAL S
  INTEGER I,J,JM1,K,KP1
!***FIRST EXECUTABLE STATEMENT  CPODI
!
!     COMPUTE DETERMINANT
!
  if (JOB/10  ==  0) go to 70
     DET(1) = 1.0E0
     DET(2) = 0.0E0
     S = 10.0E0
     DO 50 I = 1, N
        DET(1) = REAL(A(I,I))**2*DET(1)
        if (DET(1)  ==  0.0E0) go to 60
   10       if (DET(1)  >=  1.0E0) go to 20
           DET(1) = S*DET(1)
           DET(2) = DET(2) - 1.0E0
        go to 10
   20       CONTINUE
   30       if (DET(1)  <  S) go to 40
           DET(1) = DET(1)/S
           DET(2) = DET(2) + 1.0E0
        go to 30
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     COMPUTE INVERSE(R)
!
  if (MOD(JOB,10)  ==  0) go to 140
     DO 100 K = 1, N
        A(K,K) = (1.0E0,0.0E0)/A(K,K)
        T = -A(K,K)
        call CSCAL(K-1,T,A(1,K),1)
        KP1 = K + 1
        if (N  <  KP1) go to 90
        DO 80 J = KP1, N
           T = A(K,J)
           A(K,J) = (0.0E0,0.0E0)
           call CAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * CTRANS(INVERSE(R))
!
     DO 130 J = 1, N
        JM1 = J - 1
        if (JM1  <  1) go to 120
        DO 110 K = 1, JM1
           T = CONJG(A(K,J))
           call CAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
        T = CONJG(A(J,J))
        call CSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
  return
end
