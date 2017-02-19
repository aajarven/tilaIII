subroutine CPPDI (AP, N, DET, JOB)
!
!! CPPDI computes the determinant and inverse of a complex Hermitian ...
!            positive definite matrix using factors from CPPCO or CPPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2D1B, D3D1B
!***TYPE      COMPLEX (SPPDI-S, DPPDI-D, CPPDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             PACKED, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CPPDI computes the determinant and inverse
!     of a complex Hermitian positive definite matrix
!     using the factors computed by CPPCO or CPPFA .
!
!     On Entry
!
!        AP      COMPLEX (N*(N+1)/2)
!                the output from CPPCO or CPPFA.
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
!        AP      the upper triangular half of the inverse .
!                The strict lower triangle is unaltered.
!
!        DET     REAL(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  DET(1)  <  10.0
!                or  DET(1)  ==  0.0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
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
!***END PROLOGUE  CPPDI
  INTEGER N,JOB
  COMPLEX AP(*)
  REAL DET(2)
!
  COMPLEX T
  REAL S
  INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
!***FIRST EXECUTABLE STATEMENT  CPPDI
!
!     COMPUTE DETERMINANT
!
  if (JOB/10  ==  0) go to 70
     DET(1) = 1.0E0
     DET(2) = 0.0E0
     S = 10.0E0
     II = 0
     DO 50 I = 1, N
        II = II + I
        DET(1) = REAL(AP(II))**2*DET(1)
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
     KK = 0
     DO 100 K = 1, N
        K1 = KK + 1
        KK = KK + K
        AP(KK) = (1.0E0,0.0E0)/AP(KK)
        T = -AP(KK)
        call CSCAL(K-1,T,AP(K1),1)
        KP1 = K + 1
        J1 = KK + 1
        KJ = KK + K
        if (N  <  KP1) go to 90
        DO 80 J = KP1, N
           T = AP(KJ)
           AP(KJ) = (0.0E0,0.0E0)
           call CAXPY(K,T,AP(K1),1,AP(J1),1)
           J1 = J1 + J
           KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * CTRANS(INVERSE(R))
!
     JJ = 0
     DO 130 J = 1, N
        J1 = JJ + 1
        JJ = JJ + J
        JM1 = J - 1
        K1 = 1
        KJ = J1
        if (JM1  <  1) go to 120
        DO 110 K = 1, JM1
           T = CONJG(AP(KJ))
           call CAXPY(K,T,AP(J1),1,AP(K1),1)
           K1 = K1 + K
           KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
        T = CONJG(AP(JJ))
        call CSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
  return
end
