subroutine DPPDI (AP, N, DET, JOB)
!
!! DPPDI computes the determinant and inverse of a real symmetric ...
!  positive definite matrix using factors from DPPCO or DPPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B, D3B1B
!***TYPE      DOUBLE PRECISION (SPPDI-S, DPPDI-D, CPPDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             PACKED, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPPDI computes the determinant and inverse
!     of a double precision symmetric positive definite matrix
!     using the factors computed by DPPCO or DPPFA .
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the output from DPPCO or DPPFA.
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
!        DET     DOUBLE PRECISION(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                with  1.0  <=  DET(1)  <  10.0
!                or  DET(1)  ==  0.0 .
!
!     Error Condition
!
!        A division by zero will occur if the input factor contains
!        a zero on the diagonal and the inverse is requested.
!        It will not occur if the subroutines are called correctly
!        and if DPOCO or DPOFA has set INFO  ==  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPPDI
  INTEGER N,JOB
  DOUBLE PRECISION AP(*)
  DOUBLE PRECISION DET(2)
!
  DOUBLE PRECISION T
  DOUBLE PRECISION S
  INTEGER I,II,J,JJ,JM1,J1,K,KJ,KK,KP1,K1
!***FIRST EXECUTABLE STATEMENT  DPPDI
!
!     COMPUTE DETERMINANT
!
  if (JOB/10  ==  0) go to 70
     DET(1) = 1.0D0
     DET(2) = 0.0D0
     S = 10.0D0
     II = 0
     DO 50 I = 1, N
        II = II + I
        DET(1) = AP(II)**2*DET(1)
        if (DET(1)  ==  0.0D0) go to 60
   10       if (DET(1)  >=  1.0D0) go to 20
           DET(1) = S*DET(1)
           DET(2) = DET(2) - 1.0D0
        go to 10
   20       CONTINUE
   30       if (DET(1)  <  S) go to 40
           DET(1) = DET(1)/S
           DET(2) = DET(2) + 1.0D0
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
        AP(KK) = 1.0D0/AP(KK)
        T = -AP(KK)
        call DSCAL(K-1,T,AP(K1),1)
        KP1 = K + 1
        J1 = KK + 1
        KJ = KK + K
        if (N  <  KP1) go to 90
        DO 80 J = KP1, N
           T = AP(KJ)
           AP(KJ) = 0.0D0
           call DAXPY(K,T,AP(K1),1,AP(J1),1)
           J1 = J1 + J
           KJ = KJ + J
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * TRANS(INVERSE(R))
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
           T = AP(KJ)
           call DAXPY(K,T,AP(J1),1,AP(K1),1)
           K1 = K1 + K
           KJ = KJ + 1
  110       CONTINUE
  120       CONTINUE
        T = AP(JJ)
        call DSCAL(J,T,AP(J1),1)
  130    CONTINUE
  140 CONTINUE
  return
end
