subroutine SPODI (A, LDA, N, DET, JOB)
!
!! SPODI computes the determinant and inverse of a certain real symmetric ...
!  positive definite matrix using the factors
!            computed by SPOCO, SPOFA or SQRDC.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B, D3B1B
!***TYPE      SINGLE PRECISION (SPODI-S, DPODI-D, CPODI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPODI computes the determinant and inverse of a certain
!     real symmetric positive definite matrix (see below)
!     using the factors computed by SPOCO, SPOFA or SQRDC.
!
!     On Entry
!
!        A       REAL(LDA, N)
!                the output  A  from SPOCO or SPOFA
!                or the output  X  from SQRDC.
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
!        A       If SPOCO or SPOFA was used to factor  A , then
!                SPODI produces the upper half of INVERSE(A) .
!                If SQRDC was used to decompose  X , then
!                SPODI produces the upper half of INVERSE(TRANS(X)*X),
!                where TRANS(X) is the transpose.
!                Elements of  A  below the diagonal are unchanged.
!                If the units digit of JOB is zero,  A  is unchanged.
!
!        DET     REAL(2)
!                determinant of  A  or of  TRANS(X)*X  if requested.
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
!        and if SPOCO or SPOFA has set INFO  ==  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPODI
  INTEGER LDA,N,JOB
  REAL A(LDA,*)
  REAL DET(2)
!
  REAL T
  REAL S
  INTEGER I,J,JM1,K,KP1
!***FIRST EXECUTABLE STATEMENT  SPODI
!
!     COMPUTE DETERMINANT
!
  if (JOB/10  ==  0) go to 70
     DET(1) = 1.0E0
     DET(2) = 0.0E0
     S = 10.0E0
     DO 50 I = 1, N
        DET(1) = A(I,I)**2*DET(1)
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
        A(K,K) = 1.0E0/A(K,K)
        T = -A(K,K)
        call SSCAL(K-1,T,A(1,K),1)
        KP1 = K + 1
        if (N  <  KP1) go to 90
        DO 80 J = KP1, N
           T = A(K,J)
           A(K,J) = 0.0E0
           call SAXPY(K,T,A(1,K),1,A(1,J),1)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        FORM  INVERSE(R) * TRANS(INVERSE(R))
!
     DO 130 J = 1, N
        JM1 = J - 1
        if (JM1  <  1) go to 120
        DO 110 K = 1, JM1
           T = A(K,J)
           call SAXPY(K,T,A(1,J),1,A(1,K),1)
  110       CONTINUE
  120       CONTINUE
        T = A(J,J)
        call SSCAL(J,T,A(1,J),1)
  130    CONTINUE
  140 CONTINUE
  return
end
