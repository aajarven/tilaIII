subroutine CTRDI (T, LDT, N, DET, JOB, INFO)
!
!! CTRDI computes the determinant and inverse of a triangular matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C3, D3C3
!***TYPE      COMPLEX (STRDI-S, DTRDI-D, CTRDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
!             TRIANGULAR MATRIX
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     CTRDI computes the determinant and inverse of a complex
!     triangular matrix.
!
!     On Entry
!
!        T       COMPLEX(LDT,N)
!                T contains the triangular matrix.  The zero
!                elements of the matrix are not referenced, and
!                the corresponding elements of the array can be
!                used to store other information.
!
!        LDT     INTEGER
!                LDT is the leading dimension of the array T.
!
!        N       INTEGER
!                N is the order of the system.
!
!        JOB     INTEGER
!                = 010       no det, inverse of lower triangular.
!                = 011       no det, inverse of upper triangular.
!                = 100       det, no inverse.
!                = 110       det, inverse of lower triangular.
!                = 111       det, inverse of upper triangular.
!
!     On Return
!
!        T       inverse of original matrix if requested.
!                Otherwise unchanged.
!
!        DET     COMPLEX(2)
!                determinant of original matrix if requested.
!                Otherwise not referenced.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  CABS1(DET(1))  <  10.0
!                or  DET(1)  ==  0.0 .
!
!        INFO    INTEGER
!                INFO contains zero if the system is nonsingular
!                and the inverse is requested.
!                Otherwise INFO contains the index of
!                a zero diagonal element of T.
!
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
!***END PROLOGUE  CTRDI
  INTEGER LDT,N,JOB,INFO
  COMPLEX T(LDT,*),DET(2)
!
  COMPLEX TEMP
  REAL TEN
  INTEGER I,J,K,KB,KM1,KP1
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!***FIRST EXECUTABLE STATEMENT  CTRDI
!
!        COMPUTE DETERMINANT
!
     if (JOB/100  ==  0) go to 70
        DET(1) = (1.0E0,0.0E0)
        DET(2) = (0.0E0,0.0E0)
        TEN = 10.0E0
        DO 50 I = 1, N
           DET(1) = T(I,I)*DET(1)
           if (CABS1(DET(1))  ==  0.0E0) go to 60
   10          if (CABS1(DET(1))  >=  1.0E0) go to 20
              DET(1) = CMPLX(TEN,0.0E0)*DET(1)
              DET(2) = DET(2) - (1.0E0,0.0E0)
           go to 10
   20          CONTINUE
   30          if (CABS1(DET(1))  <  TEN) go to 40
              DET(1) = DET(1)/CMPLX(TEN,0.0E0)
              DET(2) = DET(2) + (1.0E0,0.0E0)
           go to 30
   40          CONTINUE
   50       CONTINUE
   60       CONTINUE
   70    CONTINUE
!
!        COMPUTE INVERSE OF UPPER TRIANGULAR
!
     if (MOD(JOB/10,10)  ==  0) go to 170
        if (MOD(JOB,10)  ==  0) go to 120
              DO 100 K = 1, N
                 INFO = K
                 if (CABS1(T(K,K))  ==  0.0E0) go to 110
                 T(K,K) = (1.0E0,0.0E0)/T(K,K)
                 TEMP = -T(K,K)
                 call CSCAL(K-1,TEMP,T(1,K),1)
                 KP1 = K + 1
                 if (N  <  KP1) go to 90
                 DO 80 J = KP1, N
                    TEMP = T(K,J)
                    T(K,J) = (0.0E0,0.0E0)
                    call CAXPY(K,TEMP,T(1,K),1,T(1,J),1)
   80                CONTINUE
   90                CONTINUE
  100             CONTINUE
              INFO = 0
  110          CONTINUE
        go to 160
  120       CONTINUE
!
!              COMPUTE INVERSE OF LOWER TRIANGULAR
!
           DO 150 KB = 1, N
              K = N + 1 - KB
              INFO = K
              if (CABS1(T(K,K))  ==  0.0E0) go to 180
              T(K,K) = (1.0E0,0.0E0)/T(K,K)
              TEMP = -T(K,K)
              if (K  /=  N) call CSCAL(N-K,TEMP,T(K+1,K),1)
              KM1 = K - 1
              if (KM1  <  1) go to 140
              DO 130 J = 1, KM1
                 TEMP = T(K,J)
                 T(K,J) = (0.0E0,0.0E0)
                 call CAXPY(N-K+1,TEMP,T(K,K),1,T(K,J),1)
  130             CONTINUE
  140             CONTINUE
  150          CONTINUE
           INFO = 0
  160       CONTINUE
  170    CONTINUE
  180 CONTINUE
  return
end
