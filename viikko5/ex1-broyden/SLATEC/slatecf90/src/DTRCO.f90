subroutine DTRCO (T, LDT, N, RCOND, Z, JOB)
!
!! DTRCO estimates the condition number of a triangular matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A3
!***TYPE      DOUBLE PRECISION (STRCO-S, DTRCO-D, CTRCO-C)
!***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             TRIANGULAR MATRIX
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DTRCO estimates the condition of a double precision triangular
!     matrix.
!
!     On Entry
!
!        T       DOUBLE PRECISION(LDT,N)
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
!                = 0         T  is lower triangular.
!                = nonzero   T  is upper triangular.
!
!     On Return
!
!        RCOND   DOUBLE PRECISION
!                an estimate of the reciprocal condition of  T .
!                For the system  T*X = B , relative perturbations
!                in  T  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND  ==  1.0
!                is true, then  T  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.
!
!        Z       DOUBLE PRECISION(N)
!                a work vector whose contents are usually unimportant.
!                If  T  is close to a singular matrix, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DASUM, DAXPY, DSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DTRCO
  INTEGER LDT,N,JOB
  DOUBLE PRECISION T(LDT,*),Z(*)
  DOUBLE PRECISION RCOND
!
  DOUBLE PRECISION W,WK,WKM,EK
  DOUBLE PRECISION TNORM,YNORM,S,SM,DASUM
  INTEGER I1,J,J1,J2,K,KK,L
  LOGICAL LOWER
!***FIRST EXECUTABLE STATEMENT  DTRCO
  LOWER = JOB  ==  0
!
!     COMPUTE 1-NORM OF T
!
  TNORM = 0.0D0
  DO 10 J = 1, N
     L = J
     if (LOWER) L = N + 1 - J
     I1 = 1
     if (LOWER) I1 = J
     TNORM = MAX(TNORM,DASUM(L,T(I1,J),1))
   10 CONTINUE
!
!     RCOND = 1/(NORM(T)*(ESTIMATE OF NORM(INVERSE(T)))) .
!     ESTIMATE = NORM(Z)/NORM(Y) WHERE  T*Z = Y  AND  TRANS(T)*Y = E .
!     TRANS(T)  IS THE TRANSPOSE OF T .
!     THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!     GROWTH IN THE ELEMENTS OF Y .
!     THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!     SOLVE TRANS(T)*Y = E
!
  EK = 1.0D0
  DO 20 J = 1, N
     Z(J) = 0.0D0
   20 CONTINUE
  DO 100 KK = 1, N
     K = KK
     if (LOWER) K = N + 1 - KK
     if (Z(K)  /=  0.0D0) EK = SIGN(EK,-Z(K))
     if (ABS(EK-Z(K))  <=  ABS(T(K,K))) go to 30
        S = ABS(T(K,K))/ABS(EK-Z(K))
        call DSCAL(N,S,Z,1)
        EK = S*EK
   30    CONTINUE
     WK = EK - Z(K)
     WKM = -EK - Z(K)
     S = ABS(WK)
     SM = ABS(WKM)
     if (T(K,K)  ==  0.0D0) go to 40
        WK = WK/T(K,K)
        WKM = WKM/T(K,K)
     go to 50
   40    CONTINUE
        WK = 1.0D0
        WKM = 1.0D0
   50    CONTINUE
     if (KK  ==  N) go to 90
        J1 = K + 1
        if (LOWER) J1 = 1
        J2 = N
        if (LOWER) J2 = K - 1
        DO 60 J = J1, J2
           SM = SM + ABS(Z(J)+WKM*T(K,J))
           Z(J) = Z(J) + WK*T(K,J)
           S = S + ABS(Z(J))
   60       CONTINUE
        if (S  >=  SM) go to 80
           W = WKM - WK
           WK = WKM
           DO 70 J = J1, J2
              Z(J) = Z(J) + W*T(K,J)
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
     Z(K) = WK
  100 CONTINUE
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
!
  YNORM = 1.0D0
!
!     SOLVE T*Z = Y
!
  DO 130 KK = 1, N
     K = N + 1 - KK
     if (LOWER) K = KK
     if (ABS(Z(K))  <=  ABS(T(K,K))) go to 110
        S = ABS(T(K,K))/ABS(Z(K))
        call DSCAL(N,S,Z,1)
        YNORM = S*YNORM
  110    CONTINUE
     if (T(K,K)  /=  0.0D0) Z(K) = Z(K)/T(K,K)
     if (T(K,K)  ==  0.0D0) Z(K) = 1.0D0
     I1 = 1
     if (LOWER) I1 = K + 1
     if (KK  >=  N) go to 120
        W = -Z(K)
        call DAXPY(N-KK,W,T(I1,K),1,Z(I1),1)
  120    CONTINUE
  130 CONTINUE
!     MAKE ZNORM = 1.0
  S = 1.0D0/DASUM(N,Z,1)
  call DSCAL(N,S,Z,1)
  YNORM = S*YNORM
!
  if (TNORM  /=  0.0D0) RCOND = YNORM/TNORM
  if (TNORM  ==  0.0D0) RCOND = 0.0D0
  return
end
