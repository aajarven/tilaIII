subroutine QZHES (NM, N, A, B, MATZ, Z)
!
!! QZHES is the first step of the QZ algorithm for generalized eigenproblems.
!  Accepts a pair of real general
!            matrices and reduces one of them to upper Hessenberg
!            and the other to upper triangular form using orthogonal
!            transformations. Usually followed by QZIT, QZVAL, QZVEC.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B3
!***TYPE      SINGLE PRECISION (QZHES-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is the first step of the QZ algorithm
!     for solving generalized matrix eigenvalue problems,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
!
!     This subroutine accepts a pair of REAL GENERAL matrices and
!     reduces one of them to upper Hessenberg form and the other
!     to upper triangular form using orthogonal transformations.
!     It is usually followed by  QZIT,  QZVAL  and, possibly,  QZVEC.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, B, and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrices A and B.  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        A contains a real general matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        B contains a real general matrix.  B is a two-dimensional
!          REAL array, dimensioned B(NM,N).
!
!        MATZ should be set to .TRUE. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
!          variable.
!
!     On Output
!
!        A has been reduced to upper Hessenberg form.  The elements
!          below the first subdiagonal have been set to zero.
!
!        B has been reduced to upper triangular form.  The elements
!          below the main diagonal have been set to zero.
!
!        Z contains the product of the right hand transformations if
!          MATZ has been set to .TRUE.  Otherwise, Z is not referenced.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  QZHES
!
  INTEGER I,J,K,L,N,LB,L1,NM,NK1,NM1,NM2
  REAL A(NM,*),B(NM,*),Z(NM,*)
  REAL R,S,T,U1,U2,V1,V2,RHO
  LOGICAL MATZ
!
!     .......... INITIALIZE Z ..........
!***FIRST EXECUTABLE STATEMENT  QZHES
  if (.NOT. MATZ) go to 10
!
  DO 3 I = 1, N
!
     DO 2 J = 1, N
        Z(I,J) = 0.0E0
    2    CONTINUE
!
     Z(I,I) = 1.0E0
    3 CONTINUE
!     .......... REDUCE B TO UPPER TRIANGULAR FORM ..........
   10 if (N  <=  1) go to 170
  NM1 = N - 1
!
  DO 100 L = 1, NM1
     L1 = L + 1
     S = 0.0E0
!
     DO 20 I = L1, N
        S = S + ABS(B(I,L))
   20    CONTINUE
!
     if (S  ==  0.0E0) go to 100
     S = S + ABS(B(L,L))
     R = 0.0E0
!
     DO 25 I = L, N
        B(I,L) = B(I,L) / S
        R = R + B(I,L)**2
   25    CONTINUE
!
     R = SIGN(SQRT(R),B(L,L))
     B(L,L) = B(L,L) + R
     RHO = R * B(L,L)
!
     DO 50 J = L1, N
        T = 0.0E0
!
        DO 30 I = L, N
           T = T + B(I,L) * B(I,J)
   30       CONTINUE
!
        T = -T / RHO
!
        DO 40 I = L, N
           B(I,J) = B(I,J) + T * B(I,L)
   40       CONTINUE
!
   50    CONTINUE
!
     DO 80 J = 1, N
        T = 0.0E0
!
        DO 60 I = L, N
           T = T + B(I,L) * A(I,J)
   60       CONTINUE
!
        T = -T / RHO
!
        DO 70 I = L, N
           A(I,J) = A(I,J) + T * B(I,L)
   70       CONTINUE
!
   80    CONTINUE
!
     B(L,L) = -S * R
!
     DO 90 I = L1, N
        B(I,L) = 0.0E0
   90    CONTINUE
!
  100 CONTINUE
!     .......... REDUCE A TO UPPER HESSENBERG FORM, WHILE
!                KEEPING B TRIANGULAR ..........
  if (N  ==  2) go to 170
  NM2 = N - 2
!
  DO 160 K = 1, NM2
     NK1 = NM1 - K
!     .......... FOR L=N-1 STEP -1 UNTIL K+1 DO -- ..........
     DO 150 LB = 1, NK1
        L = N - LB
        L1 = L + 1
!     .......... ZERO A(L+1,K) ..........
        S = ABS(A(L,K)) + ABS(A(L1,K))
        if (S  ==  0.0E0) go to 150
        U1 = A(L,K) / S
        U2 = A(L1,K) / S
        R = SIGN(SQRT(U1*U1+U2*U2),U1)
        V1 =  -(U1 + R) / R
        V2 = -U2 / R
        U2 = V2 / V1
!
        DO 110 J = K, N
           T = A(L,J) + U2 * A(L1,J)
           A(L,J) = A(L,J) + T * V1
           A(L1,J) = A(L1,J) + T * V2
  110       CONTINUE
!
        A(L1,K) = 0.0E0
!
        DO 120 J = L, N
           T = B(L,J) + U2 * B(L1,J)
           B(L,J) = B(L,J) + T * V1
           B(L1,J) = B(L1,J) + T * V2
  120       CONTINUE
!     .......... ZERO B(L+1,L) ..........
        S = ABS(B(L1,L1)) + ABS(B(L1,L))
        if (S  ==  0.0E0) go to 150
        U1 = B(L1,L1) / S
        U2 = B(L1,L) / S
        R = SIGN(SQRT(U1*U1+U2*U2),U1)
        V1 =  -(U1 + R) / R
        V2 = -U2 / R
        U2 = V2 / V1
!
        DO 130 I = 1, L1
           T = B(I,L1) + U2 * B(I,L)
           B(I,L1) = B(I,L1) + T * V1
           B(I,L) = B(I,L) + T * V2
  130       CONTINUE
!
        B(L1,L) = 0.0E0
!
        DO 140 I = 1, N
           T = A(I,L1) + U2 * A(I,L)
           A(I,L1) = A(I,L1) + T * V1
           A(I,L) = A(I,L) + T * V2
  140       CONTINUE
!
        if (.NOT. MATZ) go to 150
!
        DO 145 I = 1, N
           T = Z(I,L1) + U2 * Z(I,L)
           Z(I,L1) = Z(I,L1) + T * V1
           Z(I,L) = Z(I,L) + T * V2
  145       CONTINUE
!
  150    CONTINUE
!
  160 CONTINUE
!
  170 RETURN
end
