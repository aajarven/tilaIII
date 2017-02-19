subroutine QZVAL (NM, N, A, B, ALFR, ALFI, BETA, MATZ, Z)
!
!! QZVAL is the third step of the QZ algorithm for generalized eigenproblems.
!  Accepts a pair of real matrices, one in
!            quasi-triangular form and the other in upper triangular
!            form and computes the eigenvalues of the associated
!            eigenproblem.  Usually preceded by QZHES, QZIT, and
!            followed by QZVEC.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2C
!***TYPE      SINGLE PRECISION (QZVAL-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is the third step of the QZ algorithm
!     for solving generalized matrix eigenvalue problems,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
!
!     This subroutine accepts a pair of REAL matrices, one of them
!     in quasi-triangular form and the other in upper triangular form.
!     It reduces the quasi-triangular matrix further, so that any
!     remaining 2-by-2 blocks correspond to pairs of complex
!     eigenvalues, and returns quantities whose ratios give the
!     generalized eigenvalues.  It is usually preceded by  QZHES
!     and  QZIT  and may be followed by  QZVEC.
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
!        A contains a real upper quasi-triangular matrix.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!        B contains a real upper triangular matrix.  In addition,
!          location B(N,1) contains the tolerance quantity (EPSB)
!          computed and saved in  QZIT.  B is a two-dimensional REAL
!          array, dimensioned B(NM,N).
!
!        MATZ should be set to .TRUE. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
!          variable.
!
!        Z contains, if MATZ has been set to .TRUE., the transformation
!          matrix produced in the reductions by  QZHES  and  QZIT,  if
!          performed, or else the identity matrix.  If MATZ has been set
!          to .FALSE., Z is not referenced.  Z is a two-dimensional REAL
!          array, dimensioned Z(NM,N).
!
!     On Output
!
!        A has been reduced further to a quasi-triangular matrix in
!          which all nonzero subdiagonal elements correspond to pairs
!          of complex eigenvalues.
!
!        B is still in upper triangular form, although its elements
!          have been altered.  B(N,1) is unaltered.
!
!        ALFR and ALFI contain the real and imaginary parts of the
!          diagonal elements of the triangular matrix that would be
!          obtained if A were reduced completely to triangular form
!          by unitary transformations.  Non-zero values of ALFI occur
!          in pairs, the first member positive and the second negative.
!          ALFR and ALFI are one-dimensional REAL arrays, dimensioned
!          ALFR(N) and ALFI(N).
!
!        BETA contains the diagonal elements of the corresponding B,
!          normalized to be real and non-negative.  The generalized
!          eigenvalues are then the ratios ((ALFR+I*ALFI)/BETA).
!          BETA is a one-dimensional REAL array, dimensioned BETA(N).
!
!        Z contains the product of the right hand transformations
!          (for all three steps) if MATZ has been set to .TRUE.
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
!***END PROLOGUE  QZVAL
!
  INTEGER I,J,N,EN,NA,NM,NN,ISW
  REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
  REAL C,D,E,R,S,T,AN,A1,A2,BN,CQ,CZ,DI,DR,EI,TI,TR
  REAL U1,U2,V1,V2,A1I,A11,A12,A2I,A21,A22,B11,B12,B22
  REAL SQI,SQR,SSI,SSR,SZI,SZR,A11I,A11R,A12I,A12R
  REAL A22I,A22R,EPSB
  LOGICAL MATZ
!
!***FIRST EXECUTABLE STATEMENT  QZVAL
  EPSB = B(N,1)
  ISW = 1
!     .......... FIND EIGENVALUES OF QUASI-TRIANGULAR MATRICES.
!                FOR EN=N STEP -1 UNTIL 1 DO -- ..........
  DO 510 NN = 1, N
     EN = N + 1 - NN
     NA = EN - 1
     if (ISW  ==  2) go to 505
     if (EN  ==  1) go to 410
     if (A(EN,NA)  /=  0.0E0) go to 420
!     .......... 1-BY-1 BLOCK, ONE REAL ROOT ..........
  410    ALFR(EN) = A(EN,EN)
     if (B(EN,EN)  <  0.0E0) ALFR(EN) = -ALFR(EN)
     BETA(EN) = ABS(B(EN,EN))
     ALFI(EN) = 0.0E0
     go to 510
!     .......... 2-BY-2 BLOCK ..........
  420    if (ABS(B(NA,NA))  <=  EPSB) go to 455
     if (ABS(B(EN,EN))  >  EPSB) go to 430
     A1 = A(EN,EN)
     A2 = A(EN,NA)
     BN = 0.0E0
     go to 435
  430    AN = ABS(A(NA,NA)) + ABS(A(NA,EN)) + ABS(A(EN,NA)) &
        + ABS(A(EN,EN))
     BN = ABS(B(NA,NA)) + ABS(B(NA,EN)) + ABS(B(EN,EN))
     A11 = A(NA,NA) / AN
     A12 = A(NA,EN) / AN
     A21 = A(EN,NA) / AN
     A22 = A(EN,EN) / AN
     B11 = B(NA,NA) / BN
     B12 = B(NA,EN) / BN
     B22 = B(EN,EN) / BN
     E = A11 / B11
     EI = A22 / B22
     S = A21 / (B11 * B22)
     T = (A22 - E * B22) / B22
     if (ABS(E)  <=  ABS(EI)) go to 431
     E = EI
     T = (A11 - E * B11) / B11
  431    C = 0.5E0 * (T - S * B12)
     D = C * C + S * (A12 - E * B12)
     if (D  <  0.0E0) go to 480
!     .......... TWO REAL ROOTS.
!                ZERO BOTH A(EN,NA) AND B(EN,NA) ..........
     E = E + (C + SIGN(SQRT(D),C))
     A11 = A11 - E * B11
     A12 = A12 - E * B12
     A22 = A22 - E * B22
     if (ABS(A11) + ABS(A12)  <  &
         ABS(A21) + ABS(A22)) go to 432
     A1 = A12
     A2 = A11
     go to 435
  432    A1 = A22
     A2 = A21
!     .......... CHOOSE AND APPLY REAL Z ..........
  435    S = ABS(A1) + ABS(A2)
     U1 = A1 / S
     U2 = A2 / S
     R = SIGN(SQRT(U1*U1+U2*U2),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     U2 = V2 / V1
!
     DO 440 I = 1, EN
        T = A(I,EN) + U2 * A(I,NA)
        A(I,EN) = A(I,EN) + T * V1
        A(I,NA) = A(I,NA) + T * V2
        T = B(I,EN) + U2 * B(I,NA)
        B(I,EN) = B(I,EN) + T * V1
        B(I,NA) = B(I,NA) + T * V2
  440    CONTINUE
!
     if (.NOT. MATZ) go to 450
!
     DO 445 I = 1, N
        T = Z(I,EN) + U2 * Z(I,NA)
        Z(I,EN) = Z(I,EN) + T * V1
        Z(I,NA) = Z(I,NA) + T * V2
  445    CONTINUE
!
  450    if (BN  ==  0.0E0) go to 475
     if (AN  <  ABS(E) * BN) go to 455
     A1 = B(NA,NA)
     A2 = B(EN,NA)
     go to 460
  455    A1 = A(NA,NA)
     A2 = A(EN,NA)
!     .......... CHOOSE AND APPLY REAL Q ..........
  460    S = ABS(A1) + ABS(A2)
     if (S  ==  0.0E0) go to 475
     U1 = A1 / S
     U2 = A2 / S
     R = SIGN(SQRT(U1*U1+U2*U2),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     U2 = V2 / V1
!
     DO 470 J = NA, N
        T = A(NA,J) + U2 * A(EN,J)
        A(NA,J) = A(NA,J) + T * V1
        A(EN,J) = A(EN,J) + T * V2
        T = B(NA,J) + U2 * B(EN,J)
        B(NA,J) = B(NA,J) + T * V1
        B(EN,J) = B(EN,J) + T * V2
  470    CONTINUE
!
  475    A(EN,NA) = 0.0E0
     B(EN,NA) = 0.0E0
     ALFR(NA) = A(NA,NA)
     ALFR(EN) = A(EN,EN)
     if (B(NA,NA)  <  0.0E0) ALFR(NA) = -ALFR(NA)
     if (B(EN,EN)  <  0.0E0) ALFR(EN) = -ALFR(EN)
     BETA(NA) = ABS(B(NA,NA))
     BETA(EN) = ABS(B(EN,EN))
     ALFI(EN) = 0.0E0
     ALFI(NA) = 0.0E0
     go to 505
!     .......... TWO COMPLEX ROOTS ..........
  480    E = E + C
     EI = SQRT(-D)
     A11R = A11 - E * B11
     A11I = EI * B11
     A12R = A12 - E * B12
     A12I = EI * B12
     A22R = A22 - E * B22
     A22I = EI * B22
     if (ABS(A11R) + ABS(A11I) + ABS(A12R) + ABS(A12I)  <  &
         ABS(A21) + ABS(A22R) + ABS(A22I)) go to 482
     A1 = A12R
     A1I = A12I
     A2 = -A11R
     A2I = -A11I
     go to 485
  482    A1 = A22R
     A1I = A22I
     A2 = -A21
     A2I = 0.0E0
!     .......... CHOOSE COMPLEX Z ..........
  485    CZ = SQRT(A1*A1+A1I*A1I)
     if (CZ  ==  0.0E0) go to 487
     SZR = (A1 * A2 + A1I * A2I) / CZ
     SZI = (A1 * A2I - A1I * A2) / CZ
     R = SQRT(CZ*CZ+SZR*SZR+SZI*SZI)
     CZ = CZ / R
     SZR = SZR / R
     SZI = SZI / R
     go to 490
  487    SZR = 1.0E0
     SZI = 0.0E0
  490    if (AN  <  (ABS(E) + EI) * BN) go to 492
     A1 = CZ * B11 + SZR * B12
     A1I = SZI * B12
     A2 = SZR * B22
     A2I = SZI * B22
     go to 495
  492    A1 = CZ * A11 + SZR * A12
     A1I = SZI * A12
     A2 = CZ * A21 + SZR * A22
     A2I = SZI * A22
!     .......... CHOOSE COMPLEX Q ..........
  495    CQ = SQRT(A1*A1+A1I*A1I)
     if (CQ  ==  0.0E0) go to 497
     SQR = (A1 * A2 + A1I * A2I) / CQ
     SQI = (A1 * A2I - A1I * A2) / CQ
     R = SQRT(CQ*CQ+SQR*SQR+SQI*SQI)
     CQ = CQ / R
     SQR = SQR / R
     SQI = SQI / R
     go to 500
  497    SQR = 1.0E0
     SQI = 0.0E0
!     .......... COMPUTE DIAGONAL ELEMENTS THAT WOULD RESULT
!                if TRANSFORMATIONS WERE APPLIED ..........
  500    SSR = SQR * SZR + SQI * SZI
     SSI = SQR * SZI - SQI * SZR
     I = 1
     TR = CQ * CZ * A11 + CQ * SZR * A12 + SQR * CZ * A21 &
        + SSR * A22
     TI = CQ * SZI * A12 - SQI * CZ * A21 + SSI * A22
     DR = CQ * CZ * B11 + CQ * SZR * B12 + SSR * B22
     DI = CQ * SZI * B12 + SSI * B22
     go to 503
  502    I = 2
     TR = SSR * A11 - SQR * CZ * A12 - CQ * SZR * A21 &
        + CQ * CZ * A22
     TI = -SSI * A11 - SQI * CZ * A12 + CQ * SZI * A21
     DR = SSR * B11 - SQR * CZ * B12 + CQ * CZ * B22
     DI = -SSI * B11 - SQI * CZ * B12
  503    T = TI * DR - TR * DI
     J = NA
     if (T  <  0.0E0) J = EN
     R = SQRT(DR*DR+DI*DI)
     BETA(J) = BN * R
     ALFR(J) = AN * (TR * DR + TI * DI) / R
     ALFI(J) = AN * T / R
     if (I  ==  1) go to 502
  505    ISW = 3 - ISW
  510 CONTINUE
!
  return
end
