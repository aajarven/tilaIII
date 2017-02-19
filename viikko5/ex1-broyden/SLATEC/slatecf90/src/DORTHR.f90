subroutine DORTHR (A, N, M, NRDA, IFLAG, IRANK, ISCALE, DIAG, &
     KPIVOT, SCALES, ROWS, RS)
!
!! DORTHR is subsidiary to DBVSUP and DSUDS.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (ORTHOR-S, DORTHR-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   Reduction of the matrix A to lower triangular form by a sequence of
!   orthogonal HOUSEHOLDER transformations post-multiplying A.
!
! *********************************************************************
!   INPUT
! *********************************************************************
!
!     A -- Contains the matrix to be decomposed, must be dimensioned
!           NRDA by N.
!     N -- Number of rows in the matrix, N greater or equal to 1.
!     M -- Number of columns in the matrix, M greater or equal to N.
!     IFLAG -- Indicates the uncertainty in the matrix data.
!             = 0 when the data is to be treated as exact.
!             =-K when the data is assumed to be accurate to about
!                 K digits.
!     ISCALE -- Scaling indicator.
!               =-1 if the matrix is to be pre-scaled by
!               columns when appropriate.
!               Otherwise no scaling will be attempted.
!     NRDA -- Row dimension of A, NRDA greater or equal to N.
!     DIAG,KPIVOT,ROWS, -- Arrays of length at least N used internally
!          RS,SCALES         (except for SCALES which is M).
!
! *********************************************************************
!   OUTPUT
! *********************************************************************
!
!     IFLAG - Status indicator
!            =1 for successful decomposition.
!            =2 if improper input is detected.
!            =3 if rank of the matrix is less than N.
!     A -- Contains the reduced matrix in the strictly lower triangular
!          part and transformation information.
!     IRANK -- Contains the numerically determined matrix rank.
!     DIAG -- Contains the diagonal elements of the reduced
!             triangular matrix.
!     KPIVOT -- Contains the pivotal information, the column
!               interchanges performed on the original matrix are
!               recorded here.
!     SCALES -- Contains the column scaling parameters.
!
! *********************************************************************
!
!***SEE ALSO  DBVSUP, DSUDS
!***REFERENCES  G. Golub, Numerical methods for solving linear least
!                 squares problems, Numerische Mathematik 7, (1965),
!                 pp. 206-216.
!               P. Businger and G. Golub, Linear least squares
!                 solutions by Householder transformations, Numerische
!                 Mathematik  7, (1965), pp. 269-276.
!***ROUTINES CALLED  D1MACH, DCSCAL, DDOT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DORTHR
  DOUBLE PRECISION DDOT, D1MACH
  INTEGER IFLAG, IRANK, ISCALE, J, JROW, K, KP, KPIVOT(*), L, M, &
       MK, N, NRDA
  DOUBLE PRECISION A(NRDA,*), ACC, AKK, ANORM, AS, ASAVE, DIAG(*), &
       DIAGK, DUM, ROWS(*), RS(*), RSS, SAD, SCALES(*), SIG, SIGMA, &
       SRURO, URO
!
!     ******************************************************************
!
!          MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
!          BY THE FUNCTION D1MACH.
!
!     ******************************************************************
!
!***FIRST EXECUTABLE STATEMENT  DORTHR
  URO = D1MACH(4)
  if (M  >=  N .AND. N  >=  1 .AND. NRDA  >=  N) go to 10
     IFLAG = 2
     call XERMSG ('SLATEC', 'DORTHR', 'INVALID INPUT PARAMETERS.', &
        2, 1)
  go to 150
   10 CONTINUE
!
     ACC = 10.0D0*URO
     if (IFLAG  <  0) ACC = MAX(ACC,10.0D0**IFLAG)
     SRURO = SQRT(URO)
     IFLAG = 1
     IRANK = N
!
!        COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM
!
     ANORM = 0.0D0
     DO 20 J = 1, N
        KPIVOT(J) = J
        ROWS(J) = DDOT(M,A(J,1),NRDA,A(J,1),NRDA)
        RS(J) = ROWS(J)
        ANORM = ANORM + ROWS(J)
   20    CONTINUE
!
!        PERFORM COLUMN SCALING ON A WHEN SPECIFIED
!
     call DCSCAL(A,NRDA,N,M,SCALES,DUM,ROWS,RS,ANORM,SCALES,ISCALE, &
                 1)
!
     ANORM = SQRT(ANORM)
!
!
!        CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF
!        ORTHOGONAL TRANSFORMATIONS
!
!
     DO 130 K = 1, N
!           BEGIN BLOCK PERMITTING ...EXITS TO 80
           MK = M - K + 1
!           ...EXIT
           if (K  ==  N) go to 80
           KP = K + 1
!
!              SEARCHING FOR PIVOTAL ROW
!
           DO 60 J = K, N
!                 BEGIN BLOCK PERMITTING ...EXITS TO 50
                 if (ROWS(J)  >=  SRURO*RS(J)) go to 30
                    ROWS(J) = DDOT(MK,A(J,K),NRDA,A(J,K),NRDA)
                    RS(J) = ROWS(J)
   30                CONTINUE
                 if (J  ==  K) go to 40
!                 ......EXIT
                    if (SIGMA  >=  0.99D0*ROWS(J)) go to 50
   40                CONTINUE
                 SIGMA = ROWS(J)
                 JROW = J
   50             CONTINUE
   60          CONTINUE
!           ...EXIT
           if (JROW  ==  K) go to 80
!
!              PERFORM ROW INTERCHANGE
!
           L = KPIVOT(K)
           KPIVOT(K) = KPIVOT(JROW)
           KPIVOT(JROW) = L
           ROWS(JROW) = ROWS(K)
           ROWS(K) = SIGMA
           RSS = RS(K)
           RS(K) = RS(JROW)
           RS(JROW) = RSS
           DO 70 L = 1, M
              ASAVE = A(K,L)
              A(K,L) = A(JROW,L)
              A(JROW,L) = ASAVE
   70          CONTINUE
   80       CONTINUE
!
!           CHECK RANK OF THE MATRIX
!
        SIG = DDOT(MK,A(K,K),NRDA,A(K,K),NRDA)
        DIAGK = SQRT(SIG)
        if (DIAGK  >  ACC*ANORM) go to 90
!
!              RANK DEFICIENT PROBLEM
           IFLAG = 3
           IRANK = K - 1
           call XERMSG ('SLATEC', 'DORTHR', &
              'RANK OF MATRIX IS LESS THAN THE NUMBER OF ROWS.', 1, &
              1)
!        ......EXIT
           go to 140
   90       CONTINUE
!
!           CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
!
        AKK = A(K,K)
        if (AKK  >  0.0D0) DIAGK = -DIAGK
        DIAG(K) = DIAGK
        A(K,K) = AKK - DIAGK
        if (K  ==  N) go to 120
           SAD = DIAGK*AKK - SIG
           DO 110 J = KP, N
              AS = DDOT(MK,A(K,K),NRDA,A(J,K),NRDA)/SAD
              DO 100 L = K, M
                 A(J,L) = A(J,L) + AS*A(K,L)
  100             CONTINUE
              ROWS(J) = ROWS(J) - A(J,K)**2
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140    CONTINUE
  150 CONTINUE
!
!
  return
end
