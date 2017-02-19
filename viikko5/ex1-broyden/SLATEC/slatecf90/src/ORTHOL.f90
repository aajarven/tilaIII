subroutine ORTHOL (A, M, N, NRDA, IFLAG, IRANK, ISCALE, DIAG, &
     KPIVOT, SCALES, COLS, CS)
!
!! ORTHOL reduces a matrix to upper triangular form by Householder.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (ORTHOL-S)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   Reduction of the matrix A to upper triangular form by a sequence of
!   orthogonal HOUSEHOLDER transformations pre-multiplying A
!
!   Modeled after the ALGOL codes in the articles in the REFERENCES
!   section.
!
! **********************************************************************
!   INPUT
! **********************************************************************
!
!     A -- Contains the matrix to be decomposed, must be dimensioned
!           NRDA by N
!     M -- Number of rows in the matrix, M greater or equal to N
!     N -- Number of columns in the matrix, N greater or equal to 1
!     IFLAG -- Indicates the uncertainty in the matrix data
!             = 0 when the data is to be treated as exact
!             =-K when the data is assumed to be accurate to about
!                 K digits
!     ISCALE -- Scaling indicator
!               =-1 if the matrix A is to be pre-scaled by
!               columns when appropriate.
!               Otherwise no scaling will be attempted
!     NRDA -- Row dimension of A, NRDA greater or equal to M
!     DIAG,KPIVOT,COLS -- Arrays of length at least n used internally
!         ,CS,SCALES
!
! **********************************************************************
!   OUTPUT
! **********************************************************************
!
!     IFLAG - Status indicator
!            =1 for successful decomposition
!            =2 if improper input is detected
!            =3 if rank of the matrix is less than N
!     A -- Contains the reduced matrix in the strictly upper triangular
!          part and transformation information in the lower part
!     IRANK -- Contains the numerically determined matrix rank
!     DIAG -- Contains the diagonal elements of the reduced
!             triangular matrix
!     KPIVOT -- Contains the pivotal information, the column
!               interchanges performed on the original matrix are
!               recorded here.
!     SCALES -- Contains the column scaling parameters
!
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***REFERENCES  G. Golub, Numerical methods for solving linear least
!                 squares problems, Numerische Mathematik 7, (1965),
!                 pp. 206-216.
!               P. Businger and G. Golub, Linear least squares
!                 solutions by Householder transformations, Numerische
!                 Mathematik  7, (1965), pp. 269-276.
!***ROUTINES CALLED  CSCALE, R1MACH, SDOT, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  ORTHOL
  DIMENSION A(NRDA,*),DIAG(*),KPIVOT(*),COLS(*),CS(*),SCALES(*)
!
! **********************************************************************
!
!     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
!     BY THE FUNCTION R1MACH.
!
!***FIRST EXECUTABLE STATEMENT  ORTHOL
  URO = R1MACH(3)
!
! **********************************************************************
!
  if (M  >=  N  .AND.  N  >=  1  .AND.  NRDA  >=  M) go to 1
  IFLAG=2
  call XERMSG ('SLATEC', 'ORTHOL', 'INVALID INPUT PARAMETERS.', 2, &
     1)
  return
!
    1 ACC=10.*URO
  if (IFLAG  <  0) ACC=MAX(ACC,10.**IFLAG)
  SRURO=SQRT(URO)
  IFLAG=1
  IRANK=N
!
!     COMPUTE NORM**2 OF JTH COLUMN AND A MATRIX NORM
!
  ANORM=0.
  DO 2 J=1,N
     KPIVOT(J)=J
     COLS(J)=SDOT(M,A(1,J),1,A(1,J),1)
     CS(J)=COLS(J)
     ANORM=ANORM+COLS(J)
    2 CONTINUE
!
!     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
!
  call CSCALE(A,NRDA,M,N,COLS,CS,DUM,DUM,ANORM,SCALES,ISCALE,0)
!
  ANORM=SQRT(ANORM)
!
!
!     CONSTRUCTION OF UPPER TRIANGULAR MATRIX AND RECORDING OF
!     ORTHOGONAL TRANSFORMATIONS
!
!
  DO 50 K=1,N
     MK=M-K+1
     if (K  ==  N) go to 25
     KP=K+1
!
!        SEARCHING FOR PIVOTAL COLUMN
!
     DO 10 J=K,N
        if (COLS(J)  >=  SRURO*CS(J)) go to 5
        COLS(J)=SDOT(MK,A(K,J),1,A(K,J),1)
        CS(J)=COLS(J)
    5       if (J  ==  K) go to 7
        if (SIGMA  >=  0.99*COLS(J)) go to 10
    7       SIGMA=COLS(J)
        JCOL=J
   10    CONTINUE
     if (JCOL  ==  K) go to 25
!
!        PERFORM COLUMN INTERCHANGE
!
     L=KPIVOT(K)
     KPIVOT(K)=KPIVOT(JCOL)
     KPIVOT(JCOL)=L
     COLS(JCOL)=COLS(K)
     COLS(K)=SIGMA
     CSS=CS(K)
     CS(K)=CS(JCOL)
     CS(JCOL)=CSS
     SC=SCALES(K)
     SCALES(K)=SCALES(JCOL)
     SCALES(JCOL)=SC
     DO 20 L=1,M
        ASAVE=A(L,K)
        A(L,K)=A(L,JCOL)
   20       A(L,JCOL)=ASAVE
!
!        CHECK RANK OF THE MATRIX
!
   25    SIG=SDOT(MK,A(K,K),1,A(K,K),1)
     DIAGK=SQRT(SIG)
     if (DIAGK  >  ACC*ANORM) go to 30
!
!        RANK DEFICIENT PROBLEM
     IFLAG=3
     IRANK=K-1
     call XERMSG ('SLATEC', 'ORTHOL', &
        'RANK OF MATRIX IS LESS THAN THE NUMBER OF COLUMNS.', 1, 1)
     return
!
!        CONSTRUCT AND APPLY TRANSFORMATION TO MATRIX A
!
   30    AKK=A(K,K)
     if (AKK  >  0.) DIAGK=-DIAGK
     DIAG(K)=DIAGK
     A(K,K)=AKK-DIAGK
     if (K  ==  N) go to 50
     SAD=DIAGK*AKK-SIG
     DO 40 J=KP,N
        AS=SDOT(MK,A(K,K),1,A(K,J),1)/SAD
        DO 35 L=K,M
   35          A(L,J)=A(L,J)+AS*A(L,K)
   40       COLS(J)=COLS(J)-A(K,J)**2
   50 CONTINUE
!
!
  return
end
