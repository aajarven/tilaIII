subroutine ORTHOR (A, N, M, NRDA, IFLAG, IRANK, ISCALE, DIAG, &
     KPIVOT, SCALES, ROWS, RS)
!
!! ORTHOR reduces a matrix to lower triangular form by Householder.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (ORTHOR-S, DORTHR-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   Reduction of the matrix A to lower triangular form by a sequence of
!   orthogonal HOUSEHOLDER transformations post-multiplying A
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
!     N -- Number of rows in the matrix, N greater or equal to 1
!     M -- Number of columns in the matrix, M greater or equal to N
!     IFLAG -- Indicates the uncertainty in the matrix data
!             = 0 when the data is to be treated as exact
!             =-K when the data is assumed to be accurate to about
!                 K digits
!     ISCALE -- Scaling indicator
!               =-1 if the matrix is to be pre-scaled by
!               columns when appropriate.
!               Otherwise no scaling will be attempted
!     NRDA -- Row dimension of A, NRDA greater or equal to N
!     DIAG,KPIVOT,ROWS -- Arrays of length at least N used internally
!         ,RS,SCALES         (except for SCALES which is M)
!
! **********************************************************************
!   OUTPUT
! **********************************************************************
!
!     IFLAG - status indicator
!            =1 for successful decomposition
!            =2 if improper input is detected
!            =3 if rank of the matrix is less than N
!     A -- contains the reduced matrix in the strictly lower triangular
!          part and transformation information
!     IRANK -- contains the numerically determined matrix rank
!     DIAG -- contains the diagonal elements of the reduced
!             triangular matrix
!     KPIVOT -- Contains the pivotal information, the column
!               interchanges performed on the original matrix are
!               recorded here.
!     SCALES -- contains the column scaling parameters
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
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  ORTHOR
  DIMENSION A(NRDA,*),DIAG(*),KPIVOT(*),ROWS(*),RS(*),SCALES(*)
!
! END OF ABSTRACT
!
! **********************************************************************
!
!     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
!     BY THE FUNCTION R1MACH.
!
! **********************************************************************
!
!***FIRST EXECUTABLE STATEMENT  ORTHOR
  URO = R1MACH(4)
  if (M  >=  N  .AND.  N  >=  1  .AND.  NRDA  >=  N) go to 1
  IFLAG=2
  call XERMSG ('SLATEC', 'ORTHOR', 'INVALID INPUT PARAMETERS.', 2, &
     1)
  return
!
    1 ACC=10.*URO
  if (IFLAG  <  0) ACC=MAX(ACC,10.**IFLAG)
  SRURO=SQRT(URO)
  IFLAG=1
  IRANK=N
!
!     COMPUTE NORM**2 OF JTH ROW AND A MATRIX NORM
!
  ANORM=0.
  DO 2 J=1,N
     KPIVOT(J)=J
     ROWS(J)=SDOT(M,A(J,1),NRDA,A(J,1),NRDA)
     RS(J)=ROWS(J)
     ANORM=ANORM+ROWS(J)
    2 CONTINUE
!
!     PERFORM COLUMN SCALING ON A WHEN SPECIFIED
!
  call CSCALE(A,NRDA,N,M,SCALES,DUM,ROWS,RS,ANORM,SCALES,ISCALE,1)
!
  ANORM=SQRT(ANORM)
!
!
!     CONSTRUCTION OF LOWER TRIANGULAR MATRIX AND RECORDING OF
!     ORTHOGONAL TRANSFORMATIONS
!
!
  DO 50 K=1,N
     MK=M-K+1
     if (K  ==  N) go to 25
     KP=K+1
!
!        SEARCHING FOR PIVOTAL ROW
!
     DO 10 J=K,N
        if (ROWS(J)  >=  SRURO*RS(J)) go to 5
        ROWS(J)=SDOT(MK,A(J,K),NRDA,A(J,K),NRDA)
        RS(J)=ROWS(J)
    5       if (J  ==  K) go to 7
        if (SIGMA  >=  0.99*ROWS(J)) go to 10
    7       SIGMA=ROWS(J)
        JROW=J
   10    CONTINUE
     if (JROW  ==  K) go to 25
!
!        PERFORM ROW INTERCHANGE
!
     L=KPIVOT(K)
     KPIVOT(K)=KPIVOT(JROW)
     KPIVOT(JROW)=L
     ROWS(JROW)=ROWS(K)
     ROWS(K)=SIGMA
     RSS=RS(K)
     RS(K)=RS(JROW)
     RS(JROW)=RSS
     DO 20 L=1,M
        ASAVE=A(K,L)
        A(K,L)=A(JROW,L)
   20       A(JROW,L)=ASAVE
!
!        CHECK RANK OF THE MATRIX
!
   25    SIG=SDOT(MK,A(K,K),NRDA,A(K,K),NRDA)
     DIAGK=SQRT(SIG)
     if (DIAGK  >  ACC*ANORM) go to 30
!
!        RANK DEFICIENT PROBLEM
     IFLAG=3
     IRANK=K-1
     call XERMSG ('SLATEC', 'ORTHOR', &
        'RANK OF MATRIX IS LESS THAN THE NUMBER OF ROWS.', 1, 1)
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
        AS=SDOT(MK,A(K,K),NRDA,A(J,K),NRDA)/SAD
        DO 35 L=K,M
   35          A(J,L)=A(J,L)+AS*A(K,L)
   40       ROWS(J)=ROWS(J)-A(J,K)**2
   50 CONTINUE
!
!
  return
end
