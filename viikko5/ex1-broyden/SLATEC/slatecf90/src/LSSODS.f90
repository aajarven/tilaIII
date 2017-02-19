subroutine LSSODS (A, X, B, M, N, NRDA, IFLAG, IRANK, ISCALE, Q, &
     DIAG, KPIVOT, ITER, RESNRM, XNORM, Z, R, DIV, TD, SCALES)
!
!! LSSODS is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LSSODS-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     LSSODS solves the same problem as SODS (in fact, it is called by
!     SODS) but is somewhat more flexible in its use. In particular,
!     LSSODS allows for iterative refinement of the solution, makes the
!     transformation and triangular reduction information more
!     accessible, and enables the user to avoid destruction of the
!     original matrix A.
!
!     Modeled after the ALGOL codes in the articles in the REFERENCES
!     section.
!
! **********************************************************************
!   INPUT
! **********************************************************************
!
!     A -- Contains the matrix of M equations in N unknowns and must
!          be dimensioned NRDA by N. A remains unchanged
!     X -- Solution array of length at least N
!     B -- Given constant vector of length M, B remains unchanged
!     M -- Number of equations, M greater or equal to 1
!     N -- Number of unknowns, N not larger than M
!  NRDA -- Row dimension of A, NRDA greater or equal to M
! IFLAG -- Status indicator
!         = 0 for the first call (and for each new problem defined by
!             a new matrix A) when the matrix data is treated as exact
!         =-K for the first call (and for each new problem defined by
!             a new matrix A) when the matrix data is assumed to be
!             accurate to about K digits
!         = 1 for subsequent calls whenever the matrix A has already
!             been decomposed (problems with new vectors B but
!             same matrix a can be handled efficiently)
! ISCALE -- Scaling indicator
!         =-1 if the matrix A is to be pre-scaled by
!             columns when appropriate
!             If the scaling indicator is not equal to -1
!             no scaling will be attempted
!             For most problems scaling will probably not be necessary
!   ITER -- Maximum number of iterative improvement steps to be
!           performed,  0  <=  ITER  <=  10   (SODS uses ITER=0)
!      Q -- Matrix used for the transformation, must be dimensioned
!           NRDA by N  (SODS puts A in the Q location which conserves
!           storage but destroys A)
!           When iterative improvement of the solution is requested,
!           ITER  >  0, this additional storage for Q must be
!           made available
! DIAG,KPIVOT,Z,R, -- Arrays of length N (except for R which is M)
!   DIV,TD,SCALES     used for internal storage
!
! **********************************************************************
!   OUTPUT
! **********************************************************************
!
!  IFLAG -- Status indicator
!            =1 if solution was obtained
!            =2 if improper input is detected
!            =3 if rank of matrix is less than N
!               if the minimal length least squares solution is
!               desired, simply reset IFLAG=1 and call the code again
!
!       The next three IFLAG values can occur only when
!        the iterative improvement mode is being used.
!            =4 if the problem is ill-conditioned and maximal
!               machine accuracy is not achievable
!            =5 if the problem is very ill-conditioned and the solution
!               IS likely to have no correct digits
!            =6 if the allowable number of iterative improvement steps
!               has been completed without getting convergence
!      X -- Least squares solution of  A X = B
!  IRANK -- Contains the numerically determined matrix rank
!           the user must not alter this value on succeeding calls
!           with input values of IFLAG=1
!      Q -- Contains the strictly upper triangular part of the reduced
!           matrix and the transformation information in the lower
!           triangular part
!   DIAG -- Contains the diagonal elements of the triangular reduced
!           matrix
! KPIVOT -- Contains the pivotal information.  The column interchanges
!           performed on the original matrix are recorded here
!   ITER -- The actual number of iterative corrections used
! RESNRM -- The Euclidean norm of the residual vector  B - A X
!  XNORM -- The Euclidean norm of the solution vector
! DIV,TD -- Contains transformation information for rank
!           deficient problems
! SCALES -- Contains the column scaling parameters
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
!***ROUTINES CALLED  J4SAVE, OHTROR, ORTHOL, R1MACH, SDOT, SDSDOT,
!                    XERMAX, XERMSG, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   910408  Updated the REFERENCES section.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  LSSODS
  DIMENSION A(NRDA,*),X(*),B(*),Q(NRDA,*),DIAG(*), &
            Z(*),KPIVOT(*),R(*),DIV(*),TD(*),SCALES(*)
!
! **********************************************************************
!
!     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
!     THE FUNCTION R1MACH.
!
!***FIRST EXECUTABLE STATEMENT  LSSODS
  URO = R1MACH(3)
!
! **********************************************************************
!
  if (N  <  1  .OR.  M  <  N  .OR.  NRDA  <  M) go to 1
  if (ITER  <  0) go to 1
  if (IFLAG  <=  0) go to 5
  if (IFLAG  ==  1) go to 15
!
!     INVALID INPUT FOR LSSODS
    1 IFLAG=2
  call XERMSG ('SLATEC', 'LSSODS', 'INVALID INPUT PARAMETERS.', 2, &
     1)
  return
!
    5 call XGETF (NFATAL)
  MAXMES = J4SAVE (4,0,.FALSE.)
  if (IFLAG  ==  0) go to 7
  NFAT = -1
  if ( NFATAL  ==  0) NFAT=0
  call XSETF (NFAT)
  call XERMAX (1)
!
!     COPY MATRIX A INTO MATRIX Q
!
    7 DO 10 J=1,N
     DO 10 K=1,M
   10       Q(K,J)=A(K,J)
!
!     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO
!     UPPER TRIANGULAR FORM
!
  call ORTHOL(Q,M,N,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,SCALES,Z,TD)
!
  call XSETF (NFATAL)
  call XERMAX (MAXMES)
  if (IRANK  ==  N) go to 12
!
!     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL
!     TRANSFORMATIONS TO FURTHER REDUCE Q
!
  if (IRANK  /=  0) call OHTROR(Q,N,NRDA,DIAG,IRANK,DIV,TD)
  return
!
!     STORE DIVISORS FOR THE TRIANGULAR SOLUTION
!
   12 DO 13 K=1,N
   13    DIV(K)=DIAG(K)
!
   15 IRM=IRANK-1
  IRP=IRANK+1
  ITERP=MIN(ITER+1,11)
  ACC=10.*URO
!
!     ZERO OUT SOLUTION ARRAY
!
  DO 20 K=1,N
   20    X(K)=0.
!
  if (IRANK  >  0) go to 25
!
!     SPECIAL CASE FOR THE NULL MATRIX
  ITER=0
  XNORM=0.
  RESNRM=SQRT(SDOT(M,B(1),1,B(1),1))
  return
!
!     COPY CONSTANT VECTOR INTO R
!
   25 DO 30 K=1,M
   30    R(K)=B(K)
!
! **********************************************************************
!     SOLUTION SECTION
!     ITERATIVE REFINEMENT OF THE RESIDUAL VECTOR
! **********************************************************************
!
  DO 100 IT=1,ITERP
     ITER=IT-1
!
!        APPLY ORTHOGONAL TRANSFORMATION TO R
!
     DO 35 J=1,IRANK
        MJ=M-J+1
        GAMMA=SDOT(MJ,Q(J,J),1,R(J),1)/(DIAG(J)*Q(J,J))
        DO 35 K=J,M
   35          R(K)=R(K)+GAMMA*Q(K,J)
!
!        BACKWARD SUBSTITUTION FOR TRIANGULAR SYSTEM SOLUTION
!
     Z(IRANK)=R(IRANK)/DIV(IRANK)
     if (IRM  ==  0) go to 45
     DO 40 L=1,IRM
        K=IRANK-L
        KP=K+1
   40       Z(K)=(R(K)-SDOT(L,Q(K,KP),NRDA,Z(KP),1))/DIV(K)
!
   45    if (IRANK  ==  N) go to 60
!
!        FOR RANK DEFICIENT PROBLEMS OBTAIN THE
!        MINIMAL LENGTH SOLUTION
!
     NMIR=N-IRANK
     DO 50 K=IRP,N
   50       Z(K)=0.
     DO 55 K=1,IRANK
        GAM=((TD(K)*Z(K))+SDOT(NMIR,Q(K,IRP),NRDA,Z(IRP),1))/ &
                  (TD(K)*DIV(K))
        Z(K)=Z(K)+GAM*TD(K)
        DO 55 J=IRP,N
   55          Z(J)=Z(J)+GAM*Q(K,J)
!
!        REORDER SOLUTION COMPONENTS ACCORDING TO PIVOTAL POINTS
!        AND RESCALE ANSWERS AS DICTATED
!
   60    DO 65 K=1,N
        Z(K)=Z(K)*SCALES(K)
        L=KPIVOT(K)
   65       X(L)=X(L)+Z(K)
!
!        COMPUTE CORRECTION VECTOR NORM (SOLUTION NORM)
!
     ZNORM=SQRT(SDOT(N,Z(1),1,Z(1),1))
     if (IT  ==  1) XNORM=ZNORM
     if (ITERP  >  1) go to 80
!
!        NO ITERATIVE CORRECTIONS TO BE PERFORMED, SO COMPUTE
!        THE APPROXIMATE RESIDUAL NORM DEFINED BY THE EQUATIONS
!        WHICH ARE NOT SATISFIED BY THE SOLUTION
!        THEN WE ARE DONE
!
     MMIR=M-IRANK
     if (MMIR  ==  0) go to 70
     RESNRM=SQRT(SDOT(MMIR,R(IRP),1,R(IRP),1))
     return
   70    RESNRM=0.
     return
!
!        COMPUTE RESIDUAL VECTOR FOR THE ITERATIVE IMPROVEMENT PROCESS
!
   80    DO 85 K=1,M
   85       R(K)=-SDSDOT(N,-B(K),A(K,1),NRDA,X(1),1)
     RESNRM=SQRT(SDOT(M,R(1),1,R(1),1))
     if (IT  ==  1) go to 100
!
!        TEST FOR CONVERGENCE
!
     if (ZNORM  <=  ACC*XNORM) RETURN
!
!        COMPARE SUCCESSIVE REFINEMENT VECTOR NORMS
!        FOR LOOP TERMINATION CRITERIA
!
     if (ZNORM  <=  0.25*ZNRM0) go to 100
     if (IT  ==  2) go to 90
!
     IFLAG=4
     call XERMSG ('SLATEC', 'LSSODS', &
     'PROBLEM MAY BE ILL-CONDITIONED.  MAXIMAL MACHINE ACCURACY ' // &
     'IS NOT ACHIEVABLE.', 3, 1)
     return
!
   90    IFLAG=5
     call XERMSG ('SLATEC', 'LSSODS', &
        'PROBLEM IS VERY ILL-CONDITIONED.  ITERATIVE ' // &
        'IMPROVEMENT IS INEFFECTIVE.', 8, 1)
     return
!
  100    ZNRM0=ZNORM
! **********************************************************************
!
! **********************************************************************
  IFLAG=6
     call XERMSG ('SLATEC', 'LSSODS', &
        'CONVERGENCE HAS NOT BEEN OBTAINED WITH ALLOWABLE ' // &
        'NUMBER OF ITERATIVE IMPROVEMENT STEPS.', 8, 1)
!
  return
end
