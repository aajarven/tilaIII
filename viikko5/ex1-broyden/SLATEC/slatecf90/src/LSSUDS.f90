subroutine LSSUDS (A, X, B, N, M, NRDA, U, NRDU, IFLAG, MLSO, &
     IRANK, ISCALE, Q, DIAG, KPIVOT, S, DIV, TD, ISFLG, SCALES)
!
!! LSSUDS is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LSSUDS-S, DLSSUD-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!    LSSUDS solves the underdetermined system of equations  A Z = B,
!    where A is N by M and N  <=  M.  In particular, if rank A equals
!    IRA, a vector X and a matrix U are determined such that X is the
!    UNIQUE solution of smallest length, satisfying A X = B, and the
!    columns of U form an orthonormal basis for the null space of A,
!    satisfying A U = 0 .  Then all solutions Z are given by
!              Z = X + C(1)*U(1) + ..... + C(M-IRA)*U(M-IRA)
!    where U(J) represents the J-th column of U and the C(J) are
!    arbitrary constants.
!    If the system of equations are not compatible, only the least
!    squares solution of minimal length is computed.
!
! *********************************************************************
!   INPUT
! *********************************************************************
!
!     A -- Contains the matrix of N equations in M unknowns, A remains
!          unchanged, must be dimensioned NRDA by M.
!     X -- Solution array of length at least M.
!     B -- Given constant vector of length N, B remains unchanged.
!     N -- Number of equations, N greater or equal to 1.
!     M -- Number of unknowns, M greater or equal to N.
!     NRDA -- Row dimension of A, NRDA greater or equal to N.
!     U -- Matrix used for solution, must be dimensioned NRDU by
!          (M - rank of A).
!          (storage for U may be ignored when only the minimal length
!           solution X is desired)
!     NRDU -- Row dimension of U, NRDU greater or equal to M.
!             (if only the minimal length solution is wanted,
!              NRDU=0 is acceptable)
!     IFLAG -- Status indicator
!           =0  for the first call (and for each new problem defined by
!               a new matrix A) when the matrix data is treated as exact
!           =-K for the first call (and for each new problem defined by
!               a new matrix A) when the matrix data is assumed to be
!               accurate to about K digits.
!           =1  for subsequent calls whenever the matrix A has already
!               been decomposed (problems with new vectors B but
!               same matrix A can be handled efficiently).
!     MLSO -- =0 if only the minimal length solution is wanted.
!             =1 if the complete solution is wanted, includes the
!                linear space defined by the matrix U.
!     IRANK -- Variable used for the rank of A, set by the code.
!     ISCALE -- Scaling indicator
!               =-1 if the matrix A is to be pre-scaled by
!               columns when appropriate.
!               If the scaling indicator is not equal to -1
!               no scaling will be attempted.
!            For most problems scaling will probably not be necessary.
!     Q -- Matrix used for the transformation, must be dimensioned
!            NRDA by M.
!     DIAG,KPIVOT,S, -- Arrays of length at least N used for internal
!      DIV,TD,SCALES    storage (except for SCALES which is M).
!     ISFLG -- Storage for an internal variable.
!
! *********************************************************************
!   OUTPUT
! *********************************************************************
!
!     IFLAG -- Status indicator
!            =1 if solution was obtained.
!            =2 if improper input is detected.
!            =3 if rank of matrix is less than N.
!               To continue, simply reset IFLAG=1 and call LSSUDS again.
!            =4 if the system of equations appears to be inconsistent.
!               However, the least squares solution of minimal length
!               was obtained.
!     X -- Minimal length least squares solution of A Z = B
!     IRANK -- Numerically determined rank of A, must not be altered
!              on succeeding calls with input values of IFLAG=1.
!     U -- Matrix whose M-IRANK columns are mutually orthogonal unit
!          vectors which span the null space of A. This is to be ignored
!          when MLSO was set to zero or IFLAG=4 on output.
!     Q -- Contains the strictly upper triangular part of the reduced
!           matrix and transformation information.
!     DIAG -- Contains the diagonal elements of the triangular reduced
!             matrix.
!     KPIVOT -- Contains the pivotal information.  The row interchanges
!               performed on the original matrix are recorded here.
!     S -- Contains the solution of the lower triangular system.
!     DIV,TD -- Contains transformation information for rank
!               deficient problems.
!     SCALES -- Contains the column scaling parameters.
!
! *********************************************************************
!
!***SEE ALSO  BVSUP
!***REFERENCES  H. A. Watts, Solving linear least squares problems
!                 using SODS/SUDS/CODS, Sandia Report SAND77-0683,
!                 Sandia Laboratories, 1977.
!***ROUTINES CALLED  J4SAVE, OHTROL, ORTHOR, R1MACH, SDOT, XERMAX,
!                    XERMSG, XGETF, XSETF
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Fixed an error message.  (RWC)
!   910408  Updated the AUTHOR and REFERENCES sections.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  LSSUDS
  DIMENSION A(NRDA,*),X(*),B(*),U(NRDU,*),Q(NRDA,*), &
            DIAG(*),KPIVOT(*),S(*),DIV(*),TD(*),SCALES(*)
!
! **********************************************************************
!
!     MACHINE PRECISION (COMPUTER UNIT ROUNDOFF VALUE) IS DEFINED
!     BY THE FUNCTION R1MACH.
!
!***FIRST EXECUTABLE STATEMENT  LSSUDS
  URO = R1MACH(4)
!
! **********************************************************************
!
  if (N  <  1  .OR.  M  <  N  .OR.  NRDA  <  N) go to 1
  if (NRDU  /=  0  .AND.  NRDU  <  M) go to 1
  if (IFLAG  <=  0) go to 5
  if (IFLAG  ==  1) go to 25
!
!     INVALID INPUT FOR LSSUDS
    1 IFLAG=2
  call XERMSG ('SLATEC', 'LSSUDS', 'INVALID INPUT PARAMETERS.', 2, &
     1)
  return
!
    5 call XGETF(NFATAL)
  MAXMES = J4SAVE (4,0,.FALSE.)
  ISFLG=-15
  if (IFLAG  ==  0) go to 7
  ISFLG=IFLAG
  NFAT = -1
  if (NFATAL  ==  0) NFAT=0
  call XSETF(NFAT)
  call XERMAX(1)
!
!     COPY MATRIX A INTO MATRIX Q
!
    7 DO 10 K=1,M
     DO 10 J=1,N
   10       Q(J,K)=A(J,K)
!
!     USE ORTHOGONAL TRANSFORMATIONS TO REDUCE Q TO LOWER
!     TRIANGULAR FORM
!
  call ORTHOR(Q,N,M,NRDA,IFLAG,IRANK,ISCALE,DIAG,KPIVOT,SCALES, &
              DIV,TD)
!
  call XSETF(NFATAL)
  call XERMAX(MAXMES)
  if (IRANK  ==  N) go to 15
!
!     FOR RANK DEFICIENT PROBLEMS USE ADDITIONAL ORTHOGONAL
!     TRANSFORMATIONS TO FURTHER REDUCE Q
!
  if (IRANK  /=  0) call OHTROL(Q,N,NRDA,DIAG,IRANK,DIV,TD)
  return
!
!     STORE DIVISORS FOR THE TRIANGULAR SOLUTION
!
   15 DO 20 K=1,N
   20    DIV(K)=DIAG(K)
!
!
   25 if (IRANK  >  0) go to 40
!
!     SPECIAL CASE FOR THE NULL MATRIX
  DO 35 K=1,M
     X(K)=0.
     if (MLSO  ==  0) go to 35
     U(K,K)=1.
     DO 30 J=1,M
        if (J  ==  K) go to 30
        U(J,K)=0.
   30    CONTINUE
   35 CONTINUE
  DO 37 K=1,N
     if (B(K)  >  0.) IFLAG=4
   37 CONTINUE
  return
!
!     COPY CONSTANT VECTOR INTO S AFTER FIRST INTERCHANGING
!     THE ELEMENTS ACCORDING TO THE PIVOTAL SEQUENCE
!
   40 DO 45 K=1,N
     KP=KPIVOT(K)
   45    X(K)=B(KP)
  DO 50 K=1,N
   50    S(K)=X(K)
!
  IRP=IRANK+1
  NU=1
  if (MLSO  ==  0) NU=0
  if (IRANK  ==  N) go to 60
!
!     FOR RANK DEFICIENT PROBLEMS WE MUST APPLY THE
!     ORTHOGONAL TRANSFORMATION TO S
!     WE ALSO CHECK TO SEE if THE SYSTEM APPEARS TO BE INCONSISTENT
!
  NMIR=N-IRANK
  SS=SDOT(N,S(1),1,S(1),1)
  DO 55 L=1,IRANK
     K=IRP-L
     GAM=((TD(K)*S(K))+SDOT(NMIR,Q(IRP,K),1,S(IRP),1))/ &
               (TD(K)*DIV(K))
     S(K)=S(K)+GAM*TD(K)
     DO 55 J=IRP,N
   55       S(J)=S(J)+GAM*Q(J,K)
  RES=SDOT(NMIR,S(IRP),1,S(IRP),1)
  if (RES  <=  SS*(10.*MAX(10.**ISFLG,10.*URO))**2) go to 60
!
!     INCONSISTENT SYSTEM
  IFLAG=4
  NU=0
!
!     APPLY FORWARD SUBSTITUTION TO SOLVE LOWER TRIANGULAR SYSTEM
!
   60 S(1)=S(1)/DIV(1)
  if (IRANK  ==  1) go to 70
  DO 65 K=2,IRANK
   65    S(K)=(S(K)-SDOT(K-1,Q(K,1),NRDA,S(1),1))/DIV(K)
!
!     INITIALIZE X VECTOR AND THEN APPLY ORTHOGONAL TRANSFORMATION
!
   70 DO 75 K=1,M
     X(K)=0.
     if (K  <=  IRANK) X(K)=S(K)
   75 CONTINUE
!
  DO 80 JR=1,IRANK
     J=IRP-JR
     MJ=M-J+1
     GAMMA=SDOT(MJ,Q(J,J),NRDA,X(J),1)/(DIAG(J)*Q(J,J))
     DO 80 K=J,M
   80       X(K)=X(K)+GAMMA*Q(J,K)
!
!     RESCALE ANSWERS AS DICTATED
!
  DO 85 K=1,M
   85    X(K)=X(K)*SCALES(K)
!
  if ((NU  ==  0) .OR. (M  ==  IRANK)) RETURN
!
!     INITIALIZE U MATRIX AND THEN APPLY ORTHOGONAL TRANSFORMATION
!
  L=M-IRANK

  DO K=1,L

     DO I=1,M
        U(I,K)=0.
        if (I  ==  IRANK+K) U(I,K)=1.
     end do

     DO JR=1,IRANK
        J=IRP-JR
        MJ=M-J+1
        GAMMA=SDOT(MJ,Q(J,J),NRDA,U(J,K),1)/(DIAG(J)*Q(J,J))
        DO I=J,M
           U(I,K)=U(I,K)+GAMMA*Q(J,I)
        end do
    end do

  end do

  return
end
