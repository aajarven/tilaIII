subroutine DQELG (N, EPSTAB, RESULT, ABSERR, RES3LA, NRES)
!
!! DQELG applies the Epsilon algorithm.
!
!***SUBSIDIARY
!***PURPOSE  The routine determines the limit of a given sequence of
!            approximations, by means of the Epsilon algorithm of
!            P.Wynn. An estimate of the absolute error is also given.
!            The condensed Epsilon table is computed. Only those
!            elements needed for the computation of the next diagonal
!            are preserved.
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QELG-S, DQELG-D)
!***KEYWORDS  CONVERGENCE ACCELERATION, EPSILON ALGORITHM, EXTRAPOLATION
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!           Epsilon algorithm
!           Standard fortran subroutine
!           Double precision version
!
!           PARAMETERS
!              N      - Integer
!                       EPSTAB(N) contains the new element in the
!                       first column of the epsilon table.
!
!              EPSTAB - Double precision
!                       Vector of dimension 52 containing the elements
!                       of the two lower diagonals of the triangular
!                       epsilon table. The elements are numbered
!                       starting at the right-hand corner of the
!                       triangle.
!
!              RESULT - Double precision
!                       Resulting approximation to the integral
!
!              ABSERR - Double precision
!                       Estimate of the absolute error computed from
!                       RESULT and the 3 previous results
!
!              RES3LA - Double precision
!                       Vector of dimension 3 containing the last 3
!                       results
!
!              NRES   - Integer
!                       Number of calls to the routine
!                       (should be zero at first call)
!
!***SEE ALSO  DQAGIE, DQAGOE, DQAGPE, DQAGSE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQELG
!
  DOUBLE PRECISION ABSERR,DELTA1,DELTA2,DELTA3,D1MACH, &
    EPMACH,EPSINF,EPSTAB,ERROR,ERR1,ERR2,ERR3,E0,E1,E1ABS,E2,E3, &
    OFLOW,RES,RESULT,RES3LA,SS,TOL1,TOL2,TOL3
  INTEGER I,IB,IB2,IE,INDX,K1,K2,K3,LIMEXP,N,NEWELM,NRES,NUM
  DIMENSION EPSTAB(52),RES3LA(3)
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!          E0     - THE 4 ELEMENTS ON WHICH THE COMPUTATION OF A NEW
!           E1       ELEMENT IN THE EPSILON TABLE IS BASED
!           E2
!           E3                 E0
!                        E3    E1    NEW
!                              E2
!           NEWELM - NUMBER OF ELEMENTS TO BE COMPUTED IN THE NEW
!                    DIAGONAL
!           ERROR  - ERROR = ABS(E1-E0)+ABS(E2-E1)+ABS(NEW-E2)
!           RESULT - THE ELEMENT IN THE NEW DIAGONAL WITH LEAST VALUE
!                    OF ERROR
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!           LIMEXP IS THE MAXIMUM NUMBER OF ELEMENTS THE EPSILON
!           TABLE CAN CONTAIN. if THIS NUMBER IS REACHED, THE UPPER
!           DIAGONAL OF THE EPSILON TABLE IS DELETED.
!
!***FIRST EXECUTABLE STATEMENT  DQELG
  EPMACH = D1MACH(4)
  OFLOW = D1MACH(2)
  NRES = NRES+1
  ABSERR = OFLOW
  RESULT = EPSTAB(N)
  if ( N < 3) go to 100
  LIMEXP = 50
  EPSTAB(N+2) = EPSTAB(N)
  NEWELM = (N-1)/2
  EPSTAB(N) = OFLOW
  NUM = N
  K1 = N
  DO 40 I = 1,NEWELM
    K2 = K1-1
    K3 = K1-2
    RES = EPSTAB(K1+2)
   E0 = EPSTAB(K3)
    E1 = EPSTAB(K2)
    E2 = RES
    E1ABS = ABS(E1)
    DELTA2 = E2-E1
    ERR2 = ABS(DELTA2)
    TOL2 = MAX(ABS(E2),E1ABS)*EPMACH
    DELTA3 = E1-E0
    ERR3 = ABS(DELTA3)
    TOL3 = MAX(E1ABS,ABS(E0))*EPMACH
    if ( ERR2 > TOL2.OR.ERR3 > TOL3) go to 10
!
!           if E0, E1 AND E2 ARE EQUAL TO WITHIN MACHINE
!           ACCURACY, CONVERGENCE IS ASSUMED.
!           RESULT = E2
!           ABSERR = ABS(E1-E0)+ABS(E2-E1)
!
    RESULT = RES
    ABSERR = ERR2+ERR3
! ***JUMP OUT OF DO-LOOP
    go to 100
   10   E3 = EPSTAB(K1)
    EPSTAB(K1) = E1
    DELTA1 = E1-E3
    ERR1 = ABS(DELTA1)
    TOL1 = MAX(E1ABS,ABS(E3))*EPMACH
!
!           if TWO ELEMENTS ARE VERY CLOSE TO EACH OTHER, OMIT
!           A PART OF THE TABLE BY ADJUSTING THE VALUE OF N
!
    if ( ERR1 <= TOL1.OR.ERR2 <= TOL2.OR.ERR3 <= TOL3) go to 20
    SS = 0.1D+01/DELTA1+0.1D+01/DELTA2-0.1D+01/DELTA3
    EPSINF = ABS(SS*E1)
!
!           TEST TO DETECT IRREGULAR BEHAVIOUR IN THE TABLE, AND
!           EVENTUALLY OMIT A PART OF THE TABLE ADJUSTING THE VALUE
!           OF N.
!
    if ( EPSINF > 0.1D-03) go to 30
   20   N = I+I-1
! ***JUMP OUT OF DO-LOOP
    go to 50
!
!           COMPUTE A NEW ELEMENT AND EVENTUALLY ADJUST
!           THE VALUE OF RESULT.
!
   30   RES = E1+0.1D+01/SS
    EPSTAB(K1) = RES
    K1 = K1-2
    ERROR = ERR2+ABS(RES-E2)+ERR3
    if ( ERROR > ABSERR) go to 40
    ABSERR = ERROR
    RESULT = RES
   40 CONTINUE
!
!           SHIFT THE TABLE.
!
   50 if ( N == LIMEXP) N = 2*(LIMEXP/2)-1
  IB = 1
  if ( (NUM/2)*2 == NUM) IB = 2
  IE = NEWELM+1
  DO 60 I=1,IE
    IB2 = IB+2
    EPSTAB(IB) = EPSTAB(IB2)
    IB = IB2
   60 CONTINUE
  if ( NUM == N) go to 80
  INDX = NUM-N+1
  DO 70 I = 1,N
    EPSTAB(I)= EPSTAB(INDX)
    INDX = INDX+1
   70 CONTINUE
   80 if ( NRES >= 4) go to 90
  RES3LA(NRES) = RESULT
  ABSERR = OFLOW
  go to 100
!
!           COMPUTE ERROR ESTIMATE
!
   90 ABSERR = ABS(RESULT-RES3LA(3))+ABS(RESULT-RES3LA(2)) &
    +ABS(RESULT-RES3LA(1))
  RES3LA(1) = RES3LA(2)
  RES3LA(2) = RES3LA(3)
  RES3LA(3) = RESULT
  100 ABSERR = MAX(ABSERR,0.5D+01*EPMACH*ABS(RESULT))
  return
end
