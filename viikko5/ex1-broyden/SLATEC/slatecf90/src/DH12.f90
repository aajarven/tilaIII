subroutine DH12 (MODE, LPIVOT, L1, M, U, IUE, UP, C, ICE, ICV, NCV)
!
!! DH12 constructs or applies a Householder transformation.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DHFTI, DLSEI and DWNNLS
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (H12-S, DH12-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!      *** DOUBLE PRECISION VERSION OF H12 ******
!
!     C.L.Lawson and R.J.Hanson, Jet Propulsion Laboratory, 1973 Jun 12
!     to appear in 'Solving Least Squares Problems', Prentice-Hall, 1974
!
!     Construction and/or application of a single
!     Householder transformation..     Q = I + U*(U**T)/B
!
!     MODE    = 1 or 2   to select algorithm  H1  or  H2 .
!     LPIVOT is the index of the pivot element.
!     L1,M   If L1  <=  M   the transformation will be constructed to
!            zero elements indexed from L1 through M.   If L1 GT. M
!            THE SUBROUTINE DOES AN IDENTITY TRANSFORMATION.
!     U(),IUE,UP    On entry to H1 U() contains the pivot vector.
!                   IUE is the storage increment between elements.
!                                       On exit from H1 U() and UP
!                   contain quantities defining the vector U of the
!                   Householder transformation.   On entry to H2 U()
!                   and UP should contain quantities previously computed
!                   by H1.  These will not be modified by H2.
!     C()    On entry to H1 or H2 C() contains a matrix which will be
!            regarded as a set of vectors to which the Householder
!            transformation is to be applied.  On exit C() contains the
!            set of transformed vectors.
!     ICE    Storage increment between elements of vectors in C().
!     ICV    Storage increment between vectors in C().
!     NCV    Number of vectors in C() to be transformed. If NCV  <=  0
!            no operations will be done on C().
!
!***SEE ALSO  DHFTI, DLSEI, DWNNLS
!***ROUTINES CALLED  DAXPY, DDOT, DSWAP
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900911  Added DDOT to DOUBLE PRECISION statement.  (WRB)
!***END PROLOGUE  DH12
  INTEGER I, I2, I3, I4, ICE, ICV, INCR, IUE, J, KL1, KL2, KLP, &
       L1, L1M1, LPIVOT, M, MML1P2, MODE, NCV
  DOUBLE PRECISION B, C, CL, CLINV, ONE, UL1M1, SM, U, UP, DDOT
  DIMENSION U(IUE,*), C(*)
!     BEGIN BLOCK PERMITTING ...EXITS TO 140
!***FIRST EXECUTABLE STATEMENT  DH12
     ONE = 1.0D0
!
!     ...EXIT
     if (0  >=  LPIVOT .OR. LPIVOT  >=  L1 .OR. L1  >  M) go to 140
     CL = ABS(U(1,LPIVOT))
     if (MODE  ==  2) go to 40
!           ****** CONSTRUCT THE TRANSFORMATION. ******
        DO 10 J = L1, M
           CL = MAX(ABS(U(1,J)),CL)
   10       CONTINUE
        if (CL  >  0.0D0) go to 20
!     .........EXIT
           go to 140
   20       CONTINUE
        CLINV = ONE/CL
        SM = (U(1,LPIVOT)*CLINV)**2
        DO 30 J = L1, M
           SM = SM + (U(1,J)*CLINV)**2
   30       CONTINUE
        CL = CL*SQRT(SM)
        if (U(1,LPIVOT)  >  0.0D0) CL = -CL
        UP = U(1,LPIVOT) - CL
        U(1,LPIVOT) = CL
     go to 50
   40    CONTINUE
!        ****** APPLY THE TRANSFORMATION  I+U*(U**T)/B  TO C. ******
!
     if (CL  >  0.0D0) go to 50
!     ......EXIT
        go to 140
   50    CONTINUE
!     ...EXIT
     if (NCV  <=  0) go to 140
     B = UP*U(1,LPIVOT)
!        B  MUST BE NONPOSITIVE HERE.  if B = 0., RETURN.
!
     if (B  <  0.0D0) go to 60
!     ......EXIT
        go to 140
   60    CONTINUE
     B = ONE/B
     MML1P2 = M - L1 + 2
     if (MML1P2  <=  20) go to 80
        L1M1 = L1 - 1
        KL1 = 1 + (L1M1 - 1)*ICE
        KL2 = KL1
        KLP = 1 + (LPIVOT - 1)*ICE
        UL1M1 = U(1,L1M1)
        U(1,L1M1) = UP
        if (LPIVOT  /=  L1M1) call DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
        DO 70 J = 1, NCV
           SM = DDOT(MML1P2,U(1,L1M1),IUE,C(KL1),ICE)
           SM = SM*B
           call DAXPY(MML1P2,SM,U(1,L1M1),IUE,C(KL1),ICE)
           KL1 = KL1 + ICV
   70       CONTINUE
        U(1,L1M1) = UL1M1
!     ......EXIT
        if (LPIVOT  ==  L1M1) go to 140
        KL1 = KL2
        call DSWAP(NCV,C(KL1),ICV,C(KLP),ICV)
     go to 130
   80    CONTINUE
        I2 = 1 - ICV + ICE*(LPIVOT - 1)
        INCR = ICE*(L1 - LPIVOT)
        DO 120 J = 1, NCV
           I2 = I2 + ICV
           I3 = I2 + INCR
           I4 = I3
           SM = C(I2)*UP
           DO 90 I = L1, M
              SM = SM + C(I3)*U(1,I)
              I3 = I3 + ICE
   90          CONTINUE
           if (SM  ==  0.0D0) go to 110
              SM = SM*B
              C(I2) = C(I2) + SM*UP
              DO 100 I = L1, M
                 C(I4) = C(I4) + SM*U(1,I)
                 I4 = I4 + ICE
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
  130    CONTINUE
  140 CONTINUE
  return
end