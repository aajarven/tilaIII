subroutine OHTROL (Q, N, NRDA, DIAG, IRANK, DIV, TD)
!! OHTROL
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (OHTROL-S, DOHTRL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     For a rank deficient problem, additional orthogonal
!     HOUSEHOLDER transformations are applied to the left side
!     of Q to further reduce the triangular form.
!     Thus, after application of the routines ORTHOR and OHTROL
!     to the original matrix, the result is a nonsingular
!     triangular matrix while the remainder of the matrix
!     has been zeroed out.
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  OHTROL
  DIMENSION Q(NRDA,*),DIAG(*),DIV(*),TD(*)
!***FIRST EXECUTABLE STATEMENT  OHTROL
  NMIR=N-IRANK
  IRP=IRANK+1
  DO 30 K=1,IRANK
     KIR=IRP-K
     DIAGK=DIAG(KIR)
     SIG=(DIAGK*DIAGK)+SDOT(NMIR,Q(IRP,KIR),1,Q(IRP,KIR),1)
     DD=SIGN(SQRT(SIG),-DIAGK)
     DIV(KIR)=DD
     TDV=DIAGK-DD
     TD(KIR)=TDV
     if (K  ==  IRANK) go to 30
     KIRM=KIR-1
     SQD=DD*DIAGK-SIG
     DO 20 J=1,KIRM
        QS=((TDV*Q(KIR,J))+SDOT(NMIR,Q(IRP,J),1,Q(IRP,KIR),1)) &
                 /SQD
        Q(KIR,J)=Q(KIR,J)+QS*TDV
        DO 10 L=IRP,N
   10          Q(L,J)=Q(L,J)+QS*Q(L,KIR)
   20    CONTINUE
   30 CONTINUE
  return
end
