subroutine OHTROR (Q, N, NRDA, DIAG, IRANK, DIV, TD)
!
!! OHTROR further reduces a triangular form after ORTHOL has been applied.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (OHTROR-S)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     For a rank deficient problem, additional orthogonal
!     HOUSEHOLDER transformations are applied to the right side
!     of Q to further reduce the triangular form.
!     Thus, after application of the routines ORTHOL and OHTROR
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
!   900402  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  OHTROR
  DIMENSION Q(NRDA,*),DIAG(*),DIV(*),TD(*)
!***FIRST EXECUTABLE STATEMENT  OHTROR
  NMIR=N-IRANK
  IRP=IRANK+1
  DO 30 K=1,IRANK
     KIR=IRP-K
     DIAGK=DIAG(KIR)
     SIG=(DIAGK*DIAGK)+SDOT(NMIR,Q(KIR,IRP),NRDA,Q(KIR,IRP),NRDA)
     DD=SIGN(SQRT(SIG),-DIAGK)
     DIV(KIR)=DD
     TDV=DIAGK-DD
     TD(KIR)=TDV
     if (K  ==  IRANK) go to 30
     KIRM=KIR-1
     SQD=DD*DIAGK-SIG
     DO 20 J=1,KIRM
        QS=((TDV*Q(J,KIR))+SDOT(NMIR,Q(J,IRP),NRDA,Q(KIR,IRP),NRDA)) &
                 /SQD
        Q(J,KIR)=Q(J,KIR)+QS*TDV
        DO 10 L=IRP,N
   10          Q(J,L)=Q(J,L)+QS*Q(KIR,L)
   20    CONTINUE
   30 CONTINUE
  return
end
