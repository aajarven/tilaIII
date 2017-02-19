subroutine DOHTRL (Q, N, NRDA, DIAG, IRANK, DIV, TD)
!
!! DOHTRL is subsidiary to DBVSUP and DSUDS.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (OHTROL-S, DOHTRL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     For a rank deficient problem, additional orthogonal
!     HOUSEHOLDER transformations are applied to the left side
!     of Q to further reduce the triangular form.
!     Thus, after application of the routines DORTHR and DOHTRL
!     to the original matrix, the result is a nonsingular
!     triangular matrix while the remainder of the matrix
!     has been zeroed out.
!
!***SEE ALSO  DBVSUP, DSUDS
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DOHTRL
  DOUBLE PRECISION DDOT
  INTEGER IRANK, IRP, J, K, KIR, KIRM, L, N, NMIR, NRDA
  DOUBLE PRECISION DD, DIAG(*), DIAGK, DIV(*), Q(NRDA,*), QS, SIG, &
       SQD, TD(*), TDV
!***FIRST EXECUTABLE STATEMENT  DOHTRL
  NMIR = N - IRANK
  IRP = IRANK + 1
  DO 40 K = 1, IRANK
     KIR = IRP - K
     DIAGK = DIAG(KIR)
     SIG = (DIAGK*DIAGK) + DDOT(NMIR,Q(IRP,KIR),1,Q(IRP,KIR),1)
     DD = SIGN(SQRT(SIG),-DIAGK)
     DIV(KIR) = DD
     TDV = DIAGK - DD
     TD(KIR) = TDV
     if (K  ==  IRANK) go to 30
        KIRM = KIR - 1
        SQD = DD*DIAGK - SIG
        DO 20 J = 1, KIRM
           QS = ((TDV*Q(KIR,J)) &
                 + DDOT(NMIR,Q(IRP,J),1,Q(IRP,KIR),1))/SQD
           Q(KIR,J) = Q(KIR,J) + QS*TDV
           DO 10 L = IRP, N
              Q(L,J) = Q(L,J) + QS*Q(L,KIR)
   10          CONTINUE
   20       CONTINUE
   30    CONTINUE
   40 CONTINUE
  return
end
