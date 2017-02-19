subroutine INTYD (T, K, YH, NYH, DKY, IFLAG)
!
!! INTYD is subsidiary to DEBDF.

!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (INTYD-S, DINTYD-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   INTYD approximates the solution and derivatives at T by polynomial
!   interpolation. Must be used in conjunction with the integrator
!   package DEBDF.
! ----------------------------------------------------------------------
! INTYD computes interpolated values of the K-th derivative of the
! dependent variable vector Y, and stores it in DKY.
! This routine is called by DEBDF with K = 0,1 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (see detailed instructions in LSODE usage documentation.)
! ----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is..
!              Q
!  DKY(I)  =  sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
!             J=K
! where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
! communicated by common.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
! ----------------------------------------------------------------------
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  INTYD
!
!LLL. OPTIMIZE
  INTEGER K, NYH, IFLAG, I, IC, IER, IOWND, IOWNS, J, JB, JB2, &
     JJ, JJ1, JP1, JSTART, KFLAG, L, MAXORD, METH, MITER, N, NFE, &
     NJE, NQ, NQU, NST
  REAL T, YH, DKY, &
     ROWND, ROWNS, EL0, H, HMIN, HMXI, HU, TN, UROUND, &
     C, R, S, TP
  DIMENSION YH(NYH,*), DKY(*)
  COMMON /DEBDF1/ ROWND, ROWNS(210), &
     EL0, H, HMIN, HMXI, HU, TN, UROUND, IOWND(14), IOWNS(6), &
     IER, JSTART, KFLAG, L, METH, MITER, MAXORD, N, NQ, NST, NFE, &
     NJE, NQU
!
!***FIRST EXECUTABLE STATEMENT  INTYD
  IFLAG = 0
  if (K  <  0 .OR. K  >  NQ) go to 80
  TP = TN - HU*(1.0E0 + 100.0E0*UROUND)
  if ((T-TP)*(T-TN)  >  0.0E0) go to 90
!
  S = (T - TN)/H
  IC = 1
  if (K  ==  0) go to 15
  JJ1 = L - K
  DO 10 JJ = JJ1,NQ
 10     IC = IC*JJ
 15   C = IC
  DO 20 I = 1,N
 20     DKY(I) = C*YH(I,L)
  if (K  ==  NQ) go to 55
  JB2 = NQ - K
  DO 50 JB = 1,JB2
    J = NQ - JB
    JP1 = J + 1
    IC = 1
    if (K  ==  0) go to 35
    JJ1 = JP1 - K
    DO 30 JJ = JJ1,J
 30       IC = IC*JJ
 35     C = IC
    DO 40 I = 1,N
 40       DKY(I) = C*YH(I,JP1) + S*DKY(I)
 50     CONTINUE
  if (K  ==  0) RETURN
 55   R = H**(-K)
  DO 60 I = 1,N
 60     DKY(I) = R*DKY(I)
  return
!
 80   IFLAG = -1
  return
 90   IFLAG = -2
  return
!----------------------- END OF SUBROUTINE INTYD -----------------------
end
