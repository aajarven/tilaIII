subroutine DINTYD (T, K, YH, NYH, DKY, IFLAG)
!
!! DINTYD approximates the ODE solution at T by polynomial interpolation.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEBDF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (INTYD-S, DINTYD-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DINTYD approximates the solution and derivatives at T by polynomial
!   interpolation. Must be used in conjunction with the integrator
!   package DDEBDF.
! ----------------------------------------------------------------------
! DINTYD computes interpolated values of the K-th derivative of the
! dependent variable vector Y, and stores it in DKY.
! This routine is called by DDEBDF with K = 0,1 and T = TOUT, but may
! also be called by the user for any K up to the current order.
! (see detailed instructions in LSODE usage documentation.)
! ----------------------------------------------------------------------
! The computed values in DKY are gotten by interpolation using the
! Nordsieck history array YH.  This array corresponds uniquely to a
! vector-valued polynomial of degree NQCUR or less, and DKY is set
! to the K-th derivative of this polynomial at T.
! The formula for DKY is..
!              Q
!  DKY(I)  =  Sum  C(J,K) * (T - TN)**(J-K) * H**(-J) * YH(I,J+1)
!             J=K
! where  C(J,K) = J*(J-1)*...*(J-K+1), Q = NQCUR, TN = TCUR, H = HCUR.
! The quantities  NQ = NQCUR, L = NQ+1, N = NEQ, TN, and H are
! communicated by common.  The above sum is done in reverse order.
! IFLAG is returned negative if either K or T is out of bounds.
! ----------------------------------------------------------------------
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DINTYD
!
  INTEGER I, IC, IER, IFLAG, IOWND, IOWNS, J, JB, JB2, JJ, JJ1, &
        JP1, JSTART, K, KFLAG, L, MAXORD, METH, MITER, N, NFE, &
        NJE, NQ, NQU, NST, NYH
  DOUBLE PRECISION C, DKY, EL0, H, HMIN, HMXI, HU, R, ROWND, &
        ROWNS, S, T, TN, TP, UROUND, YH
  DIMENSION YH(NYH,*),DKY(*)
  COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND, &
                  IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER, &
                  MAXORD,N,NQ,NST,NFE,NJE,NQU
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 130
!***FIRST EXECUTABLE STATEMENT  DINTYD
     IFLAG = 0
     if (K  <  0 .OR. K  >  NQ) go to 110
        TP = TN - HU*(1.0D0 + 100.0D0*UROUND)
        if ((T - TP)*(T - TN)  <=  0.0D0) go to 10
           IFLAG = -2
!     .........EXIT
           go to 130
   10       CONTINUE
!
        S = (T - TN)/H
        IC = 1
        if (K  ==  0) go to 30
           JJ1 = L - K
           DO 20 JJ = JJ1, NQ
              IC = IC*JJ
   20          CONTINUE
   30       CONTINUE
        C = IC
        DO 40 I = 1, N
           DKY(I) = C*YH(I,L)
   40       CONTINUE
        if (K  ==  NQ) go to 90
           JB2 = NQ - K
           DO 80 JB = 1, JB2
              J = NQ - JB
              JP1 = J + 1
              IC = 1
              if (K  ==  0) go to 60
                 JJ1 = JP1 - K
                 DO 50 JJ = JJ1, J
                    IC = IC*JJ
   50                CONTINUE
   60             CONTINUE
              C = IC
              DO 70 I = 1, N
                 DKY(I) = C*YH(I,JP1) + S*DKY(I)
   70             CONTINUE
   80          CONTINUE
!     .........EXIT
           if (K  ==  0) go to 130
   90       CONTINUE
        R = H**(-K)
        DO 100 I = 1, N
           DKY(I) = R*DKY(I)
  100       CONTINUE
     go to 120
  110    CONTINUE
!
        IFLAG = -1
  120    CONTINUE
  130 CONTINUE
  return
!     ----------------------- END OF SUBROUTINE DINTYD
!     -----------------------
end
