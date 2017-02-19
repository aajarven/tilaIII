subroutine LA05BS (A, IND, IA, N, IP, IW, W, G, B, TRANS)
!
!! LA05BS is subsidiary to SPLP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LA05BS-S, LA05BD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THIS SUBPROGRAM IS A SLIGHT MODIFICATION OF A SUBPROGRAM
!     FROM THE C. 1979 AERE HARWELL LIBRARY.  THE NAME OF THE
!     CORRESPONDING HARWELL CODE CAN BE OBTAINED BY DELETING
!     THE FINAL LETTER =S= IN THE NAMES USED HERE.
!     REVISED SEP. 13, 1979.
!
!     ROYALTIES HAVE BEEN PAID TO AERE-UK FOR USE OF THEIR CODES
!     IN THE PACKAGE GIVEN HERE.  ANY PRIMARY USAGE OF THE HARWELL
!     SUBROUTINES REQUIRES A ROYALTY AGREEMENT AND PAYMENT BETWEEN
!     THE USER AND AERE-UK.  ANY USAGE OF THE SANDIA WRITTEN CODES
!     SPLP( ) (WHICH USES THE HARWELL SUBROUTINES) IS PERMITTED.
!
! IP(I,1),IP(I,2) POINT TO START OF ROW/COLUMN I OF U.
! IW(I,1),IW(I,2) ARE LENGTHS OF ROW/COL I OF U.
! IW(.,3),IW(.,4) HOLD ROW/COL NUMBERS IN PIVOTAL ORDER.
!
!***SEE ALSO  SPLP
!***ROUTINES CALLED  XERMSG, XSETUN
!***COMMON BLOCKS    LA05DS
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900402  Added TYPE section.  (WRB)
!   920410  Corrected second dimension on IW declaration.  (WRB)
!***END PROLOGUE  LA05BS
  REAL A(IA), B(*), AM, W(*), G, SMALL
  LOGICAL TRANS
  INTEGER IND(IA,2), IW(N,8)
  INTEGER IP(N,2)
  COMMON /LA05DS/ SMALL, LP, LENL, LENU, NCP, LROW, LCOL
!***FIRST EXECUTABLE STATEMENT  LA05BS
  if (G < 0.) go to 130
  KLL = IA - LENL + 1
  if (TRANS) go to 80
!
!     MULTIPLY VECTOR BY INVERSE OF L
  if (LENL <= 0) go to 20
  L1 = IA + 1
  DO 10 KK=1,LENL
     K = L1 - KK
     I = IND(K,1)
     if (B(I) == 0.) go to 10
     J = IND(K,2)
     B(J) = B(J) + A(K)*B(I)
   10 CONTINUE
   20 DO 30 I=1,N
     W(I) = B(I)
     B(I) = 0.
   30 CONTINUE
!
!     MULTIPLY VECTOR BY INVERSE OF U
  N1 = N + 1
  DO 70 II=1,N
     I = N1 - II
     I = IW(I,3)
     AM = W(I)
     KP = IP(I,1)
     if (KP > 0) go to 50
     KP = -KP
     IP(I,1) = KP
     NZ = IW(I,1)
     KL = KP - 1 + NZ
     K2 = KP + 1
     DO 40 K=K2,KL
        J = IND(K,2)
        AM = AM - A(K)*B(J)
   40    CONTINUE
   50    if (AM == 0.) go to 70
     J = IND(KP,2)
     B(J) = AM/A(KP)
     KPC = IP(J,2)
     KL = IW(J,2) + KPC - 1
     if (KL == KPC) go to 70
     K2 = KPC + 1
     DO 60 K=K2,KL
        I = IND(K,1)
        IP(I,1) = -ABS(IP(I,1))
   60    CONTINUE
   70 CONTINUE
  go to 140
!
!     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF U
   80 DO 90 I=1,N
     W(I) = B(I)
     B(I) = 0.
   90 CONTINUE
  DO 110 II=1,N
     I = IW(II,4)
     AM = W(I)
     if (AM == 0.) go to 110
     J = IW(II,3)
     KP = IP(J,1)
     AM = AM/A(KP)
     B(J) = AM
     KL = IW(J,1) + KP - 1
     if (KP == KL) go to 110
     K2 = KP + 1
     DO 100 K=K2,KL
        I = IND(K,2)
        W(I) = W(I) - AM*A(K)
  100    CONTINUE
  110 CONTINUE
!
!     MULTIPLY VECTOR BY INVERSE OF TRANSPOSE OF L
  if (KLL > IA) RETURN
  DO 120 K=KLL,IA
     J = IND(K,2)
     if (B(J) == 0.) go to 120
     I = IND(K,1)
     B(I) = B(I) + A(K)*B(J)
  120 CONTINUE
  go to 140
!
  130 call XSETUN(LP)
  if (LP  >  0) call XERMSG ('SLATEC', 'LA05BS', &
     'EARLIER ENTRY GAVE ERROR RETURN.', -8, 2)
  140 RETURN
end
