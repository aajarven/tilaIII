subroutine DCFOD (METH, ELCO, TESCO)
!
!! DCFOD is subsidiary to DDEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (CFOD-S, DCFOD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   DCFOD defines coefficients needed in the integrator package DDEBDF
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DCFOD
!
!
  INTEGER I, IB, METH, NQ, NQM1, NQP1
  DOUBLE PRECISION AGAMQ, ELCO, FNQ, FNQM1, PC, PINT, RAGQ, &
        RQ1FAC, RQFAC, TESCO, TSIGN, XPIN
  DIMENSION ELCO(13,12),TESCO(3,12)
!     ------------------------------------------------------------------
!      DCFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
!      NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
!      GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
!      THE MAXIMUM ORDER ASSUMED HERE IS 12 if METH = 1 AND 5 IF METH =
!      2.  (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
!      DCFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
!      AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
!
!      THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
!      THE COEFFICIENTS EL(I), 1  <=  I  <=  NQ+1, FOR THE METHOD OF
!      ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A
!      GENERATING POLYNOMIAL, I.E.,
!          L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
!      FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
!          DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) =
!      0.  FOR THE BDF METHODS, L(X) IS GIVEN BY
!          L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
!      WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
!
!      THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
!      LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
!      AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
!      SIZE AT ORDER NQ - 1 if K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
!      NQ + 1 if K = 3.
!     ------------------------------------------------------------------
  DIMENSION PC(12)
!
!***FIRST EXECUTABLE STATEMENT  DCFOD
  go to (10,60), METH
!
   10 CONTINUE
     ELCO(1,1) = 1.0D0
     ELCO(2,1) = 1.0D0
     TESCO(1,1) = 0.0D0
     TESCO(2,1) = 2.0D0
     TESCO(1,2) = 1.0D0
     TESCO(3,12) = 0.0D0
     PC(1) = 1.0D0
     RQFAC = 1.0D0
     DO 50 NQ = 2, 12
!           ------------------------------------------------------------
!            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
!                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ-1).
!            INITIALLY, P(X) = 1.
!           ------------------------------------------------------------
        RQ1FAC = RQFAC
        RQFAC = RQFAC/NQ
        NQM1 = NQ - 1
        FNQM1 = NQM1
        NQP1 = NQ + 1
!           FORM COEFFICIENTS OF P(X)*(X+NQ-1).
!           ----------------------------------
        PC(NQ) = 0.0D0
        DO 20 IB = 1, NQM1
           I = NQP1 - IB
           PC(I) = PC(I-1) + FNQM1*PC(I)
   20       CONTINUE
        PC(1) = FNQM1*PC(1)
!           COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X).
!           -----------------------
        PINT = PC(1)
        XPIN = PC(1)/2.0D0
        TSIGN = 1.0D0
        DO 30 I = 2, NQ
           TSIGN = -TSIGN
           PINT = PINT + TSIGN*PC(I)/I
           XPIN = XPIN + TSIGN*PC(I)/(I+1)
   30       CONTINUE
!           STORE COEFFICIENTS IN ELCO AND TESCO.
!           --------------------------------
        ELCO(1,NQ) = PINT*RQ1FAC
        ELCO(2,NQ) = 1.0D0
        DO 40 I = 2, NQ
           ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
   40       CONTINUE
        AGAMQ = RQFAC*XPIN
        RAGQ = 1.0D0/AGAMQ
        TESCO(2,NQ) = RAGQ
        if (NQ  <  12) TESCO(1,NQP1) = RAGQ*RQFAC/NQP1
        TESCO(3,NQM1) = RAGQ
   50    CONTINUE
  go to 100
!
   60 CONTINUE
     PC(1) = 1.0D0
     RQ1FAC = 1.0D0
     DO 90 NQ = 1, 5
!           ------------------------------------------------------------
!            THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE
!                POLYNOMIAL P(X) = (X+1)*(X+2)*...*(X+NQ).
!            INITIALLY, P(X) = 1.
!           ------------------------------------------------------------
        FNQ = NQ
        NQP1 = NQ + 1
!           FORM COEFFICIENTS OF P(X)*(X+NQ).
!           ------------------------------------
        PC(NQP1) = 0.0D0
        DO 70 IB = 1, NQ
           I = NQ + 2 - IB
           PC(I) = PC(I-1) + FNQ*PC(I)
   70       CONTINUE
        PC(1) = FNQ*PC(1)
!           STORE COEFFICIENTS IN ELCO AND TESCO.
!           --------------------------------
        DO 80 I = 1, NQP1
           ELCO(I,NQ) = PC(I)/PC(2)
   80       CONTINUE
        ELCO(2,NQ) = 1.0D0
        TESCO(1,NQ) = RQ1FAC
        TESCO(2,NQ) = NQP1/ELCO(1,NQ)
        TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
        RQ1FAC = RQ1FAC/FNQ
   90    CONTINUE
  100 CONTINUE
  return
!     ----------------------- END OF SUBROUTINE DCFOD
!     -----------------------
end
