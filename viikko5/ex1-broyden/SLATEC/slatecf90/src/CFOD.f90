subroutine CFOD (METH, ELCO, TESCO)
!
!! CFOD is subsidiary to DEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CFOD-S, DCFOD-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   CFOD defines coefficients needed in the integrator package DEBDF
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  CFOD
!
!
!LLL. OPTIMIZE
  INTEGER METH, I, IB, NQ, NQM1, NQP1
  REAL ELCO, TESCO, AGAMQ, FNQ, FNQM1, PC, PINT, RAGQ, &
     RQFAC, RQ1FAC, TSIGN, XPIN
  DIMENSION ELCO(13,12), TESCO(3,12)
!-----------------------------------------------------------------------
! CFOD  IS CALLED BY THE INTEGRATOR ROUTINE TO SET COEFFICIENTS
! NEEDED THERE.  THE COEFFICIENTS FOR THE CURRENT METHOD, AS
! GIVEN BY THE VALUE OF METH, ARE SET FOR ALL ORDERS AND SAVED.
! THE MAXIMUM ORDER ASSUMED HERE IS 12 if METH = 1 AND 5 IF METH = 2.
! (A SMALLER VALUE OF THE MAXIMUM ORDER IS ALSO ALLOWED.)
! CFOD  IS CALLED ONCE AT THE BEGINNING OF THE PROBLEM,
! AND IS NOT CALLED AGAIN UNLESS AND UNTIL METH IS CHANGED.
!
! THE ELCO ARRAY CONTAINS THE BASIC METHOD COEFFICIENTS.
! THE COEFFICIENTS EL(I), 1  <=  I  <=  NQ+1, FOR THE METHOD OF
! ORDER NQ ARE STORED IN ELCO(I,NQ).  THEY ARE GIVEN BY A GENERATING
! POLYNOMIAL, I.E.,
!     L(X) = EL(1) + EL(2)*X + ... + EL(NQ+1)*X**NQ.
! FOR THE IMPLICIT ADAMS METHODS, L(X) IS GIVEN BY
!     DL/DX = (X+1)*(X+2)*...*(X+NQ-1)/FACTORIAL(NQ-1),    L(-1) = 0.
! FOR THE BDF METHODS, L(X) IS GIVEN BY
!     L(X) = (X+1)*(X+2)* ... *(X+NQ)/K,
! WHERE         K = FACTORIAL(NQ)*(1 + 1/2 + ... + 1/NQ).
!
! THE TESCO ARRAY CONTAINS TEST CONSTANTS USED FOR THE
! LOCAL ERROR TEST AND THE SELECTION OF STEP SIZE AND/OR ORDER.
! AT ORDER NQ, TESCO(K,NQ) IS USED FOR THE SELECTION OF STEP
! SIZE AT ORDER NQ - 1 if K = 1, AT ORDER NQ IF K = 2, AND AT ORDER
! NQ + 1 if K = 3.
!-----------------------------------------------------------------------
  DIMENSION PC(12)
!
!***FIRST EXECUTABLE STATEMENT  CFOD
  go to (100, 200), METH
!
 100  ELCO(1,1) = 1.0E0
  ELCO(2,1) = 1.0E0
  TESCO(1,1) = 0.0E0
  TESCO(2,1) = 2.0E0
  TESCO(1,2) = 1.0E0
  TESCO(3,12) = 0.0E0
  PC(1) = 1.0E0
  RQFAC = 1.0E0
  DO 140 NQ = 2,12
!-----------------------------------------------------------------------
! THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
!     P(X) = (X+1)*(X+2)*...*(X+NQ-1).
! INITIALLY, P(X) = 1.
!-----------------------------------------------------------------------
    RQ1FAC = RQFAC
    RQFAC = RQFAC/NQ
    NQM1 = NQ - 1
    FNQM1 = NQM1
    NQP1 = NQ + 1
! FORM COEFFICIENTS OF P(X)*(X+NQ-1). ----------------------------------
    PC(NQ) = 0.0E0
    DO 110 IB = 1,NQM1
      I = NQP1 - IB
 110      PC(I) = PC(I-1) + FNQM1*PC(I)
    PC(1) = FNQM1*PC(1)
! COMPUTE INTEGRAL, -1 TO 0, OF P(X) AND X*P(X). -----------------------
    PINT = PC(1)
    XPIN = PC(1)/2.0E0
    TSIGN = 1.0E0
    DO 120 I = 2,NQ
      TSIGN = -TSIGN
      PINT = PINT + TSIGN*PC(I)/I
 120      XPIN = XPIN + TSIGN*PC(I)/(I+1)
! STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
    ELCO(1,NQ) = PINT*RQ1FAC
    ELCO(2,NQ) = 1.0E0
    DO 130 I = 2,NQ
 130      ELCO(I+1,NQ) = RQ1FAC*PC(I)/I
    AGAMQ = RQFAC*XPIN
    RAGQ = 1.0E0/AGAMQ
    TESCO(2,NQ) = RAGQ
  if ( NQ < 12)TESCO(1,NQP1)=RAGQ*RQFAC/NQP1
    TESCO(3,NQM1) = RAGQ
 140    CONTINUE
  return
!
 200  PC(1) = 1.0E0
  RQ1FAC = 1.0E0
  DO 230 NQ = 1,5
!-----------------------------------------------------------------------
! THE PC ARRAY WILL CONTAIN THE COEFFICIENTS OF THE POLYNOMIAL
!     P(X) = (X+1)*(X+2)*...*(X+NQ).
! INITIALLY, P(X) = 1.
!-----------------------------------------------------------------------
    FNQ = NQ
    NQP1 = NQ + 1
! FORM COEFFICIENTS OF P(X)*(X+NQ). ------------------------------------
    PC(NQP1) = 0.0E0
    DO 210 IB = 1,NQ
      I = NQ + 2 - IB
 210      PC(I) = PC(I-1) + FNQ*PC(I)
    PC(1) = FNQ*PC(1)
! STORE COEFFICIENTS IN ELCO AND TESCO. --------------------------------
    DO 220 I = 1,NQP1
 220      ELCO(I,NQ) = PC(I)/PC(2)
    ELCO(2,NQ) = 1.0E0
    TESCO(1,NQ) = RQ1FAC
    TESCO(2,NQ) = NQP1/ELCO(1,NQ)
    TESCO(3,NQ) = (NQ+2)/ELCO(1,NQ)
    RQ1FAC = RQ1FAC/FNQ
 230    CONTINUE
  return
!----------------------- END OF SUBROUTINE CFOD  -----------------------
end
