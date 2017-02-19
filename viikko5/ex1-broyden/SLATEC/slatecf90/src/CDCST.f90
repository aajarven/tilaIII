subroutine CDCST (MAXORD, MINT, ISWFLG, EL, TQ)
!
!! CDCST sets coefficients used by the core integrator CDSTP.
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDCST-S, DDCST-D, CDCST-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!  CDCST is called by CDNTL.  The array EL determines the basic method.
!  The array TQ is involved in adjusting the step size in relation
!  to truncation error.  EL and TQ depend upon MINT, and are calculated
!  for orders 1 to MAXORD( <=  12).  For each order NQ, the coefficients
!  EL are calculated from the generating polynomial:
!    L(T) = EL(1,NQ) + EL(2,NQ)*T + ... + EL(NQ+1,NQ)*T**NQ.
!  For the implicit Adams methods, L(T) is given by
!    dL/dT = (1+T)*(2+T)* ... *(NQ-1+T)/K,   L(-1) = 0,
!    where      K = factorial(NQ-1).
!  For the Gear methods,
!    L(T) = (1+T)*(2+T)* ... *(NQ+T)/K,
!    where      K = factorial(NQ)*(1 + 1/2 + ... + 1/NQ).
!  For each order NQ, there are three components of TQ.
!
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDCST
  REAL EL(13,12), FACTRL(12), GAMMA(14), SUM, TQ(3,12)
  INTEGER I, ISWFLG, J, MAXORD, MINT, MXRD
!***FIRST EXECUTABLE STATEMENT  CDCST
  FACTRL(1) = 1.E0
  DO 10 I = 2,MAXORD
 10     FACTRL(I) = I*FACTRL(I-1)
!                                             Compute Adams coefficients
  if (MINT  ==  1) THEN
    GAMMA(1) = 1.E0
    DO 40 I = 1,MAXORD+1
      SUM = 0.E0
      DO 30 J = 1,I
 30         SUM = SUM - GAMMA(J)/(I-J+2)
 40       GAMMA(I+1) = SUM
    EL(1,1) = 1.E0
    EL(2,1) = 1.E0
    EL(2,2) = 1.E0
    EL(3,2) = 1.E0
    DO 60 J = 3,MAXORD
      EL(2,J) = FACTRL(J-1)
      DO 50 I = 3,J
 50         EL(I,J) = (J-1)*EL(I,J-1) + EL(I-1,J-1)
 60       EL(J+1,J) = 1.E0
    DO 80 J = 2,MAXORD
      EL(1,J) = EL(1,J-1) + GAMMA(J)
      EL(2,J) = 1.E0
      DO 80 I = 3,J+1
 80         EL(I,J) = EL(I,J)/((I-1)*FACTRL(J-1))
    DO 100 J = 1,MAXORD
      TQ(1,J) = -1.E0/(FACTRL(J)*GAMMA(J))
      TQ(2,J) = -1.E0/GAMMA(J+1)
 100      TQ(3,J) = -1.E0/GAMMA(J+2)
!                                              Compute Gear coefficients
  ELSE if (MINT  ==  2) THEN
    EL(1,1) = 1.E0
    EL(2,1) = 1.E0
    DO 130 J = 2,MAXORD
      EL(1,J) = FACTRL(J)
      DO 120 I = 2,J
 120        EL(I,J) = J*EL(I,J-1) + EL(I-1,J-1)
 130      EL(J+1,J) = 1.E0
    SUM = 1.E0
    DO 150 J = 2,MAXORD
      SUM = SUM + 1.E0/J
      DO 150 I = 1,J+1
 150        EL(I,J) = EL(I,J)/(FACTRL(J)*SUM)
    DO 170 J = 1,MAXORD
      if (J  >  1) TQ(1,J) = 1.E0/FACTRL(J-1)
      TQ(2,J) = (J+1)/EL(1,J)
 170      TQ(3,J) = (J+2)/EL(1,J)
  end if
!                          Compute constants used in the stiffness test.
!                          These are the ratio of TQ(2,NQ) for the Gear
!                          methods to those for the Adams methods.
  if (ISWFLG  ==  3) THEN
    MXRD = MIN(MAXORD, 5)
    if (MINT  ==  2) THEN
      GAMMA(1) = 1.E0
      DO 190 I = 1,MXRD
        SUM = 0.E0
        DO 180 J = 1,I
 180          SUM = SUM - GAMMA(J)/(I-J+2)
 190        GAMMA(I+1) = SUM
    end if
    SUM = 1.E0
    DO 200 I = 2,MXRD
      SUM = SUM + 1.E0/I
 200      EL(1+I,1) = -(I+1)*SUM*GAMMA(I+1)
  end if
  return
end
