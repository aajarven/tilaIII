  DOUBLE PRECISION FUNCTION DXPSI (A, IPSIK, IPSIX)
!
!! DXPSI computes values of the Psi function for DXLEGF.
!
!***LIBRARY   SLATEC
!***CATEGORY  C7C
!***TYPE      DOUBLE PRECISION (XPSI-S, DXPSI-D)
!***KEYWORDS  PSI FUNCTION
!***AUTHOR  Smith, John M., (NBS and George Mason University)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820728  DATE WRITTEN
!   890126  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  DXPSI
  DOUBLE PRECISION A,B,C,CNUM,CDENOM
  DIMENSION CNUM(12),CDENOM(12)
  SAVE CNUM, CDENOM
!
!        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
!        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
!        NUMBER.
!
  DATA CNUM(1),CNUM(2),CNUM(3),CNUM(4),CNUM(5),CNUM(6),CNUM(7), &
  CNUM(8),CNUM(9),CNUM(10),CNUM(11),CNUM(12) &
      / 1.D0,     -1.D0,    1.D0,     -1.D0, 1.D0, &
     -691.D0,  1.D0,     -3617.D0, 43867.D0, -174611.D0, 77683.D0, &
     -236364091.D0/
  DATA CDENOM(1),CDENOM(2),CDENOM(3),CDENOM(4),CDENOM(5),CDENOM(6), &
   CDENOM(7),CDENOM(8),CDENOM(9),CDENOM(10),CDENOM(11),CDENOM(12) &
  /12.D0,120.D0,   252.D0,   240.D0,132.D0, &
    32760.D0, 12.D0,  8160.D0, 14364.D0, 6600.D0, 276.D0, 65520.D0/
!***FIRST EXECUTABLE STATEMENT  DXPSI
  N=MAX(0,IPSIX-INT(A))
  B=N+A
  K1=IPSIK-1
!
!        SERIES EXPANSION FOR A  >  IPSIX USING IPSIK-1 TERMS.
!
  C=0.D0
  DO 12 I=1,K1
  K=IPSIK-I
   12 C=(C+CNUM(K)/CDENOM(K))/B**2
  DXPSI=LOG(B)-(C+.5D0/B)
  if ( N == 0) go to 20
  B=0.D0
!
!        RECURRENCE FOR A  <=  IPSIX.
!
  DO 15 M=1,N
   15 B=B+1.D0/(N-M+A)
  DXPSI=DXPSI-B
   20 RETURN
end
