FUNCTION XPSI (A, IPSIK, IPSIX)
!
!! XPSI computes values of the Psi function for XLEGF.
!
!***LIBRARY   SLATEC
!***CATEGORY  C7C
!***TYPE      SINGLE PRECISION (XPSI-S, DXPSI-D)
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
!***END PROLOGUE  XPSI
  REAL A,B,C,CNUM,CDENOM
  real XPSI
  DIMENSION CNUM(12),CDENOM(12)
  SAVE CNUM, CDENOM
!
!        CNUM(I) AND CDENOM(I) ARE THE ( REDUCED ) NUMERATOR
!        AND 2*I*DENOMINATOR RESPECTIVELY OF THE 2*I TH BERNOULLI
!        NUMBER.
!
  DATA CNUM(1),CNUM(2),CNUM(3),CNUM(4),CNUM(5),CNUM(6),CNUM(7), &
  CNUM(8),CNUM(9),CNUM(10),CNUM(11),CNUM(12) &
      / 1.,     -1.,    1.,     -1., 1., &
     -691.,  1.,     -3617., 43867., -174611., 77683., &
     -236364091./
  DATA CDENOM(1),CDENOM(2),CDENOM(3),CDENOM(4),CDENOM(5),CDENOM(6), &
   CDENOM(7),CDENOM(8),CDENOM(9),CDENOM(10),CDENOM(11),CDENOM(12) &
  /12.,120.,   252.,   240.,132., &
    32760., 12.,  8160., 14364., 6600., 276., 65520./
!***FIRST EXECUTABLE STATEMENT  XPSI
  N=MAX(0,IPSIX-INT(A))
  B=N+A
  K1=IPSIK-1
!
!        SERIES EXPANSION FOR A  >  IPSIX USING IPSIK-1 TERMS.
!
  C=0.
  DO 12 I=1,K1
  K=IPSIK-I
   12 C=(C+CNUM(K)/CDENOM(K))/B**2
  XPSI=LOG(B)-(C+.5/B)
  if ( N == 0) go to 20
  B=0.
!
!        RECURRENCE FOR A  <=  IPSIX.
!
  DO 15 M=1,N
   15 B=B+1./(N-M+A)
  XPSI=XPSI-B
   20 RETURN
end
