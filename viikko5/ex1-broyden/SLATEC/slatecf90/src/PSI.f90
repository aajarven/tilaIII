function PSI (X)
!
!! PSI computes the Psi (or Digamma) function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7C
!***TYPE      SINGLE PRECISION (PSI-S, DPSI-D, CPSI-C)
!***KEYWORDS  DIGAMMA FUNCTION, FNLIB, PSI FUNCTION, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! PSI(X) calculates the psi (or digamma) function for real argument X.
! PSI(X) is the logarithmic derivative of the gamma function of X.
!
! Series for PSI        on the interval  0.          to  1.00000D+00
!                                        with weighted error   2.03E-17
!                                         log weighted error  16.69
!                               significant figures required  16.39
!                                    decimal places required  17.37
!
! Series for APSI       on the interval  0.          to  2.50000D-01
!                                        with weighted error   5.54E-17
!                                         log weighted error  16.26
!                               significant figures required  14.42
!                                    decimal places required  16.86
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  COT, CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  PSI
  DIMENSION PSICS(23), APSICS(16)
  LOGICAL FIRST
  EXTERNAL COT
  SAVE PSICS, APSICS, PI, NTPSI, NTAPSI, XBIG, DXREL, FIRST
  DATA PSICS( 1) /   -.038057080835217922E0 /
  DATA PSICS( 2) /    .49141539302938713E0 /
  DATA PSICS( 3) /   -.056815747821244730E0 /
  DATA PSICS( 4) /    .008357821225914313E0 /
  DATA PSICS( 5) /   -.001333232857994342E0 /
  DATA PSICS( 6) /    .000220313287069308E0 /
  DATA PSICS( 7) /   -.000037040238178456E0 /
  DATA PSICS( 8) /    .000006283793654854E0 /
  DATA PSICS( 9) /   -.000001071263908506E0 /
  DATA PSICS(10) /    .000000183128394654E0 /
  DATA PSICS(11) /   -.000000031353509361E0 /
  DATA PSICS(12) /    .000000005372808776E0 /
  DATA PSICS(13) /   -.000000000921168141E0 /
  DATA PSICS(14) /    .000000000157981265E0 /
  DATA PSICS(15) /   -.000000000027098646E0 /
  DATA PSICS(16) /    .000000000004648722E0 /
  DATA PSICS(17) /   -.000000000000797527E0 /
  DATA PSICS(18) /    .000000000000136827E0 /
  DATA PSICS(19) /   -.000000000000023475E0 /
  DATA PSICS(20) /    .000000000000004027E0 /
  DATA PSICS(21) /   -.000000000000000691E0 /
  DATA PSICS(22) /    .000000000000000118E0 /
  DATA PSICS(23) /   -.000000000000000020E0 /
  DATA APSICS( 1) /   -.0204749044678185E0 /
  DATA APSICS( 2) /   -.0101801271534859E0 /
  DATA APSICS( 3) /    .0000559718725387E0 /
  DATA APSICS( 4) /   -.0000012917176570E0 /
  DATA APSICS( 5) /    .0000000572858606E0 /
  DATA APSICS( 6) /   -.0000000038213539E0 /
  DATA APSICS( 7) /    .0000000003397434E0 /
  DATA APSICS( 8) /   -.0000000000374838E0 /
  DATA APSICS( 9) /    .0000000000048990E0 /
  DATA APSICS(10) /   -.0000000000007344E0 /
  DATA APSICS(11) /    .0000000000001233E0 /
  DATA APSICS(12) /   -.0000000000000228E0 /
  DATA APSICS(13) /    .0000000000000045E0 /
  DATA APSICS(14) /   -.0000000000000009E0 /
  DATA APSICS(15) /    .0000000000000002E0 /
  DATA APSICS(16) /   -.0000000000000000E0 /
  DATA PI     / 3.14159265358979324E0/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  PSI
  if (FIRST) THEN
     NTPSI = INITS (PSICS, 23, 0.1*R1MACH(3))
     NTAPSI = INITS (APSICS, 16, 0.1*R1MACH(3))
!
     XBIG = 1.0/SQRT(R1MACH(3))
     DXREL = SQRT (R1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (Y >= 2.0) go to 30
!
! PSI(X) FOR -2.  <  X  <  2.
!
  N = X
  if (X < 0.) N = N - 1
  Y = X - N
  N = N - 1
  PSI = CSEVL (2.*Y-1., PSICS, NTPSI)
  if (N == 0) RETURN
!
  N = -N
  if (X  ==  0.) call XERMSG ('SLATEC', 'PSI', 'X IS 0', 2, 2)
  if (X  <  0. .AND. X+N-2  ==  0.) call XERMSG ('SLATEC', 'PSI', &
     'X IS A NEGATIVE INTEGER', 3, 2)
  if (X  <  (-0.5) .AND. ABS((X-AINT(X-0.5))/X)  <  DXREL) &
     call XERMSG ('SLATEC', 'PSI', &
     'ANSWER LT HALF PRECISION BECAUSE X TOO NEAR NEGATIVE INTEGER', &
     1, 1)
!
  DO 20 I=1,N
    PSI = PSI - 1.0/(X+I-1)
 20   CONTINUE
  return
!
! PSI(X) FOR ABS(X)  >=  2.
!
 30   AUX = 0.
  if (Y < XBIG) AUX = CSEVL (8./Y**2-1., APSICS, NTAPSI)
  if (X < 0.) PSI = LOG(ABS(X)) - 0.5/X + AUX - PI*COT(PI*X)
  if (X > 0.) PSI = LOG(X) - 0.5/X + AUX
  return
!
end
