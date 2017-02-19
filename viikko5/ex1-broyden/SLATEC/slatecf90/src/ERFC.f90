function ERFC (X)
!
!! ERFC computes the complementary error function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C8A, L5A1E
!***TYPE      SINGLE PRECISION (ERFC-S, DERFC-D)
!***KEYWORDS  COMPLEMENTARY ERROR FUNCTION, ERFC, FNLIB,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! ERFC(X) calculates the single precision complementary error
! function for single precision argument X.
!
! Series for ERF        on the interval  0.          to  1.00000D+00
!                                        with weighted error   7.10E-18
!                                         log weighted error  17.15
!                               significant figures required  16.31
!                                    decimal places required  17.71
!
! Series for ERFC       on the interval  0.          to  2.50000D-01
!                                        with weighted error   4.81E-17
!                                         log weighted error  16.32
!                        approx significant figures required  15.0
!
!
! Series for ERC2       on the interval  2.50000D-01 to  1.00000D+00
!                                        with weighted error   5.22E-17
!                                         log weighted error  16.28
!                        approx significant figures required  15.0
!                                    decimal places required  16.96
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  ERFC
  DIMENSION ERFCS(13), ERFCCS(24), ERC2CS(23)
  LOGICAL FIRST
  SAVE ERFCS, ERC2CS, ERFCCS, SQRTPI, NTERF, NTERFC, &
   NTERC2, XSML, XMAX, SQEPS, FIRST
  DATA ERFCS( 1) /   -.049046121234691808E0 /
  DATA ERFCS( 2) /   -.14226120510371364E0 /
  DATA ERFCS( 3) /    .010035582187599796E0 /
  DATA ERFCS( 4) /   -.000576876469976748E0 /
  DATA ERFCS( 5) /    .000027419931252196E0 /
  DATA ERFCS( 6) /   -.000001104317550734E0 /
  DATA ERFCS( 7) /    .000000038488755420E0 /
  DATA ERFCS( 8) /   -.000000001180858253E0 /
  DATA ERFCS( 9) /    .000000000032334215E0 /
  DATA ERFCS(10) /   -.000000000000799101E0 /
  DATA ERFCS(11) /    .000000000000017990E0 /
  DATA ERFCS(12) /   -.000000000000000371E0 /
  DATA ERFCS(13) /    .000000000000000007E0 /
  DATA ERC2CS( 1) /   -.069601346602309501E0 /
  DATA ERC2CS( 2) /   -.041101339362620893E0 /
  DATA ERC2CS( 3) /    .003914495866689626E0 /
  DATA ERC2CS( 4) /   -.000490639565054897E0 /
  DATA ERC2CS( 5) /    .000071574790013770E0 /
  DATA ERC2CS( 6) /   -.000011530716341312E0 /
  DATA ERC2CS( 7) /    .000001994670590201E0 /
  DATA ERC2CS( 8) /   -.000000364266647159E0 /
  DATA ERC2CS( 9) /    .000000069443726100E0 /
  DATA ERC2CS(10) /   -.000000013712209021E0 /
  DATA ERC2CS(11) /    .000000002788389661E0 /
  DATA ERC2CS(12) /   -.000000000581416472E0 /
  DATA ERC2CS(13) /    .000000000123892049E0 /
  DATA ERC2CS(14) /   -.000000000026906391E0 /
  DATA ERC2CS(15) /    .000000000005942614E0 /
  DATA ERC2CS(16) /   -.000000000001332386E0 /
  DATA ERC2CS(17) /    .000000000000302804E0 /
  DATA ERC2CS(18) /   -.000000000000069666E0 /
  DATA ERC2CS(19) /    .000000000000016208E0 /
  DATA ERC2CS(20) /   -.000000000000003809E0 /
  DATA ERC2CS(21) /    .000000000000000904E0 /
  DATA ERC2CS(22) /   -.000000000000000216E0 /
  DATA ERC2CS(23) /    .000000000000000052E0 /
  DATA ERFCCS( 1) /   0.0715179310202925E0 /
  DATA ERFCCS( 2) /   -.026532434337606719E0 /
  DATA ERFCCS( 3) /    .001711153977920853E0 /
  DATA ERFCCS( 4) /   -.000163751663458512E0 /
  DATA ERFCCS( 5) /    .000019871293500549E0 /
  DATA ERFCCS( 6) /   -.000002843712412769E0 /
  DATA ERFCCS( 7) /    .000000460616130901E0 /
  DATA ERFCCS( 8) /   -.000000082277530261E0 /
  DATA ERFCCS( 9) /    .000000015921418724E0 /
  DATA ERFCCS(10) /   -.000000003295071356E0 /
  DATA ERFCCS(11) /    .000000000722343973E0 /
  DATA ERFCCS(12) /   -.000000000166485584E0 /
  DATA ERFCCS(13) /    .000000000040103931E0 /
  DATA ERFCCS(14) /   -.000000000010048164E0 /
  DATA ERFCCS(15) /    .000000000002608272E0 /
  DATA ERFCCS(16) /   -.000000000000699105E0 /
  DATA ERFCCS(17) /    .000000000000192946E0 /
  DATA ERFCCS(18) /   -.000000000000054704E0 /
  DATA ERFCCS(19) /    .000000000000015901E0 /
  DATA ERFCCS(20) /   -.000000000000004729E0 /
  DATA ERFCCS(21) /    .000000000000001432E0 /
  DATA ERFCCS(22) /   -.000000000000000439E0 /
  DATA ERFCCS(23) /    .000000000000000138E0 /
  DATA ERFCCS(24) /   -.000000000000000048E0 /
  DATA SQRTPI /1.7724538509055160E0/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  ERFC
  if (FIRST) THEN
     ETA = 0.1*R1MACH(3)
     NTERF = INITS (ERFCS, 13, ETA)
     NTERFC = INITS (ERFCCS, 24, ETA)
     NTERC2 = INITS (ERC2CS, 23, ETA)
!
     XSML = -SQRT (-LOG(SQRTPI*R1MACH(3)))
     TXMAX = SQRT (-LOG(SQRTPI*R1MACH(1)))
     XMAX = TXMAX - 0.5*LOG(TXMAX)/TXMAX - 0.01
     SQEPS = SQRT (2.0*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (X > XSML) go to 20
!
! ERFC(X) = 1.0 - ERF(X) FOR X  <  XSML
!
  ERFC = 2.
  return
!
 20   if (X > XMAX) go to 40
  Y = ABS(X)
  if (Y > 1.0) go to 30
!
! ERFC(X) = 1.0 - ERF(X) FOR -1.  <=  X  <=  1.
!
  if (Y < SQEPS) ERFC = 1.0 - 2.0*X/SQRTPI
  if (Y >= SQEPS) ERFC = 1.0 - &
    X*(1.0 + CSEVL (2.*X*X-1., ERFCS, NTERF) )
  return
!
! ERFC(X) = 1.0 - ERF(X) FOR 1.  <  ABS(X)  <=  XMAX
!
 30   Y = Y*Y
  if (Y <= 4.) ERFC = EXP(-Y)/ABS(X) * (0.5 + CSEVL ((8./Y-5.)/3., &
    ERC2CS, NTERC2) )
  if (Y > 4.) ERFC = EXP(-Y)/ABS(X) * (0.5 + CSEVL (8./Y-1., &
    ERFCCS, NTERFC) )
  if (X < 0.) ERFC = 2.0 - ERFC
  return
!
 40   call XERMSG ('SLATEC', 'ERFC', 'X SO BIG ERFC UNDERFLOWS', 1, 1)
  ERFC = 0.
  return
!
end
