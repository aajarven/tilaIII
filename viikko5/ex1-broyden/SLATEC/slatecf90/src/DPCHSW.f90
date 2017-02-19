subroutine DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!! DPCHSW limits excursion from data for DPCHCS.
!
!***LIBRARY   SLATEC (PCHIP)
!***TYPE      DOUBLE PRECISION (PCHSW-S, DPCHSW-D)
!***AUTHOR  Fritsch, F. N., (LLNL)
!***DESCRIPTION
!
!         DPCHSW:  DPCHCS Switch Excursion Limiter.
!
!     Called by  DPCHCS  to adjust D1 and D2 if necessary to insure that
!     the extremum on this interval is not further than DFMAX from the
!     extreme data value.
!
! ----------------------------------------------------------------------
!
!  Calling sequence:
!
!        INTEGER  IEXTRM, IERR
!        DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!
!        call  DPCHSW (DFMAX, IEXTRM, D1, D2, H, SLOPE, IERR)
!
!   Parameters:
!
!     DFMAX -- (input) maximum allowed difference between F(IEXTRM) and
!           the cubic determined by derivative values D1,D2.  (assumes
!           DFMAX > 0.)
!
!     IEXTRM -- (input) index of the extreme data value.  (assumes
!           IEXTRM = 1 or 2 .  Any value  /= 1 is treated as 2.)
!
!     D1,D2 -- (input) derivative values at the ends of the interval.
!           (Assumes D1*D2  <=  0.)
!          (output) may be modified if necessary to meet the restriction
!           imposed by DFMAX.
!
!     H -- (input) interval length.  (Assumes  H > 0.)
!
!     SLOPE -- (input) data slope on the interval.
!
!     IERR -- (output) error flag.  should be zero.
!           If IERR=-1, assumption on D1 and D2 is not satisfied.
!           If IERR=-2, quadratic equation locating extremum has
!                       negative discriminant (should never occur).
!
!    -------
!    WARNING:  This routine does no validity-checking of arguments.
!    -------
!
!  Fortran intrinsics used:  ABS, SIGN, SQRT.
!
!***SEE ALSO  DPCHCS
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   820218  DATE WRITTEN
!   820805  Converted to SLATEC library version.
!   870707  Corrected XERROR calls for d.p. name(s).
!   870707  Replaced DATA statement for SMALL with a use of D1MACH.
!   870813  Minor cosmetic changes.
!   890206  Corrected XERROR calls.
!   890411  1. Added SAVE statements (Vers. 3.2).
!           2. Added DOUBLE PRECISION declaration for D1MACH.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated AUTHOR and DATE WRITTEN sections in prologue.  (WRB)
!   920526  Eliminated possible divide by zero problem.  (FNF)
!   930503  Improved purpose.  (FNF)
!***END PROLOGUE  DPCHSW
!
!**End
!
!  DECLARE ARGUMENTS.
!
  INTEGER  IEXTRM, IERR
  DOUBLE PRECISION  DFMAX, D1, D2, H, SLOPE
!
!  DECLARE LOCAL VARIABLES.
!
  DOUBLE PRECISION  CP, FACT, HPHI, LAMBDA, NU, ONE, PHI, RADCAL, &
                    RHO, SIGMA, SMALL, THAT, THIRD, THREE, TWO, ZERO
  SAVE ZERO, ONE, TWO, THREE, FACT
  SAVE THIRD
  DOUBLE PRECISION  D1MACH
!
  DATA  ZERO /0.D0/,  ONE /1.D0/,  TWO /2.D0/, THREE /3.D0/, &
        FACT /100.D0/
!        THIRD SHOULD BE SLIGHTLY LESS THAN 1/3.
  DATA  THIRD /0.33333D0/
!
!  NOTATION AND GENERAL REMARKS.
!
!     RHO IS THE RATIO OF THE DATA SLOPE TO THE DERIVATIVE BEING TESTED.
!     LAMBDA IS THE RATIO OF D2 TO D1.
!     THAT = T-HAT(RHO) IS THE NORMALIZED LOCATION OF THE EXTREMUM.
!     PHI IS THE NORMALIZED VALUE OF P(X)-F1 AT X = XHAT = X-HAT(RHO),
!           WHERE  THAT = (XHAT - X1)/H .
!        THAT IS, P(XHAT)-F1 = D*H*PHI,  WHERE D=D1 OR D2.
!     SIMILARLY,  P(XHAT)-F2 = D*H*(PHI-RHO) .
!
!      SMALL SHOULD BE A FEW ORDERS OF MAGNITUDE GREATER THAN MACHEPS.
!***FIRST EXECUTABLE STATEMENT  DPCHSW
  SMALL = FACT*D1MACH(4)
!
!  DO MAIN CALCULATION.
!
  if (D1  ==  ZERO)  THEN
!
!        SPECIAL CASE -- D1 == ZERO .
!
!          if D2 IS ALSO ZERO, THIS ROUTINE SHOULD NOT HAVE BEEN CALLED.
     if (D2  ==  ZERO)  go to 5001
!
     RHO = SLOPE/D2
!          EXTREMUM IS OUTSIDE INTERVAL WHEN RHO  >=  1/3 .
     if (RHO  >=  THIRD)  go to 5000
     THAT = (TWO*(THREE*RHO-ONE)) / (THREE*(TWO*RHO-ONE))
     PHI = THAT**2 * ((THREE*RHO-ONE)/THREE)
!
!          CONVERT TO DISTANCE FROM F2 if IEXTRM /= 1 .
     if (IEXTRM  /=  1)  PHI = PHI - RHO
!
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
     HPHI = H * ABS(PHI)
     if (HPHI*ABS(D2)  >  DFMAX)  THEN
!           AT THIS POINT, HPHI > 0, SO DIVIDE IS OK.
        D2 = SIGN (DFMAX/HPHI, D2)
     ENDIF
  ELSE
!
     RHO = SLOPE/D1
     LAMBDA = -D2/D1
     if (D2  ==  ZERO)  THEN
!
!           SPECIAL CASE -- D2 == ZERO .
!
!             EXTREMUM IS OUTSIDE INTERVAL WHEN RHO  >=  1/3 .
        if (RHO  >=  THIRD)  go to 5000
        CP = TWO - THREE*RHO
        NU = ONE - TWO*RHO
        THAT = ONE / (THREE*NU)
     ELSE
        if (LAMBDA  <=  ZERO)  go to 5001
!
!           NORMAL CASE -- D1 AND D2 BOTH NONZERO, OPPOSITE SIGNS.
!
        NU = ONE - LAMBDA - TWO*RHO
        SIGMA = ONE - RHO
        CP = NU + SIGMA
        if (ABS(NU)  >  SMALL)  THEN
           RADCAL = (NU - (TWO*RHO+ONE))*NU + SIGMA**2
           if (RADCAL  <  ZERO)  go to 5002
           THAT = (CP - SQRT(RADCAL)) / (THREE*NU)
        ELSE
           THAT = ONE/(TWO*SIGMA)
        ENDIF
     ENDIF
     PHI = THAT*((NU*THAT - CP)*THAT + ONE)
!
!          CONVERT TO DISTANCE FROM F2 if IEXTRM /= 1 .
     if (IEXTRM  /=  1)  PHI = PHI - RHO
!
!          TEST FOR EXCEEDING LIMIT, AND ADJUST ACCORDINGLY.
     HPHI = H * ABS(PHI)
     if (HPHI*ABS(D1)  >  DFMAX)  THEN
!           AT THIS POINT, HPHI > 0, SO DIVIDE IS OK.
        D1 = SIGN (DFMAX/HPHI, D1)
        D2 = -LAMBDA*D1
     ENDIF
  end if
!
!  NORMAL RETURN.
!
 5000 CONTINUE
  IERR = 0
  return
!
!  ERROR RETURNS.
!
 5001 CONTINUE
!     D1 AND D2 BOTH ZERO, OR BOTH NONZERO AND SAME SIGN.
  IERR = -1
  call XERMSG ('SLATEC', 'DPCHSW', 'D1 AND/OR D2 INVALID', IERR, 1)
  return
!
 5002 CONTINUE
!     NEGATIVE VALUE OF RADICAL (SHOULD NEVER OCCUR).
  IERR = -2
  call XERMSG ('SLATEC', 'DPCHSW', 'NEGATIVE RADICAL', IERR, 1)
  return
!------------- LAST LINE OF DPCHSW FOLLOWS -----------------------------
end
