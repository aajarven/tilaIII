FUNCTION D1MACH (I)
!
!! D1MACH returns floating point machine dependent constants.
!
!***LIBRARY   SLATEC
!***CATEGORY  R1
!***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
!***KEYWORDS  MACHINE CONSTANTS
!***AUTHOR  Fox, P. A., (Bell Labs)
!           Hall, A. D., (Bell Labs)
!           Schryer, N. L., (Bell Labs)
!***DESCRIPTION
!
!   D1MACH can be used to obtain machine-dependent parameters for the
!   local machine environment.  It is a function subprogram with one
!   (input) argument, and can be referenced as follows:
!
!        D = D1MACH(I)
!
!   where I=1,...,5.  The (output) value of D above is determined by
!   the (input) value of I.  The results for various values of I are
!   discussed below.
!
!   D1MACH( 1) = B**(EMIN-1), the smallest positive magnitude.
!   D1MACH( 2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
!   D1MACH( 3) = B**(-T), the smallest relative spacing.
!   D1MACH( 4) = B**(1-T), the largest relative spacing.
!   D1MACH( 5) = LOG10(B)
!
!   Assume double precision numbers are represented in the T-digit,
!   base-B form
!
!              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
!
!   where 0  <=  X(I)  <  B for I=1,...,T, 0  <  X(1), and
!   EMIN  <=  E  <=  EMAX.
!
!   The values of B, T, EMIN and EMAX are provided in I1MACH as
!   follows:
!   I1MACH(10) = B, the base.
!   I1MACH(14) = T, the number of base-B digits.
!   I1MACH(15) = EMIN, the smallest exponent E.
!   I1MACH(16) = EMAX, the largest exponent E.
!
!   To alter this function for a particular environment, the desired
!   set of DATA statements should be activated by removing the C from
!   column 1.  Also, the values of D1MACH(1) - D1MACH(4) should be
!   checked for consistency with the local operating system.
!
!***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
!                 a portable library, ACM Transactions on Mathematical
!                 Software 4, 2 (June 1978), pp. 177-188.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   750101  DATE WRITTEN
!   890213  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900618  Added DEC RISC constants.  (WRB)
!   900723  Added IBM RS 6000 constants.  (WRB)
!   900911  Added SUN 386i constants.  (WRB)
!   910710  Added HP 730 constants.  (SMR)
!   911114  Added Convex IEEE constants.  (WRB)
!   920121  Added SUN -r8 compiler option constants.  (WRB)
!   920229  Added Touchstone Delta i860 constants.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920625  Added CONVEX -p8 and -pd8 compiler option constants.
!           (BKS, WRB)
!   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
!   090212  A. Kuronen, antti.kuronen@helsinki.fi: Use Fortran90 numeric inquiry functions
!   130205  A. Kuronen, antti.kuronen@helsinki.fi: Removed the old data and equivalence statements
!                                                  and unused variables (particularly the log10 
!                                                  which messed up d1mach(5))
!***END PROLOGUE  D1MACH
!
  double precision d1mach
!
  DOUBLE PRECISION DMACH
!
!***FIRST EXECUTABLE STATEMENT  D1MACH
!
  if ( I < 1 .OR. I > 5 ) then
    call XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
  end if

!-----------------------------------------------
! A.Kuronen, 2007, antti.kuronen@helsinki.fi
! Use Fortran90 numeric inquiry functions
!-----------------------------------------------

  select case (i)
     case (1)
        d1mach=tiny(dmach)
     case (2)
        d1mach=huge(dmach)
     case (3)
        d1mach=2.0*epsilon(dmach)
     case (4)
        d1mach=epsilon(dmach)
     case (5)
        d1mach=log10(1.0*radix(dmach))
     case default
        call XERMSG ('SLATEC', 'D1MACH', 'I OUT OF BOUNDS', 1, 2)
  end select


  return
end FUNCTION D1MACH
