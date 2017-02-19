subroutine FZERO (F, B, C, R, RE, AE, IFLAG)
!
!! FZERO searches for a zero of a function F(X) in a given interval ...
!  (B,C).  It is designed primarily for problems where F(B)
!  and F(C) have opposite signs.
!
!***LIBRARY   SLATEC
!***CATEGORY  F1B
!***TYPE      SINGLE PRECISION (FZERO-S, DFZERO-D)
!***KEYWORDS  BISECTION, NONLINEAR EQUATIONS, ROOTS, ZEROS
!***AUTHOR  Shampine, L. F., (SNLA)
!           Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     FZERO searches for a zero of a REAL function F(X) between the
!     given REAL values B and C until the width of the interval (B,C)
!     has collapsed to within a tolerance specified by the stopping
!     criterion,
!        ABS(B-C)  <=  2.*(RW*ABS(B)+AE).
!     The method used is an efficient combination of bisection and the
!     secant rule and is due to T. J. Dekker.
!
!     Description Of Arguments
!
!   F     :EXT   - Name of the REAL external function.  This name must
!                  be in an EXTERNAL statement in the calling program.
!                  F must be a function of one REAL argument.
!
!   B     :INOUT - One end of the REAL interval (B,C).  The value
!                  returned for B usually is the better approximation
!                  to a zero of F.
!
!   C     :INOUT - The other end of the REAL interval (B,C)
!
!   R     :OUT   - A (better) REAL guess of a zero of F which could help
!                  in speeding up convergence.  If F(B) and F(R) have
!                  opposite signs, a root will be found in the interval
!                  (B,R); if not, but F(R) and F(C) have opposite signs,
!                  a root will be found in the interval (R,C);
!                  otherwise, the interval (B,C) will be searched for a
!                  possible root.  When no better guess is known, it is
!                  recommended that r be set to B or C, since if R is
!                  not interior to the interval (B,C), it will be
!                  ignored.
!
!   RE    :IN    - Relative error used for RW in the stopping criterion.
!                  If the requested RE is less than machine precision,
!                  then RW is set to approximately machine precision.
!
!   AE    :IN    - Absolute error used in the stopping criterion.  If
!                  the given interval (B,C) contains the origin, then a
!                  nonzero value should be chosen for AE.
!
!   IFLAG :OUT   - A status code.  User must check IFLAG after each
!                  call.  Control returns to the user from FZERO in all
!                  cases.
!
!                1  B is within the requested tolerance of a zero.
!                   The interval (B,C) collapsed to the requested
!                   tolerance, the function changes sign in (B,C), and
!                   F(X) decreased in magnitude as (B,C) collapsed.
!
!                2  F(B) = 0.  However, the interval (B,C) may not have
!                   collapsed to the requested tolerance.
!
!                3  B may be near a singular point of F(X).
!                   The interval (B,C) collapsed to the requested tol-
!                   erance and the function changes sign in (B,C), but
!                   F(X) increased in magnitude as (B,C) collapsed, i.e.
!                     ABS(F(B out))  >  MAX(ABS(F(B in)),ABS(F(C in)))
!
!                4  No change in sign of F(X) was found although the
!                   interval (B,C) collapsed to the requested tolerance.
!                   The user must examine this case and decide whether
!                   B is near a local minimum of F(X), or B is near a
!                   zero of even multiplicity, or neither of these.
!
!                5  Too many ( >  500) function evaluations used.
!
!***REFERENCES  L. F. Shampine and H. A. Watts, FZERO, a root-solving
!                 code, Report SC-TM-70-631, Sandia Laboratories,
!                 September 1970.
!               T. J. Dekker, Finding a zero by means of successive
!                 linear interpolation, Constructive Aspects of the
!                 Fundamental Theorem of Algebra, edited by B. Dejon
!                 and P. Henrici, Wiley-Interscience, 1969.
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   700901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  FZERO
  REAL A,ACBS,ACMB,AE,AW,B,C,CMB,ER,FA,FB,FC,FX,FZ,P,Q,R, &
       RE,RW,T,TOL,Z
  INTEGER IC,IFLAG,KOUNT
!***FIRST EXECUTABLE STATEMENT  FZERO
!
!   ER is two times the computer unit roundoff value which is defined
!   here by the function R1MACH.
!
  ER = 2.0E0 * R1MACH(4)
!
!   Initialize.
!
  Z = R
  if (R  <=  MIN(B,C)  .OR.  R  >=  MAX(B,C)) Z = C
  RW = MAX(RE,ER)
  AW = MAX(AE,0.E0)
  IC = 0
  T = Z
  FZ = F(T)
  FC = FZ
  T = B
  FB = F(T)
  KOUNT = 2
  if (SIGN(1.0E0,FZ)  ==  SIGN(1.0E0,FB)) go to 1
  C = Z
  go to 2
    1 if (Z  ==  C) go to 2
  T = C
  FC = F(T)
  KOUNT = 3
  if (SIGN(1.0E0,FZ)  ==  SIGN(1.0E0,FC)) go to 2
  B = Z
  FB = FZ
    2 A = C
  FA = FC
  ACBS = ABS(B-C)
  FX = MAX(ABS(FB),ABS(FC))
!
    3 if (ABS(FC)  >=  ABS(FB)) go to 4
!
!   Perform interchange.
!
  A = B
  FA = FB
  B = C
  FB = FC
  C = A
  FC = FA
!
    4 CMB = 0.5E0*(C-B)
  ACMB = ABS(CMB)
  TOL = RW*ABS(B) + AW
!
!   Test stopping criterion and function count.
!
  if (ACMB  <=  TOL) go to 10
  if (FB  ==  0.E0) go to 11
  if (KOUNT  >=  500) go to 14
!
!   Calculate new iterate implicitly as B+P/Q, where we arrange
!   P  >=  0.  The implicit form is used to prevent overflow.
!
  P = (B-A)*FB
  Q = FA - FB
  if (P  >=  0.E0) go to 5
  P = -P
  Q = -Q
!
!   Update A and check for satisfactory reduction in the size of the
!   bracketing interval.  If not, perform bisection.
!
    5 A = B
  FA = FB
  IC = IC + 1
  if (IC  <  4) go to 6
  if (8.0E0*ACMB  >=  ACBS) go to 8
  IC = 0
  ACBS = ACMB
!
!   Test for too small a change.
!
    6 if (P  >  ABS(Q)*TOL) go to 7
!
!   Increment by TOLerance.
!
  B = B + SIGN(TOL,CMB)
  go to 9
!
!   Root ought to be between B and (C+B)/2.
!
    7 if (P  >=  CMB*Q) go to 8
!
!   Use secant rule.
!
  B = B + P/Q
  go to 9
!
!   Use bisection (C+B)/2.
!
    8 B = B + CMB
!
!   Have completed computation for new iterate B.
!
    9 T = B
  FB = F(T)
  KOUNT = KOUNT + 1
!
!   Decide whether next step is interpolation or extrapolation.
!
  if (SIGN(1.0E0,FB)  /=  SIGN(1.0E0,FC)) go to 3
  C = A
  FC = FA
  go to 3
!
!   Finished.  Process results for proper setting of IFLAG.
!
   10 if (SIGN(1.0E0,FB)  ==  SIGN(1.0E0,FC)) go to 13
  if (ABS(FB)  >  FX) go to 12
  IFLAG = 1
  return
   11 IFLAG = 2
  return
   12 IFLAG = 3
  return
   13 IFLAG = 4
  return
   14 IFLAG = 5
  return
end
