subroutine CDZRO (AE, F, H, N, NQ, IROOT, RE, T, YH, UROUND, B, C, &
     FB, FC, Y)
!
!! CDZRO searches for a zero of a function F(N, T, Y, IROOT) ...
!            between the given values B and C until the width of the ...
!            interval (B, C) has collapsed to within a tolerance ...
!            specified by the stopping criterion, ...
!              ABS(B - C)  <=  2.*(RW*ABS(B) + AE).
!
!***LIBRARY   SLATEC (SDRIVE)
!***TYPE      COMPLEX (SDZRO-S, DDZRO-D, CDZRO-C)
!***AUTHOR  Kahaner, D. K., (NIST)
!             National Institute of Standards and Technology
!             Gaithersburg, MD  20899
!           Sutherland, C. D., (LANL)
!             Mail Stop D466
!             Los Alamos National Laboratory
!             Los Alamos, NM  87545
!***DESCRIPTION
!
!     This is a special purpose version of ZEROIN, modified for use with
!     the CDRIV package.
!
!     Sandia Mathematical Program Library
!     Mathematical Computing Services Division 5422
!     Sandia Laboratories
!     P. O. Box 5800
!     Albuquerque, New Mexico  87115
!     Control Data 6600 Version 4.5, 1 November 1971
!
!     PARAMETERS
!        F     - Name of the external function, which returns a
!                real result.  This name must be in an
!                EXTERNAL statement in the calling program.
!        B     - One end of the interval (B, C).  The value returned for
!                B usually is the better approximation to a zero of F.
!        C     - The other end of the interval (B, C).
!        RE    - Relative error used for RW in the stopping criterion.
!                If the requested RE is less than machine precision,
!                then RW is set to approximately machine precision.
!        AE    - Absolute error used in the stopping criterion.  If the
!                given interval (B, C) contains the origin, then a
!                nonzero value should be chosen for AE.
!
!***REFERENCES  L. F. Shampine and H. A. Watts, ZEROIN, a root-solving
!                 routine, SC-TM-70-631, Sept 1970.
!               T. J. Dekker, Finding a zero by means of successive
!                 linear interpolation, Constructive Aspects of the
!                 Fundamental Theorem of Algebra, edited by B. Dejon
!                 and P. Henrici, 1969.
!***ROUTINES CALLED  CDNTP
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   900329  Initial submission to SLATEC.
!***END PROLOGUE  CDZRO
  INTEGER IC, IROOT, KOUNT, N, NQ
  COMPLEX Y(*), YH(N,*)
  REAL A, ACBS, ACMB, AE, B, C, CMB, ER, F, FA, FB, FC, &
       H, P, Q, RE, RW, T, TOL, UROUND
!***FIRST EXECUTABLE STATEMENT  CDZRO
  ER = 4.E0*UROUND
  RW = MAX(RE, ER)
  IC = 0
  ACBS = ABS(B - C)
  A = C
  FA = FC
  KOUNT = 0
!                                                    Perform interchange
 10   if (ABS(FC)  <  ABS(FB)) THEN
    A = B
    FA = FB
    B = C
    FB = FC
    C = A
    FC = FA
  end if
  CMB = 0.5E0*(C - B)
  ACMB = ABS(CMB)
  TOL = RW*ABS(B) + AE
!                                                Test stopping criterion
  if (ACMB  <=  TOL) RETURN
  if (KOUNT  >  50) RETURN
!                                    Calculate new iterate implicitly as
!                                    B + P/Q, where we arrange P  >=  0.
!                         The implicit form is used to prevent overflow.
  P = (B - A)*FB
  Q = FA - FB
  if (P  <  0.E0) THEN
    P = -P
    Q = -Q
  end if
!                          Update A and check for satisfactory reduction
!                          in the size of our bounding interval.
  A = B
  FA = FB
  IC = IC + 1
  if (IC  >=  4) THEN
    if (8.E0*ACMB  >=  ACBS) THEN
!                                                                 Bisect
      B = 0.5E0*(C + B)
      go to 20
    end if
    IC = 0
  end if
  ACBS = ACMB
!                                            Test for too small a change
  if (P  <=  ABS(Q)*TOL) THEN
!                                                 Increment by tolerance
    B = B + SIGN(TOL, CMB)
!                                               Root ought to be between
!                                               B and (C + B)/2.
  ELSE if (P  <  CMB*Q) THEN
!                                                            Interpolate
    B = B + P/Q
  ELSE
!                                                                 Bisect
    B = 0.5E0*(C + B)
  end if
!                                             Have completed computation
!                                             for new iterate B.
 20   call CDNTP (H, 0, N, NQ, T, B, YH,  Y)
  FB = F(N, B, Y, IROOT)
  if (N  ==  0) RETURN
  if (FB  ==  0.E0) RETURN
  KOUNT = KOUNT + 1
!
!             Decide whether next step is interpolation or extrapolation
!
  if (SIGN(1.0E0, FB)  ==  SIGN(1.0E0, FC)) THEN
    C = A
    FC = FA
  end if
  go to 10
end
