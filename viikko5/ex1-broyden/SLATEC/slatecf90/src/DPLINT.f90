subroutine DPLINT (N, X, Y, C)
!
!! DPLINT produces the polynomial which interpolates a set of discrete ...
!  data points.
!
!***LIBRARY   SLATEC
!***CATEGORY  E1B
!***TYPE      DOUBLE PRECISION (POLINT-S, DPLINT-D)
!***KEYWORDS  POLYNOMIAL INTERPOLATION
!***AUTHOR  Huddleston, R. E., (SNLL)
!***DESCRIPTION
!
!     Abstract
!        Subroutine DPLINT is designed to produce the polynomial which
!     interpolates the data  (X(I),Y(I)), I=1,...,N.  DPLINT sets up
!     information in the array C which can be used by subroutine DPOLVL
!     to evaluate the polynomial and its derivatives and by subroutine
!     DPOLCF to produce the coefficients.
!
!     Formal Parameters
!     *** All TYPE REAL variables are DOUBLE PRECISION ***
!     N  - the number of data points  (N  >=  1)
!     X  - the array of abscissas (all of which must be distinct)
!     Y  - the array of ordinates
!     C  - an array of information used by subroutines
!     *******  Dimensioning Information  *******
!     Arrays X,Y, and C must be dimensioned at least N in the calling
!     program.
!
!***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!                 Curve fitting by polynomials in one variable, Report
!                 SLA-74-0270, Sandia Laboratories, June 1974.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   740601  DATE WRITTEN
!   891006  Cosmetic changes to prologue.  (WRB)
!   891006  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPLINT
  INTEGER I,K,KM1,N
  DOUBLE PRECISION DIF,C(*),X(*),Y(*)
!***FIRST EXECUTABLE STATEMENT  DPLINT
  if (N  <=  0) go to 91
  C(1)=Y(1)
  if ( N  ==  1) RETURN
  DO 10010 K=2,N
  C(K)=Y(K)
  KM1=K-1
  DO 10010 I=1,KM1
!     CHECK FOR DISTINCT X VALUES
  DIF = X(I)-X(K)
  if (DIF  ==  0.0) go to 92
  C(K) = (C(I)-C(K))/DIF
10010 CONTINUE
  return
   91 call XERMSG ('SLATEC', 'DPLINT', 'N IS ZERO OR NEGATIVE.', 2, 1)
  return
   92 call XERMSG ('SLATEC', 'DPLINT', &
     'THE ABSCISSAS ARE NOT DISTINCT.', 2, 1)
  return
end
