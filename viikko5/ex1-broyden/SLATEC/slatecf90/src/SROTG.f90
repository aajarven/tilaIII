subroutine SROTG (SA, SB, SC, SS)
!
!! SROTG constructs a plane Givens rotation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      SINGLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
!***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
!             LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!       SA  single precision scalar
!       SB  single precision scalar
!
!     --Output--
!       SA  single precision result R
!       SB  single precision result Z
!       SC  single precision result
!       SS  single precision result
!
!     Construct the Givens transformation
!
!         ( SC  SS )
!     G = (        ) ,    SC**2 + SS**2 = 1 ,
!         (-SS  SC )
!
!     which zeros the second entry of the 2-vector  (SA,SB)**T.
!
!     The quantity R = (+/-)SQRT(SA**2 + SB**2) overwrites SA in
!     storage.  The value of SB is overwritten by a value Z which
!     allows SC and SS to be recovered by the following algorithm:
!
!           If Z=1  set  SC=0.0  and  SS=1.0
!           If ABS(Z)  <  1  set  SC=SQRT(1-Z**2)  and  SS=Z
!           If ABS(Z)  >  1  set  SC=1/Z  and  SS=SQRT(1-SC**2)
!
!     Normally, the subprogram SROT(N,SX,INCX,SY,INCY,SC,SS) will
!     next be called to apply the transformation to a 2 by N matrix.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SROTG
!***FIRST EXECUTABLE STATEMENT  SROTG
  if (ABS(SA)  <=  ABS(SB)) go to 10
!
! *** HERE ABS(SA)  >  ABS(SB) ***
!
  U = SA + SA
  V = SB / U
!
!     NOTE THAT U AND R HAVE THE SIGN OF SA
!
  R = SQRT(0.25E0 + V**2) * U
!
!     NOTE THAT SC IS POSITIVE
!
  SC = SA / R
  SS = V * (SC + SC)
  SB = SS
  SA = R
  return
!
! *** HERE ABS(SA)  <=  ABS(SB) ***
!
   10 if (SB  ==  0.0E0) go to 20
  U = SB + SB
  V = SA / U
!
!     NOTE THAT U AND R HAVE THE SIGN OF SB
!     (R IS IMMEDIATELY STORED IN SA)
!
  SA = SQRT(0.25E0 + V**2) * U
!
!     NOTE THAT SS IS POSITIVE
!
  SS = SB / SA
  SC = V * (SS + SS)
  if (SC  ==  0.0E0) go to 15
  SB = 1.0E0 / SC
  return
   15 SB = 1.0E0
  return
!
! *** HERE SA = SB = 0.0 ***
!
   20 SC = 1.0E0
  SS = 0.0E0
  return
!
end
