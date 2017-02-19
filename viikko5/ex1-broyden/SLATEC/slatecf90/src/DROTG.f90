subroutine DROTG (DA, DB, DC, DS)
!
!! DROTG constructs a plane Givens rotation.
!
!***PURPOSE  Construct a plane Givens rotation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      DOUBLE PRECISION (SROTG-S, DROTG-D, CROTG-C)
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
!       DA  double precision scalar
!       DB  double precision scalar
!
!     --Output--
!       DA  double precision result R
!       DB  double precision result Z
!       DC  double precision result
!       DS  double precision result
!
!     Construct the Givens transformation
!
!         ( DC  DS )
!     G = (        ) ,    DC**2 + DS**2 = 1 ,
!         (-DS  DC )
!
!     which zeros the second entry of the 2-vector  (DA,DB)**T .
!
!     The quantity R = (+/-)SQRT(DA**2 + DB**2) overwrites DA in
!     storage.  The value of DB is overwritten by a value Z which
!     allows DC and DS to be recovered by the following algorithm.
!
!           If Z=1  set  DC=0.0  and  DS=1.0
!           If ABS(Z)  <  1  set  DC=SQRT(1-Z**2)  and  DS=Z
!           If ABS(Z)  >  1  set  DC=1/Z  and  DS=SQRT(1-DC**2)
!
!     Normally, the subprogram DROT(N,DX,INCX,DY,INCY,DC,DS) will
!     next be called to apply the transformation to a 2 by N matrix.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DROTG
  DOUBLE PRECISION  DA, DB, DC, DS, U, V, R
!***FIRST EXECUTABLE STATEMENT  DROTG
  if (ABS(DA)  <=  ABS(DB)) go to 10
!
! *** HERE ABS(DA)  >  ABS(DB) ***
!
  U = DA + DA
  V = DB / U
!
!     NOTE THAT U AND R HAVE THE SIGN OF DA
!
  R = SQRT(0.25D0 + V**2) * U
!
!     NOTE THAT DC IS POSITIVE
!
  DC = DA / R
  DS = V * (DC + DC)
  DB = DS
  DA = R
  return
!
! *** HERE ABS(DA)  <=  ABS(DB) ***
!
   10 if (DB  ==  0.0D0) go to 20
  U = DB + DB
  V = DA / U
!
!     NOTE THAT U AND R HAVE THE SIGN OF DB
!     (R IS IMMEDIATELY STORED IN DA)
!
  DA = SQRT(0.25D0 + V**2) * U
!
!     NOTE THAT DS IS POSITIVE
!
  DS = DB / DA
  DC = V * (DS + DS)
  if (DC  ==  0.0D0) go to 15
  DB = 1.0D0 / DC
  return
   15 DB = 1.0D0
  return
!
! *** HERE DA = DB = 0.0 ***
!
   20 DC = 1.0D0
  DS = 0.0D0
  return
!
end
