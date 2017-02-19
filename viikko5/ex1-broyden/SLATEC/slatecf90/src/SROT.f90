subroutine SROT (N, SX, INCX, SY, INCY, SC, SS)
!
!! SROT applies a plane Givens rotation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A8
!***TYPE      SINGLE PRECISION (SROT-S, DROT-D, CSROT-C)
!***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
!             LINEAR ALGEBRA, PLANE ROTATION, VECTOR
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
!        N  number of elements in input vector(s)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!       SC  element of rotation matrix
!       SS  element of rotation matrix
!
!     --Output--
!       SX  rotated vector SX (unchanged if N  <=  0)
!       SY  rotated vector SY (unchanged if N  <=  0)
!
!     Multiply the 2 x 2 matrix  ( SC SS) times the 2 x N matrix (SX**T)
!                                (-SS SC)                        (SY**T)
!     where **T indicates transpose.  The elements of SX are in
!     SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX  >=  0, else
!     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
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
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SROT
  REAL SX, SY, SC, SS, ZERO, ONE, W, Z
  DIMENSION SX(*), SY(*)
  SAVE ZERO, ONE
  DATA ZERO, ONE /0.0E0, 1.0E0/
!***FIRST EXECUTABLE STATEMENT  SROT
  if (N  <=  0 .OR. (SS  ==  ZERO .AND. SC  ==  ONE)) go to 40
  if (.NOT. (INCX  ==  INCY .AND. INCX  >  0)) go to 20
!
!          Code for equal and positive increments.
!
       NSTEPS=INCX*N
       DO 10 I = 1,NSTEPS,INCX
            W=SX(I)
            Z=SY(I)
            SX(I)=SC*W+SS*Z
            SY(I)=-SS*W+SC*Z
   10           CONTINUE
       go to 40
!
!     Code for unequal or nonpositive increments.
!
   20 CONTINUE
       KX=1
       KY=1
!
       if (INCX  <  0) KX = 1-(N-1)*INCX
       if (INCY  <  0) KY = 1-(N-1)*INCY
!
       DO 30 I = 1,N
            W=SX(KX)
            Z=SY(KY)
            SX(KX)=SC*W+SS*Z
            SY(KY)=-SS*W+SC*Z
            KX=KX+INCX
            KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
!
  return
end
