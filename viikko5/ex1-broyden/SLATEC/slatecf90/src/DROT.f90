subroutine DROT (N, DX, INCX, DY, INCY, DC, DS)
!
!! DROT applies a plane Givens rotation.
!
!***PURPOSE  Apply a plane Givens rotation.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A8
!***TYPE      DOUBLE PRECISION (SROT-S, DROT-D, CSROT-C)
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
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!       DC  D.P. element of rotation matrix
!       DS  D.P. element of rotation matrix
!
!     --Output--
!       DX  rotated vector DX (unchanged if N  <=  0)
!       DY  rotated vector DY (unchanged if N  <=  0)
!
!     Multiply the 2 x 2 matrix  ( DC DS) times the 2 x N matrix (DX**T)
!                                (-DS DC)                        (DY**T)
!     where **T indicates transpose.  The elements of DX are in
!     DX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX  >=  0, else
!     LX = 1+(1-N)*INCX, and similarly for DY using LY and INCY.
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
!***END PROLOGUE  DROT
  DOUBLE PRECISION DX, DY, DC, DS, ZERO, ONE, W, Z
  DIMENSION DX(*), DY(*)
  SAVE ZERO, ONE
  DATA ZERO, ONE /0.0D0, 1.0D0/
!***FIRST EXECUTABLE STATEMENT  DROT
  if (N  <=  0 .OR. (DS  ==  ZERO .AND. DC  ==  ONE)) go to 40
  if (.NOT. (INCX  ==  INCY .AND. INCX  >  0)) go to 20
!
!          Code for equal and positive increments.
!
       NSTEPS=INCX*N
       DO 10 I = 1,NSTEPS,INCX
            W=DX(I)
            Z=DY(I)
            DX(I)=DC*W+DS*Z
            DY(I)=-DS*W+DC*Z
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
            W=DX(KX)
            Z=DY(KY)
            DX(KX)=DC*W+DS*Z
            DY(KY)=-DS*W+DC*Z
            KX=KX+INCX
            KY=KY+INCY
   30           CONTINUE
   40 CONTINUE
!
  return
end
