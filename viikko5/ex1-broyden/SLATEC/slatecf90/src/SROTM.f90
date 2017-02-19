subroutine SROTM (N, SX, INCX, SY, INCY, SPARAM)
!
!! SROTM applies a modified Givens transformation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A8
!***TYPE      SINGLE PRECISION (SROTM-S, DROTM-D)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MODIFIED GIVENS ROTATION, VECTOR
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
!   SPARAM  5-element vector. SPARAM(1) is SFLAG described below.
!           Locations 2-5 of SPARAM contain elements of the
!           transformation matrix H described below.
!
!     --Output--
!       SX  rotated vector (unchanged if N  <=  0)
!       SY  rotated vector (unchanged if N  <=  0)
!
!     Apply the modified Givens transformation, H, to the 2 by N matrix
!     (SX**T)
!     (SY**T) , where **T indicates transpose.  The elements of SX are
!     in SX(LX+I*INCX), I = 0 to N-1, where LX = 1 if INCX  >=  0, else
!     LX = 1+(1-N)*INCX, and similarly for SY using LY and INCY.
!
!     With SPARAM(1)=SFLAG, H has one of the following forms:
!
!     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!
!       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!     H=(          )    (          )    (          )    (          )
!       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!
!     See SROTMG for a description of data storage in SPARAM.
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
!***END PROLOGUE  SROTM
  DIMENSION SX(*), SY(*), SPARAM(5)
  SAVE ZERO, TWO
  DATA ZERO, TWO /0.0E0, 2.0E0/
!***FIRST EXECUTABLE STATEMENT  SROTM
  SFLAG=SPARAM(1)
  if (N <= 0 .OR. (SFLAG+TWO == ZERO)) go to 140
      if (.NOT.(INCX == INCY.AND. INCX  > 0)) go to 70
!
           NSTEPS=N*INCX
           if (SFLAG) 50,10,30
   10          CONTINUE
           SH12=SPARAM(4)
           SH21=SPARAM(3)
                DO 20 I = 1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=W+Z*SH12
                SY(I)=W*SH21+Z
   20               CONTINUE
           go to 140
   30          CONTINUE
           SH11=SPARAM(2)
           SH22=SPARAM(5)
                DO 40 I = 1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=W*SH11+Z
                SY(I)=-W+SH22*Z
   40               CONTINUE
           go to 140
   50          CONTINUE
           SH11=SPARAM(2)
           SH12=SPARAM(4)
           SH21=SPARAM(3)
           SH22=SPARAM(5)
                DO 60 I = 1,NSTEPS,INCX
                W=SX(I)
                Z=SY(I)
                SX(I)=W*SH11+Z*SH12
                SY(I)=W*SH21+Z*SH22
   60               CONTINUE
           go to 140
   70     CONTINUE
      KX=1
      KY=1
      if (INCX  <  0) KX = 1+(1-N)*INCX
      if (INCY  <  0) KY = 1+(1-N)*INCY
!
      if (SFLAG) 120,80,100
   80     CONTINUE
      SH12=SPARAM(4)
      SH21=SPARAM(3)
           DO 90 I = 1,N
           W=SX(KX)
           Z=SY(KY)
           SX(KX)=W+Z*SH12
           SY(KY)=W*SH21+Z
           KX=KX+INCX
           KY=KY+INCY
   90          CONTINUE
      go to 140
  100     CONTINUE
      SH11=SPARAM(2)
      SH22=SPARAM(5)
           DO 110 I = 1,N
           W=SX(KX)
           Z=SY(KY)
           SX(KX)=W*SH11+Z
           SY(KY)=-W+SH22*Z
           KX=KX+INCX
           KY=KY+INCY
  110          CONTINUE
      go to 140
  120     CONTINUE
      SH11=SPARAM(2)
      SH12=SPARAM(4)
      SH21=SPARAM(3)
      SH22=SPARAM(5)
           DO 130 I = 1,N
           W=SX(KX)
           Z=SY(KY)
           SX(KX)=W*SH11+Z*SH12
           SY(KY)=W*SH21+Z*SH22
           KX=KX+INCX
           KY=KY+INCY
  130          CONTINUE
  140     CONTINUE
      return
end
