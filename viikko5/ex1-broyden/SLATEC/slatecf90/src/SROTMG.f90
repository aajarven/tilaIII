subroutine SROTMG (SD1, SD2, SX1, SY1, SPARAM)
!
!! SROTMG constructs a modified Givens transformation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      SINGLE PRECISION (SROTMG-S, DROTMG-D)
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
!      SD1  single precision scalar
!      SD2  single precision scalar
!      SX1  single precision scalar
!      SY2  single precision scalar
!   SPARAM  S.P. 5-vector. SPARAM(1)=SFLAG defined below.
!           Locations 2-5 contain the rotation matrix.
!
!     --Output--
!      SD1  changed to represent the effect of the transformation
!      SD2  changed to represent the effect of the transformation
!      SX1  changed to represent the effect of the transformation
!      SY2  unchanged
!
!     Construct the modified Givens transformation matrix H which zeros
!     the second component of the 2-vector  (SQRT(SD1)*SX1,SQRT(SD2)*
!     SY2)**T.
!     With SPARAM(1)=SFLAG, H has one of the following forms:
!
!     SFLAG=-1.E0     SFLAG=0.E0        SFLAG=1.E0     SFLAG=-2.E0
!
!       (SH11  SH12)    (1.E0  SH12)    (SH11  1.E0)    (1.E0  0.E0)
!     H=(          )    (          )    (          )    (          )
!       (SH21  SH22),   (SH21  1.E0),   (-1.E0 SH22),   (0.E0  1.E0).
!
!     Locations 2-5 of SPARAM contain SH11, SH21, SH12, and SH22,
!     respectively.  (Values of 1.E0, -1.E0, or 0.E0 implied by the
!     value of SPARAM(1) are not stored in SPARAM.)
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780301  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920316  Prologue corrected.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SROTMG
  DIMENSION SPARAM(5)
  SAVE ZERO, ONE, TWO, GAM, GAMSQ, RGAMSQ
  DATA ZERO, ONE, TWO /0.0E0, 1.0E0, 2.0E0/
  DATA GAM, GAMSQ, RGAMSQ /4096.0E0, 1.67772E7, 5.96046E-8/
!***FIRST EXECUTABLE STATEMENT  SROTMG
  if (.NOT. SD1  <  ZERO) go to 10
!       GO ZERO-H-D-AND-SX1..
      go to 60
   10 CONTINUE
!     CASE-SD1-NONNEGATIVE
  SP2=SD2*SY1
  if (.NOT. SP2  ==  ZERO) go to 20
      SFLAG=-TWO
      go to 260
!     REGULAR-CASE..
   20 CONTINUE
  SP1=SD1*SX1
  SQ2=SP2*SY1
  SQ1=SP1*SX1
!
  if (.NOT. ABS(SQ1)  >  ABS(SQ2)) go to 40
      SH21=-SY1/SX1
      SH12=SP2/SP1
!
      SU=ONE-SH12*SH21
!
      if (.NOT. SU  <=  ZERO) go to 30
!         GO ZERO-H-D-AND-SX1..
           go to 60
   30     CONTINUE
           SFLAG=ZERO
           SD1=SD1/SU
           SD2=SD2/SU
           SX1=SX1*SU
!         GO SCALE-CHECK..
           go to 100
   40 CONTINUE
      if (.NOT. SQ2  <  ZERO) go to 50
!         GO ZERO-H-D-AND-SX1..
           go to 60
   50     CONTINUE
           SFLAG=ONE
           SH11=SP1/SP2
           SH22=SX1/SY1
           SU=ONE+SH11*SH22
           STEMP=SD2/SU
           SD2=SD1/SU
           SD1=STEMP
           SX1=SY1*SU
!         GO SCALE-CHECK
           go to 100
!     PROCEDURE..ZERO-H-D-AND-SX1..
   60 CONTINUE
      SFLAG=-ONE
      SH11=ZERO
      SH12=ZERO
      SH21=ZERO
      SH22=ZERO
!
      SD1=ZERO
      SD2=ZERO
      SX1=ZERO
!         return..
      go to 220
!     PROCEDURE..FIX-H..
   70 CONTINUE
  if (.NOT. SFLAG  >=  ZERO) go to 90
!
      if (.NOT. SFLAG  ==  ZERO) go to 80
      SH11=ONE
      SH22=ONE
      SFLAG=-ONE
      go to 90
   80     CONTINUE
      SH21=-ONE
      SH12=ONE
      SFLAG=-ONE
   90 CONTINUE
  go to IGO,(120,150,180,210)
!     PROCEDURE..SCALE-CHECK
  100 CONTINUE
  110     CONTINUE
      if (.NOT. SD1  <=  RGAMSQ) go to 130
           if (SD1  ==  ZERO) go to 160
           ASSIGN 120 TO IGO
!              FIX-H..
           go to 70
  120          CONTINUE
           SD1=SD1*GAM**2
           SX1=SX1/GAM
           SH11=SH11/GAM
           SH12=SH12/GAM
      go to 110
  130 CONTINUE
  140     CONTINUE
      if (.NOT. SD1  >=  GAMSQ) go to 160
           ASSIGN 150 TO IGO
!              FIX-H..
           go to 70
  150          CONTINUE
           SD1=SD1/GAM**2
           SX1=SX1*GAM
           SH11=SH11*GAM
           SH12=SH12*GAM
      go to 140
  160 CONTINUE
  170     CONTINUE
      if (.NOT. ABS(SD2)  <=  RGAMSQ) go to 190
           if (SD2  ==  ZERO) go to 220
           ASSIGN 180 TO IGO
!              FIX-H..
           go to 70
  180          CONTINUE
           SD2=SD2*GAM**2
           SH21=SH21/GAM
           SH22=SH22/GAM
      go to 170
  190 CONTINUE
  200     CONTINUE
      if (.NOT. ABS(SD2)  >=  GAMSQ) go to 220
           ASSIGN 210 TO IGO
!              FIX-H..
           go to 70
  210          CONTINUE
           SD2=SD2/GAM**2
           SH21=SH21*GAM
           SH22=SH22*GAM
      go to 200
  220 CONTINUE
      if (SFLAG) 250,230,240
  230     CONTINUE
           SPARAM(3)=SH21
           SPARAM(4)=SH12
           go to 260
  240     CONTINUE
           SPARAM(2)=SH11
           SPARAM(5)=SH22
           go to 260
  250     CONTINUE
           SPARAM(2)=SH11
           SPARAM(3)=SH21
           SPARAM(4)=SH12
           SPARAM(5)=SH22
  260 CONTINUE
      SPARAM(1)=SFLAG
      return
end
