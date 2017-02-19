subroutine BVDER (X, Y, YP, G, IPAR)
!
!! BVDER is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BVDER-S, DBVDER-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!     NFC = Number of base solution vectors
!
!     NCOMP = Number of components per solution vector
!
!              1 -- Nonzero particular solution
!     INHOMO =
!              2 or 3 -- Zero particular solution
!
!             0 -- Inhomogeneous vector term G(X) identically zero
!     IGOFX =
!             1 -- Inhomogeneous vector term G(X) not identically zero
!
!     G = Inhomogeneous vector term G(X)
!
!     XSAV = Previous value of X
!
!     C = Normalization factor for the particular solution
!
!           0   ( if  NEQIVP = 0 )
!     IVP =
!           Number of differential equations integrated due to
!           the original boundary value problem   ( if  NEQIVP  >  0 )
!
!     NOFST - For problems with auxiliary initial value equations,
!             NOFST communicates to the routine FMAT how to access
!             the dependent variables corresponding to this initial
!             value problem.  For example, during any call to FMAT,
!             the first dependent variable for the initial value
!             problem is in position  Y(NOFST + 1).
!             See example in SAND77-1328.
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    ML8SZ, MLIVP
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910701  Corrected ROUTINES CALLED section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920618  Minor restructuring of code.  (RWC, WRB)
!***END PROLOGUE  BVDER
  DIMENSION Y(*),YP(*),G(*)
!
! **********************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
!
! **********************************************************************
!     The COMMON block below is used to communicate with the user
!     supplied subroutine FMAT.  The user should not alter this
!     COMMON block.
!
  COMMON /MLIVP/ NOFST
! **********************************************************************
!
!***FIRST EXECUTABLE STATEMENT  BVDER
  if (IVP  >  0) call UIVP(X,Y(IVP+1),YP(IVP+1))
  NOFST = IVP
  NA = 1
  DO 10 K=1,NFC
     call FMAT(X,Y(NA),YP(NA))
     NOFST = NOFST - NCOMP
     NA = NA + NCOMP
   10 CONTINUE
!
  if (INHOMO  /=  1) RETURN
  call FMAT(X,Y(NA),YP(NA))
!
  if (IGOFX  ==  0) RETURN
  if (X  /=  XSAV) THEN
     if (IVP  ==  0) call GVEC(X,G)
     if (IVP  >  0) call UVEC(X,Y(IVP+1),G)
     XSAV = X
  end if
!
!     If the user has chosen not to normalize the particular
!     solution, then C is defined in BVPOR to be 1.0
!
!     The following loop is just
!     call SAXPY (NCOMP, 1.0E0/C, G, 1, YP(NA), 1)
!
  DO 20 J=1,NCOMP
     L = NA + J - 1
     YP(L) = YP(L) + G(J)/C
   20 CONTINUE
  return
end
