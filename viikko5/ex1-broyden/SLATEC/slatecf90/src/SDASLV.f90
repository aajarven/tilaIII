subroutine SDASLV (NEQ, DELTA, WM, IWM)
!
!! SDASLV is the linear system solver for SDASSL.
!
!***LIBRARY   SLATEC (DASSL)
!***TYPE      SINGLE PRECISION (SDASLV-S, DDASLV-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
!     SYSTEM ARISING IN THE NEWTON ITERATION.
!     MATRICES AND REAL TEMPORARY STORAGE AND
!     REAL INFORMATION ARE STORED IN THE ARRAY WM.
!     INTEGER MATRIX INFORMATION IS STORED IN
!     THE ARRAY IWM.
!     FOR A DENSE MATRIX, THE LINPACK ROUTINE
!     SGESL IS CALLED.
!     FOR A BANDED MATRIX,THE LINPACK ROUTINE
!     SGBSL IS CALLED.
!-----------------------------------------------------------------------
!***ROUTINES CALLED  SGBSL, SGESL
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  SDASLV
!
  INTEGER  NEQ, IWM(*)
  REAL  DELTA(*), WM(*)
!
  EXTERNAL  SGBSL, SGESL
!
  INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
  PARAMETER (NPD=1)
  PARAMETER (LML=1)
  PARAMETER (LMU=2)
  PARAMETER (LMTYPE=4)
  PARAMETER (LIPVT=21)
!
!***FIRST EXECUTABLE STATEMENT  SDASLV
  MTYPE=IWM(LMTYPE)
  go to(100,100,300,400,400),MTYPE
!
!     DENSE MATRIX
100   call SGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
  return
!
!     DUMMY SECTION FOR MTYPE=3
300   CONTINUE
  return
!
!     BANDED MATRIX
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
  call SGBSL(WM(NPD),MEBAND,NEQ,IWM(LML), &
    IWM(LMU),IWM(LIPVT),DELTA,0)
  return
!------END OF SUBROUTINE SDASLV------
end
