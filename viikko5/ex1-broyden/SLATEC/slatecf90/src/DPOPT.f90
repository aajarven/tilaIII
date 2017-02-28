subroutine DPOPT (PRGOPT, MRELAS, NVARS, INFO, CSC, IBASIS, ROPT, &
     INTOPT, LOPT)
!
!! DPOPT is subsidiary to DSPLP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SPOPT-S, DPOPT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     THE EDITING REQUIRED TO CONVERT THIS SUBROUTINE FROM SINGLE TO
!     DOUBLE PRECISION INVOLVES THE FOLLOWING CHARACTER STRING CHANGES.
!
!     USE AN EDITING COMMAND (CHANGE) /STRING-1/(TO)STRING-2/.
!     /REAL (12 BLANKS)/DOUBLE PRECISION/,/R1MACH/D1MACH/,/E0/D0/
!
!     REVISED 821122-1045
!     REVISED YYMMDD-HHMM
!
!     THIS SUBROUTINE PROCESSES THE OPTION VECTOR, PRGOPT(*),
!     AND VALIDATES ANY MODIFIED DATA.
!
!***SEE ALSO  DSPLP
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   811215  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890605  Removed unreferenced labels.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900328  Added TYPE section.  (WRB)
!   900510  Fixed an error message.  (RWC)
!***END PROLOGUE  DPOPT
  DOUBLE PRECISION ABIG,ASMALL,COSTSC,CSC(*),EPS,ONE,PRGOPT(*), &
   ROPT(07),TOLLS,TUNE,ZERO,D1MACH,TOLABS
  INTEGER IBASIS(*),INTOPT(08)
  LOGICAL CONTIN,USRBAS,SIZEUP,SAVEDT,COLSCP,CSTSCP,MINPRB, &
   STPEDG,LOPT(8)
!
!***FIRST EXECUTABLE STATEMENT  DPOPT
  IOPT=1
  ZERO=0.D0
  ONE=1.D0
  go to 30001
20002 CONTINUE
  go to 30002
!
20003 LOPT(1)=CONTIN
  LOPT(2)=USRBAS
  LOPT(3)=SIZEUP
  LOPT(4)=SAVEDT
  LOPT(5)=COLSCP
  LOPT(6)=CSTSCP
  LOPT(7)=MINPRB
  LOPT(8)=STPEDG
!
  INTOPT(1)=IDG
  INTOPT(2)=IPAGEF
  INTOPT(3)=ISAVE
  INTOPT(4)=MXITLP
  INTOPT(5)=KPRINT
  INTOPT(6)=ITBRC
  INTOPT(7)=NPP
  INTOPT(8)=LPRG
!
  ROPT(1)=EPS
  ROPT(2)=ASMALL
  ROPT(3)=ABIG
  ROPT(4)=COSTSC
  ROPT(5)=TOLLS
  ROPT(6)=TUNE
  ROPT(7)=TOLABS
  return
!
!
!     PROCEDURE (INITIALIZE PARAMETERS AND PROCESS USER OPTIONS)
30001 CONTIN = .FALSE.
  USRBAS = .FALSE.
  SIZEUP = .FALSE.
  SAVEDT = .FALSE.
  COLSCP = .FALSE.
  CSTSCP = .FALSE.
  MINPRB = .TRUE.
  STPEDG = .TRUE.
!
!     GET THE MACHINE REL. FLOATING POINT ACCURACY VALUE FROM THE
!     LIBRARY SUBPROGRAM, D1MACH( ).
  EPS=D1MACH(4)
  TOLLS=D1MACH(4)
  TUNE=ONE
  TOLABS=ZERO
!
!     DEFINE NOMINAL FILE NUMBERS FOR MATRIX PAGES AND DATA SAVING.
  IPAGEF=1
  ISAVE=2
  ITBRC=10
  MXITLP=3*(NVARS+MRELAS)
  KPRINT=0
  IDG=-4
  NPP=NVARS
  LPRG=0
!
  LAST = 1
  IADBIG=10000
  ICTMAX=1000
  ICTOPT= 0
20004 NEXT=PRGOPT(LAST)
  if (.NOT.(NEXT <= 0 .OR. NEXT > IADBIG)) go to 20006
!
!     THE CHECKS FOR SMALL OR LARGE VALUES OF NEXT ARE TO PREVENT
!     WORKING WITH UNDEFINED DATA.
  NERR=14
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, THE USER OPTION ARRAY HAS UNDEFINED DATA.', NERR, &
     IOPT)
  INFO=-NERR
  return
20006 if (.NOT.(NEXT == 1)) go to 10001
  go to 20005
10001 if (.NOT.(ICTOPT > ICTMAX)) go to 10002
  NERR=15
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, OPTION ARRAY PROCESSING IS CYCLING.', NERR, IOPT)
  INFO=-NERR
  return
10002 CONTINUE
  KEY = PRGOPT(LAST+1)
!
!     if KEY = 50, THIS IS TO BE A MAXIMIZATION PROBLEM
!     INSTEAD OF A MINIMIZATION PROBLEM.
  if (.NOT.(KEY == 50)) go to 20010
  MINPRB = PRGOPT(LAST+2) == ZERO
  LDS=3
  go to 20009
20010 CONTINUE
!
!     if KEY = 51, THE LEVEL OF OUTPUT IS BEING MODIFIED.
!     KPRINT = 0, NO OUTPUT
!            = 1, SUMMARY OUTPUT
!            = 2, LOTS OF OUTPUT
!            = 3, EVEN MORE OUTPUT
  if (.NOT.(KEY == 51)) go to 20013
  KPRINT=PRGOPT(LAST+2)
  LDS=3
  go to 20009
20013 CONTINUE
!
!     if KEY = 52, REDEFINE THE FORMAT AND PRECISION USED
!     IN THE OUTPUT.
  if (.NOT.(KEY == 52)) go to 20016
  if (PRGOPT(LAST+2) /= ZERO) IDG=PRGOPT(LAST+3)
  LDS=4
  go to 20009
20016 CONTINUE
!
!     if KEY = 53, THE ALLOTTED SPACE FOR THE SPARSE MATRIX
!     STORAGE AND/OR SPARSE EQUATION SOLVING HAS BEEN CHANGED.
!    (PROCESSED IN DSPLP(). THIS IS TO COMPUTE THE LENGTH OF PRGOPT(*).)
  if (.NOT.(KEY == 53)) go to 20019
  LDS=5
  go to 20009
20019 CONTINUE
!
!     if KEY = 54, REDEFINE THE FILE NUMBER WHERE THE PAGES
!     FOR THE SPARSE MATRIX ARE STORED.
  if (.NOT.(KEY == 54)) go to 20022
  if ( PRGOPT(LAST+2) /= ZERO) IPAGEF = PRGOPT(LAST+3)
  LDS=4
  go to 20009
20022 CONTINUE
!
!     if KEY = 55,  A CONTINUATION FOR A PROBLEM MAY BE REQUESTED.
  if (.NOT.(KEY  ==  55)) go to 20025
  CONTIN = PRGOPT(LAST+2) /= ZERO
  LDS=3
  go to 20009
20025 CONTINUE
!
!     if KEY = 56, REDEFINE THE FILE NUMBER WHERE THE SAVED DATA
!     WILL BE STORED.
  if (.NOT.(KEY == 56)) go to 20028
  if ( PRGOPT(LAST+2) /= ZERO) ISAVE = PRGOPT(LAST+3)
  LDS=4
  go to 20009
20028 CONTINUE
!
!     if KEY = 57, SAVE DATA (ON EXTERNAL FILE)  AT MXITLP ITERATIONS OR
!     THE OPTIMUM, WHICHEVER COMES FIRST.
  if (.NOT.(KEY == 57)) go to 20031
  SAVEDT=PRGOPT(LAST+2) /= ZERO
  LDS=3
  go to 20009
20031 CONTINUE
!
!     if KEY = 58,  SEE IF PROBLEM IS TO RUN ONLY A GIVEN
!     NUMBER OF ITERATIONS.
  if (.NOT.(KEY == 58)) go to 20034
  if (PRGOPT(LAST+2) /= ZERO) MXITLP = PRGOPT(LAST+3)
  LDS=4
  go to 20009
20034 CONTINUE
!
!     if KEY = 59,  SEE IF USER PROVIDES THE BASIS INDICES.
  if (.NOT.(KEY  ==  59)) go to 20037
  USRBAS = PRGOPT(LAST+2)  /=  ZERO
  if (.NOT.(USRBAS)) go to 20040
  I=1
  N20043=MRELAS
  go to 20044
20043 I=I+1
20044 if ((N20043-I) < 0) go to 20045
  IBASIS(I) = PRGOPT(LAST+2+I)
  go to 20043
20045 CONTINUE
20040 CONTINUE
  LDS=MRELAS+3
  go to 20009
20037 CONTINUE
!
!     if KEY = 60,  SEE IF USER HAS PROVIDED SCALING OF COLUMNS.
  if (.NOT.(KEY  ==  60)) go to 20047
  COLSCP = PRGOPT(LAST+2) /= ZERO
  if (.NOT.(COLSCP)) go to 20050
  J=1
  N20053=NVARS
  go to 20054
20053 J=J+1
20054 if ((N20053-J) < 0) go to 20055
  CSC(J)=ABS(PRGOPT(LAST+2+J))
  go to 20053
20055 CONTINUE
20050 CONTINUE
  LDS=NVARS+3
  go to 20009
20047 CONTINUE
!
!     if KEY = 61,  SEE IF USER HAS PROVIDED SCALING OF COSTS.
  if (.NOT.(KEY  ==  61)) go to 20057
  CSTSCP = PRGOPT(LAST+2) /= ZERO
  if (CSTSCP) COSTSC = PRGOPT(LAST+3)
  LDS=4
  go to 20009
20057 CONTINUE
!
!     if KEY = 62,  SEE IF SIZE PARAMETERS ARE PROVIDED WITH THE DATA.
!     THESE WILL BE CHECKED AGAINST THE MATRIX ELEMENT SIZES LATER.
  if (.NOT.(KEY  ==  62)) go to 20060
  SIZEUP = PRGOPT(LAST+2) /= ZERO
  if (.NOT.(SIZEUP)) go to 20063
  ASMALL = PRGOPT(LAST+3)
  ABIG = PRGOPT(LAST+4)
20063 CONTINUE
  LDS=5
  go to 20009
20060 CONTINUE
!
!     if KEY = 63, SEE IF TOLERANCE FOR LINEAR SYSTEM RESIDUAL ERROR IS
!     PROVIDED.
  if (.NOT.(KEY  ==  63)) go to 20066
  if (PRGOPT(LAST+2) /= ZERO) TOLLS = MAX(EPS,PRGOPT(LAST+3))
  LDS=4
  go to 20009
20066 CONTINUE
!
!     if KEY = 64,  SEE IF MINIMUM REDUCED COST OR STEEPEST EDGE
!     DESCENT IS TO BE USED FOR SELECTING VARIABLES TO ENTER BASIS.
  if (.NOT.(KEY == 64)) go to 20069
  STPEDG = PRGOPT(LAST+2) == ZERO
  LDS=3
  go to 20009
20069 CONTINUE
!
!     if KEY = 65, SET THE NUMBER OF ITERATIONS BETWEEN RECALCULATING
!     THE ERROR IN THE PRIMAL SOLUTION.
  if (.NOT.(KEY == 65)) go to 20072
  if (PRGOPT(LAST+2) /= ZERO) ITBRC=MAX(ONE,PRGOPT(LAST+3))
  LDS=4
  go to 20009
20072 CONTINUE
!
!     if KEY = 66, SET THE NUMBER OF NEGATIVE REDUCED COSTS TO BE FOUND
!     IN THE PARTIAL PRICING STRATEGY.
  if (.NOT.(KEY == 66)) go to 20075
  if (.NOT.(PRGOPT(LAST+2) /= ZERO)) go to 20078
  NPP=MAX(PRGOPT(LAST+3),ONE)
  NPP=MIN(NPP,NVARS)
20078 CONTINUE
  LDS=4
  go to 20009
20075 CONTINUE
!     if KEY = 67, CHANGE THE TUNING PARAMETER TO APPLY TO THE ERROR
!     ESTIMATES FOR THE PRIMAL AND DUAL SYSTEMS.
  if (.NOT.(KEY == 67)) go to 20081
  if (.NOT.(PRGOPT(LAST+2) /= ZERO)) go to 20084
  TUNE=ABS(PRGOPT(LAST+3))
20084 CONTINUE
  LDS=4
  go to 20009
20081 CONTINUE
  if (.NOT.(KEY == 68)) go to 20087
  LDS=6
  go to 20009
20087 CONTINUE
!
!     RESET THE ABSOLUTE TOLERANCE TO BE USED ON THE FEASIBILITY
!     DECISION PROVIDED THE RELATIVE ERROR TEST FAILED.
  if (.NOT.(KEY == 69)) go to 20090
  if ( PRGOPT(LAST+2) /= ZERO)TOLABS=PRGOPT(LAST+3)
  LDS=4
  go to 20009
20090 CONTINUE
  CONTINUE
!
20009 ICTOPT = ICTOPT+1
  LAST = NEXT
  LPRG=LPRG+LDS
  go to 20004
20005 CONTINUE
  go to 20002
!
!     PROCEDURE (VALIDATE OPTIONALLY MODIFIED DATA)
!
!     if USER HAS DEFINED THE BASIS, CHECK FOR VALIDITY OF INDICES.
30002 if (.NOT.(USRBAS)) go to 20093
  I=1
  N20096=MRELAS
  go to 20097
20096 I=I+1
20097 if ((N20096-I) < 0) go to 20098
  ITEST=IBASIS(I)
  if (.NOT.(ITEST <= 0 .OR.ITEST > (NVARS+MRELAS))) go to 20100
  NERR=16
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, AN INDEX OF USER-SUPPLIED BASIS IS OUT OF RANGE.', &
     NERR, IOPT)
  INFO=-NERR
  return
20100 CONTINUE
  go to 20096
20098 CONTINUE
20093 CONTINUE
!
!     if USER HAS PROVIDED SIZE PARAMETERS, MAKE SURE THEY ARE ORDERED
!     AND POSITIVE.
  if (.NOT.(SIZEUP)) go to 20103
  if (.NOT.(ASMALL <= ZERO .OR. ABIG < ASMALL)) go to 20106
  NERR=17
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, SIZE PARAMETERS FOR MATRIX MUST BE SMALLEST AND ' // &
     'LARGEST MAGNITUDES OF NONZERO ENTRIES.', NERR, IOPT)
  INFO=-NERR
  return
20106 CONTINUE
20103 CONTINUE
!
!     THE NUMBER OF ITERATIONS OF REV. SIMPLEX STEPS MUST BE POSITIVE.
  if (.NOT.(MXITLP <= 0)) go to 20109
  NERR=18
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, THE NUMBER OF REVISED SIMPLEX STEPS BETWEEN ' // &
     'CHECK-POINTS MUST BE POSITIVE.', NERR, IOPT)
  INFO=-NERR
  return
20109 CONTINUE
!
!  CHECK THAT SAVE AND PAGE FILE NUMBERS ARE DEFINED AND NOT EQUAL.
!
  if (.NOT.(ISAVE <= 0.OR.IPAGEF <= 0.OR.(ISAVE == IPAGEF))) go to 20112
  NERR=19
  call XERMSG ('SLATEC', 'DPOPT', &
     'IN DSPLP, FILE NUMBERS FOR SAVED DATA AND MATRIX PAGES ' // &
     'MUST BE POSITIVE AND NOT EQUAL.', NERR, IOPT)
  INFO=-NERR
  return
20112 CONTINUE
  CONTINUE
  go to 20003
end