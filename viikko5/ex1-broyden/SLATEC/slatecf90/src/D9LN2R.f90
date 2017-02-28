DOUBLE PRECISION FUNCTION D9LN2R (X)
!
!! D9LN2R evaluates LOG(1+X) from second order relative accuracy ...
!  so that LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4B
!***TYPE      DOUBLE PRECISION (R9LN2R-S, D9LN2R-D, C9LN2R-C)
!***KEYWORDS  ELEMENTARY FUNCTIONS, FNLIB, LOGARITHM, SECOND ORDER
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate  LOG(1+X)  from 2-nd order with relative error accuracy so
! that    LOG(1+X) = X - X**2/2 + X**3*D9LN2R(X)
!
! Series for LN21       on the interval -6.25000E-01 to  0.
!                                        with weighted error   1.82E-32
!                                         log weighted error  31.74
!                               significant figures required  31.00
!                                    decimal places required  32.59
!
! Series for LN22       on the interval  0.          to  8.12500E-01
!                                        with weighted error   6.10E-32
!                                         log weighted error  31.21
!                               significant figures required  30.32
!                                    decimal places required  32.00
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DCSEVL, INITDS, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  D9LN2R
  DOUBLE PRECISION X, XBIG, TXBIG, XMAX, TXMAX, XMIN, LN21CS(50), &
    LN22CS(37), DCSEVL, D1MACH
  LOGICAL FIRST
  SAVE LN21CS, LN22CS, NTLN21, NTLN22, XMIN, XBIG, XMAX, FIRST
  DATA LN21CS(  1) / +.18111962513478809875894953043071D+0     /
  DATA LN21CS(  2) / -.15627123192872462669625155541078D+0     /
  DATA LN21CS(  3) / +.28676305361557275209540627102051D-1     /
  DATA LN21CS(  4) / -.55586996559481398781157725126781D-2     /
  DATA LN21CS(  5) / +.11178976652299837657335666279727D-2     /
  DATA LN21CS(  6) / -.23080508982327947182299279585705D-3     /
  DATA LN21CS(  7) / +.48598853341100175874681558068750D-4     /
  DATA LN21CS(  8) / -.10390127388903210765514242633338D-4     /
  DATA LN21CS(  9) / +.22484563707390128494621804946408D-5     /
  DATA LN21CS( 10) / -.49140592739266484875327802597091D-6     /
  DATA LN21CS( 11) / +.10828256507077483336620152971597D-6     /
  DATA LN21CS( 12) / -.24025872763420701435976675416719D-7     /
  DATA LN21CS( 13) / +.53624600472708133762984443250163D-8     /
  DATA LN21CS( 14) / -.12029951362138772264671646424377D-8     /
  DATA LN21CS( 15) / +.27107889277591860785622551632266D-9     /
  DATA LN21CS( 16) / -.61323562618319010068796728430690D-10    /
  DATA LN21CS( 17) / +.13920858369159469857436908543978D-10    /
  DATA LN21CS( 18) / -.31699300330223494015283057260883D-11    /
  DATA LN21CS( 19) / +.72383754044307505335214326197011D-12    /
  DATA LN21CS( 20) / -.16570017184764411391498805506268D-12    /
  DATA LN21CS( 21) / +.38018428663117424257364422631876D-13    /
  DATA LN21CS( 22) / -.87411189296972700259724429899137D-14    /
  DATA LN21CS( 23) / +.20135619845055748302118751028154D-14    /
  DATA LN21CS( 24) / -.46464456409033907031102008154477D-15    /
  DATA LN21CS( 25) / +.10739282147018339453453338554925D-15    /
  DATA LN21CS( 26) / -.24858534619937794755534021833960D-16    /
  DATA LN21CS( 27) / +.57620197950800189813888142628181D-17    /
  DATA LN21CS( 28) / -.13373063769804394701402199958050D-17    /
  DATA LN21CS( 29) / +.31074653227331824966533807166805D-18    /
  DATA LN21CS( 30) / -.72288104083040539906901957917627D-19    /
  DATA LN21CS( 31) / +.16833783788037385103313258186888D-19    /
  DATA LN21CS( 32) / -.39239463312069958052519372739925D-20    /
  DATA LN21CS( 33) / +.91551468387536789746385528640853D-21    /
  DATA LN21CS( 34) / -.21378895321320159520982095801002D-21    /
  DATA LN21CS( 35) / +.49964507479047864699828564568746D-22    /
  DATA LN21CS( 36) / -.11686240636080170135360806147413D-22    /
  DATA LN21CS( 37) / +.27353123470391863775628686786559D-23    /
  DATA LN21CS( 38) / -.64068025084792111965050345881599D-24    /
  DATA LN21CS( 39) / +.15016293204334124162949071940266D-24    /
  DATA LN21CS( 40) / -.35217372410398479759497145002666D-25    /
  DATA LN21CS( 41) / +.82643901014814767012482733397333D-26    /
  DATA LN21CS( 42) / -.19404930275943401918036617898666D-26    /
  DATA LN21CS( 43) / +.45587880018841283562451588437333D-27    /
  DATA LN21CS( 44) / -.10715492087545202154378625023999D-27    /
  DATA LN21CS( 45) / +.25199408007927592978096674133333D-28    /
  DATA LN21CS( 46) / -.59289088400120969341750476800000D-29    /
  DATA LN21CS( 47) / +.13955864061057513058237153279999D-29    /
  DATA LN21CS( 48) / -.32864578813478583431436697599999D-30    /
  DATA LN21CS( 49) / +.77424967950478166247254698666666D-31    /
  DATA LN21CS( 50) / -.18247735667260887638125226666666D-31    /
  DATA LN22CS(  1) / -.2224253253502046082986015223552D+0      /
  DATA LN22CS(  2) / -.6104710010807862398680104755764D-1      /
  DATA LN22CS(  3) / +.7427235009750394590519629755729D-2      /
  DATA LN22CS(  4) / -.9335018261636970565612779606397D-3      /
  DATA LN22CS(  5) / +.1200499076872601283350731287359D-3      /
  DATA LN22CS(  6) / -.1570472295282004112823352608243D-4      /
  DATA LN22CS(  7) / +.2081874781051271096050783592759D-5      /
  DATA LN22CS(  8) / -.2789195577646713654057213051375D-6      /
  DATA LN22CS(  9) / +.3769355823760132058422895135447D-7      /
  DATA LN22CS( 10) / -.5130902896527711258240589938003D-8      /
  DATA LN22CS( 11) / +.7027141178150694738206218215392D-9      /
  DATA LN22CS( 12) / -.9674859550134342389243972005137D-10     /
  DATA LN22CS( 13) / +.1338104645924887306588496449748D-10     /
  DATA LN22CS( 14) / -.1858102603534063981628453846591D-11     /
  DATA LN22CS( 15) / +.2589294422527919749308600123070D-12     /
  DATA LN22CS( 16) / -.3619568316141588674466025382172D-13     /
  DATA LN22CS( 17) / +.5074037398016623088006858917396D-14     /
  DATA LN22CS( 18) / -.7131012977031127302700938748927D-15     /
  DATA LN22CS( 19) / +.1004490328554567481853386784126D-15     /
  DATA LN22CS( 20) / -.1417906532184025791904405075285D-16     /
  DATA LN22CS( 21) / +.2005297034743326117891086396074D-17     /
  DATA LN22CS( 22) / -.2840996662339803305365396717567D-18     /
  DATA LN22CS( 23) / +.4031469883969079899599878662826D-19     /
  DATA LN22CS( 24) / -.5729325241832207320455498956799D-20     /
  DATA LN22CS( 25) / +.8153488253890010675848928733866D-21     /
  DATA LN22CS( 26) / -.1161825588549721787606027468799D-21     /
  DATA LN22CS( 27) / +.1657516611662538343659339775999D-22     /
  DATA LN22CS( 28) / -.2367336704710805190114017280000D-23     /
  DATA LN22CS( 29) / +.3384670367975521386076569599999D-24     /
  DATA LN22CS( 30) / -.4843940829215718204296396799999D-25     /
  DATA LN22CS( 31) / +.6938759162514273718676138666666D-26     /
  DATA LN22CS( 32) / -.9948142607031436571923797333333D-27     /
  DATA LN22CS( 33) / +.1427440611211698610634752000000D-27     /
  DATA LN22CS( 34) / -.2049794721898234911566506666666D-28     /
  DATA LN22CS( 35) / +.2945648756401362222885546666666D-29     /
  DATA LN22CS( 36) / -.4235973185184957027669333333333D-30     /
  DATA LN22CS( 37) / +.6095532614003832040106666666666D-31     /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  D9LN2R
  if (FIRST) THEN
     EPS = D1MACH(3)
     NTLN21 = INITDS (LN21CS, 50, 0.1*EPS)
     NTLN22 = INITDS (LN22CS, 37, 0.1*EPS)
!
     XMIN = -1.0D0 + SQRT(D1MACH(4))
     SQEPS = SQRT (EPS)
     TXMAX = 8.0/SQEPS
     XMAX = TXMAX - (EPS*TXMAX**2 - 2.D0*LOG(TXMAX)) &
       / (2.D0*EPS*TXMAX)
     TXBIG = 6.0/SQRT(SQEPS)
     XBIG = TXBIG - (SQEPS*TXBIG**2 - 2.D0*LOG(TXBIG)) &
       / (2.D0*SQEPS*TXBIG)
  end if
  FIRST = .FALSE.
!
  if (X < (-.625D0) .OR. X > 0.8125D0) go to 20
!
  if (X < 0.0D0) D9LN2R = 0.375D0 + DCSEVL (16.D0*X/5.D0+1.D0, &
    LN21CS, NTLN21)
  if (X >= 0.0D0) D9LN2R = 0.375D0 + DCSEVL (32.D0*X/13.D0-1.D0, &
    LN22CS, NTLN22)
  return
!
 20   if (X  <  XMIN) call XERMSG ('SLATEC', 'D9LN2R', &
     'ANSWER LT HALF PRECISION BECAUSE X IS TOO NEAR -1', 1, 1)
  if (X  >  XMAX) call XERMSG ('SLATEC', 'D9LN2R', &
     'NO PRECISION IN ANSWER BECAUSE X IS TOO BIG', 3, 2)
  if (X  >  XBIG) call XERMSG ('SLATEC', 'D9LN2R', &
     'ANSWER LT HALF PRECISION BECAUSE X IS TOO BIG', 2, 1)
!
  D9LN2R = (LOG(1.D0+X) - X*(1.D0 - 0.5D0*X)) / X**3
  return
!
end