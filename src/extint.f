!************************************************************************
! Initial evaluation of the two components of the function
! Here qp is an array going from 0 to 1. Corresponding to r running from
! infinity to rmax=rb.
      subroutine initext(next,qp,phibye,phiei,rb,xlambda)
      integer next
      real qp(next),phibye(next),phiei(next)
      real rb
! statement function:
      eeih(x)=alog(1.+1./x)-.56/(1.+4.1*x+0.9*x*x)

! precalculate values for an integration over qp range.

      if(rb.le.0.)then 
         write(*,*)'Error in external integration initialization.'
     $        ,' rb=',rb
         stop
      endif
      rbol=rb/xlambda
! Values multiplied by exp functions.
      if(rbol.lt.50.)then
         E1rb=exp(rbol)*Eone(rbol)
      else
         E1rb=eeih(rbol)
      endif
      Eirb=EXPEi(rbol)
      do i=1,next
! Get q in proper units.
         q=qp(i)/rb
         rol=1./(q*xlambda)
         if(rol-rbol.gt.50) then
            phibye(i)=0.
         else
            phibye(i)=exp(rbol-rol)*rb*q
         endif
         if(rol.gt.50.)then
!            ee1=eeih(rol) 
            ee1=1/rol
         else
            ee1=exp(rol)*Eone(rol)
         endif
            phiei(i)=-( 1./(2.*rol))*(
     $           exp(rbol-rol)*(E1rb+Eirb)
     $           -ee1
     $           -(EXPEi(rol)) )
      enddo
      end
!********************************************************************
! From here on is the NETLIB Exponential integral package.
      SUBROUTINE CALCEI(ARG,RESULT,INT)
!----------------------------------------------------------------------
!
! This Fortran 77 packet computes the exponential integrals Ei(x),
!  E1(x), and  exp(-x)*Ei(x)  for real arguments  x  where
!
!           integral (from t=-infinity to t=x) (exp(t)/t),  x > 0,
!  Ei(x) =
!          -integral (from t=-x to t=infinity) (exp(t)/t),  x < 0,
!
!  and where the first integral is a principal value integral.
!  The packet contains three function type subprograms: EI, EONE,
!  and EXPEI;  and one subroutine type subprogram: CALCEI.  The
!  calling statements for the primary entries are
!
!                 Y = EI(X),            where  X .NE. 0,
!
!                 Y = EONE(X),          where  X .GT. 0,
!  and
!                 Y = EXPEI(X),         where  X .NE. 0,
!
!  and where the entry points correspond to the functions Ei(x),
!  E1(x), and exp(-x)*Ei(x), respectively.  The routine CALCEI
!  is intended for internal packet use only, all computations within
!  the packet being concentrated in this routine.  The function
!  subprograms invoke CALCEI with the Fortran statement
!         CALL CALCEI(ARG,RESULT,INT)
!  where the parameter usage is as follows
!
!     Function                  Parameters for CALCEI
!       Call                 ARG             RESULT         INT
!
!      EI(X)              X .NE. 0          Ei(X)            1
!      EONE(X)            X .GT. 0         -Ei(-X)           2
!      EXPEI(X)           X .NE. 0          exp(-X)*Ei(X)    3
!
!  The main computation involves evaluation of rational Chebyshev
!  approximations published in Math. Comp. 22, 641-649 (1968), and
!  Math. Comp. 23, 289-303 (1969) by Cody and Thacher.  This
!  transportable program is patterned after the machine-dependent
!  FUNPACK packet  NATSEI,  but cannot match that version for
!  efficiency or accuracy.  This version uses rational functions
!  that theoretically approximate the exponential integrals to
!  at least 18 significant decimal digits.  The accuracy achieved
!  depends on the arithmetic system, the compiler, the intrinsic
!  functions, and proper selection of the machine-dependent
!  constants.
!
!
!*******************************************************************
!*******************************************************************
!
! Explanation of machine-dependent constants
!
!   beta = radix for the floating-point system.
!   minexp = smallest representable power of beta.
!   maxexp = smallest power of beta that overflows.
!   XBIG = largest argument acceptable to EONE; solution to
!          equation:
!                     exp(-x)/x * (1 + 1/x) = beta ** minexp.
!   XINF = largest positive machine number; approximately
!                     beta ** maxexp
!   XMAX = largest argument acceptable to EI; solution to
!          equation:  exp(x)/x * (1 + 1/x) = beta ** maxexp.
!
!     Approximate values for some important machines are:
!
!                           beta      minexp      maxexp
!
!  CRAY-1        (S.P.)       2       -8193        8191
!  Cyber 180/185 
!    under NOS   (S.P.)       2        -975        1070
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)       2        -126         128
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)       2       -1022        1024
!  IBM 3033      (D.P.)      16         -65          63
!  VAX D-Format  (D.P.)       2        -128         127
!  VAX G-Format  (D.P.)       2       -1024        1023
!
!                           XBIG       XINF       XMAX
!
!  CRAY-1        (S.P.)    5670.31  5.45E+2465   5686.21
!  Cyber 180/185 
!    under NOS   (S.P.)     669.31  1.26E+322     748.28
!  IEEE (IBM/XT,
!    SUN, etc.)  (S.P.)      82.93  3.40E+38       93.24
!  IEEE (IBM/XT,
!    SUN, etc.)  (D.P.)     701.84  1.79D+308     716.35
!  IBM 3033      (D.P.)     175.05  7.23D+75      179.85
!  VAX D-Format  (D.P.)      84.30  1.70D+38       92.54
!  VAX G-Format  (D.P.)     703.22  8.98D+307     715.66
!
!*******************************************************************
!*******************************************************************
!
! Error returns
!
!  The following table shows the types of error that may be
!  encountered in this routine and the function value supplied
!  in each case.
!
!       Error       Argument         Function values for
!                    Range         EI      EXPEI     EONE
!
!     UNDERFLOW  (-)X .GT. XBIG     0        -         0
!     OVERFLOW      X .GE. XMAX    XINF      -         -
!     ILLEGAL X       X = 0       -XINF    -XINF     XINF
!     ILLEGAL X      X .LT. 0       -        -     USE ABS(X)
!
! Intrinsic functions required are:
!
!     ABS, SQRT, EXP
!
!
!  Author: W. J. Cody
!          Mathematics abd Computer Science Division
!          Argonne National Laboratory
!          Argonne, IL 60439
!
!  Latest modification: September 9, 1988
!
!----------------------------------------------------------------------
      INTEGER I,INT
!D    DOUBLE PRECISION 
      REAL
     1       A,ARG,B,C,D,EXP40,E,EI,F,FOUR,FOURTY,FRAC,HALF,ONE,P,
     2       PLG,PX,P037,P1,P2,Q,QLG,QX,Q1,Q2,R,RESULT,S,SIX,SUMP,
     3       SUMQ,T,THREE,TWELVE,TWO,TWO4,W,X,XBIG,XINF,XMAX,XMX0,
     4       X0,X01,X02,X11,Y,YSQ,ZERO
      DIMENSION  A(7),B(6),C(9),D(9),E(10),F(10),P(10),Q(10),R(10),
     1   S(9),P1(10),Q1(9),P2(10),Q2(9),PLG(4),QLG(4),PX(10),QX(10)
!----------------------------------------------------------------------
!  Mathematical constants
!   EXP40 = exp(40)
!   X0 = zero of Ei
!   X01/X11 + X02 = zero of Ei to extra precision
!----------------------------------------------------------------------
      DATA ZERO,P037,HALF,ONE,TWO/0.0E0,0.037E0,0.5E0,1.0E0,2.0E0/,
     1     THREE,FOUR,SIX,TWELVE,TWO4/3.0E0,4.0E0,6.0E0,12.E0,24.0E0/,
     2     FOURTY,EXP40/40.0E0,2.3538526683701998541E17/,
     3     X01,X11,X02/381.5E0,1024.0E0,-5.1182968633365538008E-5/,
     4     X0/3.7250741078136663466E-1/
!D    DATA ZERO,P037,HALF,ONE,TWO/0.0D0,0.037D0,0.5D0,1.0D0,2.0D0/,
!D   1     THREE,FOUR,SIX,TWELVE,TWO4/3.0D0,4.0D0,6.0D0,12.D0,24.0D0/,
!D   2     FOURTY,EXP40/40.0D0,2.3538526683701998541D17/,
!D   3     X01,X11,X02/381.5D0,1024.0D0,-5.1182968633365538008D-5/,
!D   4     X0/3.7250741078136663466D-1/
!----------------------------------------------------------------------
! Machine-dependent constants
!----------------------------------------------------------------------
      DATA XINF/3.40E+38/,XMAX/93.246E0/,XBIG/82.93E0/
!D    DATA XINF/1.79D+308/,XMAX/716.351D0/,XBIG/701.84D0/
!----------------------------------------------------------------------
! Coefficients  for -1.0 <= X < 0.0
!----------------------------------------------------------------------
      DATA A/1.1669552669734461083368E2, 2.1500672908092918123209E3,
     1       1.5924175980637303639884E4, 8.9904972007457256553251E4,
     2       1.5026059476436982420737E5,-1.4815102102575750838086E5,
     3       5.0196785185439843791020E0/
      DATA B/4.0205465640027706061433E1, 7.5043163907103936624165E2,
     1       8.1258035174768735759855E3, 5.2440529172056355429883E4,
     2       1.8434070063353677359298E5, 2.5666493484897117319268E5/
!D    DATA A/1.1669552669734461083368D2, 2.1500672908092918123209D3,
!D   1       1.5924175980637303639884D4, 8.9904972007457256553251D4,
!D   2       1.5026059476436982420737D5,-1.4815102102575750838086D5,
!D   3       5.0196785185439843791020D0/
!D    DATA B/4.0205465640027706061433D1, 7.5043163907103936624165D2,
!D   1       8.1258035174768735759855D3, 5.2440529172056355429883D4,
!D   2       1.8434070063353677359298D5, 2.5666493484897117319268D5/
!----------------------------------------------------------------------
! Coefficients for -4.0 <= X < -1.0
!----------------------------------------------------------------------
      DATA C/3.828573121022477169108E-1, 1.107326627786831743809E+1,
     1       7.246689782858597021199E+1, 1.700632978311516129328E+2,
     2       1.698106763764238382705E+2, 7.633628843705946890896E+1,
     3       1.487967702840464066613E+1, 9.999989642347613068437E-1,
     4       1.737331760720576030932E-8/
      DATA D/8.258160008564488034698E-2, 4.344836335509282083360E+0,
     1       4.662179610356861756812E+1, 1.775728186717289799677E+2,
     2       2.953136335677908517423E+2, 2.342573504717625153053E+2,
     3       9.021658450529372642314E+1, 1.587964570758947927903E+1,
     4       1.000000000000000000000E+0/
!D    DATA C/3.828573121022477169108D-1, 1.107326627786831743809D+1,
!D   1       7.246689782858597021199D+1, 1.700632978311516129328D+2,
!D   2       1.698106763764238382705D+2, 7.633628843705946890896D+1,
!D   3       1.487967702840464066613D+1, 9.999989642347613068437D-1,
!D   4       1.737331760720576030932D-8/
!D    DATA D/8.258160008564488034698D-2, 4.344836335509282083360D+0,
!D   1       4.662179610356861756812D+1, 1.775728186717289799677D+2,
!D   2       2.953136335677908517423D+2, 2.342573504717625153053D+2,
!D   3       9.021658450529372642314D+1, 1.587964570758947927903D+1,
!D   4       1.000000000000000000000D+0/
!----------------------------------------------------------------------
! Coefficients for X < -4.0
!----------------------------------------------------------------------
      DATA E/1.3276881505637444622987E+2,3.5846198743996904308695E+4,
     1       1.7283375773777593926828E+5,2.6181454937205639647381E+5,
     2       1.7503273087497081314708E+5,5.9346841538837119172356E+4,
     3       1.0816852399095915622498E+4,1.0611777263550331766871E03,
     4       5.2199632588522572481039E+1,9.9999999999999999087819E-1/
      DATA F/3.9147856245556345627078E+4,2.5989762083608489777411E+5,
     1       5.5903756210022864003380E+5,5.4616842050691155735758E+5,
     2       2.7858134710520842139357E+5,7.9231787945279043698718E+4,
     3       1.2842808586627297365998E+4,1.1635769915320848035459E+3,
     4       5.4199632588522559414924E+1,1.0E0/
!D    DATA E/1.3276881505637444622987D+2,3.5846198743996904308695D+4,
!D   1       1.7283375773777593926828D+5,2.6181454937205639647381D+5,
!D   2       1.7503273087497081314708D+5,5.9346841538837119172356D+4,
!D   3       1.0816852399095915622498D+4,1.0611777263550331766871D03,
!D   4       5.2199632588522572481039D+1,9.9999999999999999087819D-1/
!D    DATA F/3.9147856245556345627078D+4,2.5989762083608489777411D+5,
!D   1       5.5903756210022864003380D+5,5.4616842050691155735758D+5,
!D   2       2.7858134710520842139357D+5,7.9231787945279043698718D+4,
!D   3       1.2842808586627297365998D+4,1.1635769915320848035459D+3,
!D   4       5.4199632588522559414924D+1,1.0D0/
!----------------------------------------------------------------------
!  Coefficients for rational approximation to ln(x/a), |1-x/a| < .1
!----------------------------------------------------------------------
      DATA PLG/-2.4562334077563243311E+01,2.3642701335621505212E+02,
     1         -5.4989956895857911039E+02,3.5687548468071500413E+02/
      DATA QLG/-3.5553900764052419184E+01,1.9400230218539473193E+02,
     1         -3.3442903192607538956E+02,1.7843774234035750207E+02/
!D    DATA PLG/-2.4562334077563243311D+01,2.3642701335621505212D+02,
!D   1         -5.4989956895857911039D+02,3.5687548468071500413D+02/
!D    DATA QLG/-3.5553900764052419184D+01,1.9400230218539473193D+02,
!D   1         -3.3442903192607538956D+02,1.7843774234035750207D+02/
!----------------------------------------------------------------------
! Coefficients for  0.0 < X < 6.0,
!  ratio of Chebyshev polynomials
!----------------------------------------------------------------------
      DATA P/-1.2963702602474830028590E01,-1.2831220659262000678155E03,
     1       -1.4287072500197005777376E04,-1.4299841572091610380064E06,
     2       -3.1398660864247265862050E05,-3.5377809694431133484800E08,
     3        3.1984354235237738511048E08,-2.5301823984599019348858E10,
     4        1.2177698136199594677580E10,-2.0829040666802497120940E11/
      DATA Q/ 7.6886718750000000000000E01,-5.5648470543369082846819E03,
     1        1.9418469440759880361415E05,-4.2648434812177161405483E06,
     2        6.4698830956576428587653E07,-7.0108568774215954065376E08,
     3        5.4229617984472955011862E09,-2.8986272696554495342658E10,
     4        9.8900934262481749439886E10,-8.9673749185755048616855E10/
!D    DATA P/-1.2963702602474830028590D01,-1.2831220659262000678155D03,
!D   1       -1.4287072500197005777376D04,-1.4299841572091610380064D06,
!D   2       -3.1398660864247265862050D05,-3.5377809694431133484800D08,
!D   3        3.1984354235237738511048D08,-2.5301823984599019348858D10,
!D   4        1.2177698136199594677580D10,-2.0829040666802497120940D11/
!D    DATA Q/ 7.6886718750000000000000D01,-5.5648470543369082846819D03,
!D   1        1.9418469440759880361415D05,-4.2648434812177161405483D06,
!D   2        6.4698830956576428587653D07,-7.0108568774215954065376D08,
!D   3        5.4229617984472955011862D09,-2.8986272696554495342658D10,
!D   4        9.8900934262481749439886D10,-8.9673749185755048616855D10/
!----------------------------------------------------------------------
! J-fraction coefficients for 6.0 <= X < 12.0
!----------------------------------------------------------------------
      DATA R/-2.645677793077147237806E00,-2.378372882815725244124E00,
     1       -2.421106956980653511550E01, 1.052976392459015155422E01,
     2        1.945603779539281810439E01,-3.015761863840593359165E01,
     3        1.120011024227297451523E01,-3.988850730390541057912E00,
     4        9.565134591978630774217E00, 9.981193787537396413219E-1/
      DATA S/ 1.598517957704779356479E-4, 4.644185932583286942650E00,
     1        3.697412299772985940785E02,-8.791401054875438925029E00,
     2        7.608194509086645763123E02, 2.852397548119248700147E01,
     3        4.731097187816050252967E02,-2.369210235636181001661E02,
     4        1.249884822712447891440E00/
!D    DATA R/-2.645677793077147237806D00,-2.378372882815725244124D00,
!D   1       -2.421106956980653511550D01, 1.052976392459015155422D01,
!D   2        1.945603779539281810439D01,-3.015761863840593359165D01,
!D   3        1.120011024227297451523D01,-3.988850730390541057912D00,
!D   4        9.565134591978630774217D00, 9.981193787537396413219D-1/
!D    DATA S/ 1.598517957704779356479D-4, 4.644185932583286942650D00,
!D   1        3.697412299772985940785D02,-8.791401054875438925029D00,
!D   2        7.608194509086645763123D02, 2.852397548119248700147D01,
!D   3        4.731097187816050252967D02,-2.369210235636181001661D02,
!D   4        1.249884822712447891440D00/
!----------------------------------------------------------------------
! J-fraction coefficients for 12.0 <= X < 24.0
!----------------------------------------------------------------------
      DATA P1/-1.647721172463463140042E00,-1.860092121726437582253E01,
     1        -1.000641913989284829961E01,-2.105740799548040450394E01,
     2        -9.134835699998742552432E-1,-3.323612579343962284333E01,
     3         2.495487730402059440626E01, 2.652575818452799819855E01,
     4        -1.845086232391278674524E00, 9.999933106160568739091E-1/
      DATA Q1/ 9.792403599217290296840E01, 6.403800405352415551324E01,
     1         5.994932325667407355255E01, 2.538819315630708031713E02,
     2         4.429413178337928401161E01, 1.192832423968601006985E03,
     3         1.991004470817742470726E02,-1.093556195391091143924E01,
     4         1.001533852045342697818E00/
!D    DATA P1/-1.647721172463463140042D00,-1.860092121726437582253D01,
!D   1        -1.000641913989284829961D01,-2.105740799548040450394D01,
!D   2        -9.134835699998742552432D-1,-3.323612579343962284333D01,
!D   3         2.495487730402059440626D01, 2.652575818452799819855D01,
!D   4        -1.845086232391278674524D00, 9.999933106160568739091D-1/
!D    DATA Q1/ 9.792403599217290296840D01, 6.403800405352415551324D01,
!D   1         5.994932325667407355255D01, 2.538819315630708031713D02,
!D   2         4.429413178337928401161D01, 1.192832423968601006985D03,
!D   3         1.991004470817742470726D02,-1.093556195391091143924D01,
!D   4         1.001533852045342697818D00/
!----------------------------------------------------------------------
! J-fraction coefficients for  X .GE. 24.0
!----------------------------------------------------------------------
      DATA P2/ 1.75338801265465972390E02,-2.23127670777632409550E02,
     1        -1.81949664929868906455E01,-2.79798528624305389340E01,
     2        -7.63147701620253630855E00,-1.52856623636929636839E01,
     3        -7.06810977895029358836E00,-5.00006640413131002475E00,
     4        -3.00000000320981265753E00, 1.00000000000000485503E00/
      DATA Q2/ 3.97845977167414720840E04, 3.97277109100414518365E00,
     1         1.37790390235747998793E02, 1.17179220502086455287E02,
     2         7.04831847180424675988E01,-1.20187763547154743238E01,
     3        -7.99243595776339741065E00,-2.99999894040324959612E00,
     4         1.99999999999048104167E00/
!D    DATA P2/ 1.75338801265465972390D02,-2.23127670777632409550D02,
!D   1        -1.81949664929868906455D01,-2.79798528624305389340D01,
!D   2        -7.63147701620253630855D00,-1.52856623636929636839D01,
!D   3        -7.06810977895029358836D00,-5.00006640413131002475D00,
!D   4        -3.00000000320981265753D00, 1.00000000000000485503D00/
!D    DATA Q2/ 3.97845977167414720840D04, 3.97277109100414518365D00,
!D   1         1.37790390235747998793D02, 1.17179220502086455287D02,
!D   2         7.04831847180424675988D01,-1.20187763547154743238D01,
!D   3        -7.99243595776339741065D00,-2.99999894040324959612D00,
!D   4         1.99999999999048104167D00/
!----------------------------------------------------------------------
      X = ARG
      IF (X .EQ. ZERO) THEN
            EI = -XINF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF ((X .LT. ZERO) .OR. (INT .EQ. 2)) THEN 
!----------------------------------------------------------------------
! Calculate EI for negative argument or for E1.
!----------------------------------------------------------------------
            Y = ABS(X)
            IF (Y .LE. ONE) THEN
                  SUMP = A(7) * Y + A(1)
                  SUMQ = Y + B(1)
                  DO 110 I = 2, 6
                     SUMP = SUMP * Y + A(I)
                     SUMQ = SUMQ * Y + B(I)
  110             CONTINUE
                  EI = LOG(Y) - SUMP / SUMQ
                  IF (INT .EQ. 3) EI = EI * EXP(Y)
               ELSE IF (Y .LE. FOUR) THEN
                  W = ONE / Y
                  SUMP = C(1)
                  SUMQ = D(1)
                  DO 130 I = 2, 9
                     SUMP = SUMP * W + C(I)
                     SUMQ = SUMQ * W + D(I)
  130             CONTINUE
                  EI = - SUMP / SUMQ
                  IF (INT .NE. 3) EI = EI * EXP(-Y)
               ELSE
                  IF ((Y .GT. XBIG) .AND. (INT .LT. 3)) THEN
                        EI = ZERO
                     ELSE
                        W = ONE / Y
                        SUMP = E(1) 
                        SUMQ = F(1)
                        DO 150 I = 2, 10
                           SUMP = SUMP * W + E(I)
                           SUMQ = SUMQ * W + F(I)
  150                   CONTINUE
                        EI = -W * (ONE - W * SUMP / SUMQ )
                        IF (INT .NE. 3) EI = EI * EXP(-Y)
                  END IF
            END IF
            IF (INT .EQ. 2) EI = -EI
         ELSE IF (X .LT. SIX) THEN
!----------------------------------------------------------------------
!  To improve conditioning, rational approximations are expressed
!    in terms of Chebyshev polynomials for 0 <= X < 6, and in
!    continued fraction form for larger X.
!----------------------------------------------------------------------
            T = X + X
            T = T / THREE - TWO
            PX(1) = ZERO
            QX(1) = ZERO
            PX(2) = P(1)
            QX(2) = Q(1)
            DO 210 I = 2, 9
               PX(I+1) = T * PX(I) - PX(I-1) + P(I)
               QX(I+1) = T * QX(I) - QX(I-1) + Q(I)
  210       CONTINUE
            SUMP = HALF * T * PX(10) - PX(9) + P(10)
            SUMQ = HALF * T * QX(10) - QX(9) + Q(10)
            FRAC = SUMP / SUMQ
            XMX0 = (X - X01/X11) - X02
            IF (ABS(XMX0) .GE. P037) THEN
                  EI = LOG(X/X0) + XMX0 * FRAC
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
               ELSE
!----------------------------------------------------------------------
! Special approximation to  ln(X/X0)  for X close to X0
!----------------------------------------------------------------------
                  Y = XMX0 / (X + X0)
                  YSQ = Y*Y
                  SUMP = PLG(1)
                  SUMQ = YSQ + QLG(1)
                  DO 220 I = 2, 4
                     SUMP = SUMP*YSQ + PLG(I)
                     SUMQ = SUMQ*YSQ + QLG(I)
  220             CONTINUE
                  EI = (SUMP / (SUMQ*(X+X0)) + FRAC) * XMX0
                  IF (INT .EQ. 3) EI = EXP(-X) * EI
            END IF
         ELSE IF (X .LT. TWELVE) THEN
            FRAC = ZERO
            DO 230 I = 1, 9
               FRAC = S(I) / (R(I) + X + FRAC)
  230       CONTINUE
            EI = (R(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE IF (X .LE. TWO4) THEN
            FRAC = ZERO
            DO 240 I = 1, 9
               FRAC = Q1(I) / (P1(I) + X + FRAC)
  240       CONTINUE
            EI = (P1(10) + FRAC) / X
            IF (INT .NE. 3) EI = EI * EXP(X)
         ELSE
            IF ((X .GE. XMAX) .AND. (INT .LT. 3)) THEN
                  EI = XINF
               ELSE
                  Y = ONE / X
                  FRAC = ZERO
                  DO 250 I = 1, 9
                     FRAC = Q2(I) / (P2(I) + X + FRAC)
  250             CONTINUE
                  FRAC = P2(10) + FRAC
                  EI = Y + Y * Y * FRAC
                  IF (INT .NE. 3) THEN
                        IF (X .LE. XMAX-TWO4) THEN
                              EI = EI * EXP(X)
                           ELSE
!----------------------------------------------------------------------
! Calculation reformulated to avoid premature overflow
!----------------------------------------------------------------------
                              EI = (EI * EXP(X-FOURTY)) * EXP40
                        END IF
                  END IF
            END IF
      END IF
      RESULT = EI
      RETURN
!---------- Last line of CALCEI ----------
      END
      FUNCTION EI(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   exponential integral  Ei(x), where  x  is real.
!
!  Author: W. J. Cody
!
!  Latest modification: January 12, 1988
!
!--------------------------------------------------------------------
      INTEGER INT
      REAL  EI, X, RESULT
!D    DOUBLE PRECISION  EI, X, RESULT
!--------------------------------------------------------------------
      INT = 1
      CALL CALCEI(X,RESULT,INT)
      EI = RESULT
      RETURN
!---------- Last line of EI ----------
      END
      FUNCTION EXPEI(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   function  exp(-x) * Ei(x), where  Ei(x)  is the exponential
!   integral, and  x  is real.
!
!  Author: W. J. Cody
!
!  Latest modification: January 12, 1988
!
!--------------------------------------------------------------------
      INTEGER INT
      REAL  EXPEI, X, RESULT
!D    DOUBLE PRECISION  EXPEI, X, RESULT
!--------------------------------------------------------------------
      INT = 3
      CALL CALCEI(X,RESULT,INT)
      EXPEI = RESULT
      RETURN
!---------- Last line of EXPEI ----------
      END
      FUNCTION EONE(X)
!--------------------------------------------------------------------
!
! This function program computes approximate values for the
!   exponential integral E1(x), where  x  is real.
!
!  Author: W. J. Cody
!
!  Latest modification: January 12, 1988
!
!--------------------------------------------------------------------
      INTEGER INT
      REAL  EONE, X, RESULT
!D    DOUBLE PRECISION  EONE, X, RESULT
!--------------------------------------------------------------------
      INT = 2
      CALL CALCEI(X,RESULT,INT)
      EONE = RESULT
      RETURN
!---------- Last line of EONE ----------
      END
