      SUBROUTINE RANLUX(RVEC,LENV)
!         Subtract-and-borrow random number generator proposed by
!         Marsaglia and Zaman, implemented by F. James with the name
!         RCARRY in 1991, and later improved by Martin Luescher
!         in 1993 to produce "Luxury Pseudorandom Numbers".
!     Fortran 77 coded by F. James, 1993
! Edited to remove status messages Ian Hutchinson 2014
!          
!       references: 
!  M. Luscher, Computer Physics Communications  79 (1994) 100
!  F. James, Computer Physics Communications 79 (1994) 111
!
!   LUXURY LEVELS.
!   ------ ------      The available luxury levels are:
!
!  level 0  (p=24): equivalent to the original RCARRY of Marsaglia
!           and Zaman, very long period, but fails many tests.
!  level 1  (p=48): considerable improvement in quality over level 0,
!           now passes the gap test, but still fails spectral test.
!  level 2  (p=97): passes all known tests, but theoretically still
!           defective.
!  level 3  (p=223): DEFAULT VALUE.  Any theoretically possible
!           correlations have very small chance of being observed.
!  level 4  (p=389): highest possible luxury, all 24 bits chaotic.
!
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!!!!  Calling sequences for RANLUX:                                  ++
!!!!      CALL RANLUX (RVEC, LEN)   returns a vector RVEC of LEN     ++
!!!!                   32-bit random floating point numbers between  ++
!!!!                   zero (not included) and one (also not incl.). ++
!!!!      CALL RLUXGO(LUX,INT,K1,K2) initializes the generator from  ++
!!!!               one 32-bit integer INT and sets Luxury Level LUX  ++
!!!!               which is integer between zero and MAXLEV, or if   ++
!!!!               LUX .GT. 24, it sets p=LUX directly.  K1 and K2   ++
!!!!               should be set to zero unless restarting at a break++ 
!!!!               point given by output of RLUXAT (see RLUXAT).     ++
!!!!      CALL RLUXAT(LUX,INT,K1,K2) gets the values of four integers++
!!!!               which can be used to restart the RANLUX generator ++
!!!!               at the current point by calling RLUXGO.  K1 and K2++
!!!!               specify how many numbers were generated since the ++
!!!!               initialization with LUX and INT.  The restarting  ++
!!!!               skips over  K1+K2*E9   numbers, so it can be long.++
!!!!   A more efficient but less convenient way of restarting is by: ++
!!!!      CALL RLUXIN(ISVEC)    restarts the generator from vector   ++
!!!!                   ISVEC of 25 32-bit integers (see RLUXUT)      ++
!!!!      CALL RLUXUT(ISVEC)    outputs the current values of the 25 ++
!!!!                 32-bit integer seeds, to be used for restarting ++
!!!!      ISVEC must be dimensioned 25 in the calling program        ++
!!!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      DIMENSION RVEC(LENV)
      DIMENSION SEEDS(24), ISEEDS(24), ISDEXT(25)
      PARAMETER (MAXLEV=4, LXDFLT=3)
      DIMENSION NDSKIP(0:MAXLEV)
      DIMENSION NEXT(24)
      PARAMETER (TWOP12=4096., IGIGA=1000000000,JSDFLT=314159265)
      PARAMETER (ITWO24=2**24, ICONS=2147483563)
      SAVE NOTYET, I24, J24, CARRY, SEEDS, TWOM24, TWOM12, LUXLEV
      SAVE NSKIP, NDSKIP, IN24, NEXT, KOUNT, MKOUNT, INSEED
      INTEGER LUXLEV
      LOGICAL NOTYET
! IHH
      logical messages
      data messages/.false./
      DATA NOTYET, LUXLEV, IN24, KOUNT, MKOUNT /.TRUE., LXDFLT, 0,0,0/
      DATA I24,J24,CARRY/24,10,0./
!                               default
!  Luxury Level   0     1     2   *3*    4
      DATA NDSKIP/0,   24,   73,  199,  365 /
!orresponds to p=24    48    97   223   389
!     time factor 1     2     3     6    10   on slow workstation
!                 1    1.5    2     3     5   on fast mainframe
!
!  NOTYET is .TRUE. if no initialization has been performed yet.
!              Default Initialization by Multiplicative Congruential
      IF (NOTYET) THEN
         NOTYET = .FALSE.
         JSEED = JSDFLT  
         INSEED = JSEED
! IHH
         if(messages)then
         WRITE(6,'(A,I12)') ' RANLUX DEFAULT INITIALIZATION: ',JSEED
         endif
         LUXLEV = LXDFLT
         NSKIP = NDSKIP(LUXLEV)
         LP = NSKIP + 24
         IN24 = 0
         KOUNT = 0
         MKOUNT = 0
! IHH
         if(messages)then
         WRITE(6,'(A,I2,A,I4)')  ' RANLUX DEFAULT LUXURY LEVEL =  ',
     +        LUXLEV,'      p =',LP
         endif
         TWOM24 = 1.
         DO 25 I= 1, 24
            TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
   25    CONTINUE
         TWOM12 = TWOM24 * 4096.
         DO 50 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
   50    CONTINUE
         NEXT(1) = 24
         I24 = 24
         J24 = 10
         CARRY = 0.
         IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
      ENDIF
!
!          The Generator proper: "Subtract-with-borrow",
!          as proposed by Marsaglia and Zaman,
!          Florida State University, March, 1989
!
      DO 100 IVEC= 1, LENV
      UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
      IF (UNI .LT. 0.)  THEN
         UNI = UNI + 1.0
         CARRY = TWOM24
      ELSE
         CARRY = 0.
      ENDIF
      SEEDS(I24) = UNI
      I24 = NEXT(I24)
      J24 = NEXT(J24)
      RVEC(IVEC) = UNI
!  small numbers (with less than 12 "significant" bits) are "padded".
      IF (UNI .LT. TWOM12)  THEN
         RVEC(IVEC) = RVEC(IVEC) + TWOM24*SEEDS(J24)
!        and zero is forbidden in case someone takes a logarithm
         IF (RVEC(IVEC) .EQ. 0.)  RVEC(IVEC) = TWOM24*TWOM24
      ENDIF
!        Skipping to luxury.  As proposed by Martin Luscher.
      IN24 = IN24 + 1
      IF (IN24 .EQ. 24)  THEN
         IN24 = 0
         KOUNT = KOUNT + NSKIP
         DO 90 ISK= 1, NSKIP
         UNI = SEEDS(J24) - SEEDS(I24) - CARRY
         IF (UNI .LT. 0.)  THEN
            UNI = UNI + 1.0
            CARRY = TWOM24
         ELSE
            CARRY = 0.
         ENDIF
         SEEDS(I24) = UNI
         I24 = NEXT(I24)
         J24 = NEXT(J24)
   90    CONTINUE
      ENDIF
  100 CONTINUE
      KOUNT = KOUNT + LENV
      IF (KOUNT .GE. IGIGA)  THEN
         MKOUNT = MKOUNT + 1
         KOUNT = KOUNT - IGIGA
      ENDIF
      RETURN
!
!           Entry to input and float integer seeds from previous run
      ENTRY RLUXIN(ISDEXT)
      TWOM24 = 1.
      DO 195 I= 1, 24
         NEXT(I) = I-1
 195     TWOM24 = TWOM24 * 0.5
      NEXT(1) = 24
      TWOM12 = TWOM24 * 4096.
! IHH
      if(messages)then
      WRITE(6,'(A)') ' FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:'
      WRITE(6,'(5X,5I12)') ISDEXT
      endif
      DO 200 I= 1, 24
      SEEDS(I) = REAL(ISDEXT(I))*TWOM24
  200 CONTINUE
      CARRY = 0.
      IF (ISDEXT(25) .LT. 0)  CARRY = TWOM24
      ISD = IABS(ISDEXT(25))
      I24 = MOD(ISD,100)
      ISD = ISD/100
      J24 = MOD(ISD,100)
      ISD = ISD/100
      IN24 = MOD(ISD,100)
      ISD = ISD/100
      LUXLEV = ISD
        IF (LUXLEV .LE. MAXLEV) THEN
          NSKIP = NDSKIP(LUXLEV)
! IHH
          if(messages)then
          WRITE (6,'(A,I2)') ' RANLUX LUXURY LEVEL SET BY RLUXIN TO: ',
     +                         LUXLEV
          endif
        ELSE  IF (LUXLEV .GE. 24) THEN
          NSKIP = LUXLEV - 24
! IHH
          if(messages)then
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXIN TO:',LUXLEV
          endif
        ELSE
          NSKIP = NDSKIP(MAXLEV)
          WRITE (6,'(A,I5)') ' RANLUX ILLEGAL LUXURY RLUXIN: ',LUXLEV
          LUXLEV = MAXLEV
        ENDIF
      INSEED = -1
      RETURN
!
!                    Entry to ouput seeds as integers
      ENTRY RLUXUT(ISDEXT)
      DO 300 I= 1, 24
         ISDEXT(I) = INT(SEEDS(I)*TWOP12*TWOP12)
  300 CONTINUE
      ISDEXT(25) = I24 + 100*J24 + 10000*IN24 + 1000000*LUXLEV
      IF (CARRY .GT. 0.)  ISDEXT(25) = -ISDEXT(25)
      RETURN
!
!                    Entry to output the "convenient" restart point
      ENTRY RLUXAT(LOUT,INOUT,K1,K2)
      LOUT = LUXLEV
      INOUT = INSEED
      K1 = KOUNT
      K2 = MKOUNT
      RETURN
!
!                    Entry to initialize from one or three integers
      ENTRY RLUXGO(LUX,INS,K1,K2)
         IF (LUX .LT. 0) THEN
            LUXLEV = LXDFLT
         ELSE IF (LUX .LE. MAXLEV) THEN
            LUXLEV = LUX
         ELSE IF (LUX .LT. 24 .OR. LUX .GT. 2000) THEN
            LUXLEV = MAXLEV
            WRITE (6,'(A,I7)') ' RANLUX ILLEGAL LUXURY RLUXGO: ',LUX
         ELSE
            LUXLEV = LUX
            DO 310 ILX= 0, MAXLEV
              IF (LUX .EQ. NDSKIP(ILX)+24)  LUXLEV = ILX
  310       CONTINUE
         ENDIF
      IF (LUXLEV .LE. MAXLEV)  THEN
         NSKIP = NDSKIP(LUXLEV)
! IHH
         if(messages)then
         WRITE(6,'(A,I2,A,I4)') ' RANLUX LUXURY LEVEL SET BY RLUXGO :',
     +        LUXLEV,'     P=', NSKIP+24
         endif
      ELSE
          NSKIP = LUXLEV - 24
! IHH
          if(messages)then
          WRITE (6,'(A,I5)') ' RANLUX P-VALUE SET BY RLUXGO TO:',LUXLEV
          endif
      ENDIF
      IN24 = 0
      IF (INS .LT. 0)  WRITE (6,'(A)')   
     +   ' Illegal initialization by RLUXGO, negative input seed'
      IF (INS .GT. 0)  THEN
        JSEED = INS
! IHH
        if(messages)then
        WRITE(6,'(A,3I12)') ' RANLUX INITIALIZED BY RLUXGO FROM SEEDS',
     +      JSEED, K1,K2
        endif
      ELSE
        JSEED = JSDFLT
! IHH
        if(messages)then
        WRITE(6,'(A)')' RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED'
        endif
      ENDIF
      INSEED = JSEED
      NOTYET = .FALSE.
      TWOM24 = 1.
         DO 325 I= 1, 24
           TWOM24 = TWOM24 * 0.5
         K = JSEED/53668
         JSEED = 40014*(JSEED-K*53668) -K*12211
         IF (JSEED .LT. 0)  JSEED = JSEED+ICONS
         ISEEDS(I) = MOD(JSEED,ITWO24)
  325    CONTINUE
      TWOM12 = TWOM24 * 4096.
         DO 350 I= 1,24
         SEEDS(I) = REAL(ISEEDS(I))*TWOM24
         NEXT(I) = I-1
  350    CONTINUE
      NEXT(1) = 24
      I24 = 24
      J24 = 10
      CARRY = 0.
      IF (SEEDS(24) .EQ. 0.) CARRY = TWOM24
!        If restarting at a break point, skip K1 + IGIGA*K2
!        Note that this is the number of numbers delivered to
!        the user PLUS the number skipped (if luxury .GT. 0).
      KOUNT = K1
      MKOUNT = K2
      IF (K1+K2 .NE. 0)  THEN
        DO 500 IOUTER= 1, K2+1
          INNER = IGIGA
          IF (IOUTER .EQ. K2+1)  INNER = K1
          DO 450 ISK= 1, INNER
            UNI = SEEDS(J24) - SEEDS(I24) - CARRY 
            IF (UNI .LT. 0.)  THEN
               UNI = UNI + 1.0
               CARRY = TWOM24
            ELSE
               CARRY = 0.
            ENDIF
            SEEDS(I24) = UNI
            I24 = NEXT(I24)
            J24 = NEXT(J24)
  450     CONTINUE
  500   CONTINUE
!         Get the right value of IN24 by direct calculation
        IN24 = MOD(KOUNT, NSKIP+24)
        IF (MKOUNT .GT. 0)  THEN
           IZIP = MOD(IGIGA, NSKIP+24)
           IZIP2 = MKOUNT*IZIP + IN24
           IN24 = MOD(IZIP2, NSKIP+24)
        ENDIF
!       Now IN24 had better be between zero and 23 inclusive
        IF (IN24 .GT. 23) THEN
           WRITE (6,'(A/A,3I11,A,I5)')  
     +    '  Error in RESTARTING with RLUXGO:','  The values', INS,
     +     K1, K2, ' cannot occur at luxury level', LUXLEV
           IN24 = 0
        ENDIF
      ENDIF
      RETURN
      END
! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! Make this file consistent with calling ranlux in other contexts.
      subroutine luxtst
! To run as a test program comment the above and uncomment the below.
!      PROGRAM LUXTST
!         Exercise for the RANLUX Pseudorandom number generator.
!
      DIMENSION RVEC(1000)
      DIMENSION ISDEXT(25)
!
!         check that we get the right numbers (machine-indep.)
      WRITE (6,'(/A)')  '  CALL RANLUX(RVEC,100)'
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX default numbers 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')  ' CALL RLUXGO(0,0,0,0)'
      CALL RLUXGO(0,0,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury level 0, 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')  '   CALL RLUXGO(389,1,0,0)'
      CALL RLUXGO(389,1,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p=389, 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')  '  CALL RLUXGO(75,0,0,0)'
      CALL RLUXGO(75,0,0,0)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75,   1-  5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX luxury p= 75, 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')  '  test restarting from the full vector'
      CALL RLUXUT(ISDEXT)
      WRITE (6,'(/A/(1X,5I14))') '  current RANLUX status saved:',ISDEXT
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')   '   previous RANLUX status will be restored'
      CALL RLUXIN(ISDEXT)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 1- 5:',
     +    (RVEC(L),L=1,5)
      CALL RANLUX(RVEC,100)
      WRITE (6,'(A/9X,5F12.8)') ' RANLUX numbers 101-105:',
     +    (RVEC(L),L=1,5)
!
      WRITE (6,'(/A)')  '     test the restarting by skipping'
      CALL RLUXGO(4,7674985,0,0)
      CALL RLUXAT(I1,I2,I3,I4)
      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
      DO 150 LI= 1, 10
  150 CALL RANLUX(RVEC,1000)
      CALL RLUXAT(I1,I2,I3,I4)
      WRITE (6,'(A,4I10)')  '  RLUXAT values =',I1,I2,I3,I4
      CALL RANLUX(RVEC,200)
      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
     +                             RVEC(1), RVEC(200)
      CALL RLUXGO(I1,I2,I3,I4)
      CALL RANLUX(RVEC,200)
      WRITE (6,'(A,2F10.6)')  '  Next and 200th numbers are:',
     +                             RVEC(1), RVEC(200)
!
      WRITE (6,'(/A)') ' The following should provoke an error message'
      CALL RLUXGO(4,11111,31,0)
      STOP
!
!   OUTPUT FROM THE ABOVE TEST PROGRAM SHOULD BE:
!   --------------------------------------------
!  CALL RANLUX(RVEC,100)
! RANLUX DEFAULT INITIALIZATION:    314159265
! RANLUX DEFAULT LUXURY LEVEL =   3      p = 223
! RANLUX default numbers   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX default numbers 101-105:
!           0.43156743  0.03774416  0.24897110  0.00147784  0.90274453
!
!  CALL RLUXGO(0,0,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 0     P=  24
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury level 0,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury level 0, 101-105:
!           0.41538775  0.05330932  0.58195311  0.91397446  0.67034441
!
!   CALL RLUXGO(389,1,0,0)
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS           1           0           0
! RANLUX luxury p=389,   1-  5:
!           0.94589490  0.47347850  0.95152789  0.42971975  0.09127384
! RANLUX luxury p=389, 101-105:
!           0.02618265  0.03775346  0.97274780  0.13302165  0.43126065
!
!  CALL RLUXGO(75,0,0,0)
! RANLUX P-VALUE SET BY RLUXGO TO:   75
! RANLUX INITIALIZED BY RLUXGO FROM DEFAULT SEED
! RANLUX luxury p= 75,   1-  5:
!           0.53981817  0.76155043  0.06029940  0.79600263  0.30631220
! RANLUX luxury p= 75, 101-105:
!           0.25600731  0.23443210  0.59164381  0.59035838  0.07011414
!
!  test restarting from the full vector
!
!  current RANLUX status saved:
!       16156027      16534309      15243811       2751687       6002207
!        7979506       1301976       4567313       4305996       5872599
!       12003090       2146823      12606367       4111505       5979640
!       12739666      10489318      14036909      11729352       8061448
!        7832659       6069758       3197719       1832730      75080216
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
!
!   previous RANLUX status will be restored
! FULL INITIALIZATION OF RANLUX WITH 25 INTEGERS:
!         16156027    16534309    15243811     2751687     6002207
!          7979506     1301976     4567313     4305996     5872599
!         12003090     2146823    12606367     4111505     5979640
!         12739666    10489318    14036909    11729352     8061448
!          7832659     6069758     3197719     1832730    75080216
! RANLUX P-VALUE SET BY RLUXIN TO:   75
! RANLUX numbers 1- 5:
!           0.22617835  0.60655993  0.86417443  0.43920082  0.23382509
! RANLUX numbers 101-105:
!           0.08107197  0.21466845  0.84856731  0.94078046  0.85626233
!
!     test the restarting by skipping
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985           0           0
!  RLUXAT values =         4   7674985         0         0
!  RLUXAT values =         4   7674985    161840         0
!  Next and 200th numbers are:  0.019648  0.590586
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS     7674985      161840           0
!  Next and 200th numbers are:  0.019648  0.590586
!
! The following should provoke an error message
! RANLUX LUXURY LEVEL SET BY RLUXGO : 4     P= 389
! RANLUX INITIALIZED BY RLUXGO FROM SEEDS       11111          31           0
!  Error in RESTARTING with RLUXGO:
!  The values      11111         31          0 cannot occur at luxury level    4
      END
!**************************************************************************
!***********************************************************************
! Restartable gasdev based on NR. Calculates gaussian-distributed
! random numbers two at a time, costing a log and a sqrt.
      FUNCTION gasdev(ireset)
      integer ireset
      real vr(2),v1,v2
      equivalence (v1,vr(1)),(v2,vr(2))
      include 'rancom.f'
      logical lfirst
      data lfirst/.true./

      if(lfirst)gd_iset=0
      if(ireset.lt.0)then
         gd_iset=0
         write(*,*)'gasdev reset',ireset
      endif
      if(gd_iset.ne.1)then
 1       continue
         call ranlux(vr,2)
         v1=2.*v1-1.
         v2=2.*v2-1.
         r=v1**2+v2**2
        if(r.ge.1..or.r.eq.0.)go to 1
        fac=sqrt(-2.*log(r)/r)
        gd_gset=v1*fac
        gasdev=v2*fac
        gd_iset=1
      else
        gasdev=gd_gset
        gd_iset=0
      endif
      end
!**********************************************************************
