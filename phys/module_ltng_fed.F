!WRF:MODEL_LAYER:PHYSICS
!
MODULE MODULE_LTNG_FED

! USE WRF PHYSICS CONSTANTS
      USE MODULE_MODEL_CONSTANTS, ONLY:                           &
                                        PICONST                   &
                                      , GRAV => G

! USE GAMMA FUNCTION FROM THOMPSON MICROPHYSICS
      USE MODULE_MP_THOMPSON,     ONLY: WGAMMA

CONTAINS
!===================================================================
!
      SUBROUTINE CALCFED(                                         &
                         QV,QC,QR,QI,QS,QG,QH                     &
                        ,W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY     & ! DZ8W is thickness (meters) of each model layer
                        ,DX,DY                                    &
                        !FLASH_EXTENT_DENSITY                     &
                        ,FED                                      &
                        ,IDS,IDE, JDS,JDE, KDS,KDE                & ! domain
                        ,IMS,IME, JMS,JME, KMS,KME                & ! memory
                        ,ITS,ITE, JTS,JTE, KTS,KTE                & ! tile
                        ,FED1, FED2, FED3, FED4, FED5, FED6       & ! OPTIONAL
                        ,FED7, FED8, FED9                         &
                                                                  )

      USE MODULE_WRF_ERROR
!-------------------------------------------------------------------
      IMPLICIT NONE
!-------------------------------------------------------------------
!-------------------------------------------------------------------

      INTEGER, INTENT(IN) ::          IDS,IDE, JDS,JDE, KDS,KDE , &
                                      IMS,IME, JMS,JME, KMS,KME , &
                                      ITS,ITE, JTS,JTE, KTS,KTE
      REAL, DIMENSION( IMS:IME , KMS:KME , JMS:JME ),             &
            INTENT(IN) ::                                         &
                                                             QV , &
                                                             QC , &
                                                             QR , &
                                                             QI , &
                                                             QS , &
                                                             QG , &
                                                             QH

      REAL, DIMENSION( IMS:IME, KMS:KME, JMS:JME ),               &
            INTENT(IN ) :: W,Z,DZ8W,PI_PHY,TH_PHY,P_PHY,RHO_PHY

      REAL, INTENT(IN ) :: DX, DY ! GRID SPACING IN METERS

      REAL, INTENT(OUT), DIMENSION(IMS:IME,JMS:JME)::             &
          !FLASH_EXTENT_DENSITY
          FED

      REAL, INTENT(OUT), DIMENSION(IMS:IME,JMS:JME), OPTIONAL ::  &
          FED1, FED2, FED3, FED4, FED5, FED6, FED7, FED8, FED9

! DEFINE LOCAL DATA
      REAL, DIMENSION(KTS:KTE):: TEMPC
      REAL :: AIR_VOL, AIR_MASS, GRAUP_MASS, ICE_MASS, SNOW_MASS
      REAL :: TMP1, TMP2, VTG, VTS
      REAL :: ICE_MASS_FLUX, PRECIP_MASS_FLUX
      REAL :: IWP
      REAL :: FLUXPROD, LIGHTFLUX, LIGHTIWP, LIGHTPRECIP_ICE
      REAL :: LIGHT_GRAUPEL_VOL, LIGHT_WMAX
      REAL :: W_VOL, LIGHTWVOL

      INTEGER         ::            I,J,K

      REAL, PARAMETER :: TMAXC    =  -5.    ! TEMPERATURE THRESHOLD IN CELSIUS
      REAL, PARAMETER :: WMIN     =   5.    ! VERTICAL VEL THRESHOLD (M/S)
      REAL, PARAMETER :: CELS2K   = 273.15
      REAL, PARAMETER :: XNOG     =   4.0E4 !LIN SNOW FALL SPEED
      REAL, PARAMETER :: CDRAG    =   0.6
      REAL, PARAMETER :: RHOGRAUL = 917. ! DENSITY OF GRAUPEL KG/M^3 ! ONLY 500 IN MODULE_MP_THOMPSON.F
      REAL, PARAMETER :: XNOS     =   3.0E6
      REAL, PARAMETER :: RHOSNOW  = 100.! KG/M3 DENSITY OF SNOW
      REAL, PARAMETER :: CONSTD   =   0.25
      REAL, PARAMETER :: EXPON    = 1-CONSTD
      REAL, PARAMETER :: CONSTC   = 152.93*0.01**EXPON
      REAL, PARAMETER :: GAMDP4   = GAMMA(CONSTD+4.)
      REAL, PARAMETER :: GAM4PT5  = GAMMA(4.5)

      CHARACTER*256 :: MESSAGE
!-------------------------------------------------------------------

      CALL WRF_DEBUG(150, ' LTNG_FED: ENTERED CALCFED')

      DO J = JTS,JTE
         DO I = ITS,ITE

            ! INITIALIZE VERTICAL TOTALS TO ZERO
            AIR_MASS         = 0.
            GRAUP_MASS       = 0.
            ICE_MASS         = 0.
            SNOW_MASS        = 0.
            PRECIP_MASS_FLUX = 0.
            ICE_MASS_FLUX    = 0.
            W_VOL            = 0.

            ! CALCULATE MASS AND FLUX WITH A SINGLE LOOP THRU THE VERTICAL
            DO K = KTS,KTE
               WRITE(MESSAGE,'(" LTNG_FED: K=",I8)') K
               CALL WRF_DEBUG(600, MESSAGE)

               ! TEMPERATURE OF LAYER (CELSIUS)
               TEMPC(K)     = TH_PHY(I,K,J)*PI_PHY(I,K,J)-CELS2K

               ! CHECK TEMPERATURE
               ! ADD TO MASS OF AIR AND GRAUPEL IN GRID VOLUME IF TEMPERATURE IS AT OR BELOW TMAXC.
               ! ALSO ADD TO ICE AND PRECIP MASS FLUX IF TEMP <= TMAXC
               IF(TEMPC(K).LE.TMAXC)THEN 
                   
                  ! COULD FORCE POSITIVE MIXING RATIOS HERE, BUT NEGATIVE VALUES ARE SMALL?

                  AIR_VOL    = DX * DY * DZ8W(I,K,J)
                  AIR_MASS   = AIR_MASS   + (AIR_VOL * RHO_PHY(I,K,J)) 
                  GRAUP_MASS = GRAUP_MASS + (QG(I,K,J) * AIR_MASS)
                  ICE_MASS   = ICE_MASS   + (QI(I,K,J) * AIR_MASS)
                  SNOW_MASS  = SNOW_MASS  + (QS(I,K,J) * AIR_MASS)

                  ! ADD TO PRECIP MASS FLUX IF QG>1E-8
                  IF(QG(I,K,J).GT.1E-8)THEN

                     TMP1 = 1./(PICONST*RHOGRAUL*XNOG/(RHO_PHY(I,K,J)*QG(I,K,J)))**0.25
                     TMP2 = SQRT(4.*0.33334*GRAV*RHOGRAUL/(CDRAG*RHO_PHY(I,K,J)))
                     VTG = (1./6.)*GAM4PT5*TMP2*SQRT(TMP1)

                     PRECIP_MASS_FLUX = PRECIP_MASS_FLUX + ( GRAUP_MASS * VTG ) ! IF YOU SUBTRACT W FROM VTG, YOU MAY GET NEGATIVE FLUX
                  END IF

                  ! ADD TO ICE MASS FLUX IF QS>1E-8
                  IF(QS(I,K,J).GT.1E-8)THEN

                     TMP1 = (PICONST*RHOSNOW*XNOS/(RHO_PHY(I,K,J)*QS(I,K,J)))**0.25
                     TMP2 = SQRT(1.29/RHO_PHY(I,K,J))
                     VTS = (1./6.) * CONSTC * GAMDP4 * TMP2 / (TMP1**CONSTD)

                     ! IF VERTICAL MOTION EXCEEDS SNOW FALL SPEED AND VERTICAL MOTION > 0,
                     !  ADD SNOW AND ICE MASS FLUX TO TOTAL ICE MASS FLUX
                     IF(W(I,K,J).GT.VTS .AND. W(I,K,J).GT.0) THEN
                        ICE_MASS_FLUX = ICE_MASS_FLUX + SNOW_MASS * (W(I,K,J)-VTS) + ICE_MASS * W(I,K,J)
                     ELSE
                        ! ELSE IF FALL SPEED EXCEEDS VERTICAL MOTION OR VERTICAL MOTION < 0,
                        !   ADD SNOW MASS FLUX TO PRECIP MASS FLUX
                        PRECIP_MASS_FLUX = PRECIP_MASS_FLUX + ( SNOW_MASS * VTS )
                     END IF
                   END IF

                   ! ADD TO UPDRAFT VOLUME IF W EXCEEDS W THRESHOLD
                   IF (W(I,K,J).GT.WMIN) W_VOL = W_VOL + AIR_VOL

               END IF
            END DO

            IWP       = SNOW_MASS + ICE_MASS + GRAUP_MASS

            FLUXPROD  = ICE_MASS_FLUX * PRECIP_MASS_FLUX


            ! VARIOUS LIGHTNING FLASH EXTENT DENSITY ESTIMATES
            LIGHTFLUX       = FLUXPROD * 9E-17
            LIGHTIWP        =     IWP  * 0.20 / 2000**2 * DX * DY ! McCaul et al. 2009
            LIGHTPRECIP_ICE =       7. * GRAUP_MASS *3.4E-8
            LIGHTWVOL       =    W_VOL * 6.75E-10                ! Deierling and Petersen 2008

            !FLASH_EXTENT_DENSITY(I,J) = LIGHTFLUX
            FED(I,J)  = LIGHTIWP
            FED1(I,J) = LIGHTWVOL

         END DO
      END DO
      RETURN
      END SUBROUTINE CALCFED

  END MODULE MODULE_LTNG_FED
