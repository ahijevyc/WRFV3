!WRF:MEDIATION_LAYER:PHYSICS
!
! Contains initialization subroutine lightning_init and driver subroutine
! lightning_driver.
!
! History:
!   3.5?  - rewritten and added init, separate out flash rate
!           parameterization from emission
!   3.4.0 - Added cpm option
!   3.3.x - lightning_driver written by M. Barth called in
!           emission_driver in chem
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************

 MODULE module_lightning_driver
 CONTAINS

!**********************************************************************
!
! SUBROUTINE lightning_init
!
! Performs compatibility checks and zero out flash arrays at first timestep.
!
!**********************************************************************

 SUBROUTINE lightning_init ( &
                              itimestep, restart, dt, dx               &
                            ! Namelist control options
                             ,cu_physics,mp_physics,do_radar_ref       &
                             ,lightning_option, lightning_dt           &
                             ,lightning_start_seconds                  &
                             ,iccg_prescribed_num, iccg_prescribed_den &
                             ,cellcount_method                         &
                            ! Order dependent args for domain, mem, and tile dims
                             ,ids, ide, jds, jde, kds, kde             &
                             ,ims, ime, jms, jme, kms, kme             &
                             ,its, ite, jts, jte, kts, kte             &
                            ! IC and CG flash rates and accumulated flash count
                             ,ic_flashcount, ic_flashrate              &
                             ,cg_flashcount, cg_flashrate              &
#ifdef WRF_CHEM
                             ,lnox_opt,lnox_passive                    &
                            ! LNOx tracers (chemistry only)
                             ,lnox_total, lnox_ic, lnox_cg             &
#endif
                            )
!-----------------------------------------------------------------
 USE module_state_description
 USE module_wrf_error
 IMPLICIT NONE
!-----------------------------------------------------------------

 INTEGER,  INTENT(IN)        :: itimestep
 LOGICAL,  INTENT(IN)        :: restart
 REAL,     INTENT(IN)        :: dt,dx
 INTEGER,  INTENT(IN)        :: cu_physics,mp_physics,do_radar_ref,lightning_option
 REAL,     INTENT(IN)        :: lightning_dt, lightning_start_seconds
 REAL,     INTENT(IN)        :: iccg_prescribed_num, iccg_prescribed_den
 INTEGER,  INTENT(INOUT)     :: cellcount_method
 INTEGER , INTENT(IN)        :: ids, ide, jds, jde, kds, kde,  &
                                ims, ime, jms, jme, kms, kme,  &
                                its, ite, jts, jte, kts, kte

! Making these optional just in case qualitative lightning indices get implemented
 REAL, OPTIONAL, DIMENSION( ims:ime,jms:jme ), &
                 INTENT(OUT) :: ic_flashcount, ic_flashrate, &
                                cg_flashcount, cg_flashrate

#ifdef WRF_CHEM
 INTEGER, INTENT(IN)         :: lnox_opt
 LOGICAL, INTENT(IN)         :: lnox_passive
 REAL, OPTIONAL, DIMENSION( ims:ime,kms:kme,jms:jme ), &
                 INTENT(OUT) :: lnox_total, lnox_ic, lnox_cg
#endif

 CHARACTER (LEN=80) :: message

!-----------------------------------------------------------------

!-- do not reset unless it is the first timestep or lightning_option is on
 IF (itimestep .gt. 0 .or. lightning_option .eq. 0) return

!-- check to see if lightning_dt is a proper multiple of dt
 IF ( MOD(lightning_dt,dt) .ne. 0. ) THEN
    CALL wrf_error_fatal (' lightning_init: lightning_dt needs to be a multiple of model time step dt')
 ENDIF

!--  check to see if the prescribed IC:CG ratio is valid (0/0 and -1 are not allowed)
 IF (iccg_prescribed_den .eq. 0. .and. iccg_prescribed_num .eq. 0.) THEN
    CALL wrf_error_fatal (' lightning_init: iccg_prescribed cannot be 0.0/0.0')
 ENDIF
 IF (iccg_prescribed_den .ne. 0.) THEN
    IF (iccg_prescribed_num/iccg_prescribed_den .eq. -1.) THEN
        CALL wrf_error_fatal (' lightning_init: iccg_prescribed cannot be -1')
    ENDIF
 ENDIF

!
!-- check to see if lightning_option is valid
!
!   Add new schemes here so it is recognized and proper checks are performed
!
 ltng_select: SELECT CASE(lightning_option)

    ! Convective resolved/permitted
    !CASE (ltng_crm_PR92w,ltng_crm_PR92z)
    !+++ kac +++
    CASE (ltng_crm_PR92w,ltng_crm_PR92z,ltng_crm_DP08_up,ltng_crm_P05_iwp,ltng_crm_D08_pim,ltng_crm_D08_imfp,&
    		ltng_crm_CSU_GEV,ltng_crm_CSU_35EV,ltng_crm_CSU_PIM, &
    		ltng_crm_AL_GEV10,ltng_crm_AL_GEV5,ltng_crm_AL_UP510, &
		ltng_crm_AL_supcell_GEV10,ltng_crm_AL_supcell_GEV5,ltng_crm_AL_supcell_UP510, &
		ltng_crm_F14_iceflux,ltng_crm_CSU_29May)
    !--- kac ---
        IF ( do_radar_ref .eq. 0 .or. mp_physics .eq. 0) THEN
          CALL wrf_error_fatal( ' lightning_init: Selected lightning option requires microphysics and do_radar_ref=1' )
        ENDIF

        WRITE(message, * ) ' lightning_init: CRM lightning option used: ', lightning_option
        CALL wrf_debug ( 100 , message )

    ! Convective parameterized
    CASE (ltng_cpm_PR92z)
        IF ( cu_physics .ne. GDSCHEME .and. cu_physics .ne. G3SCHEME  ) THEN
          CALL wrf_error_fatal( ' lightning_init: Selected lightning option requires GD or G3 convective parameterization' )
        ENDIF

        WRITE(message, * ) ' lightning_init: CPM lightning option selected: ', lightning_option
        CALL wrf_debug ( 100 , message )

    ! Non-existing options
    CASE DEFAULT
        CALL wrf_error_fatal ( ' lightning_init: invalid lightning_option')
 END SELECT ltng_select

!-- do not re-initialize for restarts
 IF (restart) return

!-- zero out arrays
 IF ( PRESENT( ic_flashcount ) .and. PRESENT( ic_flashrate ) .and. &
      PRESENT( cg_flashcount ) .and. PRESENT( cg_flashrate ) ) THEN
    CALL wrf_debug ( 100 , ' lightning_init: flash initializing lightning flash arrays' )

    ic_flashrate(:,:)  = 0.
    ic_flashcount(:,:) = 0.
    cg_flashrate(:,:)  = 0.
    cg_flashcount(:,:) = 0.
 ELSE
    CALL wrf_error_fatal ( ' lightning_init: flash arrays not present' )
 ENDIF

!-- Resolve auto-cellcount method option (cellcount_method=0)
 IF ( ( cellcount_method .eq. 0 ) .and. (lightning_option .eq. ltng_crm_PR92w )) THEN
   IF ( (ime-ims+1)*dx .gt. 1E4 ) THEN ! use patch only if path size > 10 km
     cellcount_method = 1
     WRITE(message, * ) ' lightning_init: setting auto cellcount_method to patch (cellcount_method=1'
   ELSE
     cellcount_method = 2
     WRITE(message, * ) ' lightning_init: setting auto cellcount_method to domain (cellcount_method=2'
   ENDIF
   CALL wrf_debug( 100, message )
 ENDIF

#ifdef WRF_CHEM

 CALL wrf_debug( 100, ' lightning_init: initializing and validating WRF-Chem only arrays and settings')

 IF ( lnox_opt .ne. lnox_opt_none .and. lightning_option .eq. 0 ) THEN
   CALL wrf_error_fatal ( ' lightning_init: cannot set LNOx without lightning_option')
 ENDIF

 IF ( lnox_opt .eq. lnox_opt_decaria .and. ( do_radar_ref .eq. 0 .or. mp_physics .eq. 0 ) ) THEN
   CALL wrf_error_fatal ( ' lightning_init: lnox_opt_decaria requires microphysics and do_radar_ref' )
 ENDIF
 
 ! +++KAC 10/17/15 +++
 IF ( lnox_opt .eq. lnox_opt_lma .and. ( do_radar_ref .eq. 0 .or. mp_physics .eq. 0 ) ) THEN
   CALL wrf_error_fatal ( ' lightning_init: lnox_opt_lma requires microphysics and do_radar_ref' )
 ENDIF
 ! ---KAC---

 IF (PRESENT( lnox_total )) lnox_total(:,:,:) = 0.
 IF (PRESENT( lnox_cg    )) lnox_cg(:,:,:)    = 0.
 IF (PRESENT( lnox_ic    )) lnox_ic(:,:,:)    = 0.

#endif

 CALL wrf_debug( 200, ' lightning_init: finishing')

 END SUBROUTINE lightning_init


!**********************************************************************
!
! SUBROUTINE lightning_driver
!
! Redirect to the appropriate lightning subroutine.
!
!**********************************************************************

 SUBROUTINE lightning_driver ( &
                          ! Frequently used prognostics
                            itimestep, dt, dx, dy,                &
                            xlat, xlon, xland, ht,                &
                            t_phy, p_phy, rho, u, v, w,           &
                            z, moist,                             &
                          ! Scheme specific prognostics
                            ktop_deep,                            &
                            refl,                                 &
                            current_time,                         &
                          ! Mandatory namelist inputs
                            lightning_option,                     &
                            lightning_dt,                         &
                            lightning_start_seconds,              &
                            flashrate_factor,                     &
                          ! IC:CG namelist settings
                            iccg_method,                          &
                            iccg_prescribed_num,                  &
                            iccg_prescribed_den,                  &
                          ! IC:CG inputs
                            iccg_in_num, iccg_in_den,             &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                            cldtop_adjustment,                    &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
			  !+++ kac
			  ! Morrison fallspeed for snow and graupel
			    vts_morr, vtg_morr,                   &
                          !--- kac
			  ! Mandatory outputs for all quantitative schemes
                            ic_flashcount, ic_flashrate,          &
                            cg_flashcount, cg_flashrate,           &
			  !+++kac
			  !Need number concentration
			    scalar                                &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description
 USE module_utility

! Model layer
 USE module_model_constants
 USE module_wrf_error

! Parameterization options
 USE module_ltng_crmpr92       ! lightning_option == 1,   ltng_crm_PR92w
                               ! lightning_option == 2,   ltng_crm_PR92z
			       ! lightning_option == 3,   ltng_crm_DP08_up
			       ! lightning_option == 4,   ltng_crm_P05_iwp
			       ! lightning_option == 5,   ltng_crm_D08_pim
			       ! lightning_option == 6,   ltng_crm_D08_imfp
			       
			       ! lightning_option == 20,  ltng_crm_CSU_GEV
			       ! lightning_option == 21,  ltng_crm_CSU_35EV
			       ! lightning_option == 22,  ltng_crm_CSU_PIM
			       
			       ! lightning_option == 30,  ltng_crm_AL_GEV10
			       ! lightning_option == 31,  ltng_crm_AL_GEV5
			       ! lightning_option == 32,  ltng_crm_AL_UP510
			       ! lightning_option == 33,  ltng_crm_AL_supcell_GEV10
			       ! lightning_option == 34,  ltng_crm_AL_supcell_GEV5
			       ! lightning_option == 35,  ltng_crm_AL_supcell_UP510
			       
			       ! lightning_option == 40,  ltng_crm_F14_iceflux
			       
			       ! lightning_option == 50,  ltng_crm_CSU_29May
			       
 USE module_ltng_cpmpr92z      ! lightning_option == 11,  ltng_cpm_PR92z

! IC:CG methods
 USE module_ltng_iccg

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
 INTEGER, INTENT(IN   )    ::       itimestep
 REAL,    INTENT(IN   )    ::       dt, dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: xlat, xlon, xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: t_phy, p_phy, rho
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: u, v, w, z
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme, num_moist), INTENT(IN   ) :: moist
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme, num_scalar),INTENT(IN   ) :: scalar
!---kac 

! Scheme specific prognostics
 INTEGER, DIMENSION( ims:ime,          jms:jme ),           INTENT(IN   ) :: ktop_deep     ! model LNB from cu_physics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: refl          ! reflectivity from mp_physics
 TYPE(WRFU_Time),                                           INTENT(IN   ) :: current_time  ! For use of IC:CG input

! Mandatory namelist inputs
 INTEGER, INTENT(IN   )    ::       lightning_option
 REAL,    INTENT(IN   )    ::       lightning_dt, lightning_start_seconds, flashrate_factor

! IC:CG namelist settings
 INTEGER, INTENT(IN   )    ::       iccg_method
 REAL,    INTENT(IN   )    ::       iccg_prescribed_num, iccg_prescribed_den
 REAL,    DIMENSION( ims:ime, jms:jme, 12), INTENT(IN   ) :: iccg_in_num, iccg_in_den

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method                    ! used in CRM
 REAL,    INTENT(IN   )    ::       cldtop_adjustment                   ! used in CPM

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte

!+++ kac
! Morrison fallspeed snow and graupel
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ),           INTENT(IN   ) :: vts_morr, vtg_morr  
!--- kac

! Mandatory outputs for all quantitative schemes
 REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(INOUT) :: ic_flashcount , cg_flashcount
 REAL, OPTIONAL, DIMENSION( ims:ime, jms:jme ), INTENT(  OUT) :: ic_flashrate  , cg_flashrate

! Local variables
 REAL, DIMENSION( ims:ime, jms:jme ) :: total_flashrate
 CHARACTER (LEN=80) :: message

 REAL, PARAMETER            :: reflthreshold = 20. ! reflectivity threshold for CRM schemes
 REAL, DIMENSION( kms:kme ) :: cellcount

!-----------------------------------------------------------------

 IF ( lightning_option .eq. 0 ) RETURN

 IF ( itimestep * dt .lt. lightning_start_seconds ) RETURN

 IF ( MOD((itimestep * dt - lightning_start_seconds), lightning_dt ) .ne. 0 ) RETURN

!+++ KAC Added 8/13/16
 WRITE(message,*)' KAC, itimestep, dt ' , itimestep, dt
 CALL wrf_debug ( 0,message )
!--- KAC 

!-----------------------------------------------------------------
! This driver performs several steps in order to produce lightning
! flash rate and flash count diagnostics:
!
! 1. Determine cloud extents for specific CRM schemes
! 2. Total flash rate assignment to 2D array
! 3. Partitioning of total lightning into IC & CG
! 4. Scale flash rate by flashrate_factor and lightning_dt
!
!-----------------------------------------------------------------

 IF ( lightning_option .eq. ltng_crm_PR92w .or. &
      lightning_option .eq. ltng_crm_PR92z .or. &
      
      !+++ kac +++
      
      lightning_option .eq. ltng_crm_DP08_up  .or. &
      lightning_option .eq. ltng_crm_P05_iwp  .or. &
      lightning_option .eq. ltng_crm_D08_pim  .or. &
      lightning_option .eq. ltng_crm_D08_imfp .or. &
      
      lightning_option .eq. ltng_crm_CSU_GEV  .or. &
      lightning_option .eq. ltng_crm_CSU_35EV .or. &
      lightning_option .eq. ltng_crm_CSU_PIM  .or. &
      
      lightning_option .eq. ltng_crm_AL_GEV10          .or. &
      lightning_option .eq. ltng_crm_AL_GEV5           .or. &
      lightning_option .eq. ltng_crm_AL_UP510          .or. &
      lightning_option .eq. ltng_crm_AL_supcell_GEV10  .or. &
      lightning_option .eq. ltng_crm_AL_supcell_GEV5   .or. &
      lightning_option .eq. ltng_crm_AL_supcell_UP510  .or. &
      
      lightning_option .eq. ltng_crm_F14_iceflux       .or. &
      
      lightning_option .eq. ltng_crm_CSU_29May) THEN
      
      !--- kac ---
      
   CALL wrf_debug ( 100, ' lightning_driver: determining cloud extents for CRM' )
   CALL countCells( &
          ! Inputs
            refl, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            its, ite, jts, jte, kts, kte,              &
          ! Outputs
            cellcount )
   WRITE(message, * ) ' lightning_driver: Max cell count = ', maxval(cellcount)
   CALL wrf_debug ( 100, message )
 ENDIF

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: calculating flash rate' )
 flashrate_select: SELECT CASE(lightning_option)

    ! CRM lightning options
    CASE( ltng_crm_PR92w )
        CALL wrf_debug ( 100, ' lightning_driver: calling Price and Rind 1992 (w_max, CRM)' )

        CALL ltng_crmpr92w ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )
    CASE( ltng_crm_PR92z )
        CALL wrf_debug ( 100, ' lightning_driver: calling Price and Rind 1992 (z_top, CRM)' )
        CALL ltng_crmpr92z ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )
		  
!+++ kac +++

    CASE( ltng_crm_DP08_up )
        CALL wrf_debug ( 100, ' lightning_driver: calling Deierling and Petersen 2008 (up_vol, CRM)' )
        CALL ltng_crm_up ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac	
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )	
		  
    CASE( ltng_crm_P05_iwp )
        CALL wrf_debug ( 100, ' lightning_driver: calling Petersen et al. 2005 (iwp, CRM)' )
        CALL ltng_crm_iwp ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    rho, 				  &
		    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qi=moist(ims,kms,jms,P_QI),		  &
		  !+++kac
		  !New variables
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )

    CASE( ltng_crm_D08_pim )
        CALL wrf_debug ( 100, ' lightning_driver: calling Deierling et al. 2008 (pim, CRM)' )
        CALL ltng_crm_pim ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    rho, 				  &
		    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
		  ! Fallspeed for snow
		    vts_morr,                             &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )

    CASE( ltng_crm_D08_imfp )
        CALL wrf_debug ( 100, ' lightning_driver: calling Deierling et al. 2008 (imfp, CRM)' )
        CALL ltng_crm_imfp ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    rho, 				  &
		    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
		  ! Fallspeed for snow & graupel
		    vts_morr, vtg_morr,                   &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qi=moist(ims,kms,jms,P_QI),		  &
		  !+++kac
		  !New variables
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_CSU_GEV )
        CALL wrf_debug ( 100, ' lightning_driver: calling CSU Basarb 2014 (GEV, CRM)' )
        CALL ltng_crm_csugev ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )

    CASE( ltng_crm_CSU_35EV )
        CALL wrf_debug ( 100, ' lightning_driver: calling CSU Basarb 2014 (35EV, CRM)' )
        CALL ltng_crm_csu35ev ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_CSU_PIM )
        CALL wrf_debug ( 100, ' lightning_driver: calling CSU Basarb 2014 (PIM, CRM)' )
        CALL ltng_crm_csupim ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )

    CASE( ltng_crm_AL_GEV10 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (GEV10, CRM)' )
        CALL ltng_crm_algev10 ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_AL_GEV5 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (GEV5, CRM)' )
        CALL ltng_crm_algev5 ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_AL_UP510 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (UP510, CRM)' )
        CALL ltng_crm_alup510 ( &
                  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac	
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )
		  
    CASE( ltng_crm_AL_supcell_GEV10 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (supcell_GEV10, CRM)' )
        CALL ltng_crm_alsupcellgev10 ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_AL_supcell_GEV5 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (supcell_GEV5, CRM)' )
        CALL ltng_crm_alsupcellgev5 ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		  !+++kac
		  !New variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )
		  
    CASE( ltng_crm_AL_supcell_UP510 )
        CALL wrf_debug ( 100, ' lightning_driver: calling AL Carey 2014 (supcell_UP510, CRM)' )
        CALL ltng_crm_alsupcellup510 ( &
                  ! Frequently used prognostics	
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )

    CASE( ltng_crm_F14_iceflux )
        CALL wrf_debug ( 100, ' lightning_driver: calling Finney 2014 (iceflux, CRM)' )
        CALL ltng_crm_iceflux ( &
                  ! Frequently used prognostics 
		  !+++kac
		    itimestep, dt,                        &
		  !---kac	
                    dx, dy, xland, ht, z, t_phy, p_phy,   &
                  ! Scheme specific prognostics
		    rho,				  &
                    w, refl, reflthreshold, cellcount,    &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Fallspeed for snow and graupel
                    vts_morr, vtg_morr,                   &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  ! Moisture variables
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qi=moist(ims,kms,jms,P_QI),		  &
		    qc=moist(ims,kms,jms,P_QC),	          &
		    qr=moist(ims,kms,jms,P_QR),           &
		  !+++kac
		  !New variables
		    qv=moist(ims,kms,jms,P_QV),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
                  )

    CASE( ltng_crm_CSU_29May )
        CALL wrf_debug ( 100, ' lightning_driver: calling CSU Davis 29May (35EV, CRM)' )
        CALL ltng_crm_csu29May ( &
		  ! Frequently used prognostics
		  !+++kac
		    itimestep, dt,                        &
		  !---kac
                    dx, dy, xland, ht, z, t_phy,          &
		  !+++kac
		    p_phy,                                &
		  !---kac
                  ! Scheme specific prognostics
                    refl, reflthreshold, cellcount,       &
                  ! Scheme specific namelist inputs
                    cellcount_method,                     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate,                      &
		  !+++kac
		  !New variables
		    qg=moist(ims,kms,jms,P_QG),		  &
		    qs=moist(ims,kms,jms,P_QS),		  &
		    qv=moist(ims,kms,jms,P_QV),           &
		    qr=moist(ims,kms,jms,P_QR),           &
		    nr=scalar(ims,kms,jms,P_QNR),         &
		    ns=scalar(ims,kms,jms,P_QNS),         &
		    ng=scalar(ims,kms,jms,P_QNG)          &
		  !---kac
		  )		  

!--- kac ---		  
		
		  
!   CASE ( another_crm_option)
!       CALL ...

    ! CPM lightning options
    CASE( ltng_cpm_PR92z )
        CALL wrf_debug ( 100, ' lightning_driver: calling Price and Rind 1992 (z_top, CPM)' )

        CALL ltng_cpmpr92z ( &
                  ! Frequently used prognostics
                    dx, dy, xland, ht, z, t_phy,      &
                    ktop_deep, cldtop_adjustment,     &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,     &
                    ims, ime, jms, jme, kms, kme,     &
                    its, ite, jts, jte, kts, kte,     &
                  ! Mandatory output for all quantitative schemes
                    total_flashrate                   &
                  )

!   CASE ( another_cpm_option)
!       CALL ...

    ! Invalid lightning options
    CASE DEFAULT
        WRITE(wrf_err_message, * ) ' lightning_driver: The lightning option does not exist: lightning_opt = ', lightning_option
        CALL wrf_error_fatal ( wrf_err_message )

 END SELECT flashrate_select

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: partitioning IC:CG')
 iccg_select: SELECT CASE(iccg_method)
    ! Flash rate option defaults
    CASE( 0 ) iccg_select
        CALL wrf_debug( 100, ' lightning_driver: using option-default IC:CG method' )
        iccg_method_default: SELECT CASE(lightning_option)

            CASE( ltng_crm_PR92w, ltng_crm_PR92z, ltng_cpm_PR92z ) iccg_method_default
                CALL iccg_boccippio( &
                            xlat, xlon,                                &
                            iccg_prescribed_num, iccg_prescribed_den,  &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,              &
                            ims, ime, jms, jme, kms, kme,              &
                            its, ite, jts, jte, kts, kte,              &
                          ! Input
                            total_flashrate,                           &
                          ! Output
                            ic_flashrate, cg_flashrate                 &
                          )

            CASE DEFAULT iccg_method_default
                CALL wrf_debug ( 100, ' lightning_driver: no method-default IC:CG implemented, using user-prescribed constant')
                CALL iccg_user_prescribed( &
                            iccg_prescribed_num,                  &
                            iccg_prescribed_den,                  &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            its, ite, jts, jte, kts, kte,         &
                          ! Input
                            total_flashrate,                      &
                          ! Output
                            ic_flashrate, cg_flashrate            &
                          )

        END SELECT iccg_method_default

    ! Used-prescribed constant
    CASE( 1 ) iccg_select
        WRITE(message, * ) ' lightning_driver: using user-prescribed IC:CG ratio = ', iccg_prescribed_num, iccg_prescribed_den
        CALL wrf_debug ( 100, message )
        CALL iccg_user_prescribed( &
                    iccg_prescribed_num,                  &
                    iccg_prescribed_den,                  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,         &
                    ims, ime, jms, jme, kms, kme,         &
                    its, ite, jts, jte, kts, kte,         &
                  ! Input
                    total_flashrate,                      &
                  ! Output
                    ic_flashrate, cg_flashrate            &
                  )

    ! Boccippio et al, 2001
    CASE( 2 ) iccg_select
        CALL wrf_debug ( 100, ' lightning_driver: using Boccippio 2001 IC:CG climatology')
        CALL iccg_boccippio( &
                    xlat, xlon,                                &
                    iccg_prescribed_num, iccg_prescribed_den,  &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                  )

    ! Price and Rind, 1993
    CASE( 3 ) iccg_select
        iccg_pr93_select: SELECT CASE(lightning_option)
        CASE( ltng_crm_PR92w, ltng_crm_PR92z ) iccg_pr93_select
            CALL wrf_debug ( 100, ' lightning_driver: using Price and Rind 1993 IC:CG ratio (CRM)')
            CALL iccg_crm_pr93( &
                    refl, reflthreshold, t_phy, z,             &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                )

        CASE DEFAULT iccg_pr93_select
            CALL wrf_debug ( 100, ' lightning_driver: using Price and Rind 1993 IC:CG ratio (CPM)')
            CALL iccg_pr93( &
                    ktop_deep, cldtop_adjustment, t_phy, z,    &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                )
        END SELECT iccg_pr93_select

    CASE( 4 ) iccg_select
        CALL wrf_debug ( 100, ' lightning_driver: using input IC:CG ratio from iccg_in_(num|den)' )
        CALL iccg_input( &
                    iccg_prescribed_num, iccg_prescribed_den,  &
                    iccg_in_num, iccg_in_den, current_time,    &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                  )

!+++ kac +++
    ! Custom IC:CG based on total LMA & NLDN CG data for a specific storm
    CASE( 5 ) iccg_select
        CALL wrf_debug ( 0, ' lightning_driver: using custom IC:CG from LMA and NLDN for specific storm')
        CALL iccg_custom( &
                    xlat, xlon,                                &
                    iccg_prescribed_num, iccg_prescribed_den,  &
                    itimestep, dt,                             &
                  ! Order dependent args for domain, mem, and tile dims
                    ids, ide, jds, jde, kds, kde,              &
                    ims, ime, jms, jme, kms, kme,              &
                    its, ite, jts, jte, kts, kte,              &
                  ! Input
                    total_flashrate,                           &
                  ! Output
                    ic_flashrate, cg_flashrate                 &
                  )
!--- kac ---


    ! Invalid IC:CG method
    CASE DEFAULT iccg_select
        WRITE(wrf_err_message, * ) ' lightning_driver: Invalid IC:CG method (iccg_method) = ', lightning_option
        CALL wrf_error_fatal ( wrf_err_message )

 END SELECT iccg_select

!-----------------------------------------------------------------

 CALL wrf_debug( 200, ' lightning_driver: converting flash rates to flash counts')

 ic_flashrate(its:ite,jts:jte) = ic_flashrate(its:ite,jts:jte) * flashrate_factor
 cg_flashrate(its:ite,jts:jte) = cg_flashrate(its:ite,jts:jte) * flashrate_factor

 ic_flashcount(its:ite,jts:jte) = ic_flashcount(its:ite,jts:jte) + ic_flashrate(its:ite,jts:jte) * lightning_dt
 cg_flashcount(its:ite,jts:jte) = cg_flashcount(its:ite,jts:jte) + cg_flashrate(its:ite,jts:jte) * lightning_dt

!-----------------------------------------------------------------

 CALL wrf_debug ( 100, ' lightning_driver: returning from')

 END SUBROUTINE lightning_driver


!**********************************************************************
!
! SUBROUTINE countCells
!
! For counting number of cells where reflectivity exceeds a certain
! threshold. Typically used by CRM schemes to redistribute lightning
! within convective cores.
!
! Departure from original implementation:
! Output includes domain-wide cellcounts if cellcount_method = 2
!
!**********************************************************************

 SUBROUTINE countCells( &
          ! Inputs
            refl, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            its, ite, jts, jte, kts, kte,              &
          ! Outputs
            cellcount )

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Inputs
 REAL,    DIMENSION( ims:ime,kms:kme,jms:jme ), INTENT(IN   ) :: refl
 REAL,    INTENT(IN   ) :: reflthreshold
 INTEGER, INTENT(IN   ) :: cellcount_method

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       its,ite, jts,jte, kts,kte


! Outputs
 REAL,    DIMENSION( kms:kme ), INTENT(  OUT) :: cellcount

! Local vars
 INTEGER :: i,k,j

!-----------------------------------------------------------------

 cellcount(kts:kte) = 0.
 DO j=jts,jte
   DO k=kts,kte
     DO i=its,ite
       IF ( refl(i,k,j) .gt. reflthreshold ) THEN
         cellcount(k) = cellcount(k) + 1
       ENDIF
     ENDDO
   ENDDO
 ENDDO

 IF ( cellcount_method .eq. 2 ) THEN
   DO k=kts,kte
     cellcount(k) = wrf_dm_sum_real(cellcount(k))
   ENDDO
 ENDIF

 END SUBROUTINE

 END MODULE module_lightning_driver
