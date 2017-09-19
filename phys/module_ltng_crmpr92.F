! WRF:MODEL_LAYER:PHYSICS
!
! Lightning flash rate prediction based on max vert. verlocity. Implemented
! for resolutions permitting resolved deep convection.
!
! Price, C., and D. Rind (1992), A Simple Lightning Parameterization for Calculating
!   Global Lightning Distributions, J. Geophys. Res., 97(D9), 9919–9933, doi:10.1029/92JD00719.
!
! Wong, J., M. Barth, and D. Noone (2012), Evaluating a Lightning Parameterization
!   at Resolutions with Partially-Resolved Convection, GMDD, in preparation.
!
! Unlike previous implementation, this version will produce slightly inconsistent
! IC and CG grid-flash rates against NO emission after production via calling
! lightning_nox_decaria.
!
! Contact: J. Wong <johnwong@ucar.edu>
!
!**********************************************************************

 MODULE module_ltng_crmpr92
 CONTAINS

 SUBROUTINE ltng_crmpr92w ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            w, refl, reflthreshold, cellcount,    &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_max_real
 
 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg,qs,qv,qr,nr,ns,ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ 
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: wmax            ! max w in patch or domain
 REAL :: total_fr,ave_fr ! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do
 
 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )
 
!---kac


 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN

! Compute flash rate across cell
 wmax = maxval(w(ips:ipe,kps:kpe,jps:jpe))
 IF ( cellcount_method .eq. 2 ) THEN
   wmax = wrf_dm_max_real(wmax)
 ENDIF

 total_fr = 5.0e-6 * wmax**4.5


! Write total_fr and associated wmax to rsl.out files
WRITE(message, * )' wmax, total_fr ' , wmax,total_fr
CALL wrf_debug ( 0,message )


! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO

 ! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN			!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE									!KAC Added
   ave_fr = 0.								!KAC Added
 ENDIF									!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crmpr92w

 SUBROUTINE ltng_crmpr92z ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_max_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg,qs,qv,qr,nr,ns,ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ 
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: zmax            ! max w in patch or domain
 REAL :: total_fr,ave_fr ! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount, count
 REAL :: maxcount, mostlyLand
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac


 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN

! Compute flash rate across cell
 k = kpe
 do while ( cellcount_kac(k) .eq. 0 .and. k .gt. kps)
   k = k-1
 ENDDO
 zmax = 0.
 mostlyland = 0.
 count = 0
 DO i=ips,ipe
   DO j=jps,jpe
     IF ( (refl10cm_kac(i,k,j) .gt. reflthreshold) .and. (t(i,k,j) .lt. 273.15) ) THEN
       IF (z(i,k,j)-ht(i,j) .gt. zmax) THEN
         zmax = z(i,k,j)-ht(i,j)
       ENDIF
       count = count + 1
       mostlyland = mostlyland + xland(i,j)
     ENDIF
   ENDDO
 ENDDO
 mostlyland = mostlyland/count

 zmax = zmax * 1.e-3
 WRITE(message, * ) ' ltng_crmpr92z: reflectivity cloud top height: ', zmax
 CALL wrf_debug ( 15, message )

 if ( cellcount_method .eq. 2 ) THEN
   zmax = wrf_dm_max_real(zmax)
 endif

 if ( mostlyLand .lt. 1.5 ) then
    total_fr = 3.44E-5 * (zmax**4.9)  ! PR 92 continental eq
 else
    total_fr = 6.57E-6 * (zmax**4.9)  ! Michalon 99 marine eq
 ENDIF
 
 
 ! Write total_fr and associated zmax to rsl.out files
WRITE(message, * )' zmax, total_fr ' , zmax,total_fr
CALL wrf_debug ( 0,message )
 

! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN			!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE									!KAC Added
   ave_fr = 0.								!KAC Added
 ENDIF	

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold  )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crmpr92z
 
 !+++ kac +++
 
 ! UPDRAFT VOLUME FRPS
 SUBROUTINE ltng_crm_up ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            w, refl, reflthreshold, cellcount,    &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg,qs,qv,qr,nr,ns,ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ 
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: up_5ms          	! upvol (m3)
 REAL :: total_fr,ave_fr 	! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 INTEGER :: cnt_grd		! grid cell count
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )
    
!---kac

 cnt_grd = 0								!KAC Added
 ave_fr  = 0.								!KAC Added
 total_fr= 0.								!KAC Added
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN

! Calculate dz here so that it can be used for updraft volume
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate across cell
 up_5ms = 0.
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
          if (t(i,k,j) .lt. 273.15-5.) then		!temp lt -5C=268.15K
	     if (w(i,k,j) .gt. 5.) then
	        up_5ms = up_5ms + dx * dy * dz(k)
		
		! Grid cell count meeting criteria			!KAC Added
		cnt_grd = cnt_grd + 1
		
	     endif
	  endif
       END DO
    END DO
 END DO


 IF ( cellcount_method .eq. 2 ) THEN
   up_5ms = wrf_dm_sum_real(up_5ms)
 ENDIF

 IF ( cnt_grd .gt. 0 ) THEN						!KAC Added
   total_fr = 6.75e-11 * up_5ms - 13.9
 ENDIF									!KAC Added
 
 ! KAC 8/8/16: Given the above equation, total_fr can be (-13.9) if there are no grid cells meeting the
 ! criteria for updraft vol, which means up_5ms=0 in above equation. If up_5ms=0, need total_fr=0. The
 ! total_fr can also be a negative value between -13.9 and 0 if the product of up_5ms and 6.75e-11 does
 ! not exceed (13.9), so need to ignore theses negative flashes as well.
  total_fr = max(total_fr,0.)
 
 
 ! Write total_fr and associated up_5ms to rsl.out files
WRITE(message, * )' up_5ms, total_fr ' , up_5ms,total_fr
CALL wrf_debug ( 0,message )
 

! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO

! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN			!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE									!KAC Added
   ave_fr = 0.								!KAC Added
 ENDIF									!KAC Added
 
 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 ! +++KAC
 WRITE(message, * )' up_5ms, ave_fr, maxcount ' , ave_fr, maxcount
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' up_5ms, cnt_grd ' , cnt_grd
 CALL wrf_debug ( 0,message )
! ---KAC 

 END SUBROUTINE ltng_crm_up

 
 !ICE WATER PATH
 SUBROUTINE ltng_crm_iwp ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            rho,				  &
			    refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qs, qg, qi,			          &
			  !+++kac
			  ! New variables
			    qv, qr, nr, ns, ng                    &
			  !---kac
			  )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    :: dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: rho
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qg, qi	!snow, grauple, ice mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qv, qr, nr, ns, ng
 REAL, DIMENSION( ims:ime, kms:kme, jms:jme )                   :: qs_k, qg_k, qi_k, qr_k, nr_k, ns_k, ng_k
 REAL, DIMENSION(          kps:kpe          )                   :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL, DIMENSION(          kps:kpe          )                   :: nr1d,ns1d,ng1d
 REAL, DIMENSION(          kps:kpe          )                   :: dBZ
!---kac 
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: qice							
 REAL :: iwp          						! ice water path (kg m-2)
 REAL :: total_fd,total_fr,ave_fr 				! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 INTEGER, DIMENSION( ips:ipe,jps:jpe ) :: flag
 REAL :: count
 CHARACTER (LEN=250) :: message
 REAL :: grid_res_P05, patch_area, grid_area			! variables for converting flash density
 !+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
 !---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
	   qi_k(i,:,j) = qi(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
	   qi_k(i,:,j) = qi(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
   do k=kps,kpe
   
   qv1d(k) = qv(i,k,j)
   qg1d(k) = qg_k(i,k,j)
   qr1d(k) = qr_k(i,k,j)
   qs1d(k) = qs_k(i,k,j)
   t1d(k)  = t(i,k,j)
   p1d(k)  = pres(i,k,j)
   nr1d(k) = nr(i,k,j)
   ns1d(k) = ns(i,k,j)
   ng1d(k) = ng(i,k,j)
   
   end do
   
   !Call lyy_lda
   call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
               kms,kme,dt,i,j)
   

	 
   !Call refl10cm_hm		 
   call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
               t1d,p1d,dBZ,kps,kpe,i,j)
	       
   do k = kps, kpe
      refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
   enddo	  

 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_fd = 0.							!KAC Added
 total_fr = 0.							!KAC Added
 ave_fr   = 0.							!KAC Added
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute variable masses and flash rate
 qice = 0.
 iwp = 0.
 flag = 0
 count = 0.


 DO i = ips, ipe
  DO j = jps, jpe
   DO k = kps, kpe
    
        ! Compute total ice mass
         !IF (t(i,k,j) .lt. 263.15 .and. refl10cm_kac(i,k,j) .gt. 18.) THEN
	 IF (t(i,k,j) .lt. 263.15 .and. refl10cm_kac(i,k,j) .gt. 25.) THEN
          qice = qice + (qs_k(i,k,j) + qi_k(i,k,j) + qg_k(i,k,j)) * rho(i,k,j) * dx * dy * dz(k)
	  flag(i,j) = flag(i,j) + 1
         ENDIF
 
   END DO
   
    IF (flag(i,j) .ge. 1) THEN
        count = count + 1.
    ENDIF
   
  END DO
 END DO
 
 ! Compute flash rate across cell
  IF (count .gt. 0.) THEN
      iwp = qice / (count * dx * dy)
  ELSE
      iwp = 0.
  ENDIF
 

 IF ( cellcount_method .eq. 2 ) THEN
   iwp = wrf_dm_sum_real(iwp)
 ENDIF

! Units are in flash density (flashes/km2*day)
 total_fd = 33.33 * iwp - 0.17
 
! Convert units to flashes/minute (flash density to flash rate)
 !grid_res_P05 = (50.*50.)						! km2; reps est of grid res used by P05 (0.5x0.5deg)
 patch_area = (dx*((ipe-ips)+1)) * (dy*((jpe-jps)+1)) * 1.e-6		! convert from m2 to km2; area of processor/tile
 
 total_fr = total_fd * patch_area / 1440.


! KAC Added 8/11/16: You don't want negative flashes
 total_fr = max(total_fr,0.)

 
! Write total_fr and associated iwp and qice to rsl.out files
WRITE(message, * )' IWP, qice, count ' , qice,count
CALL wrf_debug ( 0,message )

WRITE(message, * )' IWP, total_fd ' , iwp,total_fd
CALL wrf_debug ( 0,message )

WRITE(message, * )' IWP, total_fr ' , iwp,total_fr
CALL wrf_debug ( 0,message )


! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF


! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added this line
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added this line
   ave_fr = 0.							!KAC Added this line
 ENDIF								!KAC Added this line
 
 
 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crm_iwp
 
 
!PRECIPITATION ICE MASS
 SUBROUTINE ltng_crm_pim ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            rho,				  &
			    w, refl, reflthreshold, cellcount,    &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
			  ! Fallspeed for snow
			    vts_morr,                             &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qs, qg,				  &
			  !+++kac
			  ! New variables
			    qv, qr, nr, ns, ng                    &
			  !---kac
			  )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: rho
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qg	!snow, grauple mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Fallspeed for snow
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: vts_morr

! Local variables
 REAL :: qs_p_mass, qg_mass					! precip snow mass and graupel							
 REAL :: precip_ice_mass					! precip ice mass (kg)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_fr = 0.							!KAC Added
 ave_fr   = 0.							!KAC Added
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute variable masses and flash rate
 qs_p_mass = 0.
 qg_mass = 0.
 precip_ice_mass = 0.
 
 
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    
        ! Compute variable masses
         IF (t(i,k,j) .lt. 268.15) THEN

            ! Compute precipitation snow mass
	    ! +++kac - Changed snow threshold to Conrad's 0.00025 kg kg-1
            ! IF (qs(i,k,j) .gt. 1.0e-12) THEN
	     IF (qs_k(i,k,j) .gt. 0.00025) THEN
              IF (w(i,k,j) .lt. vts_morr(i,k,j)) THEN
	          qs_p_mass = qs_p_mass + qs_k(i,k,j) * rho(i,k,j) * dx * dy * dz(k)
              ENDIF
             ENDIF
	
	    ! Compute graupel mass
	    ! +++kac - changed qg threshold from 1.e-12 to Conrad's 0.0005 kg kg-1
            ! IF (qg(i,k,j) .gt. 1.0e-12) THEN
	     IF (qg_k(i,k,j) .gt. 0.00050) THEN
	         qg_mass = qg_mass + qg_k(i,k,j) * rho(i,k,j) * dx * dy * dz(k)
             ENDIF
	     
	 ENDIF

   END DO
  END DO
 END DO
 
 ! Compute flash rate across cell
         IF (qg_mass .gt. 0. .or. qs_p_mass .gt. 0.) THEN
              precip_ice_mass = precip_ice_mass + qs_p_mass + qg_mass
         ENDIF


IF ( cellcount_method .eq. 2 ) THEN
   precip_ice_mass = wrf_dm_sum_real(precip_ice_mass)
 ENDIF

 total_fr = precip_ice_mass * 3.4e-8 - 18.1

! KAC Added 8/11/16: You don't want negative flashes. You can have them if PIM=0, or if the
! product of PIM and 3.4e-8 is not larger than 18.1
 total_fr = max(total_fr,0.)


! Write total_fr and associated precip_ice_mass to rsl.out files
WRITE(message, * )' precip_ice_mass, total_fr ' , precip_ice_mass,total_fr
CALL wrf_debug ( 0,message )


! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added
 
 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE


! Write variable values for validation
!WRITE(message, * )' Max vts_morr, qs_p_mass, qg_mass, precip_ice_mass ' , maxval(vts_morr),qs_p_mass,qg_mass,precip_ice_mass
!CALL wrf_debug ( 0,message )

!WRITE(message, * )' total_fr, k_maxcount, maxcount, ave_fr ' , total_fr,k_maxcount,maxcount,ave_fr
!CALL wrf_debug ( 0,message )

 END SUBROUTINE ltng_crm_pim


!ICE MASS FLUX PRODUCT
 SUBROUTINE ltng_crm_imfp ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            rho,				  &
			    w, refl, reflthreshold, cellcount,    &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
			  ! Fallspeed for snow & graupel
			    vts_morr, vtg_morr,                   &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qs, qg, qi,				  &
			  !+++kac
			    qv, qr, nr, ns, ng                    &
			  !---kac
			  )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: rho			!kg/m3
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qg, qi		!snow, graupel, ice mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qi_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Fallspeed for snow and graupel
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: vts_morr, vtg_morr

! Local variables
 REAL :: wqi, wqs_np, wqs_p, wqg				
 REAL :: flux_prod          					! flux product (kg2 m s-2)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
	   qi_k(i,:,j) = qi(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
	   qi_k(i,:,j) = qi(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

total_fr = 0.							!KAC Added
ave_fr   = 0.							!KAC Added
total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute variable masses and flash rate
 wqi = 0.
 wqs_np = 0.
 wqs_p = 0.
 wqg = 0.
 flux_prod = 0.
 
 
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    
    ! IMFP is defined as the product of the precip and the non-precip ice mass fluxes for T < -5C
     IF (t(i,k,j) .lt. 268.15) THEN
    
        ! Compute ice crystal mass (units: kg/s)
         IF (w(i,k,j) .gt. 0.) THEN
	          wqi = wqi + w(i,k,j) * qi_k(i,k,j) * rho(i,k,j) * dx * dy
         ENDIF
	
	! Compute variable mass fluxes needed for FRPS
	 
	 	! Compute precipitation (units: kg m/s) and non-precipitation (units: kg/s) snow mass fluxes
		! +++kac - Changed snow threshold from 1.0e-12 to Conrad's 0.00025 kg kg-1
		 IF (qs_k(i,k,j) .gt. 0.00025) THEN
	      	  IF (w(i,k,j) .gt. vts_morr(i,k,j)) THEN
		      wqs_np = wqs_np + (w(i,k,j) - vts_morr(i,k,j)) * qs_k(i,k,j) * rho(i,k,j) * dx * dy
		  ENDIF
		  IF (w(i,k,j) .lt. vts_morr(i,k,j)) THEN
		      wqs_p = wqs_p + (w(i,k,j) - vts_morr(i,k,j)) * qs_k(i,k,j) * rho(i,k,j) * dx * dy * dz(k)
		  ENDIF
         	 ENDIF
		 
		 ! Compute graupel mass flux (units: kg m/s)
		 ! +++kac - Changed graupel threshold from 1.e-12 to Conrad's 0.00050 kg kg-1
		 IF (qg_k(i,k,j) .gt. 0.00050) THEN
		     wqg = wqg + (w(i,k,j) - vtg_morr(i,k,j)) * qg_k(i,k,j) * rho(i,k,j) * dx * dy * dz(k)
		 ENDIF
 
    ENDIF
 
   END DO
  END DO
 END DO
 
 ! Compute flash rate across cell
  IF (wqg .lt. 0.) THEN
   IF (wqi .gt. 0.) THEN
       flux_prod = flux_prod - (wqi + wqs_np) * (wqs_p + wqg)
   ENDIF
  ENDIF
 

 IF ( cellcount_method .eq. 2 ) THEN
   flux_prod = wrf_dm_sum_real(flux_prod)
 ENDIF

! KAC Added 8/11/16: You can't have flashes when flux_prod=0
 IF ( flux_prod .gt. 0. ) THEN
    total_fr = flux_prod * 9.0e-15 + 13.4
 ELSE
    total_fr = 0.
 ENDIF

! Write total_fr and associated flux_prod to rsl.out files
WRITE(message, * )' flux_prod, total_fr ' , flux_prod,total_fr
CALL wrf_debug ( 0,message )


! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added
 
 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crm_imfp

 
!CSU GRAUPEL ECHO VOLUME
 SUBROUTINE ltng_crm_csugev ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qg,					  &
			  !+++kac
			  ! New variables
			    qs,                                   &
			    qv,qr,nr,ns,ng                        &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg		!grauple mixing ratios
!++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qg_k,qs_k,qi_k,qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: graup_echo    						! graupel echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )
    
!---kac
 
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 graup_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp lt -40C=233.15K, refl gt 35dBZ; Convert m3 to km3
	  ! +++kac - changed qg from (gt 0kg/kg) to (ge 0.0005kg/kg)
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. refl10cm_kac(i,k,j) .gt. 35. .and. qg_k(i,k,j) .gt. 0.0005) then
	     graup_echo = graup_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   graup_echo = wrf_dm_sum_real(graup_echo)
 ENDIF

 total_fr = 7.0e-2 * graup_echo


 ! Write total_fr and associated gev to rsl.out files
WRITE(message, * )' gev, total_fr ' , graup_echo,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_csugev
 
 
 !CSU 35-dBZ ECHO VOLUME
 SUBROUTINE ltng_crm_csu35ev ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg,qs,qv,qr,nr,ns,ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ 
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: tf_echo    						! 35-dBZ echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )
    
!---kac
 
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 tf_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp lt -40C=233.15K; Convert m3 to km3
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. refl10cm_kac(i,k,j) .gt. 35.) then
	     tf_echo = tf_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   tf_echo = wrf_dm_sum_real(tf_echo)
 ENDIF

 total_fr = 7.2e-2 * tf_echo


 ! Write total_fr and associated 35ev to rsl.out files
WRITE(message, * )' 35ev, total_fr ' , tf_echo,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_csu35ev
 
 
!CSU PRECIPITATION ICE MASS
 SUBROUTINE ltng_crm_csupim ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg, qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate	!flashes/sec

! Local variables
 REAL :: dbz_to_z, hail_mass, graup_mass, total_mass		! refl (mm^6/m^3), hail, graupel and total mass (g/m^3)
 REAL :: precip_ice_mass					! precipitation ice mass (kg)
 REAL :: total_fr,ave_fr 					! cloud flash rate (flashes/min)
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 precip_ice_mass = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp lt -40C=233.15K, refl gt 35dBZ
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. refl10cm_kac(i,k,j) .gt. 35. .and. qg_k(i,k,j) .gt. 0.0005) then
	  
	     ! Convert reflectivity from dBZ to mm^6 / m^3
	     dbz_to_z = 10.**(refl10cm_kac(i,k,j) / 10.)
	     
	     ! Convert (mm^6/m^3) to g/m^3
	     ! This simulation uses Morrison 2-moment MP and selects graupel over hail, so ignoring hail here.
	     !hail_mass  = (4.4 * (10.**(-5.))) * dbz_to_z**0.71
	     graup_mass = (5.2 * (10.**(-3.))) * dbz_to_z**0.5
	     !total_mass = hail_mass + graup_mass
	     total_mass = graup_mass
	     
	     ! Multiply by grid volume and convert from g to kg	     
	     precip_ice_mass = precip_ice_mass + (total_mass * dx * dy * dz(k) * 1.e-3)
	     
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   precip_ice_mass = wrf_dm_sum_real(precip_ice_mass)
 ENDIF

 total_fr = 1.2e-7 * precip_ice_mass


 ! Write total_fr and associated pim to rsl.out files
WRITE(message, * )' pim, total_fr ' , precip_ice_mass,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_csupim
 
 
! AL GRAUPEL ECHO VOLUME (T=-10C)
SUBROUTINE ltng_crm_algev10 ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qg,					  &
			  !+++kac
			  ! New variables
			    qs,                                   &
			    qv,qr,nr,ns,ng                        &
			  !---kac
			   )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg		!grauple mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: graup_echo    						! graupel echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac
 
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 graup_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -10C=263.15K, temp gt -40C=233.15K; Convert m3 to km3
	  ! +++kac - changed qg from (gt 0kg/kg) to (ge 0.0005kg/kg)
          if (t(i,k,j) .lt. 263.15 .and. t(i,k,j) .gt. 233.15 .and. qg_k(i,k,j) .gt. 0.0005) then
	     graup_echo = graup_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   graup_echo = wrf_dm_sum_real(graup_echo)
 ENDIF

 total_fr = 0.0534 * graup_echo


 ! Write total_fr and associated gev to rsl.out files
WRITE(message, * )' gev, total_fr ' , graup_echo,total_fr
CALL wrf_debug ( 0,message )

    
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_algev10
 
 
! AL GRAUPEL ECHO VOLUME (T=-5C)
SUBROUTINE ltng_crm_algev5 ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qg,					  &
			  !+++kac
			  ! New variables
			    qs,                                   &
			    qv,qr,nr,ns,ng                        &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg		!grauple mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: graup_echo    						! graupel echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac
 
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 graup_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp gt -40C=233.15K; Convert m3 to km3
	  ! +++kac - changed qg from (gt 0kg/kg) to (ge 0.0005kg/kg)
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. qg_k(i,k,j) .gt. 0.0005) then
	     graup_echo = graup_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   graup_echo = wrf_dm_sum_real(graup_echo)
 ENDIF

 total_fr = 0.0430 * graup_echo


 ! Write total_fr and associated gev to rsl.out files
WRITE(message, * )' gev, total_fr ' , graup_echo,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_algev5
 
 
! AL UPDRAFT ECHO VOLUME (W>5m/s, T==10C)
 SUBROUTINE ltng_crm_alup510 ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            w, refl, reflthreshold, cellcount,    &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg, qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: up510          	! upvol (km3)
 REAL :: total_fr,ave_fr 	! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN

! Calculate dz here so that it can be used for updraft volume
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate across cell
 up510 = 0.
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -10C=263.15K, temp gt -40C=233.15K, updraft vel gt 5m/s
          if (t(i,k,j) .lt. 263.15 .and. t(i,k,j) .gt. 233.15 .and. w(i,k,j) .gt. 5.) then
	     
	        ! Convert m3 to km3
	        up510 = up510 + dx * dy * dz(k) * 1.e-9
	     
	  endif
	  
       END DO
    END DO
 END DO


 IF ( cellcount_method .eq. 2 ) THEN
   up510 = wrf_dm_sum_real(up510)
 ENDIF

 total_fr = 0.0337 * up510
 
 
 ! Write total_fr and associated up510 to rsl.out files
WRITE(message, * )' up510, total_fr ' , up510,total_fr
CALL wrf_debug ( 0,message )
 

! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crm_alup510
 
 
! AL SUPERCELL GRAUPEL ECHO VOLUME (T=-10C)
SUBROUTINE ltng_crm_alsupcellgev10 ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qg,					  &
			  !+++kac
			  ! New variables
			    qs,qv,qr,nr,ns,ng                     &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg		!grauple mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: graup_echo    						! graupel echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 graup_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -10C=263.15K, temp gt -40C=233.15K; Convert m3 to km3
	  ! +++kac - changed qg from (gt 0kg/kg) to (ge 0.0005kg/kg)
          if (t(i,k,j) .lt. 263.15 .and. t(i,k,j) .gt. 233.15 .and. qg_k(i,k,j) .gt. 0.0005) then
	     graup_echo = graup_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   graup_echo = wrf_dm_sum_real(graup_echo)
 ENDIF

 total_fr = 0.0674 * graup_echo


 ! Write total_fr and associated gev to rsl.out files
WRITE(message, * )' gev, total_fr ' , graup_echo,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_alsupcellgev10
 
 
! AL SUPERCELL GRAUPEL ECHO VOLUME (T=-5C)
SUBROUTINE ltng_crm_alsupcellgev5 ( &
 			  ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qg,					  &
			  !+++kac
			  ! New variables
			    qs,qv,qr,nr,ns,ng                     &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg		!grauple mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: graup_echo    						! graupel echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 graup_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp gt -40C=233.15K; Convert m3 to km3
	  ! +++kac - changed qg threshold from (gt 0kg/kg) to (ge 0.0005kg/kg)
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. qg_k(i,k,j) .gt. 0.0005) then
	     graup_echo = graup_echo + dx * dy * dz(k) * 1.e-9
	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   graup_echo = wrf_dm_sum_real(graup_echo)
 ENDIF

 total_fr = 0.0569 * graup_echo


 ! Write total_fr and associated gev to rsl.out files
WRITE(message, * )' gev, total_fr ' , graup_echo,total_fr
CALL wrf_debug ( 0,message )


 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE
 
 END SUBROUTINE ltng_crm_alsupcellgev5


! AL SUPERCELL UPDRAFT ECHO VOLUME (W>5m/s, T==10C)
 SUBROUTINE ltng_crm_alsupcellup510 ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            w, refl, reflthreshold, cellcount,    &
                          ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
                          )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg, qs, qv, qr, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: up510          	! upvol (km3)
 REAL :: total_fr,ave_fr 	! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN

! Calculate dz here so that it can be used for updraft volume
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate across cell
 up510 = 0.
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -10C=263.15K, temp gt -40C=233.15K, updraft vel gt 5m/s
          if (t(i,k,j) .lt. 263.15 .and. t(i,k,j) .gt. 233.15 .and. w(i,k,j) .gt. 5.) then
	     
	        ! Convert m3 to km3
	        up510 = up510 + dx * dy * dz(k) * 1.e-9
		
	  endif
	  
       END DO
    END DO
 END DO


 IF ( cellcount_method .eq. 2 ) THEN
   up510 = wrf_dm_sum_real(up510)
 ENDIF

 total_fr = 0.0403 * up510
 
 
 ! Write total_fr and associated up510 to rsl.out files
WRITE(message, * )' up510, total_fr ' , up510,total_fr
CALL wrf_debug ( 0,message )
 

! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added
 
 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 END SUBROUTINE ltng_crm_alsupcellup510


! FINNEY UPWARD CLOUD ICE FLUX
 SUBROUTINE ltng_crm_iceflux ( &
                          ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t, pres,        &
                          ! Scheme specific prognostics
                            rho,				  &
			    w, refl, reflthreshold, cellcount,    &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
			  ! Fallspeed for snow and graupel
			    vts_morr, vtg_morr,			  &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  ! Moisture variables
			    qs, qg, qi, qc, qr,			  &
			  !+++kac
			  ! New variables
			    qv,nr,ns,ng                           &
			  !---kac
			  )
!-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t, pres		!pressure (Pa)

! Scheme specific prognostics
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: w
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: rho			!kg/m3
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qs, qg, qi, qc, qr	!snow, graupel, ice, cloud water, rain water mixing ratios
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qv, nr, ns, ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qi_k, qr_k, nr_k, ns_k, ng_k, qc_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ      
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Fallspeed for snow and graupel
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: vts_morr, vtg_morr

! Local variables
 REAL, DIMENSION (kps:kpe) :: tmp_plevel
 REAL :: plevel
 INTEGER :: k_plevel
 REAL :: frac_cld						! fractional cloud cover at 44hPa (m2 [cloud] m-2 [cell])			
 REAL :: up_cld_iceflux          				! upward cloud ice flux at 440hPa (kg [ice] m-2 [cloud] s-1)
 REAL :: total_fr						! total flash (flash s-1)
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac
 
 REAL, PARAMETER :: cldfra_thresh = 1.e-5			! This threshold is equivalent to 0.01 g/kg and is used with qtot

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
	   qi_k(i,:,j) = qi(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
	   qi_k(i,:,j) = qi(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

!---kac

 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO

 
! Loop through each grid cell within processor in the horizontal
 DO j = jps, jpe
  DO i = ips, ipe
  
  ! Compute variables
   tmp_plevel( kps:kpe ) = 0.
   frac_cld = 0.
   up_cld_iceflux = 0.
   total_fr = 0.
  
  
  ! Loop through each vertical level of the grid to find the model level that is closest to 440 hPa (44000 Pa) at the grid
   DO k = kps, kpe
    !tmp_plevel(k) = abs(44000. - pres(i,k,j))
    tmp_plevel(k) = abs(39000. - pres(i,k,j))
   ENDDO
   
   k_plevel = kps
   plevel = tmp_plevel(kps)
   DO k=kps+1,kpe
     IF ( tmp_plevel(k) .lt. plevel ) THEN
       k_plevel = k
       plevel = tmp_plevel(k)
     ENDIF
   ENDDO

  
  ! Calculate the fractional cloud cover (Qtot) 				! units: m2 [cloud] m-2 [cell])
   frac_cld = qc_k(i,k_plevel,j) + qi_k(i,k_plevel,j) + qs_k(i,k_plevel,j) + qg_k(i,k_plevel,j) + qr_k(i,k_plevel,j)

        ! CLOUD IS PRESENT
         IF (frac_cld .gt. cldfra_thresh) THEN
   
                ! Calculate the specific cloud ice water content at 440 hPa		! units: kg [ice] kg-2 [air]
		
		! Non-precipitating snow (qs threshold 0.25g/kg came from Conrad)
		 IF ( w(i,k_plevel,j) .gt. vts_morr(i,k_plevel,j) .and. qs_k(i,k_plevel,j) .ge. 0.00025 ) THEN
		   up_cld_iceflux = up_cld_iceflux + qs_k(i,k_plevel,j) * rho(i,k_plevel,j) * ( w(i,k_plevel,j) - vts_morr(i,k_plevel,j) )
		 ENDIF
		 
		! Non-precipitating ice
		 IF ( w(i,k_plevel,j) .gt. 0. .and. qi_k(i,k_plevel,j) .gt. 0. ) THEN
   	   	   up_cld_iceflux = up_cld_iceflux + qi_k(i,k_plevel,j) * rho(i,k_plevel,j) * w(i,k_plevel,j)
		 ENDIF
	
	! NO CLOUD AT 440hPa (No cloud=no ice flux)
  	 ELSE IF (frac_cld .le. cldfra_thresh) THEN
	 
	  	! Calculate the upward cloud ice flux				! units: kg [ice] m-2 [cloud] s-1
	   	 up_cld_iceflux = 0.
	   
   	 ENDIF
   
        ! Calculate the flash density over land (units: fl m-2 [cell] s-1) & convert to flash rate with (dx * dy) (units: fl s-1)
         total_fr = 6.58e-7 * up_cld_iceflux * dx * dy
  

  	! Restrict total flashes to within refl > reflthreshold
   	 IF ( refl10cm_kac(i,k_maxcount,j) .gt. reflthreshold ) THEN
     	   total_flashrate(i,j) = total_fr
   	 ENDIF
	
  ENDDO
 ENDDO

  ! Write total_fr and associated up510 to rsl.out files
  WRITE(message, * )' iceflux, maxval(total_flashrate) ' , maxval(total_flashrate(ips:ipe,jps:jpe))
  CALL wrf_debug ( 0,message )
 
 END SUBROUTINE ltng_crm_iceflux


 !CSU 35-dBZ ECHO VOLUME tailored for 29 May OK storm region
 SUBROUTINE ltng_crm_csu29May ( &
 ! Frequently used prognostics
			  !+++kac
			    itimestep, dt,                        &
			  !---kac
                            dx, dy, xland, ht, z, t,              &
			  !+++kac
			    pres,                                 &
			  !---kac
                          ! Scheme specific prognostics
                            refl, reflthreshold, cellcount,       &
			  ! Scheme specific namelist inputs
                            cellcount_method,                     &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,         &
                            ims, ime, jms, jme, kms, kme,         &
                            ips, ipe, jps, jpe, kps, kpe,         &
                          ! Mandatory output for all quantitative schemes
                            total_flashrate,                      &
			  !+++kac
			  ! New variables
			    qg,qs,qv,qr,nr,ns,ng                  &
			  !---kac
			  )
 
 !-----------------------------------------------------------------
! Framework
 USE module_state_description

! Model layer
 USE module_model_constants
 USE module_wrf_error

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Frequently used prognostics
!+++kac
 INTEGER, INTENT(IN   )    :: itimestep
 REAL,    INTENT(IN   )    :: dt
!---kac
 REAL,    INTENT(IN   )    ::       dx, dy

 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: xland, ht
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: z, t
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: pres
!---kac

! Scheme specific prognostics
!+++kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: qg,qs,qv,qr,nr,ns,ng
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme )                :: qs_k, qg_k, qr_k, nr_k, ns_k, ng_k
 REAL,    DIMENSION(          kps:kpe          )                :: qv1d,qg1d,qr1d,qs1d,t1d,p1d
 REAL,    DIMENSION(          kps:kpe          )                :: nr1d,ns1d,ng1d
 REAL,    DIMENSION(          kps:kpe          )                :: dBZ 
!---kac
 REAL,    DIMENSION( ims:ime, kms:kme, jms:jme ), INTENT(IN   ) :: refl
 REAL,                                            INTENT(IN   ) :: reflthreshold
 REAL,    DIMENSION(          kms:kme          ), INTENT(IN   ) :: cellcount

! Scheme specific namelist inputs
 INTEGER, INTENT(IN   )    ::       cellcount_method

! Order dependent args for domain, mem, and tile (patch) dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Mandatory outputs for all quantitative schemes
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: total_flashrate

! Local variables
 REAL :: tf_echo    						! 35-dBZ echo vol (km3)
 REAL :: total_fr,ave_fr 					! cloud flash rate
 INTEGER :: i,k,j
 INTEGER :: k_maxcount
 INTEGER :: cnt_grd						! grid cell count
 REAL :: maxcount
 REAL, DIMENSION (kps:kpe) :: sumz, countz, avgz, dz
 CHARACTER (LEN=250) :: message
!+++kac
 REAL, DIMENSION(ims:ime,kms:kme,jms:jme) :: refl10cm_kac
 REAL, DIMENSION(        kms:kme        ) :: cellcount_kac
!---kac

!-----------------------------------------------------------------

!+++kac
 do i=ips,ipe
 do j=jps,jpe
 
! Adjust the concentrataions of graupel, snow and ice to match "obs"

      !Graupel is ~55% too large above 9km, so reduce. Below 9km, leave as is.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 9000.)
           qg_k(i,:,j) = qg(i,:,j) * 0.65
           ng_k(i,:,j) = ng(i,:,j) * 0.65
      ELSEWHERE
           qg_k(i,:,j) = qg(i,:,j)
           ng_k(i,:,j) = ng(i,:,j)
      END WHERE

      !Snow and ice are ~42% too small between 11-13km, so increase. Below 11km
      !snow and ice are ~31% too large, so reduce. Applying same above 13km.
      WHERE ((z(i,:,j)-ht(i,j)) .gt. 11000. .and. (z(i,:,j)-ht(i,j)) .le. 13000.)
           qs_k(i,:,j) = qs(i,:,j) * 1.72
           ns_k(i,:,j) = ns(i,:,j) * 1.72
      ELSEWHERE
           qs_k(i,:,j) = qs(i,:,j) * 0.77
           ns_k(i,:,j) = ns(i,:,j) * 0.77
      END WHERE
      
      !Rain is ~70% too large, so reduce below 7.2km
      WHERE ((z(i,:,j)-ht(i,j)) .le. 7200.)
           qr_k(i,:,j) = qr(i,:,j) *0.6
           nr_k(i,:,j) = nr(i,:,j) *0.6
      ELSEWHERE
           qr_k(i,:,j) = qr(i,:,j)
           nr_k(i,:,j) = nr(i,:,j)
      END WHERE
      
      
! do i=ips,ipe
! do j=jps,jpe
 
    do k=kps,kpe
    
    qv1d(k) = qv(i,k,j)
    qg1d(k) = qg_k(i,k,j)
    qr1d(k) = qr_k(i,k,j)
    qs1d(k) = qs_k(i,k,j)
    t1d(k)  = t(i,k,j)
    p1d(k)  = pres(i,k,j)
    nr1d(k) = nr_k(i,k,j)
    ns1d(k) = ns_k(i,k,j)
    ng1d(k) = ng_k(i,k,j)
    
    end do
    
    !Call lyy_lda
    call lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                kms,kme,dt,i,j)

    !Call refl10cm_hm
    call refl_kac(qv1d,qr1d,nr1d,qs1d,ns1d,qg1d,ng1d,&
                 t1d,p1d,dBZ,kps,kpe,i,j)
    
    do k = kps,kpe
        refl10cm_kac(i,k,j) = MAX(-35., dBZ(k))
    enddo
 
 end do
 end do

 ! Determine new cellcount (cellcount_kac) based on impact of scaled hydrometeors on refl
 CALL countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )
    
!---kac

 cnt_grd = 0
 ave_fr = 0.
 total_fr = 0. 
 total_flashrate( ips:ipe,jps:jpe ) = 0.

 IF ( maxval(cellcount_kac(kps:kpe)) .eq. 0 ) RETURN


! Calculate dz here so that it can be used for calculating mass
 sumz(kps:kpe)=0
 countz(kps:kpe)=0
 DO j = jps, jpe
  DO k = kps, kpe
   DO i = ips, ipe
    sumz(k) = sumz(k) + (z(i,k,j) - ht(i,j))		!units: m
    countz(k) = countz(k) + 1
   END DO
  END DO
 END DO
 
 avgz(kps:kpe) = sumz(kps:kpe)/countz(kps:kpe)
 
 dz(kps) = 0.5*(avgz(kps+1) - avgz(kps))
 dz(kpe) = 0.5*(avgz(kpe) - avgz(kpe-1))
 
 DO k = kps+1, kpe-1
    dz(k) = 0.5*(avgz(k+1) + avgz(k)) - 0.5*(avgz(k) + avgz(k-1))
 END DO

! Compute flash rate
 tf_echo = 0.
  
 DO j = jps,jpe
    DO k = kps,kpe
       DO i = ips,ipe
       
          !temp lt -5C=268.15K, temp lt -40C=233.15K; Convert m3 to km3
          if (t(i,k,j) .lt. 268.15 .and. t(i,k,j) .gt. 233.15 .and. refl10cm_kac(i,k,j) .gt. 35.) then
	     tf_echo = tf_echo + dx * dy * dz(k) * 1.e-9

	     ! Grid cell count meeting criteria
	     cnt_grd = cnt_grd + 1	     

	  endif
	  
       END DO
    END DO
 END DO
 
 IF ( cellcount_method .eq. 2 ) THEN
   tf_echo = wrf_dm_sum_real(tf_echo)
 ENDIF


! KAC Added 8/11/16: You can't have flashes when tf_echo=0
 IF ( tf_echo .gt. 0. ) THEN
    total_fr = (-1.5242e-6)*tf_echo*tf_echo + (0.0328*tf_echo) + 12.977
 ELSE
    total_fr = 0.
 ENDIF


 ! Write total_fr and associated 35ev to rsl.out files
WRITE(message, * )' 35ev, total_fr ' , tf_echo,total_fr
CALL wrf_debug ( 0,message )

  
 ! Locating widest part of convective core
 k_maxcount = kps
 maxcount = cellcount_kac(kps)
 DO k=kps+1,kpe
   IF ( cellcount_kac(k) .gt. maxcount ) THEN
     k_maxcount = k
     maxcount = cellcount_kac(k)
   ENDIF
 ENDDO


! KAC Added 8/11/16 
 IF (maxcount .eq. 0.) THEN
   WRITE(message, * )' BEWARE maxcount equals zero '
   CALL wrf_debug ( 0,message )
 ENDIF

! Distributing across convective core
 IF ( total_fr .gt. 0. .and. maxcount .gt. 0.) THEN		!KAC Added
   ave_fr = total_fr/maxcount/60.
 ELSE								!KAC Added
   ave_fr = 0.							!KAC Added
 ENDIF								!KAC Added

 WHERE( refl10cm_kac(ips:ipe,k_maxcount,jps:jpe) .gt. reflthreshold )
   total_flashrate(ips:ipe,jps:jpe) = ave_fr
 ENDWHERE

 ! +++KAC
 WRITE(message, * )' 35ev, ave_fr, maxcount ' , ave_fr, maxcount
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, kps,kpe,kms,kme ' , kps,kpe,kms,kme
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, avgz(0:2) ' , avgz(0), avgz(1), avgz(2)
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, avgz(87:90) ' , avgz(87), avgz(88), avgz(89), avgz(90)
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, dz(0:1), dz(30) ' , dz(0), dz(1), dz(30)
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, dz(87:88), dz(89) ' , dz(87), dz(88), dz(89)
 CALL wrf_debug ( 0,message )
 
 WRITE(message, * )' 35ev, cnt_grd ' , cnt_grd
 CALL wrf_debug ( 0,message )
 ! ---KAC
 
 END SUBROUTINE ltng_crm_csu29May 
 
 !--- kac ---
 
 
 !+++ kac +++
 SUBROUTINE lda_kac(qv1d,qg1d,kps,kpe,itimestep,t1d,p1d,&
                          kms,kme,dt,i,j)
 
 ! Framework
 USE module_wrf_error
 USE module_utility
 USE module_domain
 
 IMPLICIT NONE
 INTEGER, INTENT( IN) :: KPS,KPE
 INTEGER, INTENT( IN) :: ITIMESTEP,I,J
 REAL, INTENT(IN) :: dt
 INTEGER, INTENT( IN) :: KMS,KME
!+++KAC
 REAL, PARAMETER :: ldarhtd = 0.95
 REAL, PARAMETER :: ldatmax = 289.15
 REAL, PARAMETER :: ldatmin = 268.00
 REAL, PARAMETER :: ldaa = 0.93
 REAL, PARAMETER :: ldab = 0.20
 REAL, PARAMETER :: ldac = 0.01
 REAL, PARAMETER :: ldad = 0.25
 REAL, PARAMETER :: ldarhmax = 1.03
 REAL, PARAMETER :: ldarhtd_damp = 0.75
!---KAC
 INTEGER, PARAMETER :: INUM = 12
 INTEGER, PARAMETER :: ILMM1 = 420
 INTEGER, PARAMETER :: IJMM1 = 480
 REAL, DIMENSION(kms:kme), INTENT(INOUT):: &
                          qv1d
 REAL, DIMENSION(kms:kme), INTENT(IN):: &
                          qg1d,t1d,p1d
 INTEGER :: II,JJ,IIT,KK,IM,DTIME
 INTEGER, PARAMETER :: id = 1
 !REAL DT
!===================lyyx time level:96
 REAL, DIMENSION(kps:kpe) :: QSAT
 REAL, DIMENSION(kps:kpe) :: es,angle
 REAL, DIMENSION(1:ilmm1,1:ijmm1,1:inum),SAVE ::LYYX1
 
 !dt = DT_IN
 
 dtime = int(dt)
 im = int(itimestep*dtime/600+1)-18
 
 !CALL lyy_lda
 if (id .eq. 1) then
   !read
        if (itimestep < 3602) then
          open(702,file=&
      './flashes_d01.dat',&
         form='formatted')
          read(702,'(420F8.3)')(((lyyx1(ii,jj,iit),ii=1,ilmm1),jj=1,ijmm1),iit=1,inum)
          close(702)
        endif
   !loop start
        ii=i
        jj=j
	
	 
	 do kk=kps,kpe
	 
           es(kk) = 611.2*exp(17.67*(t1d(kk)-273.15)/((t1d(kk)-273.15)+243.5))	!t1d-273.15=C; es=Pa
	   qsat(kk) = 0.622*es(kk)/p1d(kk)			!es=Pa; p1d=Pa; qsat=kg/kg
	 
	  ! Increase model Qv (assuming Qg and Qv aren't too high) because there is a flash
	   if (lyyx1(ii,jj,im) .gt. 1) then

          ! If Qg >= 3 g/kg no action
	  ! qg1d=kg/kg
             if (qg1d(kk).lt.0.003) then
	  ! If ambient Qvapor already >= 95% of env qsat no action
	  ! qv1d=kg/kg; qsat=kg/kg
	       if (qv1d(kk).lt.ldarhtd*qsat(kk)) then		
          ! Increase Qv between ldatmin and ldatmax only
	  ! t1d=K
	         if (t1d(kk).lt.ldatmax.and.t1d(kk).gt.ldatmin) then
	  ! multiply qg1d by 1000 so units are g/kg (not kg/kg) like in Fierro et al. (2012) 
		     angle(kk)=(1000.*qg1d(kk))**2.2
		     
		     ! multiply qsat by 1000. so units are in g/kg like in Fierro et al. (2012); then divide by 1000. to put qv1d into kg/kg
		     qv1d(kk)=ldaa*qsat(kk)+ldab*qsat(kk)*tanh(ldac*lyyx1(ii,jj,im))*(1-tanh(ldad*angle(kk))) 
!		     !write(*,*)'qv1d ok1'
		   if (qv1d(kk)/qsat(kk).gt.ldarhmax) then
		       qv1d(kk)=ldarhmax*qsat(kk)
		   endif
		  ! write(*,*)'qv1d ok2'
		 endif
	       endif
	     endif
	   ! Dampen the model Qv because there is no observed flash  
           else if (lyyx1(ii,jj,im).lt.0) then
             if (qv1d(kk).gt.ldarhtd_damp*qsat(kk)) then
                 qv1d(kk)=ldarhtd_damp*qsat(kk)
                ! write(*,*)'damp ok'
             endif
           else if (lyyx1(ii,jj,im).eq.0) then
               qv1d(kk)=qv1d(kk)
           endif	!lyyx1
	     
	 enddo	!kk loop
		 
 endif
 
 !CALL refl10cm_hm
 
 END SUBROUTINE lda_kac
 
 
 SUBROUTINE refl_kac(qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d, &
                      t1d, p1d, dBZ, kps, kpe, ii, jj)

! Framework
 USE module_wrf_error
 USE module_utility
 USE module_domain
 USE module_mp_radar
 USE module_model_constants, ONLY: CP, G, R => r_d, RV => r_v, EP_2

 IMPLICIT NONE

 !..Sub arguments
 INTEGER, INTENT(IN):: kps, kpe, ii, jj
 REAL, DIMENSION(kps:kpe), INTENT(IN)::                            &
                      qv1d, qr1d, nr1d, qs1d, ns1d, qg1d, ng1d, t1d, p1d		      
 REAL, DIMENSION(kps:kpe), INTENT(INOUT):: dBZ

 !..Local variables
 REAL, DIMENSION(kps:kpe):: temp, pres, qv, rho
 REAL, DIMENSION(kps:kpe):: rr, nr, rs, ns, rg, ng

 DOUBLE PRECISION, DIMENSION(kps:kpe):: ilamr, ilamg, ilams
 DOUBLE PRECISION, DIMENSION(kps:kpe):: N0_r, N0_g, N0_s
 DOUBLE PRECISION:: lamr, lamg, lams
 LOGICAL, DIMENSION(kps:kpe):: L_qr, L_qs, L_qg

 REAL, DIMENSION(kps:kpe):: ze_rain, ze_snow, ze_graupel
 DOUBLE PRECISION:: fmelt_s, fmelt_g
 DOUBLE PRECISION:: cback, x, eta, f_d

 INTEGER:: k, k_0, n
 LOGICAL:: melti
 
 !+++kac
 REAL,PARAMETER :: PI = 3.1415926535897932384626434
 !xam_r
 !xcrg(3)
 !xorg2
 !xobmr
 !xcre(2)
 !xam_s
 !xcsg(3)
 !xosg2
 !xobms
 !xcse(2)
 !xam_g
 !xcgg(3)
 !xogg2
 !xobmg
 !xcge(2)
 !xcrg(4)
 !ilamr
 !xcre(4)
 !PI
 !xcsg(4)
 !xcse(4)
 !xcgg(4)
 !xcge(4)

 !+---+
 
      do k = kps, kpe
         dBZ(k) = -35.0
      enddo

!+---+-----------------------------------------------------------------+
!..Put column of data into local arrays.
!+---+-----------------------------------------------------------------+
      do k = kps, kpe
         temp(k) = t1d(k)
         qv(k) = MAX(1.E-10, qv1d(k))
         pres(k) = p1d(k)
         rho(k) = 0.622*pres(k)/(R*temp(k)*(qv(k)+0.622))

         if (qr1d(k) .gt. 1.E-9) then
            rr(k) = qr1d(k)*rho(k)
            nr(k) = nr1d(k)*rho(k)
            lamr = (xam_r*xcrg(3)*xorg2*nr(k)/rr(k))**xobmr
            ilamr(k) = 1./lamr
            N0_r(k) = nr(k)*xorg2*lamr**xcre(2)
            L_qr(k) = .true.
         else
            rr(k) = 1.E-12
            nr(k) = 1.E-12
            L_qr(k) = .false.
         endif

         if (qs1d(k) .gt. 1.E-9) then
            rs(k) = qs1d(k)*rho(k)
            ns(k) = ns1d(k)*rho(k)
            lams = (xam_s*xcsg(3)*xosg2*ns(k)/rs(k))**xobms
            ilams(k) = 1./lams
            N0_s(k) = ns(k)*xosg2*lams**xcse(2)
            L_qs(k) = .true.
         else
            rs(k) = 1.E-12
            ns(k) = 1.E-12
            L_qs(k) = .false.
         endif

         if (qg1d(k) .gt. 1.E-9) then
            rg(k) = qg1d(k)*rho(k)
            ng(k) = ng1d(k)*rho(k)
            lamg = (xam_g*xcgg(3)*xogg2*ng(k)/rg(k))**xobmg
            ilamg(k) = 1./lamg
            N0_g(k) = ng(k)*xogg2*lamg**xcge(2)
            L_qg(k) = .true.
         else
            rg(k) = 1.E-12
            ng(k) = 1.E-12
            L_qg(k) = .false.
         endif
      enddo
      
!+---+-----------------------------------------------------------------+
!..Locate K-level of start of melting (k_0 is level above).
!+---+-----------------------------------------------------------------+
      melti = .false.
      k_0 = kps
      do k = kpe-1, kps, -1
         if ( (temp(k).gt.273.15) .and. L_qr(k)                         &
                                  .and. (L_qs(k+1).or.L_qg(k+1)) ) then
            k_0 = MAX(k+1, k_0)
            melti=.true.
            goto 195
         endif
      enddo
 195  continue     
      
 !+---+-----------------------------------------------------------------+
!..Assume Rayleigh approximation at 10 cm wavelength. Rain (all temps)
!.. and non-water-coated snow and graupel when below freezing are
!.. simple. Integrations of m(D)*m(D)*N(D)*dD.
!+---+-----------------------------------------------------------------+

      do k = kps, kpe
         ze_rain(k) = 1.e-22
         ze_snow(k) = 1.e-22
         ze_graupel(k) = 1.e-22
         if (L_qr(k)) ze_rain(k) = N0_r(k)*xcrg(4)*ilamr(k)**xcre(4)
         if (L_qs(k)) ze_snow(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)     &
                                 * (xam_s/900.0)*(xam_s/900.0)          &
                                 * N0_s(k)*xcsg(4)*ilams(k)**xcse(4)
         if (L_qg(k)) ze_graupel(k) = (0.176/0.93) * (6.0/PI)*(6.0/PI)  &
                                    * (xam_g/900.0)*(xam_g/900.0)       &
                                    * N0_g(k)*xcgg(4)*ilamg(k)**xcge(4)
      enddo     
      
 !+---+-----------------------------------------------------------------+
!..Special case of melting ice (snow/graupel) particles.  Assume the
!.. ice is surrounded by the liquid water.  Fraction of meltwater is
!.. extremely simple based on amount found above the melting level.
!.. Uses code from Uli Blahak (rayleigh_soak_wetgraupel and supporting
!.. routines).
!+---+-----------------------------------------------------------------+

      if (melti .and. k_0.ge.kps+1) then
       do k = k_0-1, kps, -1

!..Reflectivity contributed by melting snow
          if (L_qs(k) .and. L_qs(k_0) ) then
           fmelt_s = MAX(0.005d0, MIN(1.0d0-rs(k)/rs(k_0), 0.99d0))
           eta = 0.d0
           lams = 1./ilams(k)
           do n = 1, nrbins
              x = xam_s * xxDs(n)**xbm_s
              call rayleigh_soak_wetgraupel (x,DBLE(xocms),DBLE(xobms), &
                    fmelt_s, melt_outside_s, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_s, matrixstring_s,          &
                    inclusionstring_s, hoststring_s,                    &
                    hostmatrixstring_s, hostinclusionstring_s)
              f_d = N0_s(k)*xxDs(n)**xmu_s * DEXP(-lams*xxDs(n))
              eta = eta + f_d * CBACK * simpson(n) * xdts(n)
           enddo
           ze_snow(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif


!..Reflectivity contributed by melting graupel

          if (L_qg(k) .and. L_qg(k_0) ) then
           fmelt_g = MAX(0.005d0, MIN(1.0d0-rg(k)/rg(k_0), 0.99d0))
           eta = 0.d0
           lamg = 1./ilamg(k)
           do n = 1, nrbins
              x = xam_g * xxDg(n)**xbm_g
              call rayleigh_soak_wetgraupel (x,DBLE(xocmg),DBLE(xobmg), &
                    fmelt_g, melt_outside_g, m_w_0, m_i_0, lamda_radar, &
                    CBACK, mixingrulestring_g, matrixstring_g,          &
                    inclusionstring_g, hoststring_g,                    &
                    hostmatrixstring_g, hostinclusionstring_g)
              f_d = N0_g(k)*xxDg(n)**xmu_g * DEXP(-lamg*xxDg(n))
              eta = eta + f_d * CBACK * simpson(n) * xdtg(n)
           enddo
           ze_graupel(k) = SNGL(lamda4 / (pi5 * K_w) * eta)
          endif

       enddo
      endif

      do k = kpe, kps, -1
         dBZ(k) = 10.*log10((ze_rain(k)+ze_snow(k)+ze_graupel(k))*1.d18)
      enddo     
 
 END SUBROUTINE refl_kac
 
 
 SUBROUTINE countCells_kac( &
          ! Inputs
            refl10cm_kac, reflthreshold, cellcount_method,     &
          ! Order dependent args for domain, mem, and tile dims
            ids, ide, jds, jde, kds, kde,              &
            ims, ime, jms, jme, kms, kme,              &
            ips, ipe, jps, jpe, kps, kpe,              &
          ! Outputs
            cellcount_kac )

 USE module_dm, only: wrf_dm_sum_real

 IMPLICIT NONE
!-----------------------------------------------------------------

! Inputs
 REAL,    DIMENSION( ims:ime,kms:kme,jms:jme ), INTENT(IN   ) :: refl10cm_kac
 REAL,    INTENT(IN   ) :: reflthreshold
 INTEGER, INTENT(IN   ) :: cellcount_method

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe


! Outputs
 REAL,    DIMENSION( kms:kme ), INTENT(  OUT) :: cellcount_kac

! Local vars
 INTEGER :: i,k,j

!-----------------------------------------------------------------

 cellcount_kac(kps:kpe) = 0.
 DO j=jps,jpe
   DO k=kps,kpe
     DO i=ips,ipe
       IF ( refl10cm_kac(i,k,j) .gt. reflthreshold ) THEN
         cellcount_kac(k) = cellcount_kac(k) + 1
       ENDIF
     ENDDO
   ENDDO
 ENDDO

 IF ( cellcount_method .eq. 2 ) THEN
   DO k=kps,kpe
     cellcount_kac(k) = wrf_dm_sum_real(cellcount_kac(k))
   ENDDO
 ENDIF

 END SUBROUTINE countCells_kac
 
 
 !--- kac ---

!**********************************************************************
!
! Price and Rind 1993 base on cold cloud depth (CCD)
!
! Price, C. and D. Rind (1993), What determines the cloud‐to‐ground lightning
! fraction in thunderstorms?, Geophys. Res. Lett., 20(6), 463–466, doi:10.1029/93GL00226.
!
! Valid range of CCD is set to 5.5-14 km. Beyond this range CCD is assumed
! to be 5.5 or 14 for continuity.
!
!**********************************************************************
 SUBROUTINE iccg_crm_pr93( &
                            refl, reflthreshold, t, z,                 &
                          ! Order dependent args for domain, mem, and tile dims
                            ids, ide, jds, jde, kds, kde,              &
                            ims, ime, jms, jme, kms, kme,              &
                            ips, ipe, jps, jpe, kps, kpe,              &
                          ! Input
                            total_flashrate,                           &
                          ! Output
                            ic_flashrate, cg_flashrate                 &
                        )
!-----------------------------------------------------------------
 IMPLICIT NONE
!-----------------------------------------------------------------
! Inputs
 REAL,    DIMENSION( ims:ims, kms:kme, jms:jme ), INTENT(IN   ) :: refl, t, z
 REAL,                                            INTENT(IN   ) :: reflthreshold

! Order dependent args for domain, mem, and tile dims
 INTEGER, INTENT(IN   )    ::       ids,ide, jds,jde, kds,kde
 INTEGER, INTENT(IN   )    ::       ims,ime, jms,jme, kms,kme
 INTEGER, INTENT(IN   )    ::       ips,ipe, jps,jpe, kps,kpe

! Primary inputs and outpus
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(IN   ) :: total_flashrate   
 REAL,    DIMENSION( ims:ime,          jms:jme ), INTENT(  OUT) :: ic_flashrate, cg_flashrate

! Local variables
 INTEGER :: kfreeze, ktop

 INTEGER :: i,j,k
 REAL    :: ratio, cgfrac, depth

 REAL, PARAMETER :: dH_min = 5.5
 REAL, PARAMETER :: dH_max = 14.

 REAL, PARAMETER :: coef_A = 0.021
 REAL, PARAMETER :: coef_B = -0.648
 REAL, PARAMETER :: coef_C = 7.493
 REAL, PARAMETER :: coef_D = -36.54
 REAL, PARAMETER :: coef_E = 63.09
!-----------------------------------------------------------------

 ic_flashrate(ips:ipe,jps:jpe) = 0.
 cg_flashrate(ips:ipe,jps:jpe) = 0.

 jloop: DO j=jps,jpe
    iloop: DO i=ips,ipe
    IF ( total_flashrate(i,j) .gt. 0.) THEN
        ktop = kpe
        do while ( refl(i,ktop,j) .lt. reflthreshold .and. ktop .gt. kps)
          ktop = ktop-1
        enddo

        kfreeze = ktop
        DO WHILE ( t(i,kfreeze,j) .lt. 273.15 .and. ktop .gt. kps )
            kfreeze = kfreeze - 1
        ENDDO

        depth = ( z(i,ktop,j) - z(i,kfreeze,j) ) * 1E-3
        IF (depth .le. 0.) CONTINUE
        depth = max( dH_min, min( dH_max, depth ))

        ratio = (((coef_A*depth+coef_B )*depth+coef_C)*depth+coef_D)*depth+coef_E
        cgfrac = 1./(ratio+1.)

        cg_flashrate(i,j) = total_flashrate(i,j) * cgfrac
        ic_flashrate(i,j) = total_flashrate(i,j) - cg_flashrate(i,j)
    ENDIF
    ENDDO iloop
 ENDDO jloop

 END SUBROUTINE iccg_crm_pr93

 END MODULE module_ltng_crmpr92
