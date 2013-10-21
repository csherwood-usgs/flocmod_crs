!&E==========================================================================
!&E                   ***  MODULE  FLOCMOD  ***
!&E
!&E
!&E ** Purpose : concerns all subroutines related to flocculation processes 
!&E    originally based on the FLOCMOD routine developed from SiAM3D, Verney et al., 2011
!&E
!&E ** Description :
!&E     subroutine flocmod_alloc            ! allocates variables/arrays for flocmod
!&E     subroutine flocmod_init_param       !  defines parameterization of flocculation 
!&E                                         ! processes
!&E
!&E     subroutine flocmod_kernels          ! computes flocculation kernels (agregation/fragmentation)
!&E
!&E     subroutine flocmod_agregation_statistics   ! computes shear and differential settling statistics (ASH/ADS)
!&E                                                ! used by flocmod_kernels, only at initialisation
!&E
!&E     subroutine flocmod_collfrag         ! computes fragmentation by collision term, at every sub time step as direct function of turbulence
!&E
!&E     subroutine flocmod_comp_fsd         ! computes the new floc size distribution at each sub time step. called by flocmod_main
!&E
!&E     subroutine flocmod_main             ! routine called by the main subroutine in MARS3D/ROMS before advection, at every time step
!&E
!&E     subroutine flocmod_mass_control     ! check that after flocculation, no negative masses are produced or not over f_mneg_param. called by flocmod_main
!&E
!&E     subroutine flocmod_mass_redistribute ! redistribute negative masses (below f_mneg_param) to all remaining positive classes
!&E
!&E     subroutine flocmod_comp_g           ! computes shear rates values. Actually based on Nezu and Nakagawa from u*...
!&E                                         ! but can be calculated directly from the turbulent dissipation energy if calculated by the hydrodynamic model
!&E
!&E
!&E ** History :
!&E     ! 2013-09 (Romaric Verney)
!&E
!&E==========================================================================

!!===========================================================================
PROGRAM flocmod_main

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_main  ***
  !&E
  !&E ** Purpose : main routine used to compute aggregation /fragmentation processes at every i,j,k cell
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : 
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !
  ! Modfied by csherwood@usgs.gov for testing prior to putting into ROMS

  !! * Modules used
  USE comvars

  !! * Local declarations
  IMPLICIT NONE
  LOGICAL      :: f_ld50,f_ld10,f_ld90
  INTEGER      :: iv1,i,j,k
  integer      :: fid = 60
  integer      :: nt
  REAL(KIND=rsh)     :: cvtotmud,dttemp,f_clim,Gval,mneg,dt1,f_csum,dtmin
  REAL(KIND=rsh),DIMENSION(1:nv_mud)     :: cv_tmp,NNin,NNout
  real(kind=rsh) :: tstart, tend

  !!--------------------------------------------------------------------------
  !! * Executable part

  !! CRS - Hardwire some values for now
  eps = tiny(1.0)
  dt = 1.0_rsh
  tstart = 0.0_rsh
  tend   = 750.0_rsh * 60.0_rsh
  nt = 0

  call flocmod_init_param

  call kernal_stats

  ! attempt to initialize w/ test case values
  cv_wat = 0.0_rsh
  cv_wat(5,1,1,1)=0.093_rsh
  open(UNIT=fid,file='flocmod.dat',status='UNKNOWN')

  f_clim=0.001 ! min concentration below which flocculation processes are not calculated

  t = tstart
  do while (t.lt.tend)
     DO i=limin,limax
        DO j=ljmin,ljmax
           DO k=1,kmax
              nt = nt+1
              dtmin=dt

              f_dt=dt
              dttemp=0.0_rsh

              cv_tmp(1:nv_mud)=cv_wat(imud1:nvpc,k,i,j) ! concentration of all mud classes in one grid cell
              cvtotmud=sum(cv_tmp(1:nv_mud))

              NNin(1:nv_mud)=cv_tmp(1:nv_mud)/f_mass(1:nv_mud)

              Do iv1=1,nv_mud
                 if (NNin(iv1).lt.0.0_rsh) then
                    write(*,*) '***********************************************'
                    write(*,*) 'CAUTION, negative mass at cell i,j,k : ',i,j,k
                    write(*,*) '***********************************************'        
                 endif
              ENDDO

              if (cvtotmud .gt. f_clim) then

                 call flocmod_comp_g(k,i,j,Gval) ! This used to be here

                 f_gval(k,i,j)=Gval

                 if (mod(nt,600).eq.0) write(*,*) 't, G, cvtotmud',t,Gval,cvtotmud

                 do while (dttemp .le. dt)
                    !     print*, 'f_dt:',f_dt

                    call flocmod_comp_fsd(NNin,NNout,Gval)
                    call flocmod_mass_control(NNout,mneg)

                    !     write(*,*) 'mneg',mneg

                    if (mneg .gt. f_mneg_param) then
                       DO while (mneg .gt. f_mneg_param)
                          f_dt=MIN(f_dt/2._rsh,dt-dttemp)
                          call flocmod_comp_fsd(NNin,NNout,Gval)
                          call flocmod_mass_control(NNout,mneg)   
                       ENDDO
                       !         if (f_dt.lt.1._rsh) then
                       !           write(*,*) 'apres : Gval,f_dt',Gval, f_dt,dttemp
                       !  endif
                    else

                       if (f_dt.lt.dt) then
                          DO while (mneg .lt.f_mneg_param)

                             if (dttemp+f_dt .eq. dt) then
                                call flocmod_comp_fsd(NNin,NNout,Gval)
                                exit
                             else
                                dt1=f_dt
                                f_dt=MIN(2._rsh*f_dt,dt-dttemp)
                                call flocmod_comp_fsd(NNin,NNout,Gval)
                                call flocmod_mass_control(NNout,mneg)
                                if (mneg .gt. f_mneg_param) then 
                                   f_dt=dt1
                                   call flocmod_comp_fsd(NNin,NNout,Gval)
                                   exit
                                endif
                             endif
                          ENDDO
                       endif
                    endif
                    dtmin=MIN(dtmin,f_dt)   
                    dttemp=dttemp+f_dt
                    NNin(:)=NNout(:) ! update new Floc size distribution

                    CALL flocmod_mass_redistribute(NNin) ! redistribute negative masses if any on positive classes, 
                    ! depends on f_mneg_param

                    if (abs(sum(NNin(1:nv_mud)*f_mass(1:nv_mud))-cvtotmud).gt.epsilon*100._rsh) then

                       write(*,*) 'CAUTION flocculation routine not conservative!!!'
                       write(*,*) 'time,i,j,k= ',t,i,j,k
                       write(*,*) 'f_dt= ',f_dt   
                       write(*,*) 'before : cvtotmud= ',cvtotmud
                       write(*,*) 'after  : cvtotmud= ',sum(NNin(1:nv_mud)*f_mass(1:nv_mud))
                       write(*,*) 'absolute difference  : cvtotmud= ',abs(cvtotmud-sum(NNin(1:nv_mud)*f_mass(1:nv_mud)))
                       write(*,*) 'absolute difference reference  : espilon= ',epsilon
                       write(*,*) 'before redistribution', sum(NNout(1:nv_mud)*f_mass(1:nv_mud))
                       write(*,*) 'after redistribution', sum(NNin(1:nv_mud)*f_mass(1:nv_mud)) 
                       write(*,*) 'Simultation stopped'

                       STOP
                    endif

                    if (dttemp.eq.dt) exit

                 ENDDO ! loop on full dt

              endif ! only if cvtotmud > f_clim

              if (abs(sum(NNin(1:nv_mud)*f_mass(1:nv_mud))-cvtotmud).gt.epsilon*10._rsh) then

                 write(*,*) 'CAUTION flocculation routine not conservative!!!'
                 write(*,*) 'time,i,j,k= ',t,i,j,k
                 write(*,*) 'before : cvtotmud= ',cvtotmud
                 write(*,*) 'after  : cvtotmud= ',sum(NNin(1:nv_mud)*f_mass(1:nv_mud))
                 write(*,*) 'absolute difference  : cvtotmud= ',abs(cvtotmud-sum(NNin(1:nv_mud)*f_mass(1:nv_mud)))
                 write(*,*) 'absolute difference reference  : espilon= ',epsilon
                 write(*,*) 'Simultation stopped'

                 STOP
              endif

              ! update mass concentration for all mud classes
              cv_wat(imud1:nvpc,k,i,j)=NNin(1:nv_mud)*f_mass(1:nv_mud)

              ! compute floc distribution statistics before output
              f_csum=0.0_rsh
              f_ld50=.true.
              f_ld10=.true.
              f_ld90=.true.

              f_davg(k,i,j)=sum(NNin(1:nv_mud)*f_mass(1:nv_mud)*f_diam(1:nv_mud))/(sum(NNin(1:nv_mud)*f_mass(1:nv_mud))+epsilon)
              f_dtmin(k,i,j)= dtmin

              DO iv1=1,nv_mud

                 f_csum=f_csum+NNin(iv1)*f_mass(iv1)/((sum(NNin(1:nv_mud)*f_mass(1:nv_mud)))+epsilon)
                 if (f_csum.gt.0.1_rsh .and. f_ld10) then
                    f_d10(k,i,j)=f_diam(iv1)    
                    f_ld10=.false.
                 endif

                 if (f_csum.gt.0.5_rsh .and. f_ld50) then
                    f_d50(k,i,j)=f_diam(iv1)
                    f_ld50=.false.
                 endif

                 if (f_csum.gt.0.9_rsh .and. f_ld90) then
                    f_d90(k,i,j)=f_diam(iv1)
                    f_ld90=.false.
                 endif
              ENDDO

           ENDDO
        ENDDO
     ENDDO
     write(fid,*) t, Gval, f_d50(1,1,1)*1e6
     t = t+dt
  enddo
  close(fid)
  write(*,*) 'END flocmod_main'
END program flocmod_main

!!===========================================================================
SUBROUTINE flocmod_init_param
  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_init_param  ***
  !&E
  !&E ** Purpose : Initialization of parameterizations used for flocculation model
  !&E              
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used

  !USE comvars2d,     ONLY : l_chrono,lscreen,dt
  USE comvars
  !,    ONLY : limin,limax,ljmin,ljmax,kmax,limin,ljmin, &
  !       l_testcase,grav,rhoref,diam_sed,ros,f_ws,dt,nvpc

  !! * Local declarations
  INTEGER   :: iv
  CHARACTER(len=80) :: filepc
  NAMELIST/namflocmod/ l_ADS,l_ASH,l_COLLFRAG,f_dp0,f_nf,f_dmax, &
       f_nb_frag,f_alpha,f_beta,f_ater, &
       f_ero_frac,f_ero_nbfrag,f_ero_iv,f_mneg_param, &
       f_collfragparam

  !!--------------------------------------------------------------------------
  !! * Executable part
  filepc='./paraflocmod.txt'
  OPEN(50,file=filepc,status='old',form='formatted',access='sequential')
  READ(50,namflocmod)
  CLOSE(50)
  ! fin de lecture de la namelist

  ! CRS - initialize for test case
  diam_sed = &
       (/4.0, 6.1, 9.3, 14.2, 21.8, 33.2, 50.7, 77.5, 118.3, 180.6, 275.8, 421.2, 643.2, 982.3, 1500.0/)
  diam_sed = diam_sed*1.e-6
  do iv=1,nv_mud
     ros(iv)=2650.0_rsh
  enddo

  !!--------------------------------------------------
  !! floc characteristics
  f_diam(1:nv_mud)=diam_sed(imud1:nvpc) ! floc size from variable.dat
  f_vol(1:nv_mud)=pi/6._rsh*(f_diam(1:nv_mud))**3
  f_rho(1:nv_mud)=rhoref+(ros(1:nv_mud)-rhoref)*(f_dp0/f_diam(1:nv_mud))**(3._rsh-f_nf)
  f_mass(1:nv_mud)=f_vol(1:nv_mud)*(f_rho(1:nv_mud)-rhoref)
  f_mass(nv_mud+1)=f_mass(nv_mud)*2_rsh+1._rsh  
  if (f_diam(1).eq.f_dp0)  then
     f_mass(1)=f_vol(1)*ros(1)
  endif
  f_ws(1:nv_mud)=grav*(f_rho(1:nv_mud)-rhoref)*f_diam(1:nv_mud)**2._rsh/(18._rsh*0.001_rsh)
  f_dt=dt

  !!--------------------------------------------------
  WRITE(*,*) ' '
  WRITE(*,*) '***********************'
  WRITE(*,*) '    FLOCMOD'
  WRITE(*,*) '***********************'
  WRITE(*,*) 'class  diameter (um)  volume (m3)  density (kg/m3)  mass (kg) Ws (m/s)'
  DO iv=1,nv_mud
     WRITE(*,*) iv,f_diam(iv)*1e6,f_vol(iv),f_rho(iv),f_mass(iv),f_ws(iv)
  ENDDO
  WRITE(*,*) ' '
  WRITE(*,*) ' *** PARAMETERS ***'
  WRITE(*,*) 'Primary particle size (f_dp0)                                : ',f_dp0
  WRITE(*,*) 'Fractal dimension (f_nf)                                     : ',f_nf
  WRITE(*,*) 'Flocculation efficiency (f_alpha)                            : ',f_alpha
  WRITE(*,*) 'Floc break up parameter (f_beta)                             : ',f_beta
  WRITE(*,*) 'Nb of fragments (f_nb_frag)                                  : ',f_nb_frag
  WRITE(*,*) 'Ternary fragmentation (f_ater)                               : ',f_ater
  WRITE(*,*) 'Floc erosion (% of mass) (f_ero_frac)                        : ',f_ero_frac
  WRITE(*,*) 'Nb of fragments by erosion (f_ero_nbfrag)                    : ',f_ero_nbfrag
  WRITE(*,*) 'fragment class (f_ero_iv)                                    : ',f_ero_iv
  WRITE(*,*) 'negative mass tolerated before redistribution (f_mneg_param) : ',f_mneg_param
  WRITE(*,*) 'Boolean for differential settling aggregation (L_ADS)        : ',l_ADS
  WRITE(*,*) 'Boolean for shear aggregation (L_ASH)                        : ',l_ASH
  WRITE(*,*) 'Boolean for collision fragmenation (L_COLLFRAG)              : ',l_COLLFRAG
  WRITE(*,*) 'Collision fragmentation parameter (f_collfragparam)          : ',f_collfragparam
  WRITE(*,*) ' '
  WRITE(*,*) 'Value of eps                                                 : ',eps
  WRITE(*,*) ' '
  WRITE(*,*) '*** END FLOCMOD INIT *** '    

  if (.not.l_ADS .and. .not.l_ASH) then
     write(*,*) 'CAUTION : incompatible flocculation kernel options : '
     write(*,*) '*****************************************************'
     write(*,*) 'l_ADS=',l_ADS
     write(*,*) 'l_ASH=',l_ASH
     write(*,*) 'simulation stopped'
     stop
  endif

  ! kernels copmutation
  CALL flocmod_kernels

  write(*,*) 'END FLOCMOD_INIT_PARAM'
END SUBROUTINE flocmod_init_param

!!===========================================================================
SUBROUTINE flocmod_kernels

  !&E--------------------------------------------------------------------------
  !&E                 ***  flocmod_kernels  ***
  !&E
  !&E ** Purpose : computations of agregation/fragmentation kernels for FLOCMOD
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : dyn3dxtlscoef, dyn3dytlscoef, turbk, turbpsi
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars

  !! * Local declarations
  LOGICAL        :: f_test
  REAL(KIND=rsh) :: f_weight,mult,dfragmax
  INTEGER        :: iv1,iv2,iv3

  !!--------------------------------------------------------------------------
  !! * Executable part
  f_test=.true.
  dfragmax=0.00003
  ! compute collision probability
  CALL flocmod_agregation_statistics

  !********************************************************************************
  ! agregation : GAIN : f_g1_sh and f_g1_ds
  !********************************************************************************
  DO iv1=1,nv_mud
     DO iv2=1,nv_mud
        DO iv3=iv2,nv_mud
           if((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1-1) &
                .and. ((f_mass(iv2)+f_mass(iv3)) .le. f_mass(iv1))) then

              f_weight=(f_mass(iv2)+f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

           elseif ((f_mass(iv2)+f_mass(iv3)) .gt. f_mass(iv1) &
                .and. ((f_mass(iv2)+f_mass(iv3)) .lt. f_mass(iv1+1))) then

              if (iv1 .eq. nv_mud) then
                 f_weight=1._rsh
              else
                 f_weight=1._rsh-(f_mass(iv2)+f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
              endif

           else
              f_weight=0.0_rsh
           endif

           f_g1_sh(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_sh(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)
           f_g1_ds(iv2,iv3,iv1)=f_weight*f_alpha*f_coll_prob_ds(iv2,iv3)*(f_mass(iv2)+f_mass(iv3))/f_mass(iv1)

        ENDDO
     ENDDO
  ENDDO

  !********************************************************************************
  ! Shear fragmentation : GAIN : f_g3
  !********************************************************************************
  write(*,*) 'Calc f_g3: dfragmax, f_nb_frag:',dfragmax,f_nb_frag
  write(*,*) 'f_mass(0), f_mass(nv_mud+1):',f_mass(0),f_mass(nv_mud+1)
  DO iv1=1,nv_mud
     DO iv2=iv1,nv_mud

        if (f_diam(iv2).gt.dfragmax) then
           ! binary fragmentation
           if (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1-1) &
                .and. f_mass(iv2)/f_nb_frag .le. f_mass(iv1)) then
              if (iv1 .eq. 1) then 
                 f_weight=1._rsh
              else
                 f_weight=(f_mass(iv2)/f_nb_frag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
              endif
           elseif (f_mass(iv2)/f_nb_frag .gt. f_mass(iv1) &
                .and. f_mass(iv2)/f_nb_frag .lt. f_mass(iv1+1)) then
              f_weight=1._rsh-(f_mass(iv2)/f_nb_frag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
           else
              f_weight=0._rsh
           endif
        else
           f_weight=0.0_rsh
        endif
        if(iv2.eq.1)then
        write(*,*) 'iv2, iv1, f_weight: ',iv2,iv1,f_weight
        write(*,*) f_diam(iv2)-f_dp0
        write(*,*) (f_diam(iv2)-f_dp0)/f_dp0
        write(*,*) (max(eps,(f_diam(iv2)-f_dp0)/f_dp0))**(3._rsh-f_nf)
        endif
        f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1._rsh-f_ero_frac)*(1._rsh-f_ater)*f_weight*f_beta &
             *f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**(3._rsh-f_nf)           &
             *f_mass(iv2)/f_mass(iv1)
        !write(*,*) 'f_g3(iv2,iv1)     : ',f_g3(iv2,iv1)

        ! ternary fragmentation
        if (f_diam(iv2).gt.dfragmax) then
           if (f_mass(iv2)/(2._rsh*f_nb_frag) .gt. f_mass(iv1-1) &
                .and. f_mass(iv2)/(2._rsh*f_nb_frag) .le. f_mass(iv1)) then
              if (iv1 .eq. 1) then 
                 f_weight=1._rsh
              else
                 f_weight=(f_mass(iv2)/(2._rsh*f_nb_frag)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
              endif
           elseif (f_mass(iv2)/(2._rsh*f_nb_frag) .gt. f_mass(iv1) &
                .and. f_mass(iv2)/(2._rsh*f_nb_frag) .lt. f_mass(iv1+1)) then
              f_weight=1._rsh-(f_mass(iv2)/(2._rsh*f_nb_frag)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
           else
              f_weight=0._rsh
           endif

           ! update for ternary fragments
           f_g3(iv2,iv1)=f_g3(iv2,iv1)+(1._rsh-f_ero_frac)*(f_ater)*f_weight*f_beta &
                *f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**(3._rsh-f_nf)             &
                *f_mass(iv2)/f_mass(iv1)   

           ! Floc erosion
           if ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(f_ero_iv)) then
              if (((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt.f_mass(iv1-1)) &
                   .and. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .le. f_mass(iv1)) then
                 if (iv1 .eq. 1) then
                    f_weight=1._rsh
                 else
                    f_weight=(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))
                 endif
              elseif ((f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(iv1) &
                   .and. (f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag) .lt. f_mass(iv1+1)) then
                 f_weight=1._rsh-(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1))
              else
                 f_weight=0._rsh
              endif
 
              ! update for eroded floc masses 
              f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_weight*f_beta                 &
                   *f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**(3._rsh-f_nf)           &
                                !    *f_mass(iv2)/f_mass(iv1)                        &
                   *(f_mass(iv2)-f_mass(f_ero_iv)*f_ero_nbfrag)/f_mass(iv1)

              if (iv1 .eq. f_ero_iv) then
                 f_g3(iv2,iv1)=f_g3(iv2,iv1)+f_ero_frac*f_beta                       &
                      *f_diam(iv2)*(max(eps,(f_diam(iv2)-f_dp0))/f_dp0)**(3._rsh-f_nf)        &
                                !    *f_mass(iv2)/f_mass(iv1)                        &
                      *f_ero_nbfrag*f_mass(f_ero_iv)/f_mass(iv1)
              endif
           endif
        endif ! condition on dfragmax
     ENDDO
  ENDDO

  !********************************************************************************
  !  Shear agregation : LOSS : f_l1
  !********************************************************************************

  DO iv1=1,nv_mud
     DO iv2=1,nv_mud

        if(iv2 .eq. iv1) then
           mult=2._rsh
        else
           mult=1._rsh
        endif

        f_l1_sh(iv2,iv1)=mult*f_alpha*f_coll_prob_sh(iv2,iv1) 
        f_l1_ds(iv2,iv1)=mult*f_alpha*f_coll_prob_ds(iv2,iv1) 

     ENDDO
  ENDDO

  !********************************************************************************
  !  Shear fragmentation : LOSS : f_l2
  !********************************************************************************

  DO iv1=1,nv_mud
     if (f_diam(iv1).gt.dfragmax) then
        ! shear fragmentation
        f_l3(iv1)=f_l3(iv1)+(1._rsh-f_ero_frac)*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3._rsh-f_nf)

        ! shear erosion
        if ((f_mass(iv1)-f_mass(f_ero_iv)*f_ero_nbfrag) .gt. f_mass(f_ero_iv)) then
           f_l3(iv1)=f_l3(iv1)+f_ero_frac*f_beta*f_diam(iv1)*((f_diam(iv1)-f_dp0)/f_dp0)**(3._rsh-f_nf)
        endif
     endif
  ENDDO

  write(*,*) 'END FLOCMOD_KERNELS'

END SUBROUTINE flocmod_kernels

SUBROUTINE kernal_stats
  USE comvars
  write(*,*) 'Sum of kernal coefficients:'
  write(*,*) 'f_coll_prob_sh',sum(f_coll_prob_sh)
  write(*,*) 'f_coll_prob_ds',sum(f_coll_prob_ds)
  write(*,*) 'f_g1_sh',sum(f_g1_sh)
  write(*,*) 'f_g1_ds',sum(f_g1_ds)
  write(*,*) 'f_l1_sh',sum(f_l1_sh)
  write(*,*) 'f_l1_ds',sum(f_l1_ds)
  write(*,*) 'f_g3',sum(f_g3)
  write(*,*) 'f_l3',sum(f_l3)
  write(*,*) 'f_g4',sum(f_g4)
  write(*,*) 'f_l4',sum(f_l4)
  write(*,*) ' '
END SUBROUTINE kernal_stats

!!===========================================================================
SUBROUTINE flocmod_agregation_statistics

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_agregation_statistics  ***
  !&E
  !&E ** Purpose : computation of shear / differential settling statistics
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_kernels
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars
  !,    ONLY : limin,limax,ljmin,ljmax,kmax,grav,rhoref

  !! * Local declarations
  INTEGER      :: iv1,iv2
  REAL(KIND=rsh) :: mu
  !!--------------------------------------------------------------------------
  !! * Executable part

  mu=0.001_rsh

  DO iv1=1,nv_mud
     DO iv2=1,nv_mud

        f_coll_prob_sh(iv1,iv2)=1._rsh/6._rsh*(f_diam(iv1)+f_diam(iv2))**3._rsh

        f_coll_prob_ds(iv1,iv2)=0.25_rsh*pi*(f_diam(iv1)+f_diam(iv2))**2._rsh &
             *grav/mu*abs((f_rho(iv1)-rhoref)*f_diam(iv1)**2._rsh &
             -(f_rho(iv2)-rhoref)*f_diam(iv2)**2._rsh)

     ENDDO
  ENDDO

  write(*,*) 'END flocmod_agregation_statistics'

END SUBROUTINE flocmod_agregation_statistics

!!===========================================================================
SUBROUTINE flocmod_collfrag(Gval)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_collfrag  ***
  !&E
  !&E ** Purpose : computation of collision fragmentation term, based on McAnally and Mehta, 2001
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_comp_fsd
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars
  ! ,    ONLY : limin,limax,ljmin,ljmax,kmax,grav,rhoref

  REAL(KIND=rsh),INTENT(in) :: Gval

  !! * Local declarations
  INTEGER      :: iv1,iv2,iv3
  REAL(KIND=rsh) :: f_fp,f_fy,f_cfcst,gcolfragmin,gcolfragmax,gcolfragiv1,gcolfragiv2 &
       ,f_weight,mult
  !!--------------------------------------------------------------------------
  !! * Executable part

  f_fp=0.1_rsh
  f_fy=1e-10
  f_cfcst=3._rsh/16._rsh
  f_g4(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
  f_l4(1:nv_mud,1:nv_mud)=0.0_rsh

  DO iv1=1,nv_mud
     DO iv2=1,nv_mud
        DO iv3=iv2,nv_mud
           ! fragmentation after collision probability based on Gval for particles iv2 and iv3
           ! gcolfrag=(collision induced shear) / (floc strength)

           gcolfragmin=2._rsh*(Gval*(f_diam(iv2)+f_diam(iv3)))**2._rsh*f_mass(iv2)*f_mass(iv3)  &
                /(pi*f_fy*f_fp*f_diam(iv3)**2._rsh*(f_mass(iv2)+f_mass(iv3))         &
                *((f_rho(iv3)-rhoref)/rhoref)**(2._rsh/(3._rsh-f_nf)))

           gcolfragmax=2._rsh*(Gval*(f_diam(iv2)+f_diam(iv3)))**2._rsh*f_mass(iv2)*f_mass(iv3)  &
                /(pi*f_fy*f_fp*f_diam(iv2)**2._rsh*(f_mass(iv2)+f_mass(iv3))         &
                *((f_rho(iv2)-rhoref)/rhoref)**(2._rsh/(3._rsh-f_nf)))


           ! first case : iv3 not eroded, iv2 eroded forming 2 particles : iv3+f_cfcst*iv2 / iv2-f_cfcst*iv2
           if (gcolfragmin.lt.1._rsh .and. gcolfragmax.ge.1_rsh) then  

              if (((f_mass(iv3)+f_cfcst*f_mass(iv2)).gt.f_mass(iv1-1)) .and.  &
                   ((f_mass(iv3)+f_cfcst*f_mass(iv2)).le.f_mass(iv1))) then

                 f_weight=((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))

              elseif (f_mass(iv3)+f_cfcst*f_mass(iv2).gt.f_mass(iv1)  .and. &
                   f_mass(iv3)+f_cfcst*f_mass(iv2).lt.f_mass(iv1+1)) then

                 if (iv1.eq.nv_mud) then 
                    f_weight=1._rsh
                 else

                    f_weight=1._rsh-((f_mass(iv3)+f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
                 endif

              else
                 f_weight=0.0_rsh
              endif

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
                   *(f_mass(iv3)+f_cfcst*f_mass(iv2))/f_mass(iv1)

              if (f_mass(iv2)-f_cfcst*f_mass(iv2).gt.f_mass(iv1-1)   .and. &
                   f_mass(iv2)-f_cfcst*f_mass(iv2).le.f_mass(iv1)) then

                 f_weight=((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))

              elseif (f_mass(iv2)-f_cfcst*f_mass(iv2).gt.f_mass(iv1)  .and.  &
                   f_mass(iv2)-f_cfcst*f_mass(iv2).lt.f_mass(iv1+1)) then

                 if (iv1.eq.nv_mud) then 
                    f_weight=1._rsh
                 else

                    f_weight=1._rsh-((f_mass(iv2)-f_cfcst*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
                 endif

              else
                 f_weight=0.0_rsh
              endif

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
                   *(f_mass(iv2)-f_cfcst*f_mass(iv2))/f_mass(iv1)


              ! second case : iv3 eroded and iv2 eroded forming 3 particles : iv3-f_cfcst*iv3 / iv2-f_cfcst*iv2 / f_cfcst*iv3+f_cfcst*iv2
           elseif (gcolfragmin.ge.1._rsh .and. gcolfragmax.ge.1_rsh) then  ! iv2 and iv3 eroded forming new (third) particle

              if (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).gt.f_mass(iv1-1) .and.  &
                   f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).le.f_mass(iv1)) then

                 f_weight=((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1)))

              elseif (f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).gt.f_mass(iv1) .and.  &
                   f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3).lt.f_mass(iv1+1)) then

                 if (iv1.eq.nv_mud) then 
                    f_weight=1._rsh
                 else

                    f_weight=1._rsh-((f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
                 endif

              else
                 f_weight=0.0_rsh
              endif

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
                   *(f_cfcst*f_mass(iv2)+f_cfcst*f_mass(iv3))/f_mass(iv1)

              if ((1._rsh-f_cfcst)*f_mass(iv2).gt.f_mass(iv1-1) .and.  &
                   (1._rsh-f_cfcst)*f_mass(iv2).le.f_mass(iv1)) then

                 f_weight=((1._rsh-f_cfcst)*f_mass(iv2)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

              elseif ((1._rsh-f_cfcst)*f_mass(iv2).gt.f_mass(iv1) .and.  &
                   (1._rsh-f_cfcst)*f_mass(iv2).lt.f_mass(iv1+1)) then

                 if (iv1.eq.nv_mud) then 
                    f_weight=1._rsh
                 else

                    f_weight=1._rsh-(((1._rsh-f_cfcst)*f_mass(iv2)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
                 endif

              else
                 f_weight=0.0_rsh
              endif

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
                   *((1._rsh-f_cfcst)*f_mass(iv2))/f_mass(iv1) 


              if ((1._rsh-f_cfcst)*f_mass(iv3).gt.f_mass(iv1-1) .and.  &
                   (1._rsh-f_cfcst)*f_mass(iv3).le.f_mass(iv1)) then

                 f_weight=((1._rsh-f_cfcst)*f_mass(iv3)-f_mass(iv1-1))/(f_mass(iv1)-f_mass(iv1-1))

              elseif ((1._rsh-f_cfcst)*f_mass(iv3).gt.f_mass(iv1) .and.  &
                   (1._rsh-f_cfcst)*f_mass(iv3).lt.f_mass(iv1+1)) then

                 if (iv1.eq.nv_mud) then 
                    f_weight=1._rsh
                 else

                    f_weight=1._rsh-(((1._rsh-f_cfcst)*f_mass(iv3)-f_mass(iv1))/(f_mass(iv1+1)-f_mass(iv1)))
                 endif

              else
                 f_weight=0.0_rsh
              endif

              f_g4(iv2,iv3,iv1)=f_g4(iv2,iv3,iv1)+f_weight*(f_coll_prob_sh(iv2,iv3))   &
                   *((1._rsh-f_cfcst)*f_mass(iv3))/f_mass(iv1)


           endif ! end collision test case
        ENDDO
     ENDDO
  ENDDO

  DO iv1=1,nv_mud
     DO iv2=1,nv_mud

        gcolfragiv1=2._rsh*(Gval*(f_diam(iv1)+f_diam(iv2)))**2._rsh*f_mass(iv1)*f_mass(iv2)  &
             /(pi*f_fy*f_fp*f_diam(iv1)**2._rsh*(f_mass(iv1)+f_mass(iv2))         &
             *((f_rho(iv1)-rhoref)/rhoref)**(2._rsh/(3._rsh-f_nf)))

        gcolfragiv2=2._rsh*(Gval*(f_diam(iv1)+f_diam(iv2)))**2._rsh*f_mass(iv1)*f_mass(iv2)  &
             /(pi*f_fy*f_fp*f_diam(iv2)**2._rsh*(f_mass(iv1)+f_mass(iv2))         &
             *((f_rho(iv2)-rhoref)/rhoref)**(2._rsh/(3._rsh-f_nf)))    

        mult=1._rsh
        if (iv1.eq.iv2) mult=2._rsh

        if (iv1.eq.MAX(iv1,iv2) .and. gcolfragiv1.ge.1._rsh) then
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2))
        elseif (iv1.eq.MIN(iv1,iv2) .and. gcolfragiv2.ge.1._rsh) then
           f_l4(iv2,iv1)=f_l4(iv2,iv1)+mult*(f_coll_prob_sh(iv1,iv2))
        endif

     ENDDO
  ENDDO

  f_g4(1:nv_mud,1:nv_mud,1:nv_mud)=f_g4(1:nv_mud,1:nv_mud,1:nv_mud)*f_collfragparam
  f_l4(1:nv_mud,1:nv_mud)=f_l4(1:nv_mud,1:nv_mud)*f_collfragparam

  write(*,*) 'END flocmod_collfrag'

END SUBROUTINE flocmod_collfrag

!!===========================================================================
SUBROUTINE flocmod_comp_fsd(NNin,NNout,Gval)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_fsd  ***
  !&E
  !&E ** Purpose : computation of floc size distribution
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars
  ! ,    ONLY : limin,limax,ljmin,ljmax,kmax,rhoref

  !! * Arguments
  REAL(KIND=rsh),INTENT(in) :: Gval
  REAL(KIND=rsh),DIMENSION(1:nv_mud),INTENT(in)  :: NNin
  REAL(KIND=rsh),DIMENSION(1:nv_mud),INTENT(out) :: NNout

  !! * Local declarations
  INTEGER      :: iv1,iv2,iv3
  REAL(KIND=rsh) :: tmp_g1,tmp_g3,tmp_l1,tmp_l3,tmp_l4,tmp_g4
  REAL(KIND=rsh),DIMENSION(1:nv_mud,1:nv_mud,1:nv_mud)     :: f_g1_tmp
  REAL(KIND=rsh),DIMENSION(1:nv_mud,1:nv_mud)              :: f_l1_tmp

  !!--------------------------------------------------------------------------
  !! * Executable part

  tmp_g1=0.0_rsh
  tmp_g3=0.0_rsh
  tmp_g4=0.0_rsh
  tmp_l1=0.0_rsh
  tmp_l3=0.0_rsh
  tmp_l4=0.0_rsh    
  f_g1_tmp(1:nv_mud,1:nv_mud,1:nv_mud)=0.0_rsh
  f_l1_tmp(1:nv_mud,1:nv_mud)=0.0_rsh

  if (l_COLLFRAG) CALL flocmod_collfrag(Gval)

  DO iv1=1,nv_mud
     DO iv2=1,nv_mud
        DO iv3=1,nv_mud
           if (l_ASH) then
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_sh(iv2,iv3,iv1)*Gval   
           endif
           if (l_ADS) then
              f_g1_tmp(iv2,iv3,iv1)=f_g1_tmp(iv2,iv3,iv1)+f_g1_ds(iv2,iv3,iv1)   
           endif

           tmp_g1=tmp_g1+(NNin(iv3)*(f_g1_tmp(iv2,iv3,iv1))*NNin(iv2))

           if (l_COLLFRAG) then
              tmp_g4=tmp_g4+(NNin(iv3)*(f_g4(iv2,iv3,iv1)*Gval)*NNin(iv2))
           endif
        ENDDO

        tmp_g3=tmp_g3+f_g3(iv2,iv1)*NNin(iv2)*Gval**1.5_rsh

        if (l_ASH) then
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_sh(iv2,iv1)*Gval
        endif
        if (l_ADS) then
           f_l1_tmp(iv2,iv1)=f_l1_tmp(iv2,iv1)+f_l1_ds(iv2,iv1)*Gval
        endif

        tmp_l1=tmp_l1+(f_l1_tmp(iv2,iv1))*NNin(iv2)

        if (l_COLLFRAG) then
           tmp_l4=tmp_l4+(f_l4(iv2,iv1)*Gval)*NNin(iv2)
        endif
     ENDDO

     tmp_l1=tmp_l1*NNin(iv1)
     tmp_l4=tmp_l4*NNin(iv1)

     tmp_l3=f_l3(iv1)*Gval**1.5_rsh*NNin(iv1)

     !    write(*,*) 'iv1',(iv1)
     !    write(*,*) 'NNin',NNin(iv1)
     !    write(*,*) 'tmpg1',tmp_g1
     !    write(*,*) 'tmpg3',tmp_g3
     !    write(*,*) 'tmpl1',tmp_l1
     !    write(*,*) 'tmpl3',tmp_l3
     !    write(*,*) 'tmpl3',f_l3(iv1)*Gval**1.5_rsh

     NNout(iv1)=NNin(iv1)+f_dt*(tmp_g1+tmp_g3+tmp_g4-(tmp_l1+tmp_l3+tmp_l4))

     tmp_g1=0.0_rsh
     tmp_g3=0.0_rsh
     tmp_g4=0.0_rsh
     tmp_l1=0.0_rsh
     tmp_l3=0.0_rsh
     tmp_l4=0.0_rsh    
  ENDDO
  !write(*,*) 'END flocmod_comp_fsd'

END SUBROUTINE flocmod_comp_fsd



!!===========================================================================
SUBROUTINE flocmod_mass_control(NN,mneg)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_control  ***
  !&E
  !&E ** Purpose : Compute mass in every class after flocculation and returns negative mass if any
  !&E
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars

  !! * Arguments


  !! * Local declarations
  INTEGER      :: iv1
  REAL(KIND=rsh),intent(out)     :: mneg
  REAL(KIND=rsh),DIMENSION(1:nv_mud),intent(in)     :: NN
  !!--------------------------------------------------------------------------
  !! * Executable part

  mneg=0.0_rsh

  DO iv1=1,nv_mud
     if (NN(iv1).lt.0.0_rsh) then
        mneg=mneg-NN(iv1)*f_mass(iv1)
     endif
  ENDDO
  !write(*,*) 'END flocmod_mass_control'

END SUBROUTINE flocmod_mass_control

!!===========================================================================
SUBROUTINE flocmod_mass_redistribute(NN)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_mass_redistribute  ***
  !&E
  !&E ** Purpose : based on a tolerated negative mass parameter, negative masses  
  !&E              are redistributed linearly towards remaining postive masses 
  !&E              and negative masses are set to 0
  !&E                   
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars, ONLY : rsh, nv_mud, f_mass
  !! * Arguments


  !! * Local declarations
  INTEGER      :: iv
  REAL(KIND=rsh)     :: npos
  REAL(KIND=rsh)     :: mneg
  REAL(KIND=rsh),DIMENSION(1:nv_mud),intent(inout)     :: NN
  REAL(KIND=rsh),DIMENSION(1:nv_mud)                   :: NNtmp
  !!--------------------------------------------------------------------------
  !! * Executable part

  mneg=0.0_rsh
  npos=0.0_rsh
  NNtmp(:)=NN(:)

  DO iv=1,nv_mud
     if (NN(iv).lt.0.0_rsh) then
        mneg=mneg-NN(iv)*f_mass(iv)
        NNtmp(iv)=0.0_rsh
     else
        npos=npos+1._rsh
     endif
  ENDDO

  if (mneg.gt.0.0_rsh) then 
     if (npos.eq.0._rsh) then
        write(*,*) 'CAUTION : all floc sizes have negative mass!'
        write(*,*) 'SIMULATION STOPPED'
        STOP    
     else
        DO iv=1,nv_mud
           if (NN(iv).gt.0.0_rsh) then
              NN(iv)=NN(iv)-mneg/sum(NNtmp)*NN(iv)/f_mass(iv)
           else
              NN(iv)=0.0_rsh
           endif

        ENDDO

     endif
  endif
  !write(*,*) 'END flocmod_mass_redistribute'

END SUBROUTINE flocmod_mass_redistribute

!!===========================================================================    
SUBROUTINE flocmod_comp_g(k,i,j,Gval)

  !&E--------------------------------------------------------------------------
  !&E                 ***  ROUTINE flocmod_comp_g  ***
  !&E
  !&E ** Purpose : compute shear rate to estimate shear aggregation and erosion  
  !&E 
  !&E ** Description :
  !&E
  !&E ** Called by : flocmod_main
  !&E
  !&E ** External calls : 
  !&E
  !&E ** Reference :
  !&E
  !&E ** History :
  !&E     ! 2013-09 (Romaric Verney)
  !&E
  !&E--------------------------------------------------------------------------
  !! * Modules used
  USE comvars,    ONLY : rsh, t, l_testcase
  !! * Arguments


  !! * Local declarations
  INTEGER,  intent(in)      :: k,i,j
  REAL(KIND=rsh),intent(out)     :: Gval
  REAL(KIND=rsh)     :: htn,ustar,z,diss,nueau
  !!--------------------------------------------------------------------------
  !! * Executable part
  !    nueau=1.06e-6_rsh
  !    ustar=sqrt(tenfon(i,j)/rhoref)
  !    htn=h0(i,j)+ssh(i,j)
  !    z=(1.0_rsh+sig(k))*htn

  if (l_testcase) then
     ! reproducing flocculation epxeriment Verney et al., 2011
     Gval=0.0_rsh
     if (t .lt. 7201._rsh) then
        Gval=1._rsh
     elseif (t .lt. 8401._rsh) then
        Gval=2._rsh
     elseif (t .lt. 9601._rsh) then
        Gval=3._rsh  
     elseif (t .lt. 10801._rsh) then
        Gval=4._rsh
     elseif (t .lt. 12601._rsh) then
        Gval=12._rsh
     elseif (t .lt. 13801._rsh) then
        Gval=4._rsh
     elseif (t .lt. 15001._rsh) then
        Gval=3._rsh
     elseif (t .lt. 16201._rsh) then
        Gval=2._rsh
     elseif (t .lt. 21601._rsh) then
        Gval=1._rsh
     elseif (t .lt. 25201._rsh) then
        Gval=0._rsh
     elseif (t .lt. 30601._rsh) then
        Gval=1._rsh
     elseif (t .lt. 31801._rsh) then
        Gval=2._rsh                     
     elseif (t .lt. 33001._rsh) then
        Gval=3._rsh       
     elseif (t .lt. 34201._rsh) then
        Gval=4._rsh
     elseif (t .lt. 36001._rsh) then
        Gval=12._rsh
     elseif (t .lt. 37201._rsh) then
        Gval=4._rsh
     elseif (t .lt. 38401._rsh) then
        Gval=3._rsh
     elseif (t .lt. 39601._rsh) then
        Gval=2._rsh
     elseif (t .lt. 45001._rsh) then
        Gval=1._rsh
     elseif (t .lt. 48601._rsh) then
        Gval=0._rsh
     elseif (t .lt. 54001._rsh) then
        Gval=1._rsh
     elseif (t .lt. 55201._rsh) then
        Gval=2._rsh                     
     elseif (t .lt. 56401._rsh) then
        Gval=3._rsh       
     elseif (t .lt. 57601._rsh) then
        Gval=4._rsh
     else 
        Gval=12._rsh
     endif
  else
     write(*,*) 'not calculating Gval'
     Gval = 4._rsh !CRS
     !       if (htn.lt.0.05_rsh) then
     !          Gval=0.0_rsh
     !       else
     !          diss=ustar**3._rsh/(0.41_rsh*htn)*(1._rsh-z)/z
     !          Gval=sqrt(diss/nueau)
     !       endif
  endif
  ! write(*,*) 'END flocmod_comp_g'
END SUBROUTINE flocmod_comp_g
