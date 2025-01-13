!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef OASIS

module mod_oasis_interface

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_stdio
  use mod_message
  use mod_service
  use mod_mppparam
  use mod_runparams ! , only : isocean , dtsec , alarm_out_sav , rcmtimer
  use mod_bats_common , only : rdnnsg
!  use mod_atm_interface , only : atms , sfs , flwd , flw , fsw , sinc , mddom
 ! use mod_lm_interface , only : lms , lm
  use mod_date , only : lfdomonth , lmidnight
  use mod_memutil  
  use mod_regcm_types

  use mod_oasis
  use mod_oasis_params
  use mod_oasis_signature
  use mod_oasis_generic

!! IM ADDIND EVERYTHING RIGHTNOW
  use mod_date
  use mod_stdio
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_kdinterp
  use mod_runparams
  use mod_ipcc_scenario , only : ghgval , igh_co2 , igh_ch4
  use mod_memutil
  use mod_mppparam
  use mod_ensemble
  use mod_mpmessage
  use mod_constants
  use mod_sunorbit
  use mod_regcm_types
  use mod_service
  use mod_lm_interface 
  use mod_atm_interface
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_outvars
  use mod_constants
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_bats_common
  use mod_ocn_common
  use mod_stdio
  use mod_slabocean
  use mod_heatindex
  use mod_dynparam
  use mod_stdio
  use mod_date
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_memutil
  use mod_regcm_types
  use mod_stdatm
  use mod_zita
  use mod_ncio


  implicit none

  private

  public :: comp_name , comp_id ! -> mod_oasis_params
  public :: oasis_lag           ! -> mod_oasis_params
  real(rkx) , dimension(:,:) , pointer :: temps , temps_rf , temps_snw, temp, temp2
  real(rkx) , dimension(:,:,:) , pointer :: temp3d, temp3d_2
  real(rkx) , dimension(:,:,:,:) , pointer :: temp4d
!  real(rkx), dimension(:,:,:), allocatable :: lon3d, lat3d
!  real(rkx), dimension(:,:,:), pointer :: vocemiss_slice
  integer :: k, kdim,n  !--------------------------------------------------------------------
  ! for grid and partition activation
  ! (dot/cross , external/internal)
  logical :: l_make_grdde , l_make_grddi , &
             l_make_grdce , l_make_grdci , l_make_grdtest , l_make_grd3d
  type(infogrd) , target , allocatable :: grdde , grddi , &
                                          grdce , grdci , grdtest 
  type(info3dgrd) , target , allocatable :: grd3d
  !--------------------------------------------------------------------
  ! below the namelist parameters (namelist oasisparam)
  public :: l_write_grids , write_restart_option ! -> mod_oasis_params
  public :: oasis_sync_lag                       ! -> mod_oasis_params
  logical , public :: l_cpl_im_sst , &
!                      l_cpl_im_sit , & ! not coded
                      l_cpl_im_wz0 , & ! not tested
                      l_cpl_im_wust , & ! not tested
                      l_cpl_ex_u10m , & ! not tested
                      l_cpl_ex_v10m , & ! not tested
                      l_cpl_ex_wspd , & ! not tested
                      l_cpl_ex_wdir , & ! not tested
                      l_cpl_ex_t2m , & ! not tested
!                      l_cpl_ex_t10m , & ! not coded
                      l_cpl_ex_q2m , & ! not tested
!                      l_cpl_ex_q10m , & ! not coded
                      l_cpl_ex_slp , &
                      l_cpl_ex_taux , &
                      l_cpl_ex_tauy , &
                      l_cpl_ex_z0 , & ! not tested
                      l_cpl_ex_ustr , & ! not tested
                      l_cpl_ex_evap , &
                      l_cpl_ex_prec , &
                      l_cpl_ex_nuwa , &
                      l_cpl_ex_ulhf , &
                      l_cpl_ex_ushf , &
                      l_cpl_ex_uwlw , &
                      l_cpl_ex_dwlw , &
                      l_cpl_ex_nulw , &
                      l_cpl_ex_uwsw , &
                      l_cpl_ex_dwsw , &
                      l_cpl_ex_ndsw , &
                  !    l_cpl_ex_pr
!!!! ECLM
                      l_cpl_ex_tatm, &
                      l_cpl_ex_uatm, &
                      l_cpl_ex_vatm, &
                      l_cpl_ex_qatm, &
                      l_cpl_ex_hgt, &
                      l_cpl_ex_patm, &
                      l_cpl_ex_dwrlwf, &
                      l_cpl_ex_swdir, &
                      l_cpl_ex_swdif, &
                      l_cpl_im_tgbb, &
                      l_cpl_im_t2m, &
                      l_cpl_im_q2m, &
                      l_cpl_im_u10m, &
                      l_cpl_im_sent, &
                      l_cpl_im_evpr, &
                      l_cpl_im_ram1, &
                      l_cpl_im_rah1, &
                      l_cpl_im_tauy, &
                      l_cpl_im_taux, &
                      l_cpl_im_sncv, &
                      l_cpl_im_wt, &
                      l_cpl_im_zo, &
                      l_cpl_im_tlef, &
                      l_cpl_ex_rainf, &
                      l_cpl_ex_snow, &
                      l_cpl_im_albei, &
                      l_cpl_im_albed, &
                      l_cpl_im_flxvoc, &
                      l_cpl_im_flxdst, &
                      l_cpl_im_fluxch4, &
                      l_cpl_im_ddvel, &
                      l_cpl_im_tlai, &
                      l_cpl_im_roff, &
                      l_cpl_im_srfroff, &
                      l_cpl_im_snwmelt, &
                      l_cpl_im_tsoi, &
                      l_cpl_im_h2ovol, &
                      l_cpl_im_h2oliq, &
                      l_cpl_im_h2oice, &
                      l_cpl_im_h2o10cm, &
                      l_cpl_ex_rhoa 
                                        
  ! OASIS field +++
! l_cpl_im_wust1
  !--------------------------------------------------------------------
  ! below the variables for the fields' information
  type(infofld) , allocatable :: im_sst , &
!                                 im_sit , &
                                 im_wz0 , &
                                 im_wust , &
                                 ex_u10m , &
                                 ex_v10m , &
                                 ex_wspd , &
                                 ex_wdir , &
                                 ex_t2m , &
!                                 ex_t10m , &
                                 ex_q2m , &
!                                 ex_q10m , &
                                 ex_slp , &
                                 ex_taux , &
                                 ex_tauy , &
                                 ex_z0 , &
                                 ex_ustr , &
                                 ex_evap , &
                                 ex_prec , &
                                 ex_nuwa , &
                                 ex_ulhf , &
                                 ex_ushf , &
                                 ex_uwlw , &
                                 ex_dwlw , &
                                 ex_nulw , &
                                 ex_uwsw , &
                                 ex_dwsw , &
                                 ex_ndsw , &
  ! OASIS field +++
  !!! ECLM

#ifdef ECLM
                                 ex_tatm , & 
                                 ex_uatm , & 
                                 ex_vatm , & 
                                 ex_qatm , & 
                                 ex_hgt , & 
                                 ex_patm , & 
                                 ex_dwrlwf , & 
                                 ex_swdir , & 
                                 ex_swdif , & 
                                 im_tgbb , & 
                                 im_t2m , & 
                                 im_q2m , & 
                                 im_u10m , & 
  !                               im_v10m , & 
                                 im_sent , & 
                                 im_evpr , & 
                                 im_ram1 , & 
                                 im_rah1 , & 
 !                                im_br , & 
                                 im_taux , & 
                                 im_tauy , & 
!                                 im_drag , & 
                                 im_sncv , & 
                                 im_wt , & 
                                 im_zo , & 
                                 im_tlef , & 
                                 ex_rainf , & 
                                 ex_snow , & 
!                                 im_albei , & 
 !                                im_albed , & 
                                 im_fluxch4 , & 
                                 im_tlai , & 
                                 im_roff , & 
                                 im_srfroff , & 
                                 im_snwmelt , & 
                                 im_h2o10cm , &
#endif
                                 ex_rhoa 
                               
#ifdef ECLM                            
  type(info3dfld) , allocatable :: im_albei , &
                                 im_flxdst , & 
                                 im_flxvoc , & 
                                 im_ddvel , & 
                                 im_h2oliq , & 
                                 im_h2oice , & 
                                 im_h2ovol , & 
                                 im_tsoi , & 
                                 im_albed 
#endif
                               ! OASIS needs double precision arrays. These variables are optional.
  ! Required for the fields to import; recommended for the fields to
  ! export which need a bit of work from the model variables.
  real(rkx) , dimension(:,:) , allocatable :: cpl_sst , &
!                                              cpl_sit , &
                                              cpl_wz0 , &
                                              cpl_wust , &
                                              cpl_wdir , &
#ifdef ECLM
                                              cpl_t2m , &
                                              cpl_q2m , &
                                              cpl_u10m , &
 !                                             cpl_v10m , &
                                              cpl_sent , &
                                              cpl_evpr , &
                                              cpl_ram1 , &
                                              cpl_rah1 , &
  !                                            cpl_br , &
                                              cpl_taux , &
                                              cpl_tauy , &
!                                              cpl_drag , &
                                              cpl_sncv , &
                                              cpl_wt , &
                                              cpl_zo , &
                                              cpl_tlef , &
                                              !cpl_albei , &
                                              !cpl_albed , &
                                              cpl_fluxch4 , &
                                              cpl_tlai , &
                                              cpl_roff , &
                                              cpl_srfroff , &
                                              cpl_snwmelt , &
                                              cpl_h2o10cm, & 
#endif
                                              cpl_tgbb 
#ifdef ECLM

                                            ! OASIS field +++
  real(rkx) , dimension(:,:,:) , allocatable :: cpl_albei , &
                                              cpl_flxdst , &
                                              cpl_flxvoc , &
                                              cpl_h2oliq , &
                                              cpl_ddvel , &
                                              cpl_h2oice , &
                                              cpl_h2ovol , &
                                              cpl_tsoi , &
                                              cpl_albed 
#endif

  ! OASIS needs double precision arrays. These variables are optional.
  ! They are used when calling the subroutine oasisxregcm_setup_field_array(),
  ! when an array is to be initialized. If not given, the array will be
  ! initialized with zeros.
  real(rkx) , parameter :: init_sst = tzero + 25.0
  ! OASIS field +++

  !--------------------------------------------------------------------
  ! accessibility from oasislib:
  public :: oasisxregcm_init , oasisxregcm_finalize ! -> mod_oasis_generic
  public :: oasisxregcm_header                      ! -> mod_oasis_signature

  !--------------------------------------------------------------------
  ! from this module:
  public :: oasisxregcm_params , oasisxregcm_release
  public :: oasisxregcm_def
  public :: oasisxregcm_rcv_all , oasisxregcm_snd_all
  public :: oasisxregcm_sync_wait

  contains

  subroutine oasisxregcm_params
    implicit none
    !--------------------------------------------------------------------------
    oasis_lag = 0
    ! state which grids have to be defined,
    ! depending on what field is activated.
    l_make_grdde   = .false.
    ! OASIS field +++
    l_make_grddi   = .false.
    ! OASIS field +++
    l_make_grdce   = l_cpl_ex_slp
    ! OASIS field +++
#ifdef ECLM
    l_make_grdtest = l_cpl_ex_uatm .or. &
                      l_cpl_ex_vatm .or. &
                      l_cpl_ex_qatm .or. &
                      l_cpl_ex_hgt .or. &
                      l_cpl_ex_patm .or. &
                      l_cpl_ex_dwrlwf .or. &
                      l_cpl_ex_swdir .or. &
                      l_cpl_ex_tatm .or. & 
                      l_cpl_ex_swdif .or. &
                      l_cpl_im_tgbb .or. & 
                      l_cpl_im_t2m .or. & 
                      l_cpl_im_q2m .or. & 
                      l_cpl_im_u10m .or. & 
                      l_cpl_im_sent .or. & 
                      l_cpl_im_evpr .or. & 
                      l_cpl_im_ram1 .or. & 
                      l_cpl_im_rah1 .or. & 
                      l_cpl_im_tauy .or. & 
                      l_cpl_im_taux .or. & 
                      l_cpl_im_sncv .or. & 
                      l_cpl_im_wt .or. & 
                      l_cpl_im_zo .or. & 
                      l_cpl_im_tlef .or. &  
                      l_cpl_ex_rainf .or. &
                      l_cpl_ex_snow .or. & 
                   !   l_cpl_im_albei .or. & 
                     ! l_cpl_im_albed .or. & 
                      l_cpl_im_tlai .or. & 
                      l_cpl_im_roff .or. & 
                      l_cpl_im_srfroff .or. & 
                      l_cpl_im_snwmelt .or. & 
                      l_cpl_im_h2o10cm 
    l_make_grd3d =  l_cpl_im_albei .or. & 
                      l_cpl_im_flxdst .or. & 
                      l_cpl_im_flxvoc .or. & 
                      l_cpl_im_ddvel .or. & 
                      l_cpl_im_h2oliq .or. & 
                      l_cpl_im_h2oice .or. & 
                      l_cpl_im_tsoi .or. & 
                      l_cpl_im_h2ovol .or. & 
                    l_cpl_im_albed 
#endif
    l_make_grdci   = l_cpl_im_sst  .or. &
!                     l_cpl_im_sit  .or. &
                     l_cpl_im_wz0  .or. &
                     l_cpl_im_wust .or. &
                     l_cpl_ex_u10m .or. &
                     l_cpl_ex_v10m .or. &
                     l_cpl_ex_wspd .or. &
                     l_cpl_ex_wdir .or. &
                     l_cpl_ex_t2m  .or. &
!                     l_cpl_ex_t10m .or. &
                     l_cpl_ex_q2m  .or. &
!                     l_cpl_ex_q10m .or. &
                     l_cpl_ex_taux .or. &
                     l_cpl_ex_tauy .or. &
                     l_cpl_ex_z0   .or. &
                     l_cpl_ex_ustr .or. &
                     l_cpl_ex_evap .or. &
                     l_cpl_ex_prec .or. &
                     l_cpl_ex_nuwa .or. &
                     l_cpl_ex_ulhf .or. &
                     l_cpl_ex_ushf .or. &
                     l_cpl_ex_uwlw .or. &
                     l_cpl_ex_dwlw .or. &
                     l_cpl_ex_nulw .or. &
                     l_cpl_ex_uwsw .or. &
                     l_cpl_ex_dwsw .or. &
                     l_cpl_ex_ndsw .or. &
                     l_cpl_ex_rhoa
    ! OASIS field +++
    !
    ! initialize grids: grid variable, name with no mask, name with mask,
    !                   local j1, j2, i1, i2, global j1, j2, i1, i2, number of
    !                   corners per cell
    if ( l_make_grdde )   call oasisxregcm_setup_grid(grdde, 'rden', 'rdem', &
                               jde1, jde2, ide1, ide2, 1, jx, 1, iy, 4)
    if ( l_make_grddi )   call oasisxregcm_setup_grid(grddi, 'rdin', 'rdim', &
                               jdi1, jdi2, idi1, idi2, 2, jx-1, 2, iy-1, 4)
    if ( l_make_grdce )   call oasisxregcm_setup_grid(grdce, 'rcen', 'rcem', &
                               jce1, jce2, ice1, ice2, 1, jx-1, 1, iy-1, 4)
    if ( l_make_grdci )   call oasisxregcm_setup_grid(grdci, 'rcin', 'rcim', &
                               jci1, jci2, ici1, ici2, 2, jx-2, 2, iy-2, 4)

#ifdef ECLM
    if ( l_make_grdtest )   call oasisxregcm_setup_grid(grdtest, 'rcin', 'rcim', &
                               jci1, jci2, ici1, ici2, 2, jx-2, 2, iy-2, 4)
    if ( l_make_grd3d )   call oasisxregcm_setup_3dgrid(grd3d, 'rcin', 'rcim', &
                               jci1, jci2, ici1, ici2, 1, nnsg, &
                               2, jx-2, 2, iy-2, 1, nnsg, 4)
#endif
                             !nnsg
                             !
    ! initialize fields: field variable, name, grid, field array (optional), initialization value
    !                                                                        (optional; 0 otherwise)
    if ( l_cpl_im_sst )  call oasisxregcm_setup_field(im_sst,  'RCM_SST',  grdci, cpl_sst, init_sst)
!    if ( l_cpl_im_sit )  call oasisxregcm_setup_field(im_sit,  'RCM_SIT',  grdci, cpl_sit)
    if ( l_cpl_im_wz0 )  call oasisxregcm_setup_field(im_wz0,  'RCM_WZ0',  grdci, cpl_wz0)
    if ( l_cpl_im_wust ) call oasisxregcm_setup_field(im_wust, 'RCM_WUST', grdci, cpl_wust)
    if ( l_cpl_ex_u10m ) call oasisxregcm_setup_field(ex_u10m, 'RCM_U10M', grdci)
    if ( l_cpl_ex_v10m ) call oasisxregcm_setup_field(ex_v10m, 'RCM_V10M', grdci)
    if ( l_cpl_ex_wspd ) call oasisxregcm_setup_field(ex_wspd, 'RCM_WSPD', grdci)
    if ( l_cpl_ex_wdir ) call oasisxregcm_setup_field(ex_wdir, 'RCM_WDIR', grdci, cpl_wdir)
    if ( l_cpl_ex_t2m )  call oasisxregcm_setup_field(ex_t2m,  'RCM_T2M',  grdci)
!    if ( l_cpl_ex_t10m ) call oasisxregcm_setup_field(ex_t10m, 'RCM_T10M', grdci)
    if ( l_cpl_ex_q2m )  call oasisxregcm_setup_field(ex_q2m,  'RCM_Q2M',  grdci)
!    if ( l_cpl_ex_q10m ) call oasisxregcm_setup_field(ex_q10m, 'RCM_Q10M', grdci)
    if ( l_cpl_ex_slp )  call oasisxregcm_setup_field(ex_slp,  'RCM_SLP',  grdce)
    if ( l_cpl_ex_taux ) call oasisxregcm_setup_field(ex_taux, 'RCM_TAUX', grdci)
    if ( l_cpl_ex_tauy ) call oasisxregcm_setup_field(ex_tauy, 'RCM_TAUY', grdci)
    if ( l_cpl_ex_z0 )   call oasisxregcm_setup_field(ex_z0,   'RCM_Z0',   grdci)
    if ( l_cpl_ex_ustr ) call oasisxregcm_setup_field(ex_ustr, 'RCM_USTR', grdci)
    if ( l_cpl_ex_evap ) call oasisxregcm_setup_field(ex_evap, 'RCM_EVAP', grdci)
    if ( l_cpl_ex_prec ) call oasisxregcm_setup_field(ex_prec, 'RCM_PREC', grdci)
    if ( l_cpl_ex_nuwa ) call oasisxregcm_setup_field(ex_nuwa, 'RCM_NUWA', grdci)
    if ( l_cpl_ex_ulhf ) call oasisxregcm_setup_field(ex_ulhf, 'RCM_ULHF', grdci)
    if ( l_cpl_ex_ushf ) call oasisxregcm_setup_field(ex_ushf, 'RCM_USHF', grdci)
    if ( l_cpl_ex_uwlw ) call oasisxregcm_setup_field(ex_uwlw, 'RCM_UWLW', grdci)
    if ( l_cpl_ex_dwlw ) call oasisxregcm_setup_field(ex_dwlw, 'RCM_DWLW', grdci)
    if ( l_cpl_ex_nulw ) call oasisxregcm_setup_field(ex_nulw, 'RCM_NULW', grdci)
    if ( l_cpl_ex_uwsw ) call oasisxregcm_setup_field(ex_uwsw, 'RCM_UWSW', grdci)
    if ( l_cpl_ex_dwsw ) call oasisxregcm_setup_field(ex_dwsw, 'RCM_DWSW', grdci)
    if ( l_cpl_ex_ndsw ) call oasisxregcm_setup_field(ex_ndsw, 'RCM_NDSW', grdci)
    if ( l_cpl_ex_rhoa ) call oasisxregcm_setup_field(ex_rhoa, 'RCM_RHOA', grdci)
#ifdef ECLM
    if ( l_cpl_ex_tatm ) call oasisxregcm_setup_field(ex_tatm, 'RCM_TATM', grdtest)
    if ( l_cpl_ex_uatm ) call oasisxregcm_setup_field(ex_uatm, 'RCM_UATM', grdtest)
    if ( l_cpl_ex_vatm ) call oasisxregcm_setup_field(ex_vatm, 'RCM_VATM', grdtest)
    if ( l_cpl_ex_qatm ) call oasisxregcm_setup_field(ex_qatm, 'RCM_QATM', grdtest)
    if ( l_cpl_ex_hgt ) call oasisxregcm_setup_field(ex_hgt, 'RCM_HGT', grdtest)
    if ( l_cpl_ex_patm ) call oasisxregcm_setup_field(ex_patm, 'RCM_PATM', grdtest)
    if ( l_cpl_ex_dwrlwf ) call oasisxregcm_setup_field(ex_dwrlwf, 'RCM_DWRLWF', grdtest)
    if ( l_cpl_ex_swdir ) call oasisxregcm_setup_field(ex_swdir, 'RCM_SWDIR', grdtest)
    if ( l_cpl_ex_swdif ) call oasisxregcm_setup_field(ex_swdif, 'RCM_SWDIF', grdtest)
    if ( l_cpl_im_tgbb )  call oasisxregcm_setup_field(im_tgbb,  'RCM_TGBB',grdtest, cpl_tgbb)
    if ( l_cpl_im_t2m )  call oasisxregcm_setup_field(im_t2m,  'RCM_T2M',grdtest, cpl_t2m)
    if ( l_cpl_im_q2m )  call oasisxregcm_setup_field(im_q2m,  'RCM_Q2M',grdtest, cpl_q2m)
    if ( l_cpl_im_u10m )  call oasisxregcm_setup_field(im_u10m,  'RCM_U10M',grdtest, cpl_u10m)
!    if ( l_cpl_im_v10m )  call oasisxregcm_setup_field(im_v10m,  'RCM_V10M',grdtest, cpl_v10m)
    if ( l_cpl_im_sent )  call oasisxregcm_setup_field(im_sent,  'RCM_SENT',grdtest, cpl_sent)
    if ( l_cpl_im_evpr )  call oasisxregcm_setup_field(im_evpr,  'RCM_EVPR',grdtest, cpl_evpr)
    if ( l_cpl_im_ram1 )  call oasisxregcm_setup_field(im_ram1,  'RCM_RAM1',grdtest, cpl_ram1)
    if ( l_cpl_im_rah1 )  call oasisxregcm_setup_field(im_rah1,  'RCM_RAH1',grdtest, cpl_rah1)
!    if ( l_cpl_im_br )  call oasisxregcm_setup_field(im_br,  'RCM_BR',grdtest, cpl_br)
    if ( l_cpl_im_taux )  call oasisxregcm_setup_field(im_taux,  'RCM_TAUX',grdtest, cpl_taux)
    if ( l_cpl_im_tauy )  call oasisxregcm_setup_field(im_tauy,  'RCM_TAUY',grdtest, cpl_tauy)
!    if ( l_cpl_im_drag )  call oasisxregcm_setup_field(im_drag,  'RCM_DRAG',grdtest, cpl_drag)
    if ( l_cpl_im_sncv )  call oasisxregcm_setup_field(im_sncv,  'RCM_SNCV',grdtest, cpl_sncv)
    if ( l_cpl_im_wt )  call oasisxregcm_setup_field(im_wt,  'RCM_WT',grdtest, cpl_wt)
    if ( l_cpl_im_zo )  call oasisxregcm_setup_field(im_zo,  'RCM_ZO',grdtest, cpl_zo)
    if ( l_cpl_im_tlef )  call oasisxregcm_setup_field(im_tlef,  'RCM_TLEF',grdtest, cpl_tlef)
    if ( l_cpl_ex_rainf ) call oasisxregcm_setup_field(ex_rainf, 'RCM_RAINF', grdtest)
    if ( l_cpl_ex_snow ) call oasisxregcm_setup_field(ex_snow, 'RCM_SNOW', grdtest)
    if ( l_cpl_im_albei )  call oasisxregcm_setup_3dfield(im_albei,  'RCM_ALBEI',grd3d, cpl_albei)
    if ( l_cpl_im_albed )  call oasisxregcm_setup_3dfield(im_albed,  'RCM_ALBED',grd3d, cpl_albed)
    if ( l_cpl_im_flxvoc )  call oasisxregcm_setup_3dfield(im_flxvoc,  'RCM_FLXVOC',grd3d, cpl_flxvoc)
    if ( l_cpl_im_flxdst )  call oasisxregcm_setup_3dfield(im_flxdst,  'RCM_FLXDST',grd3d, cpl_flxdst)
    if ( l_cpl_im_fluxch4 )  call oasisxregcm_setup_field(im_fluxch4,  'RCM_FLUXCH4',grdtest, cpl_fluxch4)
    if ( l_cpl_im_ddvel )  call oasisxregcm_setup_3dfield(im_ddvel,  'RCM_DDVEL',grd3d, cpl_ddvel)
    if ( l_cpl_im_tlai )  call oasisxregcm_setup_field(im_tlai,  'RCM_TLAI',grdtest, cpl_tlai)
    if ( l_cpl_im_roff )  call oasisxregcm_setup_field(im_roff,  'RCM_ROFF',grdtest, cpl_roff)
    if ( l_cpl_im_srfroff )  call oasisxregcm_setup_field(im_srfroff,  'RCM_SRFROFF',grdtest, cpl_srfroff)
    if ( l_cpl_im_snwmelt )  call oasisxregcm_setup_field(im_snwmelt,  'RCM_SNWMELT',grdtest, cpl_snwmelt)
    if ( l_cpl_im_tsoi )  call oasisxregcm_setup_3dfield(im_tsoi,  'RCM_TSOI',grd3d, cpl_tsoi)
    if ( l_cpl_im_h2ovol )  call oasisxregcm_setup_3dfield(im_h2ovol,  'RCM_H2OVOL',grd3d, cpl_h2ovol)
    if ( l_cpl_im_h2oliq )  call oasisxregcm_setup_3dfield(im_h2oliq,  'RCM_H2OLIQ',grd3d, cpl_h2oliq)
    if ( l_cpl_im_h2oice )  call oasisxregcm_setup_3dfield(im_h2oice,  'RCM_H2OICE',grd3d, cpl_h2oice)
    if ( l_cpl_im_h2o10cm )  call oasisxregcm_setup_field(im_h2o10cm,  'RCM_H2O10CM',grdtest, cpl_h2o10cm)
#endif
    ! OASIS field +++
  end subroutine oasisxregcm_params

  ! call all definition subroutines for setting up OASIS
  subroutine oasisxregcm_def
    implicit none
    !--------------------------------------------------------------------------
    ! partition definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: partitions'
#endif
    if ( l_make_grdde ) call oasisxregcm_def_partition(grdde)
    if ( l_make_grddi ) call oasisxregcm_def_partition(grddi)
    if ( l_make_grdce ) call oasisxregcm_def_partition(grdce)
    if ( l_make_grdci ) call oasisxregcm_def_partition(grdci)
#ifdef ECLM
    if ( l_make_grdtest ) call oasisxregcm_def_partition(grdtest)
    if ( l_make_grd3d ) call oasisxregcm_def_partition_3d(grd3d)
#endif    
! grid definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: grids'
#endif
    if ( l_write_grids ) then
      call oasisxregcm_def_grid
    else
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, '               >> skipped!'
#endif
      if ( myid == italk ) write(stdout,*) 'Note: OASIS grid files', &
                                   ' have not been written online.'
    end if
    ! variable definitions
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: fields'
#endif
    ! field variable, OASIS integer for the direction:
    !                 importing(OASIS_In)/exporting(OASIS_Out)
    if ( l_cpl_im_sst )  call oasisxregcm_def_field(im_sst,  OASIS_In)
!    if ( l_cpl_im_sit )  call oasisxregcm_def_field(im_sit,  OASIS_In)
    if ( l_cpl_im_wz0 )  call oasisxregcm_def_field(im_wz0,  OASIS_In)
    if ( l_cpl_im_wust ) call oasisxregcm_def_field(im_wust, OASIS_In)
    if ( l_cpl_ex_u10m ) call oasisxregcm_def_field(ex_u10m, OASIS_Out)
    if ( l_cpl_ex_v10m ) call oasisxregcm_def_field(ex_v10m, OASIS_Out)
    if ( l_cpl_ex_wspd ) call oasisxregcm_def_field(ex_wspd, OASIS_Out)
    if ( l_cpl_ex_wdir ) call oasisxregcm_def_field(ex_wdir, OASIS_Out)
    if ( l_cpl_ex_t2m )  call oasisxregcm_def_field(ex_t2m,  OASIS_Out)
!    if ( l_cpl_ex_t10m ) call oasisxregcm_def_field(ex_t10m, OASIS_Out)
    if ( l_cpl_ex_q2m )  call oasisxregcm_def_field(ex_q2m,  OASIS_Out)
!    if ( l_cpl_ex_q10m ) call oasisxregcm_def_field(ex_q10m, OASIS_Out)
    if ( l_cpl_ex_slp )  call oasisxregcm_def_field(ex_slp,  OASIS_Out)
    if ( l_cpl_ex_taux ) call oasisxregcm_def_field(ex_taux, OASIS_Out)
    if ( l_cpl_ex_tauy ) call oasisxregcm_def_field(ex_tauy, OASIS_Out)
    if ( l_cpl_ex_z0 )   call oasisxregcm_def_field(ex_z0,   OASIS_Out)
    if ( l_cpl_ex_ustr ) call oasisxregcm_def_field(ex_ustr, OASIS_Out)
    if ( l_cpl_ex_evap ) call oasisxregcm_def_field(ex_evap, OASIS_Out)
    if ( l_cpl_ex_prec ) call oasisxregcm_def_field(ex_prec, OASIS_Out)
    if ( l_cpl_ex_nuwa ) call oasisxregcm_def_field(ex_nuwa, OASIS_Out)
    if ( l_cpl_ex_ulhf ) call oasisxregcm_def_field(ex_ulhf, OASIS_Out)
    if ( l_cpl_ex_ushf ) call oasisxregcm_def_field(ex_ushf, OASIS_Out)
    if ( l_cpl_ex_uwlw ) call oasisxregcm_def_field(ex_uwlw, OASIS_Out)
    if ( l_cpl_ex_dwlw ) call oasisxregcm_def_field(ex_dwlw, OASIS_Out)
    if ( l_cpl_ex_nulw ) call oasisxregcm_def_field(ex_nulw, OASIS_Out)
    if ( l_cpl_ex_uwsw ) call oasisxregcm_def_field(ex_uwsw, OASIS_Out)
    if ( l_cpl_ex_dwsw ) call oasisxregcm_def_field(ex_dwsw, OASIS_Out)
    if ( l_cpl_ex_ndsw ) call oasisxregcm_def_field(ex_ndsw, OASIS_Out)
    if ( l_cpl_ex_rhoa ) call oasisxregcm_def_field(ex_rhoa, OASIS_Out)
    ! OASIS field +++
#ifdef ECLM
!!  ECLM
    if ( l_cpl_ex_tatm ) call oasisxregcm_def_field(ex_tatm, OASIS_Out)
    if ( l_cpl_ex_uatm ) call oasisxregcm_def_field(ex_uatm, OASIS_Out)
    if ( l_cpl_ex_vatm ) call oasisxregcm_def_field(ex_vatm, OASIS_Out)
    if ( l_cpl_ex_qatm ) call oasisxregcm_def_field(ex_qatm, OASIS_Out)
    if ( l_cpl_ex_hgt ) call oasisxregcm_def_field(ex_hgt, OASIS_Out)
    if ( l_cpl_ex_patm ) call oasisxregcm_def_field(ex_patm, OASIS_Out)
    if ( l_cpl_ex_swdif ) call oasisxregcm_def_field(ex_swdif, OASIS_Out)
    if ( l_cpl_ex_swdir ) call oasisxregcm_def_field(ex_swdir, OASIS_Out)
    if ( l_cpl_ex_dwrlwf ) call oasisxregcm_def_field(ex_dwrlwf, OASIS_Out)
    if ( l_cpl_im_tgbb )  call oasisxregcm_def_field(im_tgbb,  OASIS_In)
    if ( l_cpl_im_t2m )  call oasisxregcm_def_field(im_t2m,  OASIS_In)
    if ( l_cpl_im_q2m )  call oasisxregcm_def_field(im_q2m,  OASIS_In)
    if ( l_cpl_im_u10m )  call oasisxregcm_def_field(im_u10m,  OASIS_In)
!    if ( l_cpl_im_v10m )  call oasisxregcm_def_field(im_v10m,  OASIS_In)
    if ( l_cpl_im_sent )  call oasisxregcm_def_field(im_sent,  OASIS_In)
    if ( l_cpl_im_evpr )  call oasisxregcm_def_field(im_evpr,  OASIS_In)
    if ( l_cpl_im_ram1 )  call oasisxregcm_def_field(im_ram1,  OASIS_In)
    if ( l_cpl_im_rah1 )  call oasisxregcm_def_field(im_rah1,  OASIS_In)
!    if ( l_cpl_im_br )  call oasisxregcm_def_field(im_br,  OASIS_In)
    if ( l_cpl_im_taux )  call oasisxregcm_def_field(im_taux,  OASIS_In)
    if ( l_cpl_im_tauy )  call oasisxregcm_def_field(im_tauy,  OASIS_In)
!    if ( l_cpl_im_drag )  call oasisxregcm_def_field(im_drag,  OASIS_In)
    if ( l_cpl_im_sncv )  call oasisxregcm_def_field(im_sncv,  OASIS_In)
    if ( l_cpl_im_wt )  call oasisxregcm_def_field(im_wt,  OASIS_In)
    if ( l_cpl_im_zo )  call oasisxregcm_def_field(im_zo,  OASIS_In)
    if ( l_cpl_im_tlef )  call oasisxregcm_def_field(im_tlef,  OASIS_In)
    if ( l_cpl_ex_rainf ) call oasisxregcm_def_field(ex_rainf, OASIS_Out)
    if ( l_cpl_ex_snow ) call oasisxregcm_def_field(ex_snow, OASIS_Out)
    if ( l_cpl_im_albed )  call oasisxregcm_def_3dfield(im_albed,  OASIS_In)
    if ( l_cpl_im_albei )  call oasisxregcm_def_3dfield(im_albei,  OASIS_In)
    if ( l_cpl_im_flxvoc )  call oasisxregcm_def_3dfield(im_flxvoc,  OASIS_In)
    if ( l_cpl_im_flxdst )  call oasisxregcm_def_3dfield(im_flxdst,  OASIS_In)
    if ( l_cpl_im_fluxch4 )  call oasisxregcm_def_field(im_fluxch4,  OASIS_In)
    if ( l_cpl_im_ddvel )  call oasisxregcm_def_3dfield(im_ddvel,  OASIS_In)
    if ( l_cpl_im_tlai )  call oasisxregcm_def_field(im_tlai,  OASIS_In)
    if ( l_cpl_im_roff )  call oasisxregcm_def_field(im_roff,  OASIS_In)
    if ( l_cpl_im_srfroff )  call oasisxregcm_def_field(im_srfroff,  OASIS_In)
    if ( l_cpl_im_snwmelt )  call oasisxregcm_def_field(im_snwmelt,  OASIS_In)
    if ( l_cpl_im_tsoi )  call oasisxregcm_def_3dfield(im_tsoi,  OASIS_In)
    if ( l_cpl_im_h2ovol )  call oasisxregcm_def_3dfield(im_h2ovol,  OASIS_In)
    if ( l_cpl_im_h2oliq )  call oasisxregcm_def_3dfield(im_h2oliq,  OASIS_In)
    if ( l_cpl_im_h2oice )  call oasisxregcm_def_3dfield(im_h2oice,  OASIS_In)
    if ( l_cpl_im_h2o10cm )  call oasisxregcm_def_field(im_h2o10cm,  OASIS_In)
#endif
    ! termination of definition phase
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: end'
#endif
    call oasisxregcm_end_def

  contains

  ! define OASIS grids
  subroutine oasisxregcm_def_grid
    implicit none
    real(rkx) , pointer , dimension(:,:) :: dlon , dlat ! dot degree coordinates
    real(rkx) , pointer , dimension(:,:) :: xlon , xlat ! cross degree coordinates
    real(rkx) , pointer , dimension(:,:) :: lndcat ! land category (15 for ocean)
    real(rkx) , pointer , dimension(:,:) :: oasisgrid_lon , oasisgrid_lat ! lon, lat
    real(rkx) , pointer , dimension(:,:,:) :: lon3d, lat3d! lon, lat
    real(rkx) , pointer , dimension(:,:,:) :: oasisgrid_clon , & ! lon &
                                              oasisgrid_clat     ! lat of the corners
    ! note: 4 corners (always the case?) in the third dimension, counterclockwisely.
    real(rkx) , pointer , dimension(:,:) :: oasisgrid_srf ! surface of the grid meshes m2
    real(rkx) , pointer , dimension(:,:,:) :: srf3d ! surface of the grid meshes m2
    integer(ik4) , pointer , dimension(:,:) :: oasisgrid_mask ! mask, 0 = valid, 1 = mask
    integer(ik4) , pointer , dimension(:,:,:) :: mask3d! mask, 0 = valid, 1 = mask
    !                                                          (OASIS convention)
    integer(ik4) :: il_flag ! flag for grid writing by proc 0
    integer(ik4) :: ierror
#ifdef DEBUG
    !--------------------------------------------------------------------------
    write(ndebug,*) oasis_prefix, 'start collecting grid data'
    write(ndebug,*) oasis_prefix, 'grids writing is done by the in/out cpu only'
#endif
    ! collect the domain definition information in RegCM
    allocate(dlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(dlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(lndcat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating lndcat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    call grid_collect(mddom%dlon,dlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%dlat,dlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlon,xlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlat,xlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%lndcat,lndcat,jde1,jde2,ide1,ide2)
    if ( myid == iocpu ) then
      ! start the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'start grids writing'
#endif
      call oasis_start_grids_writing(il_flag)
      !
      if ( l_make_grdde .or. l_make_grddi) then ! dots
        !
        call oasisxregcm_make_oasisgrids_d( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdde ) call oasisxregcm_write_oasisgrids(grdde, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grddi ) call oasisxregcm_write_oasisgrids(grddi, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
      end if
      !
      if ( l_make_grdce .or. l_make_grdci .or. l_make_grdtest) then ! crosses
        !
        call oasisxregcm_make_oasisgrids_c( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdce ) call oasisxregcm_write_oasisgrids(grdce, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grdci ) call oasisxregcm_write_oasisgrids(grdci, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
#ifdef ECLM
        if ( l_make_grdtest ) call oasisxregcm_write_oasisgrids(grdtest, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
!        if ( l_make_grd3d ) call oasisxregcm_write_oasis3dgrids(grd3d, &
           !  oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        if (l_make_grd3d) then
    ! Extend lon and lat to 3D by repeating them along the depth dimension.
    ! The depth dimension will be added in the 3rd dimension.

    ! Create temporary 3D arrays for lon and lat by replicating the 2D arrays across depth

          allocate(lon3d(grd3d%jgl, grd3d%igl, kdim))
          allocate(lat3d(grd3d%jgl, grd3d%igl, kdim))
          allocate(srf3d(grd3d%jgl, grd3d%igl, kdim))
          allocate(mask3d(grd3d%jgl, grd3d%igl, kdim))

    ! Replicate 2D lon and lat into 3D along the depth (third) dimension
          do k = 1 , kdim
            lon3d(:,:, k) = oasisgrid_lon(:,:)
            lat3d(:,:, k) = oasisgrid_lat(:,:)
            srf3d(:,:, k) = oasisgrid_srf(:,:)
            mask3d(:,:, k) = oasisgrid_mask(:,:)
          end do
    ! Now, call the subroutine with the 3D arrays
          call oasisxregcm_write_oasis3dgrids(grd3d, &
            lon3d, lat3d, oasisgrid_clon, oasisgrid_clat, srf3d,mask3d)

    ! Deallocate temporary 3D arrays
          deallocate(lon3d)
          deallocate(lat3d)
          deallocate(srf3d)
        end if

           !
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat,oasisgrid_srf,oasisgrid_mask)
        !
#endif 
      end if
      ! terminate the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'terminate grids writing'
#endif
      call oasis_terminate_grids_writing()
    end if
    if ( associated(dlon) ) nullify(dlon)
    if ( associated(dlat) ) nullify(dlat)
    if ( associated(xlon) ) nullify(xlon)
    if ( associated(xlat) ) nullify(xlat)
    if ( associated(lndcat) ) nullify(lndcat)
#ifdef DEBUG
    write(ndebug,"(' ',A,A,I3,A)") oasis_prefix, 'please refer to cpu ', iocpu, ' for grid writing'&
                           //' debug statements'
#endif
  end subroutine oasisxregcm_def_grid

  subroutine oasisxregcm_make_oasisgrids_d(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: dlon , dlat , &
                                                         xlon , xlat , lndcat
    real(rkx) , pointer , dimension(:,:) , intent(out) :: lon , lat , srf
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: clon , clat
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: mask
    integer(ik4) :: i , j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    x(j-1,i)         x(j,i)
    !  i |
    !    |             center
    !  a |             d(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |   x(j-1,i-1)       x(j,i-1)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: concretely xlat and xlon are defined on (1:jx,1:iy) as for dlat and dlon
    !       though the cross grid is read on (1:jx-1,1:iy-1) only
    ! consequence: xlon and xlat exist at (jx,jy) -> no exception needed
    !              however, exceptions are needed for extrem west and south values
    !                                                         0         0
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx,iy,4)
    do i = 1 , iy
      do j = 1 , jx
        ! center
        lon(j,i) = dlon(j,i)
        lat(j,i) = dlat(j,i)
        ! upper right corner: no problem
        clon(j,i,1) = xlon(j,i)
        clat(j,i,1) = xlat(j,i)
        ! left corners:
        !   upper left corner: use of j-1 -> exception for j=1
        !                      lon -> use of the last lon step to the West
        !                      lat -> use of the same value as for j
        !                             (expected to be the same)
        !   lower left corner: same + use of i-1 -> exception for i=1
        !                      lon -> use of the same value as for i
        !                             (expected to be the same)
        !                      lat -> use of the last lat step to the South
        if ( j == 1 ) then
          clon(j,i,2) = 2 * xlon(j,i) - xlon(j+1,i)
          clat(j,i,2) = xlat(j,i)
          if ( i == 1 ) then
            clon(j,i,3) = 2 * xlon(j,i) - xlon(j+1,i)
            clat(j,i,3) = 2 * xlat(j,i) - xlat(j,i+1)
          else
            clon(j,i,3) = 2 * xlon(j,i-1) - xlon(j+1,i-1)
            clat(j,i,3) = xlat(j,i-1)
          end if
        else
          clon(j,i,2) = xlon(j-1,i)
          clat(j,i,2) = xlat(j-1,i)
          if ( i == 1 ) then
            clon(j,i,3) = xlon(j-1,i)
            clat(j,i,3) = 2 * xlat(j-1,i) - xlat(j-1,i+1)
          else
            clon(j,i,3) = xlon(j-1,i-1)
            clat(j,i,3) = xlat(j-1,i-1)
          end if
        end if
        ! lower right corner: use of i-1 -> exception for i=1
        if ( i == 1 ) then
          clon(j,i,4) = xlon(j,i)
          clat(j,i,4) = 2 * xlat(j,i) - xlat(j,i+1)
        else
          clon(j,i,4) = xlon(j,i-1)
          clat(j,i,4) = xlat(j,i-1)
        end if
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 1
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_d

  subroutine oasisxregcm_make_oasisgrids_c(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(in) :: dlon , dlat , &
                                                         xlon , xlat , lndcat
    real(rkx) , pointer , dimension(:,:) , intent(out) :: lon , lat , srf
    real(rkx) , pointer , dimension(:,:,:) , intent(out) :: clon , clat
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: mask
    integer(ik4) :: i , j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    d(j,i+1)       d(j+1,i+1)
    !  i |
    !    |             center
    !  a |             x(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |     d(j,i)         d(j+1,i)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: the cross case is easier as the cross grid is contained
    !       in the dot one thus all corners are already defined in dlon
    !
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx-1,iy-1,4)
    do i = 1 , iy-1
      do j = 1 , jx-1
        ! center
        lon(j,i) = xlon(j,i)
        lat(j,i) = xlat(j,i)
        ! upper right corner
        clon(j,i,1) = dlon(j+1,i+1)
        clat(j,i,1) = dlat(j+1,i+1)
        ! upper left corner
        clon(j,i,2) = dlon(j,i+1)
        clat(j,i,2) = dlat(j,i+1)
        ! lower left corner
        clon(j,i,3) = dlon(j,i)
        clat(j,i,3) = dlat(j,i)
        ! lower right corner
        clon(j,i,4) = dlon(j+1,i)
        clat(j,i,4) = dlat(j+1,i)
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 1
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_c

  end subroutine oasisxregcm_def

  ! call all subroutines consisting of receiving OASIS fields
  ! and optionally reworking them
  subroutine oasisxregcm_rcv_all(time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    integer(ik4) :: i , j, n
    logical :: l_act
    character(len=4), allocatable :: megan_specifier(:)
    integer :: shr_megan_mechcomps_n

    !--------------------------------------------------------------------------
#ifdef DEBUG
    if ( time == 0 ) then
      write(ndebug,*) oasis_prefix, 'Initialization:'
    else if ( debug_level > 0 ) then
      write(ndebug,*) oasis_prefix, 'OASIS time: ', time
    end if
#endif
    !
    if ( l_cpl_im_sst ) then ! sea surface temperature [K]
      call oasisxregcm_rcv(cpl_sst,im_sst,time,l_act)
      !if ( l_act ) then
        ! fill the ocean parts: out array (2d or 3d with nsg dimension)
        !                       in array, mask, field grid
        call fill_ocean(sfs%tg,cpl_sst,mddom%lndcat,im_sst%grd)
        call fill_ocean(sfs%tgbb,cpl_sst,mddom%lndcat,im_sst%grd)
      !end if
    end if
    ! NOT IMPLEMENTED YET, BELOW IS A COPY OF THE PROCESS IN
    ! MOD_LM_INTERFACE.F90; IMPORT_DATA_INTO_SURFACE.
!    if ( l_cpl_im_sit ) then ! sea ice thickness [m]
!      call oasisxregcm_rcv(cpl_sit,im_sit,time,l_act)
!!      if ( l_act ) then
!        grd => im_sit%grd
!        ! dans le mask
!        ! tol = 0.5 E 20
!          if ( impfie%sit(j,i) < tol .and. lm%ldmsk(j,i) /= 1 ) then
!            if ( impfie%sit(j,i) > iceminh ) then
!              flag = .false.
!              if ( lm%ldmsk(j,i) == 0 ) flag = .true.
!              ! set land-sea mask
!              lm%ldmsk(j,i) = 2
!              do n = 1, nnsg
!                lm%ldmsk1(n,j,i) = 2
!                ! set sea ice thikness (in meter)
!                lms%sfice(n,j,i) = impfie%sit(j,i)
!              end do
!              ! write debug info
!              if ( flag ) then
!                write(*,f99002) j, i, 'water', 'ice  ', &
!                   lm%ldmsk(j,i), lms%sfice(1,j,i)
!              end if
!            else
!              if ( ldmskb(j,i) == 0 .and. lm%ldmsk(j,i) == 2 ) then
!                do n = 1, nnsg
!                  ! reduce to one tenth surface ice: it should melt away
!                  lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
!                  ! check that sea ice is melted or not
!                  if ( lms%sfice(n,j,i) <= iceminh ) then
!                    if ( ldmskb(j,i) /= lm%ldmsk(j,i) ) flag = .true.
!                    ! set land-sea mask to its original value
!                    lm%ldmsk(j,i) = ldmskb(j,i)
!                    lm%ldmsk1(n,j,i) = ldmskb(j,i)
!                    ! set land-use type to its original value
!                    ! set sea ice thikness (in meter)
!                    lms%sfice(n,j,i) = d_zero
!                  else
!                    flag = .false.
!                  end if
!                end do
!                ! write debug info
!                if ( flag ) then
!                  write(*,f99003) j, i, 'ice  ', 'water',  &
!                    lm%ldmsk(j,i), lms%sfice(1,j,i)
!                end if
!              end if
!            end if
!          end if
!        ???
!        nullify(grd)
!!      end if
!    end if
    !
    if ( l_cpl_im_wz0 ) then ! surface roughness length [m]
      call oasisxregcm_rcv(cpl_wz0,im_wz0,time,l_act)
!      if ( l_act ) then
        call fill_ocean(sfs%zo,cpl_wz0,mddom%lndcat,im_wz0%grd)
!      end if
    end if
    !
    if ( l_cpl_im_wust ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_wust,im_wust,time,l_act)
!      if ( l_act ) then
        call fill_ocean(sfs%ustar ,cpl_wust,mddom%lndcat,im_wust%grd)
!      end if
    end if
    ! OASIS field +++
    !
#ifdef ECLM
    if ( l_cpl_im_tgbb ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_tgbb,im_tgbb,time,l_act)
!      if ( l_act ) then
        call fill_land(lms%tgbb ,cpl_tgbb,mddom%lndcat,im_tgbb%grd)
!      end if
    end if
    if ( l_cpl_im_t2m ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_t2m,im_t2m,time,l_act)
        call fill_land(lms%t2m ,cpl_t2m,mddom%lndcat,im_t2m%grd)
    end if
    if ( l_cpl_im_q2m ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_q2m,im_q2m,time,l_act)
        call fill_land(lms%q2m ,cpl_q2m,mddom%lndcat,im_q2m%grd)
    end if
    if ( l_cpl_im_u10m ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_u10m,im_u10m,time,l_act)
        call fill_land(lms%u10m ,cpl_u10m,mddom%lndcat,im_u10m%grd)
    end if
    if ( l_cpl_im_sent ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_sent,im_sent,time,l_act)
        call fill_land(lms%sent ,cpl_sent,mddom%lndcat,im_sent%grd)
    end if
    if ( l_cpl_im_evpr ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_evpr,im_evpr,time,l_act)
        call fill_land(lms%evpr ,cpl_evpr,mddom%lndcat,im_evpr%grd)
    end if
    if ( l_cpl_im_ram1 ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_ram1,im_ram1,time,l_act)
        call fill_land(lms%ram1 ,cpl_ram1,mddom%lndcat,im_ram1%grd)
    end if
    if ( l_cpl_im_rah1 ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_rah1,im_rah1,time,l_act)
        call fill_land(lms%rah1 ,cpl_rah1,mddom%lndcat,im_rah1%grd)
    end if
    if ( l_cpl_im_tauy .or. l_cpl_im_taux) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_tauy,im_tauy,time,l_act)
        call fill_land(lms%tauy ,cpl_tauy,mddom%lndcat,im_tauy%grd)
      call oasisxregcm_rcv(cpl_taux,im_taux,time,l_act)
        call fill_land(lms%taux ,cpl_taux,mddom%lndcat,im_taux%grd)
      call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
      !call getmem2d(temp2,jci1,jci2,ici1,ici2,'sendoasis:temp')
!      call getmem3d(temp3d_2,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      temp = lm%uatm**2 +lm%vatm**2
      do n = 1 , nnsg
       temp3d(n,:,:)= sqrt(lms%taux(n,:,:)**2+lms%tauy(n,:,:)**2) /temp
      end do
      lms%drag=temp3d
    end if
    if ( l_cpl_im_sncv ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_sncv,im_sncv,time,l_act)
        call fill_land(lms%sncv ,cpl_sncv,mddom%lndcat,im_sncv%grd)
    end if
    if ( l_cpl_im_wt ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_wt,im_wt,time,l_act)
        call fill_land(lms%wt ,cpl_wt,mddom%lndcat,im_wt%grd)
    end if
    if ( l_cpl_im_zo ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_zo,im_zo,time,l_act)
        call fill_land(lms%zo ,cpl_zo,mddom%lndcat,im_zo%grd)
    end if
    if ( l_cpl_im_tlef ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_tlef,im_tlef,time,l_act)
        call fill_land(lms%tlef ,cpl_tlef,mddom%lndcat,im_tlef%grd)
    end if
    if ( l_cpl_im_tlai ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_tlai,im_tlai,time,l_act)
        call fill_land(lms%xlai ,cpl_tlai,mddom%lndcat,im_tlai%grd)
    end if
    if ( l_cpl_im_roff ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_roff,im_roff,time,l_act)
        call fill_land(lms%srnof ,cpl_roff,mddom%lndcat,im_roff%grd)
    end if
    if ( l_cpl_im_srfroff ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_srfroff,im_srfroff,time,l_act)
        call fill_land(lms%trnof ,cpl_srfroff,mddom%lndcat,im_srfroff%grd)
    end if
    if ( l_cpl_im_snwmelt ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_snwmelt,im_snwmelt,time,l_act)
        call fill_land(lms%snwm ,cpl_snwmelt,mddom%lndcat,im_snwmelt%grd)
    end if
    if ( l_cpl_im_h2o10cm ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_h2o10cm,im_h2o10cm,time,l_act)
        call fill_land(lms%tsw ,cpl_h2o10cm,mddom%lndcat,im_h2o10cm%grd)
    end if
    if ( l_cpl_im_albed ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_albed,im_albed,time,l_act)
        call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      ! call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
        call fill_land(temp3d,cpl_albed(:,:,1),mddom%lndcat,im_albed%grd)
        do i = ici1 , ici2
          do j = jci1 , jci2
           do n = 1 , nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( temp3d(n,j,i) > 0.9999_rkx ) then
                temp3d(n,j,i)=0.16_rkx
              end if
            end if
           end do
          end do
        end do
        lms%swdiralb=temp3d
        deallocate(temp3d)
        call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      ! call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
        call fill_land(temp3d,cpl_albed(:,:,2),mddom%lndcat,im_albed%grd)
        do i = ici1 , ici2
          do j = jci1 , jci2
           do n = 1 , nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( temp3d(n,j,i) > 0.9999_rkx ) then
                temp3d(n,j,i)=0.32_rkx
              end if
            end if
           end do
          end do
        end do
        lms%lwdiralb=temp3d
        deallocate(temp3d)
        !call fill_land(lms%swdiralb ,cpl_albed(:,:,1),mddom%lndcat,im_albed%grd)
        !call fill_land(lms%lwdiralb ,cpl_albed(:,:,2),mddom%lndcat,im_albed%grd)
!      deallocate(temp3d)
    end if
    if ( l_cpl_im_albei ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_albei,im_albei,time,l_act)
        call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      ! call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
        call fill_land(temp3d,cpl_albei(:,:,1),mddom%lndcat,im_albei%grd)
        do i = ici1 , ici2
          do j = jci1 , jci2
           do n = 1 , nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( temp3d(n,j,i) > 0.9999_rkx ) then
                temp3d=0.16_rkx
              end if
            end if
           end do
          end do
        end do
        lms%swdifalb=temp3d
        deallocate(temp3d)
        call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      ! call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
        call fill_land(temp3d,cpl_albei(:,:,2),mddom%lndcat,im_albei%grd)
        do i = ici1 , ici2
          do j = jci1 , jci2
           do n = 1 , nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              if ( temp3d(n,j,i) > 0.9999_rkx ) then
                temp3d=0.32_rkx
              end if
            end if
           end do
          end do
        end do
        lms%lwdifalb=temp3d
        deallocate(temp3d)
      !        call fill_land(lms%swdifalb ,cpl_albei(:,:,1),mddom%lndcat,im_albei%grd)
!        call fill_land(lms%lwdifalb ,cpl_albei(:,:,2),mddom%lndcat,im_albei%grd)
    end if
!    if ( l_cpl_im_tsoi ) then ! surface friction velocity [s-1]
 !     call oasisxregcm_rcv_3d(cpl_tsoi,im_tsoi,time,l_act)
  !      n = ubound(tsoi, 3)
   !     do k = 1 , n
    !      call fill_land(lms%tsoi(:,:,:,k) ,cpl_tsoi(:,:,k),mddom%lndcat,im_tsoi%grd)
     !   end do 
!    end if
    if ( l_cpl_im_h2ovol ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_h2ovol,im_h2ovol,time,l_act)
      n = ubound(cpl_h2ovol, 3)
 !     call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
      call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      do k = 1 , n
        call fill_land(temp3d,cpl_h2ovol(:,:,k),mddom%lndcat,im_h2ovol%grd)
        lms%sw_vol(:,:,:,k) = temp3d(:,:,:)
      end do
      !deallocate(temp)
      deallocate(temp3d)
    end if
!if ( l_cpl_ex_snow .or. l_cpl_ex_rainf) then
    if ( l_cpl_im_h2oliq .or. l_cpl_im_h2oice ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_h2oliq,im_h2oliq,time,l_act)
      n = ubound(cpl_h2oliq, 3)
!      call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
      call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      call getmem3d(temp3d_2,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d_2')
      do k = 1 , n
        call fill_land(temp3d,cpl_h2oliq(:,:,k),mddom%lndcat,im_h2oliq%grd)
        call fill_land(temp3d_2,cpl_h2oice(:,:,k),mddom%lndcat,im_h2oice%grd)
        lms%sw(:,:,:,k) = temp3d(:,:,:) + temp3d_2(:,:,:)
      end do
     ! deallocate(temp)
      deallocate(temp3d)
      deallocate(temp3d_2)
    end if

    if ( l_cpl_im_tsoi ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_tsoi,im_tsoi,time,l_act)
      n = ubound(cpl_tsoi, 3)
  !    call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
      call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      do k = 1 , n
        call fill_land(temp3d,cpl_tsoi(:,:,k),mddom%lndcat,im_tsoi%grd)
        lms%tsoi(:,:,:,k) = temp3d(:,:,:)
      end do
    !  deallocate(temp)
      deallocate(temp3d)
    end if
    if ( ichem == 1 ) then
!     if ( l_cpl_im_ddvel ) then ! surface friction velocity [s-1]
 !     call oasisxregcm_rcv_3d(cpl_ddvel,im_ddvel,time,l_act)
  !     n = ubound(cpl_ddvel, 3)
   !    call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
    !   call fill_land(temp3d,cpl_ddvel(:,:,n),mddom%lndcat,im_ddvel%grd)
     !  lms%ddepv(:,:,:,4) = temp3d
      ! deallocate(temp3d)
!     end if
     if ( l_cpl_im_flxvoc ) then ! surface friction velocity [s-1]
      shr_megan_mechcomps_n = 9
      allocate(megan_specifier(9))
      megan_specifier = ['ISOP', 'LIMO', 'APIN', 'BPIN', '3CAR', 'OCIM', 'MYRC', 'SABI', 'OMTP']
      call oasisxregcm_rcv_3d(cpl_flxvoc,im_flxvoc,time,l_act)
      call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
      do k = 1, shr_megan_mechcomps_n
         if (megan_specifier(k) == 'ISOP'.and. iisop > 0 ) then
          call fill_land(temp3d,cpl_flxvoc(:,:,k),mddom%lndcat,im_flxvoc%grd)
         end if 
      end do
!      call fill_land(temp3d,cpl_flxvoc,mddom%lndcat,im_flxvoc%grd)
      lms%vocemiss(:,:,:,iisop) = temp3d(:,:,:)
      deallocate(temp3d)
     end if
     if ( l_cpl_im_flxdst ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv_3d(cpl_flxdst,im_flxdst,time,l_act)
      if ( ichdustemd == 3 ) then
       n = ubound(cpl_flxdst, 3)
       call getmem4d(temp4d,1,nnsg,jci1,jci2,ici1,ici2,1,n,'sendoasis:temp3d')
       call fill_land(temp4d,cpl_flxdst,mddom%lndcat,im_flxdst%grd)
!      do k = 1 , 4
        lms%dustemiss(:,:,:,:) = temp4d(:,:,:,:)
 !     end do
       deallocate(temp4d)
      end if
     end if
     if ( l_cpl_im_fluxch4 ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_fluxch4,im_fluxch4,time,l_act)
      if ( ich4 > 0 ) then
      !      call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
       call getmem3d(temp3d,1,nnsg,jci1,jci2,ici1,ici2,'sendoasis:temp3d')
       call getmem2d(temp2,jci1,jci2,ici1,ici2,'sendoasis:temp')
       temp2=cpl_fluxch4* 1.33_rk8
        call fill_land(temp3d,temp2,mddom%lndcat,im_fluxch4%grd)
 !     do k = 1 , 4
        lms%vocemiss(:,:,:,ich4)=temp3d(:,:,:)
!       end do
      end if
      deallocate(temp2)
     end if
    end if
#endif

!    if ( l_cpl_im_fluxch4 ) then ! surface friction velocity [s-1]
 !     call getmem2d(temp,jci1,jci2,ici1,ici2,'sendoasis:temp')
  !    temp = 0.0_rk8
   !   temp = cpl_fluxch4* 1.33_rk8
!      if (.not. associated(lms%vocemiss)) then
 !       print *, "Error: lms%vocemiss is not associated!"
  !      stop
   !   endif
      !vocemiss_slice => lms%vocemiss(:,:,:,1)
    !  call oasisxregcm_rcv(cpl_fluxch4,im_fluxch4,time,l_act)
     !   call fill_land(lms%vocemiss(:,:,:,ich4) ,temp,mddom%lndcat,im_fluxch4%grd)
    
!    end if
    
    !if ( l_cpl_im_ddvel ) then ! surface friction velocity [s-1]
     ! call oasisxregcm_rcv(cpl_ddvel,im_ddvel,time,l_act)
      !  call fill_land(sfs%ustar ,cpl_ddvel,mddom%lndcat,im_ddvel%grd)
   ! end if    
  end subroutine oasisxregcm_rcv_all

  ! call all subroutines consisting of sending OASIS fields
  ! with optional prior reworking
  subroutine oasisxregcm_snd_all(time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    type(infogrd) , pointer :: grd
    integer(ik4) :: i , j , ishift , jshift
    logical :: l_write_restart
    call getmem2d(temps,jci1,jci2,ici1,ici2,'sendoasis:temps')
    call getmem2d(temps_rf,jci1,jci2,ici1,ici2,'sendoasis:temps_rf')
    call getmem2d(temps_snw,jci1,jci2,ici1,ici2,'sendoasis:temps_snw')
    !--------------------------------------------------------------------------
    l_write_restart = .false.
    if ( associated(alarm_out_sav) ) then
      if ( alarm_out_sav%act( ) ) then
        l_write_restart = .true.
      else
        if ( lfdomonth(rcmtimer%idate) .and. &
             lmidnight(rcmtimer%idate) ) then
          l_write_restart = .true.
        end if
      end if
    end if
    if ( rcmtimer%reached_endtime ) l_write_restart = .true.
    if ( write_restart_option == 1 .and. time == 0 ) l_write_restart = .true.
    if ( l_cpl_ex_u10m ) then ! eatward near-surface wind (10m-height) [m.s-1]
      grd => ex_u10m%grd
      call oasisxregcm_snd( &
           sfs%u10m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_u10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_v10m ) then ! northward near-surface wind (10m-height) [m.s-1]
      grd => ex_v10m%grd
      call oasisxregcm_snd( &
           sfs%v10m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_v10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_wspd ) then ! near-surface wind speed (10m-height) [m.s-1]
      grd => ex_wspd%grd
      call oasisxregcm_snd( &
           sqrt(sfs%u10m(grd%j1:grd%j2 , grd%i1:grd%i2)**2   &
              + sfs%v10m(grd%j1:grd%j2 , grd%i1:grd%i2)**2), &
           ex_wspd, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_wdir ) then ! near-surface wind direction (10m-height) [degree]
      grd => ex_wdir%grd
      do i = grd%i1 , grd%i2
        ishift = i - grd%i1 + 1
        do j = grd%j1 , grd%j2
          jshift = j - grd%j1 + 1
          cpl_wdir(jshift,ishift) = atan2( sfs%u10m(j,i) , sfs%v10m(j,i) ) * raddeg
          if ( cpl_wdir(jshift,ishift) < d_zero ) &
               cpl_wdir(jshift,ishift) = cpl_wdir(jshift,ishift) + 360_rkx
        end do
      end do
      call oasisxregcm_snd(cpl_wdir, ex_wdir, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_t2m ) then ! near-surface air temperature (2m-height) [K]
      grd => ex_t2m%grd
      call oasisxregcm_snd( &
           sum( lms%t2m(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_t2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    if ( l_cpl_ex_t10m ) then ! near-surface air temperature (10m-height) [K]
!      grd => ex_t10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_t10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    if ( l_cpl_ex_q2m ) then ! near-surface air specific humidity (2m-height) [1]
      grd => ex_q2m%grd
      call oasisxregcm_snd( &
           sfs%q2m(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_q2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    if ( l_cpl_ex_q10m ) then ! near-surface air specific humidity (10m-height) [1]
!      grd => ex_q10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_q10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    if ( l_cpl_ex_slp ) then ! sea level pressure [Pa]
      grd => ex_slp%grd
      call oasisxregcm_snd( &
           atms%ps2d(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_slp, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_taux) then ! surface eastward wind stress [Pa]
      grd => ex_taux%grd
      call oasisxregcm_snd( &
           sum( lms%taux(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_taux, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_tauy) then ! surface northward wind stress [Pa]
      grd => ex_tauy%grd
      call oasisxregcm_snd( &
           sum( lms%tauy(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_tauy, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_z0 ) then ! surface roughness length [m]
      grd => ex_z0%grd
      call oasisxregcm_snd( &
           sfs%zo(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_z0, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ustr ) then ! surface friction velocity [s-1]
      grd => ex_ustr%grd
      call oasisxregcm_snd( &
           sfs%ustar(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ustr, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_evap ) then ! water evaporation flux [kg.m-2.s-1]
      grd => ex_evap%grd
      call oasisxregcm_snd( &
           sfs%qfx(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_evap, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_prec ) then ! precipitation flux [kg.m-2.s-1]
      grd => ex_prec%grd
      call oasisxregcm_snd( &
           sum( lms%prcp(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_prec, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_nuwa ) then ! net upward water flux [kg.m-2.s-1]
      grd => ex_nuwa%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(: , grd%j1:grd%j2 , grd%i1:grd%i2) - &
                lms%prcp(: , grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_nuwa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ulhf ) then ! surface upward latent heat flux [W.m-2]
      grd => ex_ulhf%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(: , grd%j1:grd%j2 , grd%i1:grd%i2) &
           *wlh(lms%t2m(: , grd%j1:grd%j2 , grd%i1:grd%i2)) , 1 ) * rdnnsg, &
           ex_ulhf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ushf ) then ! surface upward sensible heat flux [W.m-2]
      grd => ex_ushf%grd
      call oasisxregcm_snd( &
           sfs%hfx(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ushf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_uwlw ) then ! surface upwelling long-wave radiation flux [W.m-2]
      grd => ex_uwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2 , grd%i1:grd%i2) &
           + flw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_uwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_dwlw ) then ! surface downwelling long-wave radiation flux [W.m-2]
      grd => ex_dwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_dwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_nulw ) then ! surface net upward long-wave radiation flux [W.m-2]
      grd => ex_nulw%grd
      call oasisxregcm_snd( &
           flw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_nulw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_uwsw ) then ! surface upwelling short-wave radiation flux [W.m-2]
      grd => ex_uwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2 , grd%i1:grd%i2) &
           - fsw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_uwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_dwsw ) then ! surface downwelling short-wave radiation flux [W.m-2]
      grd => ex_dwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_dwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ndsw ) then ! surface net downward short-wave radiation flux [W.m-2]
      grd => ex_ndsw%grd
      call oasisxregcm_snd( &
           fsw(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_ndsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_rhoa ) then ! surface air density [kg.m-3]
      grd => ex_rhoa%grd
      call oasisxregcm_snd( &
           sum( lms%rhoa(:, grd%j1:grd%j2 , grd%i1:grd%i2) , 1 ) * rdnnsg, &
           ex_rhoa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    ! OASIS field +++
   !!! ECLM
#ifdef ECLM
    if ( l_cpl_ex_tatm ) then ! 
      grd => ex_tatm%grd
      call oasisxregcm_snd( &
           lm%tatm(grd%j1:grd%j2 , grd%i1:grd%i2), & !(:, grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_tatm, time, .false. .or. l_write_restart)
      nullify(grd)
    end if 
    if ( l_cpl_ex_uatm ) then ! 
      grd => ex_uatm%grd
      call oasisxregcm_snd( &
           lm%uatm( grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_uatm, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_vatm ) then ! 
      grd => ex_vatm%grd
      call oasisxregcm_snd( &
           lm%vatm( grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_vatm, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_qatm ) then ! 
      grd => ex_qatm%grd
      call oasisxregcm_snd( &
           lm%qvatm(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_qatm, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_hgt ) then ! 
      grd => ex_hgt%grd
      call oasisxregcm_snd( &
           lm%hgt(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_hgt, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_patm ) then ! 
      grd => ex_patm%grd
      call oasisxregcm_snd( &
           lm%patm( grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_patm, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_dwrlwf ) then ! 
      grd => ex_dwrlwf%grd
      call oasisxregcm_snd( &
           lm%dwrlwf(grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_dwrlwf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_swdir ) then ! 
      grd => ex_swdir%grd
      call oasisxregcm_snd( &
           lm%swdir( grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_swdir, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_swdif ) then ! 
      grd => ex_swdif%grd
      call oasisxregcm_snd( &
           lm%swdif( grd%j1:grd%j2 , grd%i1:grd%i2), &
           ex_swdif, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_snow .or. l_cpl_ex_rainf) then
     temps=(lm%cprate+lm%ncprate) * syncro_srf%rw
     do i = ici1 , ici2
      do j = jci1 , jci2
       if ( lm%tatm(j,i) -tzero <1.0_rk8 ) then
      !  print *, 2
        temps_rf(j,i)=d_zero
        temps_snw(j,i)=temps(j,i) !lm%cprate+lm%ncprate) * syncro_srf%rw
       else 
    !   ! print *, 3
        temps_snw(j,i)=d_zero
        temps_rf(j,i)=temps(j,i) !(lm%cprate+lm%ncprate) * syncro_srf%rw
       end if
      end do
     end do
    end if
    if ( l_cpl_ex_rainf ) then ! 
      grd => ex_rainf%grd
      call oasisxregcm_snd( &
           !lms%prcp(:, grd%j1:grd%j2 , grd%i1:grd%i2), &
!           (lm%cprate+lm%ncprate) * syncro_srf%rw , &  
           temps_rf, &
           ex_rainf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
!    do i = begg , endg


    if ( l_cpl_ex_snow ) then ! 
      grd => ex_snow%grd

!     if ( lm%tatm(grd%j1:grd%j2 , grd%i1:grd%i2)  <273.15 ) then
!         temps=(lm%cprate+lm%ncprate) * syncro_srf%rw
 !     else
  !       temps=0
!if ( clm_a2l%forc_t(i)-tfrz < tcrit ) then
 !       clm_a2l%forc_snow(i) = clm_a2l%rainf(i)
  !     clm_a2l%forc_rain(i) = d_zero
      call oasisxregcm_snd( &
           temps_snw, &
           !lms%prcp(:, grd%j1:grd%j2 , grd%i1:grd%i2), &
           !(lm%cprate+lm%ncprate) * syncro_srf%rw , &  
           ex_snow, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
#endif
    !!!
   !
    if ( myid == italk ) then
      if (write_restart_option == 1 .and. time == 0) then
        write(stdout,*) 'Note: OASIS restart files written at the' &
                      //' first time.'
      end if
    end if

    contains

#include <wlh.inc>

  end subroutine oasisxregcm_snd_all

  ! wait during the time indicated through oasis_sync_lag
  ! to make RegCM model time not synchronised with
  ! other coupled components
  ! update oasis_lag such that
  !   OASIS time = RegCM time + oasis_lag
  subroutine oasisxregcm_sync_wait(time)
    implicit none
    integer(ik4) , intent(in) :: time ! execution time
    !--------------------------------------------------------------------------
    if ( oasis_lag /= 0 ) then
      write(stderr,*) 'error initiating the lag with OASIS'
      call fatal(__FILE__,__LINE__,'SYNC WAIT')
    end if
    ! case number 1: oasis_sync_lag > 0
    ! RegCM actually starts oasis_sync_lag seconds after OASIS
    ! >> at the beginning of the run, OASIS starts while RegCM runs
    ! dummy loops with oasis_lag increasing like 0 -> oasis_sync_lag
    ! at the end of the run, OASIS and RegCM stop synchronously
    if ( oasis_sync_lag > 0 ) then
      do while ( oasis_lag < oasis_sync_lag )
        call oasisxregcm_rcv_all(time + oasis_lag)
        ! dummy loop
#ifdef DEBUG
        if ( oasis_lag == 0 ) then
          write(ndebug,*) oasis_prefix, 'RegCM is creating the lag ', &
                                        'with the other components'
          write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
        end if
#endif
        call oasisxregcm_snd_all(time + oasis_lag)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    ! case number 2: oasis_sync_lag < 0
    ! RegCM starts synchronously with OASIS
    ! however, some other components start with a lag
    ! at the end of the run, OASIS needs oasis_sync_lag seconds to run
    ! with the other components while RegCM has already finished
    ! >> at the end, OASIS keeps going while RegCM runs dummy loops
    ! with oasis_lag increasing like end_time + dt -> -oasis_sync_lag
    else if ( oasis_sync_lag < 0 ) then
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'RegCM is catching up ', &
                                    'with the other components'
      write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      do while ( oasis_lag < -oasis_sync_lag )
        call oasisxregcm_rcv_all(time + int(dtsec,ik4) + oasis_lag)
        ! dummy loop
        call oasisxregcm_snd_all(time + int(dtsec,ik4) + oasis_lag)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    end if
  end subroutine oasisxregcm_sync_wait

  ! call all subroutines linked with some deallocation of variables used in
  ! the OASIS coupling
  subroutine oasisxregcm_release
    implicit none
    !--------------------------------------------------------------------------
    call oasisxregcm_deallocate_field(im_sst, cpl_sst)
!    call oasisxregcm_deallocate_field(im_sit, cpl_sit)
    call oasisxregcm_deallocate_field(im_wz0, cpl_wz0)
    call oasisxregcm_deallocate_field(im_wust, cpl_wust)
    call oasisxregcm_deallocate_field(ex_u10m)
    call oasisxregcm_deallocate_field(ex_v10m)
    call oasisxregcm_deallocate_field(ex_wspd)
    call oasisxregcm_deallocate_field(ex_wdir, cpl_wdir)
    call oasisxregcm_deallocate_field(ex_t2m)
!    call oasisxregcm_deallocate_field(ex_t10m)
    call oasisxregcm_deallocate_field(ex_q2m)
!    call oasisxregcm_deallocate_field(ex_q10m)
    call oasisxregcm_deallocate_field(ex_slp)
    call oasisxregcm_deallocate_field(ex_taux)
    call oasisxregcm_deallocate_field(ex_tauy)
    call oasisxregcm_deallocate_field(ex_z0)
    call oasisxregcm_deallocate_field(ex_ustr)
    call oasisxregcm_deallocate_field(ex_evap)
    call oasisxregcm_deallocate_field(ex_prec)
    call oasisxregcm_deallocate_field(ex_nuwa)
    call oasisxregcm_deallocate_field(ex_ulhf)
    call oasisxregcm_deallocate_field(ex_ushf)
    call oasisxregcm_deallocate_field(ex_uwlw)
    call oasisxregcm_deallocate_field(ex_dwlw)
    call oasisxregcm_deallocate_field(ex_nulw)
    call oasisxregcm_deallocate_field(ex_uwsw)
    call oasisxregcm_deallocate_field(ex_dwsw)
    call oasisxregcm_deallocate_field(ex_ndsw)
#ifdef ECLM
    call oasisxregcm_deallocate_field(ex_tatm)
    call oasisxregcm_deallocate_field(ex_uatm)
    call oasisxregcm_deallocate_field(ex_vatm)
    call oasisxregcm_deallocate_field(ex_qatm)
    call oasisxregcm_deallocate_field(ex_hgt)
    call oasisxregcm_deallocate_field(ex_patm)
    call oasisxregcm_deallocate_field(ex_dwrlwf)
    call oasisxregcm_deallocate_field(ex_swdir)
    call oasisxregcm_deallocate_field(ex_swdif)
    call oasisxregcm_deallocate_field(im_tgbb, cpl_tgbb)
    call oasisxregcm_deallocate_field(im_t2m, cpl_t2m)
    call oasisxregcm_deallocate_field(im_q2m, cpl_q2m)
    call oasisxregcm_deallocate_field(im_u10m, cpl_u10m)
  !  call oasisxregcm_deallocate_field(im_v10m, cpl_v10m)
    call oasisxregcm_deallocate_field(im_sent, cpl_sent)
    call oasisxregcm_deallocate_field(im_evpr, cpl_evpr)
    call oasisxregcm_deallocate_field(im_ram1, cpl_ram1)
    call oasisxregcm_deallocate_field(im_rah1, cpl_rah1)
 !   call oasisxregcm_deallocate_field(im_br, cpl_br)
    call oasisxregcm_deallocate_field(im_taux, cpl_taux)
    call oasisxregcm_deallocate_field(im_tauy, cpl_tauy)
!    call oasisxregcm_deallocate_field(im_drag, cpl_drag)
    call oasisxregcm_deallocate_field(im_sncv, cpl_sncv)
    call oasisxregcm_deallocate_field(im_wt, cpl_wt)
    call oasisxregcm_deallocate_field(im_zo, cpl_zo)
    call oasisxregcm_deallocate_field(im_tlef, cpl_tlef)
    call oasisxregcm_deallocate_field(ex_rainf)
    call oasisxregcm_deallocate_field(ex_snow)
    call oasisxregcm_deallocate_3dfield(im_albei, cpl_albei)
    call oasisxregcm_deallocate_3dfield(im_albed, cpl_albed)
    call oasisxregcm_deallocate_3dfield(im_flxvoc, cpl_flxvoc)
    call oasisxregcm_deallocate_3dfield(im_flxdst, cpl_flxdst)
    call oasisxregcm_deallocate_field(im_fluxch4, cpl_fluxch4)
    call oasisxregcm_deallocate_3dfield(im_ddvel, cpl_ddvel)
    call oasisxregcm_deallocate_field(im_tlai, cpl_tlai)
    call oasisxregcm_deallocate_field(im_roff, cpl_roff)
    call oasisxregcm_deallocate_field(im_srfroff, cpl_srfroff)
    call oasisxregcm_deallocate_field(im_snwmelt, cpl_snwmelt)
    call oasisxregcm_deallocate_3dfield(im_tsoi, cpl_tsoi)
    call oasisxregcm_deallocate_3dfield(im_h2ovol, cpl_h2ovol)
    call oasisxregcm_deallocate_3dfield(im_h2oliq, cpl_h2oliq)
    call oasisxregcm_deallocate_3dfield(im_h2oice, cpl_h2oice)
    call oasisxregcm_deallocate_field(im_h2o10cm, cpl_h2o10cm)
#endif
    ! OASIS field +++
    if ( allocated(grdde) ) deallocate(grdde)
    if ( allocated(grddi) ) deallocate(grddi)
    if ( allocated(grdce) ) deallocate(grdce)
    if ( allocated(grdci) ) deallocate(grdci)
#ifdef ECLM
    if ( allocated(grdtest) ) deallocate(grdtest)
    if ( allocated(grd3d) ) deallocate(grd3d)
#endif
  end subroutine oasisxregcm_release

end module mod_oasis_interface
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2
