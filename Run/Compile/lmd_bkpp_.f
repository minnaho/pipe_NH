      subroutine lmd_bkpp_tile (istr,iend,jstr,jend, Kv,Kt,Ks,
     &                                 ws,wm, my_hbbl,wrk, Cr,
     &             Gm1,dGm1dS, Gt1,dGt1dS, Gs1,dGs1dS, my_kbbl)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=1024,  MMm0=1024,  N=128)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=16,  NP_ETA=16,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 stdout, Np, padd_X,padd_E
      parameter (stdout=6, Np=N+1)
      parameter (Lm=(LLm+NP_XI-1)/NP_XI, Mm=(MMm+NP_ETA-1)/NP_ETA)
      parameter (padd_X=(Lm+2)/2-(Lm+1)/2)
      parameter (padd_E=(Mm+2)/2-(Mm+1)/2)
      integer*4 NSA, N2d,N3d, size_XI,size_ETA
      integer*4 se,sse, sz,ssz
      parameter (NSA=28)
      parameter (size_XI=7+(Lm+NSUB_X-1)/NSUB_X)
      parameter (size_ETA=7+(Mm+NSUB_E-1)/NSUB_E)
      parameter (sse=size_ETA/Np, ssz=Np/size_ETA)
      parameter (se=sse/(sse+ssz), sz=1-se)
      parameter (N2d=size_XI*(se*size_ETA+sz*Np))
      parameter (N3d=size_XI*size_ETA*Np)
      real Vtransform
      parameter (Vtransform=2)
      integer*4   NT, itemp
      integer*4   ntrc_salt, ntrc_pas, ntrc_bio, ntrc_sed
      parameter (itemp=1)
      parameter (ntrc_salt=1)
      parameter (ntrc_pas=0)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_surf
     &          , isalt
      parameter (isalt=itemp+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_surf=0)
      real h(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hinv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real f(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real fomn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_angler/angler
      real latr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonu(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real latv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real lonv(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /grid_latr/latr /grid_lonr/lonr
      common /grid_latu/latu /grid_lonu/lonu
      common /grid_latv/latv /grid_lonv/lonv
      real pm(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real om_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real on_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pm_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pn_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real dndx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmon_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pnom_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real grdscl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real umask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real vmask(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real pmask2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3)
      real t(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hz_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real z_w(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Huon(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real Hvom(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Wi(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real rho1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real rho(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real sustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real svstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real sustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real svstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /smsdat_sustrg/sustrg /smsdat_svstrg/svstrg
      real    sustrp(2), svstrp(2), sms_time(2)
      real    sms_cycle, sms_scale
      integer*4 itsms, sms_ncycle, sms_rec, lsusgrd
      integer*4 lsvsgrd,sms_tid, susid, svsid
      common /smsdat1/ sustrp, svstrp, sms_time
      common /smsdat2/ sms_cycle, sms_scale
      common /smsdat3/ itsms, sms_ncycle, sms_rec, lsusgrd
      common /smsdat4/ lsvsgrd,sms_tid, susid, svsid
      real bustr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real bvstr(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real bustrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real bvstrg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /bmsdat_bustrg/bustrg /bmsdat_bvstrg/bvstrg
      real bms_tintrp(2), bustrp(2),    bvstrp(2), tbms(2)
      real bmsclen, bms_tstart, bms_tend,  tsbms, sclbms
      integer*4 itbms,      bmstid,busid, bvsid,     tbmsindx
      logical bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      common /bmsdat1/bms_tintrp, bustrp,       bvstrp,    tbms
      common /bmsdat2/bmsclen,    bms_tstart,   bms_tend,  tsbms,   
     &                                                            sclbms
      common /bmsdat3/itbms,      bmstid,busid, bvsid,     tbmsindx
      common /bmsdat4/bmscycle,   bms_onerec,   lbusgrd,   lbvsgrd
      real stflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_stflx/stflx
      real stflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2,NT)
      common /stfdat_stflxg/stflxg
      real stflxp(2,NT), stf_time(2,NT)
      real stf_cycle(NT), stf_scale(NT)
      integer*4 itstf(NT), stf_ncycle(NT), stf_rec(NT)
      integer*4 lstfgrd(NT), stf_tid(NT), stf_id(NT)
      common /stfdat1/ stflxp,  stf_time, stf_cycle, stf_scale
      common /stfdat2/ itstf, stf_ncycle, stf_rec, lstfgrd
      common /stfdat3/  stf_tid, stf_id
      real btflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /forces_btflx/btflx
      real dqdt(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real sst(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_dqdt/dqdt /forces_sst/sst
      real dqdtg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real sstg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /sstdat_dqdtg/dqdtg /sstdat_sstg/sstg
      real    sstp(2), dqdtp(2), sst_time(2)
      real    sst_cycle, scldqdt
      integer*4 itsst, sst_ncycle, sst_rec,  sst_tid,  sst_id
      integer*4 dqdt_id,     lsstgrd,   sstunused
      common /sstdat1/ sstp, dqdtp, sst_time
      common /sstdat2/ sst_cycle, scldqdt
      common /sstdat3/ itsst, sst_ncycle, sst_rec, sst_tid, sst_id
      common /sstdat4/ dqdt_id, lsstgrd, sstunused
      real sss(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_sss/sss
      real sssg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /sssdat_sssg/sssg
      real sssp(2),  sss_time(2)
      real sss_cycle
      integer*4 itsss, sss_ncycle, sss_rec,  sss_tid,  sss_id
      integer*4 lsssgrd,   sssunused
      common /sssdat1/sssp,  sss_time, sss_cycle
      common /sssdat2/itsss, sss_ncycle, sss_rec,  sss_tid, sss_id
      common /sssdat3/lsssgrd,   sssunused
      real srflx(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /forces_srflx/srflx
      real srflxg(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /srfdat_srflxg/srflxg
      real srflxp(2),srf_time(2)
      real srf_cycle, srf_scale
      integer*4 itsrf, srf_ncycle, srf_rec
      integer*4 lsrfgrd, srf_tid, srf_id
      common /srfdat1/ srflxp, srf_time, srf_cycle, srf_scale
      common /srfdat2/ itsrf, srf_ncycle, srf_rec, lsrfgrd, srf_tid, 
     &                                                            srf_id
      real visc2_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_r(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real visc2_sponge_p(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff2(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real diff4_sponge(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real diff4(-1:Lm+2+padd_X,-1:Mm+2+padd_E,NT)
      common /mixing_diff4_sponge/diff4_sponge
      common /mixing_diff4/diff4
      real diff3d_u(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real diff3d_v(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      common /mixing_diff3d_u/diff3d_u
      common /mixing_diff3d_v/diff3d_v
      real dRdx(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real dRde(-1:Lm+2+padd_X,-1:Mm+2+padd_E,N)
      real idRz(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_dRdx/dRdx
      common /mixing_dRde/dRde
      common /mixing_idRz/idRz
      real Rslope_max,Gslope_max
      parameter (Gslope_max=5.D0, Rslope_max=0.05D0)
      integer*4 ismooth
      real csmooth
      common /mixing_csmooth/ csmooth
      common /mixing_ismooth/ ismooth
      real Akv(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      real Akt(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real bvf(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /mixing_bvf/ bvf
      integer*4 kbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      integer*4 kbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real hbbl(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_kbl/ kbl
      common /lmd_kpp_hbbl/ hbbl
      common /lmd_kpp_kbbl/ kbbl
      real hbls(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /lmd_kpp_hbl/ hbls
      real ghats(-1:Lm+2+padd_X,-1:Mm+2+padd_E,0:N)
      common /lmd_kpp_ghats/ghats
      real ustar(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /lmd_kpp_ustar/ustar
      real dt, dtfast, time, time2, time_start, tdays
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &                       ndtfast, iic, kstp, krhs, knew, next_kstp,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       PREDICTOR_2D_STEP
      real time_avg, time2_avg, rho0
     &               , rdrg, rdrg2, Cdb_min, Cdb_max, Zob
     &               , xl, el, visc2, visc4, gamma2
      real  theta_s,   theta_b,   Tcline,  hc
      real  sc_w(0:N), Cs_w(0:N), sc_r(N), Cs_r(N)
      real  rx0, rx1
      real  tnu2(NT),tnu4(NT)
      real weight(6,0:NWEIGHT)
      real  x_sponge,   v_sponge
       real  tauT_in, tauT_out, tauM_in, tauM_out
      integer*4 numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
      logical ldefhis
      logical got_tini(NT)
      common /scalars_main/
     &             time_avg, time2_avg,  rho0,      rdrg,    rdrg2
     &           , Zob,       Cdb_min,   Cdb_max
     &           , xl, el,    visc2,     visc4,   gamma2
     &           , theta_s,   theta_b,   Tcline,  hc
     &           , sc_w,      Cs_w,      sc_r,    Cs_r
     &           , rx0,       rx1,       tnu2,    tnu4
     &                      , weight
     &                      , x_sponge,   v_sponge
     &                      , tauT_in, tauT_out, tauM_in, tauM_out
     &      , numthreads,     ntstart,   ntimes,  ninfo
     &      , nfast,  nrrec,     nrst,    nwrt
     &                                 , ntsavg,  navg
     &                      , got_tini
     &                      , ldefhis
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      real lonmin, lonmax, latmin, latmax
      common /communicators_lonlat/
     &     lonmin, lonmax, latmin, latmax
      real*8 volume, avgke, avgpe, avgkp, bc_crss
      common /communicators_rq/
     &          volume, avgke, avgpe, avgkp, bc_crss
      real*4 CPU_time(0:31,0:NPP)
      integer*4 proc(0:31,0:NPP),trd_count
      common /timers_roms/CPU_time,proc,trd_count
      logical EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      integer*4 mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE, p_NW,p_NE
      common /comm_setup/ mynode, ii,jj, p_W,p_E,p_S,p_N, p_SW,p_SE,
     &  p_NW,p_NE, EAST_INTER, WEST_INTER, NORTH_INTER, SOUTH_INTER
      real pi, deg2rad, rad2deg
      parameter (pi=3.14159265358979323846D0, deg2rad=pi/180.D0,
     &                                      rad2deg=180.D0/pi)
      real Eradius, g, day2sec,sec2day, jul_off,
     &     year2day,day2year
      parameter (Eradius=6371315.0D0,  day2sec=86400.D0,
     &           sec2day=1.D0/86400.D0, jul_off=2440000.D0,
     &           year2day=365.25D0, day2year=1.D0/365.25D0)
      parameter (g=9.81D0)
      real Cp
      parameter (Cp=3985.0D0)
      real vonKar
      parameter (vonKar=0.41D0)
      real spval
      parameter (spval=-9999.0D0)
      logical mask_val
      parameter (mask_val = .true.)
      real rhoA(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rhoS(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rvfrc(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real rufrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      real rvfrc_bak(-1:Lm+2+padd_X,-1:Mm+2+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DU_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DV_avg1(-1:Lm+2+padd_X,-1:Mm+2+padd_E,5)
      real DU_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      real DV_avg2(-1:Lm+2+padd_X,-1:Mm+2+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
      integer*4 istr,iend,jstr,jend, i,j,k
     &      , imin,imax,jmin,jmax
      real Kv(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     Kt(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     Ks(Istr-2:Iend+2,Jstr-2:Jend+2,0:N),
     &     ws(Istr-2:Iend+2,Jstr-2:Jend+2), zscale,
     &     wm(Istr-2:Iend+2,Jstr-2:Jend+2),
     &my_hbbl(Istr-2:Iend+2,Jstr-2:Jend+2),
     &    wrk(Istr-2:Iend+2,Jstr-2:Jend+2),
     &     Cr(Istr-2:Iend+2,0:N), z_bl, zsbl,
     &    Gm1(Istr-2:Iend+2),   Av_bl,
     & dGm1dS(Istr-2:Iend+2),  dAv_bl,
     &    Gt1(Istr-2:Iend+2),   At_bl, a1,
     & dGt1dS(Istr-2:Iend+2),  dAt_bl, a2,
     &    Gs1(Istr-2:Iend+2),   As_bl, a3,
     & dGs1dS(Istr-2:Iend+2),  dAs_bl
      integer*4 my_kbbl(Istr-2:Iend+2)
      real Kern, sigma, cff,cff1, cff_up,cff_dn, lmd_nu0c, Ricr,
     &     Ri_inv, C_Ek, eps, invrho, Kv0, Kt0, Ks0
      parameter (
     &   lmd_nu0c=0.1D0,
     &   Ricr=0.45D0,
     &   Ri_inv=1.D0/Ricr,
     &   C_Ek=258.D0,
     &   eps=1.D-20 )
      if (.not.WEST_INTER) then
        imin=istr
      else
        imin=istr-1
      endif
      if (.not.EAST_INTER) then
        imax=iend
      else
        imax=iend+1
      endif
      if (.not.SOUTH_INTER) then
        jmin=jstr
      else
        jmin=jstr-1
      endif
      if (.not.NORTH_INTER) then
        jmax=jend
      else
        jmax=jend+1
      endif
      do j=jmin,jmax
        do i=imin,imax
          ustar(i,j)=sqrt(sqrt( (0.5D0*(bustr(i,j)+bustr(i+1,j)))**2
     &                         +(0.5D0*(bvstr(i,j)+bvstr(i,j+1)))**2))
        enddo
      enddo
      do j=jmin,jmax
        do i=imin,imax
          wm(i,j)=vonKar*ustar(i,j)
          ws(i,j)=wm(i,j)
          my_kbbl(i)=N
          Cr(i,0)=0.D0
        enddo
        do k=1,N-1,+1
          do i=imin,imax
            zscale=z_r(i,j,k)-z_w(i,j,0)
            Kern=zscale/(zscale+Zob)
            Cr(i,k)=Cr(i,k-1) + Kern*(
     &                 0.5D0*( ( u(i  ,j,k+1,nstp)-u(i  ,j,k,nstp)
     &                        +u(i+1,j,k+1,nstp)-u(i+1,j,k,nstp) )**2
     &                      +( v(i,j  ,k+1,nstp)-v(i,j  ,k,nstp)
     &                        +v(i,j+1,k+1,nstp)-v(i,j+1,k,nstp) )**2
     &                      )/(Hz(i,j,k)+Hz(i,j,k+1))
     &               -0.5D0*(Hz(i,j,k)+Hz(i,j,k+1))*( Ri_inv*bvf(i,j,k)
     &                                            +C_Ek*f(i,j)*f(i,j)
     &                                                             ))
          enddo
        enddo
        do i=imin,imax
          Cr(i,N)=2.D0*Cr(i,N-1) -Cr(i,N-2)
        enddo
        do k=1,N,+1
          do i=imin,imax
            if (my_kbbl(i).eq.N .and. Cr(i,k).lt.0.D0) my_kbbl(i)=k
          enddo
        enddo
        do i=imin,imax
          my_hbbl(i,j)=z_w(i,j,N)-z_w(i,j,0)
          if (my_kbbl(i).lt.N) then
            k=my_kbbl(i)
            if (k.eq.1) then
              my_hbbl(i,j)=z_r(i,j,1)-z_w(i,j,0)
            else
              my_hbbl(i,j)=( z_r(i,j,k)*Cr(i,k-1)-z_r(i,j,k-1)*Cr(i,k)
     &                            )/(Cr(i,k-1)-Cr(i,k)) -z_w(i,j,0)
            endif
          endif
        enddo
      enddo
      if (.not.WEST_INTER) then
        do j=jmin,jmax
          my_hbbl(Istr-1,j)=my_hbbl(Istr,j)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=jmin,jmax
          my_hbbl(Iend+1,j)=my_hbbl(Iend,j)
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=imin,imax
          my_hbbl(i,Jstr-1)=my_hbbl(i,Jstr)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=imin,imax
          my_hbbl(i,Jend+1)=my_hbbl(i,Jend)
        enddo
      endif
      if (.not.WEST_INTER.and..not.SOUTH_INTER) then
        my_hbbl(Istr-1,Jstr-1)=my_hbbl(Istr,Jstr)
      endif
      if (.not.WEST_INTER.and..not.NORTH_INTER) then
        my_hbbl(Istr-1,Jend+1)=my_hbbl(Istr,Jend)
      endif
      if (.not.EAST_INTER.and..not.SOUTH_INTER) then
        my_hbbl(Iend+1,Jstr-1)=my_hbbl(Iend,Jstr)
      endif
      if (.not.EAST_INTER.and..not.NORTH_INTER) then
        my_hbbl(Iend+1,Jend+1)=my_hbbl(Iend,Jend)
      endif
      do j=Jstr,Jend+1
        do i=Istr,Iend+1
          wrk(i,j)=0.25D0*(my_hbbl(i,j)  +my_hbbl(i-1,j)
     &                  +my_hbbl(i,j-1)+my_hbbl(i-1,j-1))
     &                   *pmask2(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          cff=0.25D0*(pmask2(i,j)   +pmask2(i+1,j)
     &             +pmask2(i,j+1) +pmask2(i+1,j+1))
          my_hbbl(i,j)=(1.D0-cff)*my_hbbl(i,j)+
     &              0.25D0*(wrk(i,j)  +wrk(i+1,j)
     &                   +wrk(i,j+1)+wrk(i+1,j+1))
          my_hbbl(i,j)=my_hbbl(i,j)*rmask(i,j)
        enddo
      enddo
      do j=Jstr,Jend
        do i=Istr,Iend
          hbbl(i,j)=min(my_hbbl(i,j),z_w(i,j,N)-z_w(i,j,0))
     &                                          *rmask(i,j)
          kbbl(i,j)=N
        enddo
        do k=N-1,1,-1
          do i=Istr,Iend
            if (z_r(i,j,k)-z_w(i,j,0).gt.hbbl(i,j)) kbbl(i,j)=k
          enddo
        enddo
        do i=Istr,Iend
          k=kbbl(i,j)
          z_bl=z_w(i,j,0)+hbbl(i,j)
          if (z_bl.lt.z_w(i,j,k-1)) k=k-1
          cff=1.D0/(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)
          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbbl(i,j)*wm(i,j)+eps)
          dGm1dS(i)=min(0.D0, -dAv_bl/(wm(i,j)+eps))
          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbbl(i,j)*ws(i,j)+eps)
          dGt1dS(i)=min(0.D0, -dAt_bl/(ws(i,j)+eps))
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbbl(i,j)*ws(i,j)+eps)
          dGs1dS(i)=min(0.D0, -dAs_bl/(ws(i,j)+eps))
        enddo
        do i=Istr,Iend
          do k=1,N-1
            if (k.lt.kbbl(i,j)) then
              sigma=min((z_w(i,j,k)-z_w(i,j,0))/(hbbl(i,j)+eps),1.D0)
              a1=sigma-2.D0
              a2=3.D0-2.D0*sigma
              a3=sigma-1.D0
              Kv0 =wm(i,j)*hbbl(i,j)*( sigma*( 1.D0+sigma*(
     &                            a1+a2*Gm1(i)+a3*dGm1dS(i) )))
              Kt0 =ws(i,j)*hbbl(i,j)*( sigma*( 1.D0+sigma*(
     &                            a1+a2*Gt1(i)+a3*dGt1dS(i) )))
              Ks0 =ws(i,j)*hbbl(i,j)*( sigma*( 1.D0+sigma*(
     &                            a1+a2*Gs1(i)+a3*dGs1dS(i) )))
              zsbl=z_w(i,j,N)-hbls(i,j,3-nstp)
              if (z_w(i,j,k).gt.zsbl) then
                Kv0=max(Kv(i,j,k),Kv0)
                Kt0=max(Kt(i,j,k),Kt0)
                Ks0=max(Ks(i,j,k),Ks0)
              endif
              Kv(i,j,k)=Kv0
              Kt(i,j,k)=Kt0
              Ks(i,j,k)=Ks0
            else
               if (bvf(i,j,k).lt.0.D0) then
                zsbl=z_w(i,j,N)-hbls(i,j,3-nstp)
                if (z_w(i,j,k).lt.zsbl) then
                  Kv(i,j,k)=Kv(i,j,k)+lmd_nu0c
                  Kt(i,j,k)=Kt(i,j,k)+lmd_nu0c
                  Ks(i,j,k)=Ks(i,j,k)+lmd_nu0c
                endif
              endif
            endif
          enddo
        enddo
      enddo
      if (.not.WEST_INTER) then
        do j=jstr,jend
          hbbl(istr-1,j)=hbbl(istr,j)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=jstr,jend
          hbbl(iend+1,j)=hbbl(iend,j)
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do i=istr,iend
          hbbl(i,jstr-1)=hbbl(i,jstr)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=istr,iend
          hbbl(i,jend+1)=hbbl(i,jend)
        enddo
      endif
      if (.not.WEST_INTER .and. .not.SOUTH_INTER) then
        hbbl(istr-1,jstr-1)=hbbl(istr,jstr)
      endif
      if (.not.WEST_INTER .and. .not.NORTH_INTER) then
        hbbl(istr-1,jend+1)=hbbl(istr,jend)
      endif
      if (.not.EAST_INTER .and. .not.SOUTH_INTER) then
        hbbl(iend+1,jstr-1)=hbbl(iend,jstr)
      endif
      if (.not.EAST_INTER .and. .not.NORTH_INTER) then
        hbbl(iend+1,jend+1)=hbbl(iend,jend)
      endif
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,hbbl(-1,-1))
      return
      end
