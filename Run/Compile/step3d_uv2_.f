      subroutine step3d_uv2 (tile)
      implicit none
      integer*4 tile, trd, omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=1024,  MMm0=512,  N=64)
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
      integer*4 Msrc
      parameter (Msrc=15000)
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
      parameter (ntrc_pas=1)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_surf
     &          , isalt
     &          , itpas
      parameter (isalt=itemp+1)
      parameter (itpas=itemp+ntrc_salt+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_surf=0)
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,5,0:NPP-1)
      common /private_scratch/ A2d,A3d
      integer*4 chunk_size_X,margin_X,chunk_size_E,margin_E
      integer*4 Istr,Iend,Jstr,Jend, i_X,j_E
      chunk_size_X=(Lmmpi+NSUB_X-1)/NSUB_X
      margin_X=(NSUB_X*chunk_size_X-Lmmpi)/2
      chunk_size_E=(Mmmpi+NSUB_E-1)/NSUB_E
      margin_E=(NSUB_E*chunk_size_E-Mmmpi)/2
      j_E=tile/NSUB_X
      i_X=tile-j_E*NSUB_X
      Istr=1+i_X*chunk_size_X-margin_X
      Iend=Istr+chunk_size_X-1
      Istr=max(Istr,1)
      Iend=min(Iend,Lmmpi)
      Jstr=1+j_E*chunk_size_E-margin_E
      Jend=Jstr+chunk_size_E-1
      Jstr=max(Jstr,1)
      Jend=min(Jend,Mmmpi)
      trd=omp_get_thread_num()
      call step3d_uv2_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd),
     &                                  A2d(1,2,trd), A2d(1,3,trd),
     &                                                A2d(1,4,trd)
     &                                               ,A2d(1,5,trd)
     &                                                            )
      return
      end
      subroutine step3d_uv2_tile (Istr,Iend,Jstr,Jend, BC,CF,FC,DC
     &                                                         ,WC
     &                                                            )
      use nhmg, only : nhmg_solve
      use mg_grids
      use mg_tictoc, only : tic, toc
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=1024,  MMm0=512,  N=64)
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
      integer*4 Msrc
      parameter (Msrc=15000)
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
      parameter (ntrc_pas=1)
      parameter (ntrc_bio=0)
      parameter (ntrc_sed=0)
      parameter (NT=itemp+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed)
      integer*4   ntrc_diats, ntrc_diauv, ntrc_diabio
      integer*4   ntrc_diavrt, ntrc_diaek, ntrc_surf
     &          , isalt
     &          , itpas
      parameter (isalt=itemp+1)
      parameter (itpas=itemp+ntrc_salt+1)
      parameter (ntrc_diabio=0)
      parameter (ntrc_diats=0)
      parameter (ntrc_diauv=0)
      parameter (ntrc_diavrt=0)
      parameter (ntrc_diaek=0)
      parameter (ntrc_surf=0)
      real h(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real hinv(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real f(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real fomn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_h/h /grid_hinv/hinv /grid_f/f /grid_fomn/fomn
      real angler(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_angler/angler
      real xp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real xr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yp(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real yr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /grid_xr/xr /grid_xp/xp /grid_yp/yp /grid_yr/yr
      real pm(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real om_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real on_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pm_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pn_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pm/pm    /metrics_pn/pn
      common /metrics_omr/om_r /metrics_on_r/on_r
      common /metrics_omu/om_u /metrics_on_u/on_u
      common /metrics_omv/om_v /metrics_on_v/on_v
      common /metrics_omp/om_p /metrics_on_p/on_p
      common /metrics_pnu/pn_u /metrics_pmv/pm_v
      common /metrics_pmu/pm_u /metrics_pnv/pn_v
      real dmde(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real dndx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_dmde/dmde    /metrics_dndx/dndx
      real pmon_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmon_u(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pnom_v(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real grdscl(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /metrics_pmon_p/pmon_p /metrics_pnom_p/pnom_p
      common /metrics_pmon_r/pmon_r /metrics_pnom_r/pnom_r
      common /metrics_pmon_u/pmon_u /metrics_pnom_v/pnom_v
      common /metrics_grdscl/grdscl
      real rmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real umask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real vmask(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real pmask2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mask_r/rmask
      common /mask_p/pmask
      common /mask_u/umask
      common /mask_v/vmask
      common /mask_p2/pmask2
      real zeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real ubar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      real vbar(-2:Lm+3+padd_X,-2:Mm+3+padd_E,4)
      common /ocean_zeta/zeta
      common /ocean_ubar/ubar
      common /ocean_vbar/vbar
      real nh_ubar(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real nh_vbar(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real nh_wcor(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /nh_wcor/nh_ubar,nh_vbar,nh_wcor
      real u(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real v(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3)
      real t(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N,3,NT)
      common /ocean_u/u /ocean_v/v /ocean_t/t
      real Hz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hz_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real z_w(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Huon(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real Hvom(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /grid_Hz_bak/Hz_bak /grid_zw/z_w /grid_Huon/Huon
      common /grid_Hvom/Hvom
      real We(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Wi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      common /grid_Hz/Hz /grid_zr/z_r /grid_We/We
      common /grid_Wi/Wi
      real wz(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,3)
      real nhdu(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N,2)
      real nhdw(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      real dzdxi(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real dzdeta(-2:Lm+3+padd_X,-2:Mm+3+padd_E,1:N)
      real Hz_half(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_wz/wz
      common /ocean_nhdu/nhdu
      common /ocean_nhdv/nhdv
      common /ocean_nhdw/nhdw
      common /ocean_dzdxi/dzdxi
      common /ocean_dzdeta/dzdeta
      common /grid_Hz_half/Hz_half
      real rho1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      real rho(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_rho1/rho1 /ocean_rho/rho
      real qp1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /ocean_qp1/qp1
      real qp2
      parameter (qp2=0.0000172D0)
      real rhoA(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rhoS(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /coup_rhoA/rhoA           /coup_rhoS/rhoS
      real rufrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rvfrc(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real rufrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      real rvfrc_bak(-2:Lm+3+padd_X,-2:Mm+3+padd_E,2)
      common /coup_rufrc/rufrc
      common /coup_rvfrc/rvfrc
      common /coup_rufrc_bak/rufrc_bak
      common /coup_rvfrc_bak/rvfrc_bak
      real Zt_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DU_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DV_avg1(-2:Lm+3+padd_X,-2:Mm+3+padd_E,5)
      real DU_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real DV_avg2(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /ocean_Zt_avg1/Zt_avg1
      common /coup_DU_avg1/DU_avg1
      common /coup_DV_avg1/DV_avg1
      common /coup_DU_avg2/DU_avg2
      common /coup_DV_avg2/DV_avg2
      real sustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real svstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_sustr/sustr /forces_svstr/svstr
      real bustr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real bvstr(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_bustr/bustr /forces_bvstr/bvstr
      real stflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_stflx/stflx
      real btflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /forces_btflx/btflx
      real srflx(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /forces_srflx/srflx
      real visc2_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real visc2_sponge_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      common /mixing_visc2_r/visc2_r /mixing_visc2_p/visc2_p
      common /mixing_visc2_sponge_r/visc2_sponge_r
      common /mixing_visc2_sponge_p/visc2_sponge_p
      real diff2_sponge(-2:Lm+3+padd_X,-2:Mm+3+padd_E)
      real diff2(-2:Lm+3+padd_X,-2:Mm+3+padd_E,NT)
      common /mixing_diff2_sponge/diff2_sponge
      common /mixing_diff2/diff2
      real visc3d_r(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_visc3d_r/visc3d_r
      real visc3d_p(-2:Lm+3+padd_X,-2:Mm+3+padd_E,N)
      common /mixing_visc3d_p/visc3d_p
      real Akv(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N)
      real Akt(-2:Lm+3+padd_X,-2:Mm+3+padd_E,0:N,2)
      common /mixing_Akv/Akv /mixing_Akt/Akt
      real dt, dtfast, time, time2, time_start, tdays
      integer*4 ndtfast, iic, kstp, krhs, knew, next_kstp
     &      , iif, nstp, nrhs, nnew, nbstep3d
     &      , iprec1, iprec2
      logical PREDICTOR_2D_STEP
      common /time_indices/  dt,dtfast, time, time2,time_start, tdays,
     &                       ndtfast, iic, kstp, krhs, knew, next_kstp,
     &                       iif, nstp, nrhs, nnew, nbstep3d,
     &                       iprec1, iprec2,
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
      real Akv_bak
      real Akt_bak(NT)
      common /scalars_akt/ Akv_bak, Akt_bak
      logical synchro_flag
      common /sync_flag/ synchro_flag
      integer*4 may_day_flag
      integer*4 tile_count, first_time, bc_count
      common /communicators_i/
     &        may_day_flag, tile_count, first_time, bc_count
      real hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
      common /communicators_r/
     &     hmin, hmax, grdmin, grdmax, Cu_min, Cu_max
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
      real Qbar(Msrc)
      common /sources_Qbar/ Qbar
      real Qsrc(Msrc,N)
      common /source_Qsrc/ Qsrc
      real Qshape(Msrc,N)
      common /source_Qshape/ Qshape
      real Tsrc(Msrc,N,NT)
      common /source_Tsrc/ Tsrc
      real Tsrc0(Msrc,NT)
      common /source_Tsrc0/ Tsrc0
      real lasrc(Msrc)
      common /source_lasrc/ lasrc
      real losrc(Msrc)
      common /source_losrc/ losrc
      integer*4 Nsrc
      common /source_Nsrc/ Nsrc
      integer*4 Dsrc(Msrc)
      common /source_Dsrc/ Dsrc
      integer*4 Isrc(Msrc)
      common /source_Isrc/ Isrc
      integer*4 Jsrc(Msrc)
      common /source_Jsrc/ Jsrc
      logical Lsrc(Msrc,30)
      common /source_Lsrc/ Lsrc
      integer*4 Isrc_mpi(Msrc,0:NNODES-1)
      common /source_Isrc_mpi/ Isrc_mpi
      integer*4 Jsrc_mpi(Msrc,0:NNODES-1)
      common /source_Jsrc_mpi/ Jsrc_mpi
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
     &       ,is
      integer*4 halo
      real BC(Istr-2:Iend+2,0:N),
     &     CF(Istr-2:Iend+2,0:N),
     &     FC(Istr-2:Iend+2,0:N), cff,
     &     DC(Istr-2:Iend+2,0:N)
      real WC(Istr-2:Iend+2,0:N)
      real grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real dpth,aa,cc
      integer*4 IstrR,IendR,JstrR,JendR
      integer*4 IstrU
      integer*4 JstrV
      if (.not.WEST_INTER) then
        IstrR=Istr-1
        IstrU=Istr+1
      else
        IstrR=Istr
        IstrU=Istr
      endif
      if (.not.EAST_INTER) then
        IendR=Iend+1
      else
        IendR=Iend
      endif
      if (.not.SOUTH_INTER) then
        JstrR=Jstr-1
        JstrV=Jstr+1
      else
        JstrR=Jstr
        JstrV=Jstr
      endif
      if (.not.NORTH_INTER) then
        JendR=Jend+1
      else
        JendR=Jend
      endif
      do j=Jstr,Jend
        do i=IstrU,Iend
          FC(i,0)=0.D0
          FC(i,N)=0.D0
          DC(i,0)=0.25D0*dt*(pm(i,j)+pm(i-1,j))*(pn(i,j)+pn(i-1,j))
          WC(i,0)=0.D0
          WC(i,N)=0.D0
        enddo
        do k=1,N-1
          do i=IstrU,Iend
            FC(i,k)=-dt*(Akv(i,j,k)+Akv(i-1,j,k))
     &                 /( z_r(i,j,k+1)+z_r(i-1,j,k+1)
     &                   -z_r(i,j,k  )-z_r(i-1,j,k  ))
            WC(i,k)=DC(i,0)*0.5D0*(Wi(i,j,k)+Wi(i-1,j,k))
          enddo
        enddo
        do k=1,N
          do i=IstrU,Iend
            DC(i,k)=u(i,j,k,nnew)
            BC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))-FC(i,k)-FC(i,k-1)
     &              +max(WC(i,k  ),0.D0)-min(WC(i,k-1),0.D0)
          enddo
        enddo
        do i=IstrU,Iend
          DC(i,N)=DC(i,N)+dt*sustr(i,j)
          DC(i,1)=DC(i,1)-dt*bustr(i,j)
        enddo
        do i=IstrU,Iend
          cff=1.D0/BC(i,1)
          CF(i,1)=cff*( FC(i,1)
     &                +min( WC(i,1),0.D0)
     &                                 )
          DC(i,1)=cff*DC(i,1)
        enddo
        do k=2,N-1
          do i=IstrU,Iend
              aa = FC(i,k-1)-max(WC(i,k-1),0.D0)
              cc = FC(i,k  )+min(WC(i,k  ),0.D0)
            cff=1.D0/(BC(i,k)-aa*CF(i,k-1))
            CF(i,k)=cff*cc
            DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
          enddo
        enddo
        do i=IstrU,Iend
              aa = FC(i,N-1)-max(WC(i,N-1),0.D0)
          DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &           /(BC(i,N)-aa*CF(i,N-1))
          CF(i,0)=0.5D0*(Hz(i,j,N)+Hz(i-1,j,N))
          DC(i,0)=CF(i,0)*DC(i,N)
        enddo
        do k=N-1,1,-1
          do i=IstrU,Iend
            DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
            cff=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))
            CF(i,0)=CF(i,0)+cff
            DC(i,0)=DC(i,0)+cff*DC(i,k)
          enddo
        enddo
        do i=IstrU,Iend
          DC(i,0)=(DC(i,0)*on_u(i,j)-DU_avg1(i,j,nnew))
     &                             /(CF(i,0)*on_u(i,j))
        enddo
        do k=1,N
          do i=IstrU,Iend
            u(i,j,k,nnew)=DC(i,k) * umask(i,j)
          enddo
        enddo
        if (j.ge.JstrV) then
          do i=Istr,Iend
            FC(i,0)=0.D0
            FC(i,N)=0.D0
            DC(i,0)=0.25D0*dt*(pm(i,j)+pm(i,j-1))*(pn(i,j)+pn(i,j-1))
            WC(i,0)=0.D0
            WC(i,N)=0.D0
          enddo
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=-dt*(Akv(i,j,k)+Akv(i,j-1,k))
     &                   /( z_r(i,j,k+1)+z_r(i,j-1,k+1)
     &                     -z_r(i,j,k  )-z_r(i,j-1,k  ))
              WC(i,k)=DC(i,0)*0.5D0*(Wi(i,j,k)+Wi(i,j-1,k))
            enddo
          enddo
          do k=1,N
            do i=Istr,Iend
              DC(i,k)=v(i,j,k,nnew)
              BC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))-FC(i,k)-FC(i,k-1)
     &              +max(WC(i,k  ),0.D0)-min(WC(i,k-1),0.D0)
            enddo
          enddo
          do i=Istr,Iend
            DC(i,N)=DC(i,N)+dt*svstr(i,j)
            DC(i,1)=DC(i,1)-dt*bvstr(i,j)
          enddo
          do i=Istr,Iend
            cff=1.D0/BC(i,1)
            CF(i,1)=cff*( FC(i,1)
     &                +min( WC(i,1),0.D0)
     &                                 )
            DC(i,1)=cff*DC(i,1)
          enddo
          do k=2,N-1
            do i=Istr,Iend
              aa = FC(i,k-1)-max(WC(i,k-1),0.D0)
              cc = FC(i,k  )+min(WC(i,k  ),0.D0)
              cff=1.D0/(BC(i,k)-aa*CF(i,k-1))
              CF(i,k)=cff*cc
              DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
            enddo
          enddo
          do i=Istr,Iend
              aa = FC(i,N-1)-max(WC(i,N-1),0.D0)
            DC(i,N)=(DC(i,N)-aa*DC(i,N-1))
     &             /(BC(i,N)-aa*CF(i,N-1))
            CF(i,0)=0.5D0*(Hz(i,j,N)+Hz(i,j-1,N))
            DC(i,0)=CF(i,0)*DC(i,N)
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              DC(i,k)=DC(i,k)-CF(i,k)*DC(i,k+1)
              cff=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))
              CF(i,0)=CF(i,0)+cff
              DC(i,0)=DC(i,0)+cff*DC(i,k)
            enddo
          enddo
          do i=Istr,Iend
            DC(i,0)=(DC(i,0)*om_v(i,j)-DV_avg1(i,j,nnew))
     &                               /(CF(i,0)*om_v(i,j))
          enddo
          do k=1,N
            do i=Istr,Iend
              v(i,j,k,nnew)=DC(i,k) * vmask(i,j)
            enddo
          enddo
        endif
        call tic(1,'croco_step3d_uv2_1')
        do k=1,N
          do i=Istr,Iend
            FC(i,k)=-0.5D0*dt*( Akv(i,j,k)+Akv(i,j,k-1) )
     &                 /    ( z_w(i,j,k)-z_w(i,j,k-1) )
            DC(i,0)=pm(i,j)*pn(i,j)
            WC(i,k)=dt*DC(i,0)*0.5D0*(Wi(i,j,k)+Wi(i,j,k-1))
          enddo
        enddo
        do i=Istr,Iend
          FC(i,0)=0.D0
          WC(i,0)=0.D0
        enddo
        do k=1,N-1
          do i=Istr,Iend
            DC(i,k)=wz(i,j,k,nnew)
            BC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i,j,k+1)) - FC(i,k) - FC(i,k+1)
     &              +max(WC(i,k+1),0.D0)-min(WC(i,k),0.D0)
          enddo
        enddo
        do i=Istr,Iend
          DC(i,N)=wz(i,j,N,nnew)
          BC(i,N)=0.5D0*Hz(i,j,N) - FC(i,N)
     &              -min(WC(i,N),0.D0)
        enddo
        do i=Istr,Iend
          cff=1.D0/BC(i,1)
          CF(i,2)=cff*(FC(i,2)
     &              +min(WC(i,2),0.D0)
     &                                 )
          DC(i,1)=cff*DC(i,1)
        enddo
        do k=2,N-1
          do i=Istr,Iend
              aa = FC(i,k)-max(WC(i,k),0.D0)
              cc = FC(i,k+1)+min(WC(i,k+1),0.D0)
            cff=1.D0/(BC(i,k)-aa*CF(i,k))
            CF(i,k+1)=cff*cc
            DC(i,k)=cff*(DC(i,k)-aa*DC(i,k-1))
          enddo
        enddo
        do i=Istr,Iend
              aa = FC(i,N)-max(WC(i,N),0.D0)
          wz(i,j,N,nnew)=(DC(i,N)-aa*DC(i,N-1))
     &                  /(BC(i,N)-aa*CF(i,N))
        enddo
        do k=N-1,1,-1
          do i=Istr,Iend
             wz(i,j,k,nnew)=DC(i,k)-CF(i,k+1)*wz(i,j,k+1,nnew)
          enddo
        enddo
        call toc(1,'croco_step3d_uv2_1')
      enddo
      call u3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call w3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call tic(1,'croco_swap corrector')
         do k=1,N
            do j=Jstr,Jend
               do i=IstrU,Iend
                  u(i,j,k,nnew) = u(i,j,k,nnew)*umask(i,j)
               enddo
            enddo
         enddo
        do k=1,N
            do j=JstrV,Jend
               do i=Istr,Iend
                  v(i,j,k,nnew) = v(i,j,k,nnew)*vmask(i,j)
               enddo
            enddo
         enddo
         do k=1,N
            do j=Jstr,Jend
               do i=Istr,Iend
                  wz(i,j,k,nnew) = wz(i,j,k,nnew)*rmask(i,j)
               enddo
            enddo
         enddo
      do k=1,N
         do j=Jstr,Jend
            do i=Istr,Iend+1
               u(i,j,k,nnew) = u(i,j,k,nnew)
     &              *(0.5D0*(Hz(i,j,k) + Hz(i-1,j,k))) / pn_u(i,j)
            enddo
         enddo
         do j=Jstr,Jend+1
            do i=Istr,Iend
               v(i,j,k,nnew) = v(i,j,k,nnew)
     &              *(0.5D0*(Hz(i,j,k) + Hz(i,j-1,k))) / pm_v(i,j)
            enddo
         enddo
         do j=Jstr,Jend
            do i=Istr,Iend
               wz(i,j,k,nnew) = wz(i,j,k,nnew)
     &              / (pm(i,j)*pn(i,j))
            enddo
         enddo
      enddo
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                        u(-2,-2,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                        v(-2,-2,1,nnew))
       nh_ubar = 0.D0
       nh_vbar = 0.D0
       do k=1,N
         do j=Jstr,Jend
            do i=Istr,Iend+1
                nh_ubar(i,j) = nh_ubar(i,j) + u(i,j,k,nnew)
            enddo
         enddo
         do j=Jstr,Jend+1
           do i=Istr,Iend
              nh_vbar(i,j) = nh_vbar(i,j) + v(i,j,k,nnew)
           enddo
         enddo
       enddo
       do j=Jstr,Jend
          do i=Istr,Iend
             nh_wcor(i,j) = wz(i,j,N,nnew) +
     &         (nh_ubar(i+1,j)-nh_ubar(i,j)+nh_vbar(i,j+1)-nh_vbar(i,j))
          enddo
       enddo
       do k=1,N
          do j=Jstr,Jend
             do i=Istr,Iend
                wz(i,j,k,nnew) = wz(i,j,k,nnew) - nh_wcor(i,j)
     &                         * (z_w(i,j,k)-z_w(i,j,0))
     &                         / (z_w(i,j,N)-z_w(i,j,0))
             enddo
          enddo
       enddo
      halo = 3
         call nhmg_solve(Lm,Mm,N,halo,padd_X,padd_E,
     &          u(:,:,:,nnew),v(:,:,:,nnew),wz(:,:,:,nnew) )
      do k=1,N
         do j=Jstr,Jend
            do i=Istr,Iend+1
               nhdu(i,j,k,iprec1) = grid(1)%du(k,j,i)/dt
               u(i,j,k,nnew) = (u(i,j,k,nnew) + grid(1)%du(k,j,i) )
     &              /(0.5D0*(Hz(i,j,k) + Hz(i-1,j,k))) * pn_u(i,j)
            enddo
         enddo
         do j=Jstr,Jend+1
            do i=Istr,Iend
               nhdv(i,j,k,iprec1) = grid(1)%dv(k,j,i)/dt
               v(i,j,k,nnew) = (v(i,j,k,nnew) + grid(1)%dv(k,j,i) )
     &              /(0.5D0*(Hz(i,j,k) + Hz(i,j-1,k))) * pm_v(i,j)
            enddo
         enddo
         do j=Jstr,Jend
            do i=Istr,Iend
               nhdw(i,j,k,iprec1) = grid(1)%dw(k+1,j,i)/dt
               wz(i,j,k,nnew) = (wz(i,j,k,nnew) + grid(1)%dw(k+1,j,i) )
     &              * (pm(i,j)*pn(i,j))
            enddo
         enddo
      enddo
      call u3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call v3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call w3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call toc(1,'croco_swap corrector')
      do j=JstrR,JendR
        do i=Istr,IendR
          DC(i,0)=0.D0
          CF(i,0)=0.D0
          FC(i,0)=0.D0
        enddo
        do k=1,N,+1
          do i=Istr,IendR
            DC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)
            DC(i,0)=DC(i,0)+DC(i,k)
            CF(i,0)=CF(i,0)+DC(i,k)*u(i,j,k,nnew)
          enddo
        enddo
        do i=Istr,IendR
          DC(i,0)=1.D0/DC(i,0)
          CF(i,0)=DC(i,0)*(CF(i,0)-DU_avg1(i,j,nnew))
          ubar(i,j,knew)=DC(i,0)*DU_avg1(i,j,nnew)
        enddo
        do k=N,1,-1
          do i=Istr,IendR
            u(i,j,k,nnew)=(u(i,j,k,nnew)-CF(i,0))
     &                                *umask(i,j)
            FC(i,k)=0.75D0*Huon(i,j,k) +
     &              0.125D0*DC(i,k)*(u(i,j,k,nstp)+u(i,j,k,nnew))
            FC(i,0)=FC(i,0)+FC(i,k)
          enddo
        enddo
        do i=Istr,IendR
          FC(i,0)=DC(i,0)*(FC(i,0)-DU_avg2(i,j))
        enddo
        do k=1,N,+1
          do i=Istr,IendR
            Huon(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
          enddo
        enddo
        if (j.ge.Jstr) then
          do i=IstrR,IendR
            DC(i,0)=0.D0
            CF(i,0)=0.D0
            FC(i,0)=0.D0
          enddo
          do k=1,N,+1
            do i=IstrR,IendR
              DC(i,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)
              DC(i,0)=DC(i,0)+DC(i,k)
              CF(i,0)=CF(i,0)+DC(i,k)*v(i,j,k,nnew)
            enddo
          enddo
          do i=IstrR,IendR
            DC(i,0)=1.D0/DC(i,0)
            CF(i,0)=DC(i,0)*(CF(i,0)-DV_avg1(i,j,nnew))
            vbar(i,j,knew)=DC(i,0)*DV_avg1(i,j,nnew)
          enddo
          do k=N,1,-1
            do i=IstrR,IendR
              v(i,j,k,nnew)=(v(i,j,k,nnew)-CF(i,0))
     &                                  *vmask(i,j)
              FC(i,k)=0.75D0*Hvom(i,j,k)+
     &                0.125D0*DC(i,k)*(v(i,j,k,nstp)+v(i,j,k,nnew))
              FC(i,0)=FC(i,0)+FC(i,k)
            enddo
          enddo
          do i=IstrR,IendR
            FC(i,0)=DC(i,0)*(FC(i,0)-DV_avg2(i,j))
          enddo
          do k=1,N,+1
            do i=IstrR,IendR
              Hvom(i,j,k)=FC(i,k)-DC(i,k)*FC(i,0)
            enddo
          enddo
        endif
      enddo
      do is=1,Nsrc
           i=Isrc_mpi(is,mynode)
           j=Jsrc_mpi(is,mynode)
        if (IstrR.le.i .and. i.le.IendR .and.
     &      JstrR.le.j .and. j.le.JendR) then
          if (Dsrc(is).eq.0) then
            do k=1,N
              u(i,j,k,nnew)=2.D0*Qsrc(is,k)/( on_u(i,j)*(
     &                       z_w(i-1,j,k)-z_w(i-1,j,k-1)
     &                      +z_w(i  ,j,k)-z_w(i  ,j,k-1)
     &                                                ))
              Huon(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i-1,j,k))*on_u(i,j)*
     &                                            u(i,j,k,nnew)
            enddo
          elseif (Dsrc(is).eq.1) then
            do k=1,N
              v(i,j,k,nnew)=2.D0*Qsrc(is,k)/( om_v(i,j)*(
     &                       z_w(i,j-1,k)-z_w(i,j-1,k-1)
     &                      +z_w(i,j  ,k)-z_w(i,j  ,k-1)
     &                                                ))
              Hvom(i,j,k)=0.5D0*(Hz(i,j,k)+Hz(i,j-1,k))*om_v(i,j)*
     &                                            v(i,j,k,nnew)
            enddo
          endif
        endif
      enddo
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                                 u(-2,-2,1,nnew))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                                 v(-2,-2,1,nnew))
      call exchange_u3d_tile (Istr,Iend,Jstr,Jend,
     &                                   Huon(-2,-2,1))
      call exchange_v3d_tile (Istr,Iend,Jstr,Jend,
     &                                   Hvom(-2,-2,1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,
     &                                   ubar(-2,-2,knew))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,
     &                                   vbar(-2,-2,knew))
      do j=Jstr,Jend
        do k=1,N
          do i=Istr,Iend
            DC(i,0)=pm(i,j)*pn(i,j)
            wz(i,j,k,nnew)=wz(i,j,k-1,nnew)
     &      - DC(i,0)*(
     &       +0.5D0*(Hz(i,j,k)+Hz(i+1,j,k))*on_u(i+1,j)*u(i+1,j,k,nnew)
     &       -0.5D0*(Hz(i-1,j,k)+Hz(i,j,k))*on_u(i,j)*u(i,j,k,nnew)
     &       +0.5D0*(Hz(i,j,k)+Hz(i,j+1,k))*om_v(i,j+1)*v(i,j+1,k,nnew)
     &       -0.5D0*(Hz(i,j-1,k)+Hz(i,j,k))*om_v(i,j)*v(i,j,k,nnew) )
          enddo
        enddo
      enddo
      call w3dbc_tile (Istr,Iend,Jstr,Jend, grad)
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend,
     &                                 wz(-2,-2,0,nnew))
      return
      end
