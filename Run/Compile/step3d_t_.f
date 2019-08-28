      subroutine step3d_t (tile)
      implicit none
      integer*4 tile, trd, omp_get_thread_num
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=256,  MMm0=128,  N=64)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=16,  NP_ETA=8,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Msrc
      parameter (Msrc=270)
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
      call step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                    A2d(1,1,trd), A2d(1,2,trd), A2d(1,3,trd),
     &                    A2d(1,4,trd), A2d(1,5,trd), A2d(1,6,trd),
     &                    A2d(1,7,trd), A2d(1,8,trd), A2d(1,9,trd),
     &                                                A3d(1,1,trd))
      return
      end
      subroutine step3d_t_tile (Istr,Iend,Jstr,Jend,
     &                          FX,FE, WORK, FC,CF,BC,DC,EC,GC, swdk)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=256,  MMm0=128,  N=64)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=16,  NP_ETA=8,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Msrc
      parameter (Msrc=270)
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
      integer*4 Istr,Iend,Jstr,Jend, itrc, i,j,k, indx, kmld
     &       ,imin,imax,jmin,jmax,iAkt,nadv
     &       ,is,iii,jjj
      real FX(Istr-2:Iend+2,Jstr-2:Jend+2),
     &     FE(Istr-2:Iend+2,Jstr-2:Jend+2),   cff,
     &     WORK(Istr-2:Iend+2,Jstr-2:Jend+2), epsil,
     &     FC(Istr-2:Iend+2,0:N),
     &     CF(Istr-2:Iend+2,0:N),
     &     BC(Istr-2:Iend+2,0:N),
     &     DC(Istr-2:Iend+2,0:N),
     &     EC(Istr-2:Iend+2,0:N),
     &     GC(Istr-2:Iend+2,0:N),
     &     swdk(Istr-2:Iend+2,Jstr-2:Jend+2,0:N)
      real cff1,cff2,gama,dRz,hbltmp,sig,dXmax,dEmax,
     &     dpth,smax,amax,amaxx
      parameter (epsil=1.D-16)
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2
      REAL    :: ua, vel, cdiff, cdif
      REAL    :: flux1, flux2, flux3, flux4, flux5, flux6
      REAL    :: flx2, flx3, flx4, flx5
      REAL    :: mask0, mask1, mask2, mask3
      flux2(q_im1, q_i, ua, cdiff) = 0.5D0*( q_i + q_im1 )
      flux1(q_im1, q_i, ua, cdiff) = flux2(q_im1, q_i, ua, cdiff) -
     &      0.5D0*cdiff*sign(1.D0,ua)*(q_i-q_im1)
      flux4(q_im2, q_im1, q_i, q_ip1, ua) =
     &      ( 7.D0*(q_i + q_im1) - (q_ip1 + q_im2) )/12.0D0
      flux3(q_im2, q_im1, q_i, q_ip1, ua) =
     &      flux4(q_im2, q_im1, q_i, q_ip1, ua) +
     &      sign(1.D0,ua)*((q_ip1 -
     &      q_im2)-3.D0*(q_i-q_im1))/12.0D0
      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      ( 37.D0*(q_i+q_im1) - 8.D0*(q_ip1+q_im2)
     &      +(q_ip2+q_im3) )/60.0D0
      flux5(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua) =
     &      flux6(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
     &      -sign(1.D0,ua)*(
     &      (q_ip2-q_im3)-5.D0*(q_ip1-q_im2)+10.D0*(q_i-q_im1) )/60.0D0
      REAL    :: flux3_weno, flux5_weno
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
      nadv = 3
      do itrc=1,NT
        do k=1,N
          if (itrc .le. 2) then
          if (WEST_INTER) then
            imin=Istr-1
          else
            imin=max(Istr-1,1)
          endif
          if (EAST_INTER) then
            imax=Iend+2
          else
            imax=min(Iend+2,Lmmpi+1)
          endif
          if (SOUTH_INTER) then
            jmin=Jstr-1
          else
            jmin=max(Jstr-1,1)
          endif
          if (NORTH_INTER) then
            jmax=Jend+2
          else
            jmax=min(Jend+2,Mmmpi+1)
          endif
          do j=Jstr,Jend
            do i=imin,imax
              FX(i,j)=(t(i,j,k,nadv,itrc)-t(i-1,j,k,nadv,itrc))
     &                                               *umask(i,j)
            enddo
          enddo
          if (.not.WEST_INTER) then
            do j=Jstr,Jend
              FX(0,j)=FX(1,j)
            enddo
          endif
          if (.not.EAST_INTER) then
            do j=Jstr,Jend
              FX(Lmmpi+2,j)=FX(Lmmpi+1,j)
            enddo
          endif
          do j=Jstr,Jend
            do i=Istr-1,Iend+1
              WORK(i,j)=FX(i+1,j)-FX(i,j)
            enddo
          enddo
          do j=Jstr,Jend
            do i=Istr,Iend+1
              if (Huon(i,j,k) .gt. 0.D0) then
                cff=WORK(i-1,j)
              else
                cff=WORK(i,j)
              endif
              FX(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i-1,j,k,nadv,itrc)
     &                           -0.333333333333D0*cff )*Huon(i,j,k)
            enddo
          enddo
          do j=jmin,jmax
            do i=Istr,Iend
              FE(i,j)=(t(i,j,k,nadv,itrc)-t(i,j-1,k,nadv,itrc))
     &                                               *vmask(i,j)
            enddo
          enddo
          if (.not.SOUTH_INTER) then
            do i=Istr,Iend
              FE(i,0)=FE(i,1)
            enddo
          endif
          if (.not.NORTH_INTER) then
            do i=Istr,Iend
              FE(i,Mmmpi+2)=FE(i,Mmmpi+1)
            enddo
          endif
          do j=Jstr-1,Jend+1
            do i=Istr,Iend
              WORK(i,j)=FE(i,j+1)-FE(i,j)
            enddo
          enddo
          do j=Jstr,Jend+1
            do i=Istr,Iend
              if (Hvom(i,j,k) .gt. 0.D0) then
                cff=WORK(i,j-1)
              else
                cff=WORK(i,j)
              endif
              FE(i,j)=0.5D0*( t(i,j,k,nadv,itrc)+t(i,j-1,k,nadv,itrc)
     &                          -0.333333333333D0*cff )*Hvom(i,j,k)
            enddo
          enddo
          else
            cdif=1.D0
          if (SOUTH_INTER) then
            jmin=1
          else
            jmin=3
          endif
          if (NORTH_INTER) then
            jmax=Mmmpi+1
          else
            jmax=Mmmpi-1
          endif
          if (WEST_INTER) then
            imin=1
          else
            imin=3
          endif
          if (EAST_INTER) then
            imax=Lmmpi+1
          else
            imax=Lmmpi-1
          endif
          DO j = Jstr,Jend+1
            IF ( j.ge.jmin .and. j.le.jmax ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx5 = vel*flux5_weno(
     &             t(i,j-3,k,nrhs,itrc), t(i,j-2,k,nrhs,itrc),
     &             t(i,j-1,k,nrhs,itrc), t(i,j  ,k,nrhs,itrc),
     &             t(i,j+1,k,nrhs,itrc), t(i,j+2,k,nrhs,itrc),  vel )
                flx3 = vel*flux3_weno(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask0=rmask(i,j-1)*rmask(i,j)
                mask2=rmask(i,j-2)*mask0*rmask(i,j+1)
                IF (vel.gt.0) THEN
                  mask1=rmask(i,j-2)*mask0
                  mask3=rmask(i,j-3)*mask2
                ELSE
                  mask1=rmask(i,j+1)*mask0
                  mask3=rmask(i,j+2)*mask2
                ENDIF
                FE(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
              ENDDO
            ELSE IF ( j.eq.jmin-2 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                FE(i,j) = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( j.eq.jmin-1 .and. jmax.ge.jmin ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx3 = vel*flux3_weno(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                FE(i,j)=mask1*flx3+(1-mask1)*flx2
              ENDDO
            ELSE IF ( j.eq.jmax+2 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                FE(i,j) = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( j.eq.jmax+1 ) THEN
              DO i = Istr,Iend
                vel = Hvom(i,j,k)
                flx3 = vel*flux3_weno(
     &             t(i,j-2,k,nrhs,itrc), t(i,j-1,k,nrhs,itrc),
     &             t(i,j  ,k,nrhs,itrc), t(i,j+1,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i,j-1,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask1=rmask(i,j-2)*rmask(i,j+1)
                FE(i,j)=mask1*flx3+(1-mask1)*flx2
              ENDDO
            ENDIF
          ENDDO
          DO i = Istr,Iend+1
            IF ( i.ge.imin .and. i.le.imax ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx5 = vel*flux5_weno(
     &             t(i-3,j,k,nrhs,itrc), t(i-2,j,k,nrhs,itrc),
     &             t(i-1,j,k,nrhs,itrc), t(i  ,j,k,nrhs,itrc),
     &             t(i+1,j,k,nrhs,itrc), t(i+2,j,k,nrhs,itrc),  vel )
                flx3 = vel*flux3_weno(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask0=rmask(i-1,j)*rmask(i,j)
                mask2=rmask(i-2,j)*mask0*rmask(i+1,j)
                IF (vel.gt.0) THEN
                  mask1=rmask(i-2,j)*mask0
                  mask3=rmask(i-3,j)*mask2
                ELSE
                  mask1=rmask(i+1,j)*mask0
                  mask3=rmask(i+2,j)*mask2
                ENDIF
                FX(i,j)=mask3*flx5+(1-mask3)*mask1*flx3+
     &                             (1-mask3)*(1-mask1)*mask0*flx2
              ENDDO
            ELSE IF ( i.eq.imin-2 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                FX(i,j) = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( i.eq.imin-1 .and. imax.ge.imin ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx3 = vel*flux3_weno(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                FX(i,j)=mask1*flx3+(1-mask1)*flx2
              ENDDO
            ELSE IF ( i.eq.imax+2 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                FX(i,j) = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
              ENDDO
            ELSE IF ( i.eq.imax+1 ) THEN
              DO j = Jstr,Jend
                vel = Huon(i,j,k)
                flx3 = vel*flux3_weno(
     &             t(i-2,j,k,nrhs,itrc), t(i-1,j,k,nrhs,itrc),
     &             t(i  ,j,k,nrhs,itrc), t(i+1,j,k,nrhs,itrc),  vel )
                flx2 = vel*flux1(
     &             t(i-1,j,k,nrhs,itrc), t(i,j,k,nrhs,itrc), vel, cdif)
                mask1=rmask(i-2,j)*rmask(i+1,j)
                FX(i,j)=mask1*flx3+(1-mask1)*flx2
              ENDDO
            ENDIF
          ENDDO
          endif
         do is=1,Nsrc
           i=Isrc_mpi(is,mynode)
           j=Jsrc_mpi(is,mynode)
           if (Istr.le.i .and. i.le.Iend+1
     &                   .and. Jstr.le.j .and. j.le.Jend+1) then
             if (Dsrc(is).eq.0) then
               if (Lsrc(is,itrc)) then
                  FX(i,j)=Huon(i,j,k)*Tsrc(is,k,itrc)
               else
                 if (rmask(i,j).eq.0 .and. rmask(i-1,j).eq.1) then
                    FX(i,j)=Huon(i,j,k)*t(i-1,j,k,3,itrc)
                 elseif (rmask(i,j).eq.1.D0 .and. rmask(i-1,j).eq.0) 
     &                                                              then
                    FX(i,j)=Huon(i,j,k)*t(i  ,j,k,3,itrc)
                 endif
               endif
             elseif (Dsrc(is).eq.1) then
               if (Lsrc(is,itrc)) then
                 FE(i,j)=Hvom(i,j,k)*Tsrc(is,k,itrc)
               else
                 if (rmask(i,j).eq.0 .and. rmask(i,j-1).eq.1) then
                   FE(i,j)=Hvom(i,j,k)*t(i,j-1,k,3,itrc)
                 elseif (rmask(i,j).eq.1 .and. rmask(i,j-1).eq.0) then
                   FE(i,j)=Hvom(i,j,k)*t(i,j  ,k,3,itrc)
                 endif
               endif
             endif
           endif
         enddo
          do j=Jstr,Jend
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=Hz_bak(i,j,k)*t(i,j,k,nstp,itrc)
     &                     -dt*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                          +FE(i,j+1)-FE(i,j)
     &                                                           )
            enddo
          enddo
        enddo
      enddo
      do j=Jstr,Jend
        do itrc=1,NT
          do i=Istr,Iend
            FC(i,0)=0.D0
            CF(i,0)=0.D0
          enddo
          do k=1,N-1,+1
            do i=Istr,Iend
              cff    = 1.D0
     &                    /(2.D0*Hz(i,j,k+1)+Hz(i,j,k)*(2.D0-FC(i,k-1)))
              FC(i,k)= cff*Hz(i,j,k+1)
              CF(i,k)= cff*( 6.D0*( t(i,j,k+1,nadv,itrc)
     &                           -t(i,j,k  ,nadv,itrc) 
     &                                             )-Hz(i,j,k)*CF(i,k-1)
     &                                                                  
     &                                                                 )
            enddo
          enddo
          do i=Istr,Iend
            CF(i,N)=0.D0
          enddo
          do k=N-1,1,-1
            do i=Istr,Iend
              CF(i,k)=CF(i,k)-FC(i,k)*CF(i,k+1)
            enddo
          enddo
          cff=1.D0/3.D0
          do k=1,N-1
            do i=Istr,Iend
              FC(i,k)=We(i,j,k)*( t(i,j,k,nadv,itrc)+cff*Hz(i,j,k)
     &                                  *(CF(i,k)+0.5D0*CF(i,k-1))
     &                                                            )
            enddo
          enddo
          do i=Istr,Iend
            FC(i,N)=0.D0
            FC(i,0)=0.D0
            CF(i,0)=dt*pm(i,j)*pn(i,j)
          enddo
          do k=1,N
            do i=Istr,Iend
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)-CF(i,0)*(FC(i,k)
     &                                                      -FC(i,k-1))
            enddo
          enddo
          do is=1,Nsrc
            iii=Isrc_mpi(is,mynode)
            jjj=Jsrc_mpi(is,mynode)
            if (Istr.le.iii .and. iii.le.Iend+1
     &         .and. Jstr.le.jjj .and. jjj.le.Jend+1.and.jjj.eq.j) then
              if (Dsrc(is).eq.2) then
                if (Lsrc(is,itrc)) then
                  cff =dt*pn(iii,jjj)*pm(iii,jjj)
                  do k=1,N
                    t(iii,jjj,k,nnew,itrc)=(t(iii,jjj,k,nnew,itrc)
     &                                 +cff*Qsrc(is,k)*Tsrc(is,k,itrc))
                  enddo
                endif
              endif
            endif
          enddo
          do i=Istr,Iend
            FC(i,N)=dt*stflx(i,j,itrc)
            FC(i,0)=-dt*btflx(i,j,itrc)
          enddo
          if (itrc.eq.itemp) then
            do k=1,N-1
              do i=Istr,Iend
                FC(i,k)=0.D0
              enddo
            enddo
          elseif (itrc.eq.isalt) then
            do k=1,N-1
              do i=Istr,Iend
                FC(i,k)=0.D0
              enddo
            enddo
          endif
          if (itrc.eq.itemp .or. itrc.eq.isalt) then
            do k=1,N
              do i=Istr,Iend
                t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)+FC(i,k )
     &                                               -FC(i,k-1)
              enddo
            enddo
          endif
       do i=Istr,Iend
            DC(i,0)=dt*pn(i,j)*pm(i,j)
       enddo
          indx=min(itrc,isalt)
          do i=istr,iend
            FC(i,1)=dt* Akt(i,j,1,indx)
     &                               /( z_r(i,j,2)-z_r(i,j,1) )
            BC(i,1)=DC(i,0)*Wi(i,j,1)
            cff=1.D0/(Hz(i,j,1)      +FC(i,1)+max(BC(i,1),0.D0))
            CF(i,1)=cff*(           FC(i,1)-min(BC(i,1),0.D0))
            DC(i,1)= cff* t(i,j,1,nnew,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=dt* Akt(i,j,k,indx)
     &                              /( z_r(i,j,k+1)-z_r(i,j,k) )
              BC(i,k)=DC(i,0)*Wi(i,j,k)
              cff=1.D0/(      Hz(i,j,k) +FC(i,k)+max(BC(i,k),0.D0)
     &                              +FC(i,k-1)-min(BC(i,k-1),0.D0)
     &                   -CF(i,k-1)*(FC(i,k-1)+max(BC(i,k-1),0.D0))
     &                                                          )
              CF(i,k)=cff*(FC(i,k)-min(BC(i,k),0.D0))
              DC(i,k)=cff*( t(i,j,k,nnew,itrc) +DC(i,k-1)*(
     &                          FC(i,k-1)+max(BC(i,k-1),0.D0) )
     &                                                          )
            enddo
          enddo
          do i=istr,iend
            t(i,j,N,nnew,itrc)=( t(i,j,N,nnew,itrc)
     &                                           +DC(i,N-1)*(
     &                                FC(i,N-1)+max(BC(i,N-1),0.D0) )
     &               )/( Hz(i,j,N) +FC(i,N-1)-min(BC(i,N-1),0.D0)
     &                      -CF(i,N-1)*(FC(i,N-1)+max(BC(i,N-1),0.D0))
     &                                                            )
        enddo
          do k=N-1,1,-1
            do i=istr,iend
              t(i,j,k,nnew,itrc)=DC(i,k)+CF(i,k)*t(i,j,k+1,nnew,itrc)
            enddo
          enddo
        enddo
      enddo
      do itrc=1,NT
        call t3dbc_tile (Istr,Iend,Jstr,Jend, nnew,itrc, WORK)
      enddo
      do itrc=1,NT
        do k=1,N
          do j=JstrR,JendR
            do i=IstrR,IendR
              t(i,j,k,nnew,itrc)=t(i,j,k,nnew,itrc)*rmask(i,j)
            enddo
          enddo
        enddo
        call exchange_r3d_3pts_tile (Istr,Iend,Jstr,Jend,
     &                          t(-2,-2,1,nnew,itrc))
      enddo
      return
      end
      function flux5_weno(q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua)
      implicit none
      REAL    :: flux5_weno
      REAL    :: q_im3, q_im2, q_im1, q_i, q_ip1, q_ip2, ua
      REAL    :: IS0, IS1, IS2
      REAl    :: d0, d1, d2
      REAl    :: a0, a1, a2
      REAL    :: w0, w1, w2
      REAL    :: p0, p1, p2
      REAL    :: Eps, cff1, cff2, T5
      Eps = 1.D-40
      d0=1.D0/10.D0
      d1=6.D0/10.D0
      d2=3.D0/10.D0
      cff1=13.D0/12.D0
      cff2=1.D0/6.D0
      if (ua .ge. 0.D0) then
        IS0 = cff1*(q_im3 - 2.D0*q_im2 + q_im1)**2
     &                        + 0.25D0*(q_im3 - 4.D0*q_im2 + 3*q_im1)**2
        IS1 = cff1*(q_im2 - 2.D0*q_im1 + q_i)**2
     &                        + 0.25D0*(q_im2 - q_i)**2
        IS2 = cff1*(q_im1 - 2.D0*q_i + q_ip1)**2
     &                        + 0.25D0*(3.D0*q_im1 - 4.D0*q_i + 
     &                                                         q_ip1)**2
        T5 = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_im3 - 7.D0*q_im2 + 11.D0*q_im1)
        p1  = cff2*(-q_im2   + 5.D0*q_im1 + 2.D0*q_i)
        p2  = cff2*(2.D0*q_im1 + 5.D0*q_i   - q_ip1)
        flux5_weno = w0*p0 + w1*p1 +w2*p2
      else
        IS0 = cff1*(q_ip2 - 2.D0*q_ip1 + q_i)**2
     &                         + 0.25D0*(q_ip2 -4.D0*q_ip1 + 3*q_i)**2
        IS1 = cff1*(q_ip1 - 2.D0*q_i + q_im1)**2
     &                         + 0.25D0*(q_ip1-q_im1)**2
        IS2 = cff1*(q_i   - 2.D0*q_im1 + q_im2)**2
     &                         + 0.25D0*(3.D0*q_i -4.D0*q_im1 + 
     &                                                         q_im2)**2
        T5 = abs(IS2-IS0)
        a0  = d0*(1+T5/(Eps+IS0))
        a1  = d1*(1+T5/(Eps+IS1))
        a2  = d2*(1+T5/(Eps+IS2))
        w0  = a0/(a0+a1+a2)
        w1  = a1/(a0+a1+a2)
        w2  = a2/(a0+a1+a2)
        p0  = cff2*(2.D0*q_ip2 - 7.D0*q_ip1 + 11.D0*q_i)
        p1  = cff2*(-q_ip1   + 5.D0*q_i   + 2.D0*q_im1)
        p2  = cff2*(2.D0*q_i   + 5.D0*q_im1 - q_im2)
        flux5_weno = w0*p0 + w1*p1 +w2*p2
      endif
      return
      end
      function flux3_weno( q_im2, q_im1, q_i, q_ip1, ua)
      implicit none
      REAL    :: flux3_weno
      REAL    :: q_im2, q_im1, q_i, q_ip1, ua
      REAL    :: IS0, IS1
      REAl    :: a0, a1
      REAL    :: w0, w1
      REAL    :: p0, p1
      REAL    :: Eps, d0,d1, T3
      Eps = 1.D-40
      d0=1.D0/3.D0
      d1=2.D0/3.D0
      if (ua .ge. 0.D0) then
        IS0 = (q_im1-q_im2)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_im1-q_im2)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      else
        IS0 = (q_i-q_ip1)**2
        IS1 = (q_im1-q_i)**2
        T3 = abs(IS1-IS0)
        a0  = d0*(1+T3/(Eps+IS0))
        a1  = d1*(1+T3/(Eps+IS1))
        w0  = a0/(a0+a1)
        w1  = a1/(a0+a1)
        p0  = 1.D0/2.D0*(3.D0*q_i-q_ip1)
        p1  = 1.D0/2.D0*(q_im1+q_i)
        flux3_weno = w0*p0 + w1*p1
      endif
      return
      end
