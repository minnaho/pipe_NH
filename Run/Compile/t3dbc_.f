      subroutine t3dbc_tile(Istr,Iend,Jstr,Jend,indx,itrc,grad)
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
      real bry_time(2)
      common /bry_indices_array/ bry_time
      real bry_cycle
      common /bry_indices_real/ bry_cycle
      integer*4 bry_id, bry_time_id, bry_ncycle, bry_rec, itbry, ntbry
      common /bry_indices_integer/ bry_id, bry_time_id, bry_ncycle,
     &                             bry_rec, itbry, ntbry
      integer*4 zetabry_west_id
      common /zeta_west_id/ zetabry_west_id
      integer*4 ubarbry_west_id, vbarbry_west_id
      common /ubar_west_id/ ubarbry_west_id, vbarbry_west_id
      integer*4 ubry_west_id, vbry_west_id
      common /u_west_id/ ubry_west_id, vbry_west_id
      integer*4 tbry_west_id(NT)
      common /t_west_id/ tbry_west_id
      integer*4 zetabry_east_id
      common /zeta_east_id/ zetabry_east_id
      integer*4 ubarbry_east_id, vbarbry_east_id
      common /ubar_east_id/ ubarbry_east_id, vbarbry_east_id
      integer*4 ubry_east_id, vbry_east_id
      common /u_east_id/ ubry_east_id, vbry_east_id
      integer*4 tbry_east_id(NT)
      common /t_east_id/ tbry_east_id
      integer*4 zetabry_south_id
      common /zeta_south_id/ zetabry_south_id
      integer*4 ubarbry_south_id, vbarbry_south_id
      common /ubar_south_id/ ubarbry_south_id, vbarbry_south_id
      integer*4 ubry_south_id, vbry_south_id
      common /u_south_id/ ubry_south_id, vbry_south_id
      integer*4 tbry_south_id(NT)
      common /t_south_id/ tbry_south_id
      integer*4 zetabry_north_id
      common /zeta_north_id/ zetabry_north_id
      integer*4 ubarbry_north_id, vbarbry_north_id
      common /ubar_north_id/ ubarbry_north_id, vbarbry_north_id
      integer*4 ubry_north_id, vbry_north_id
      common /u_north_id/ ubry_north_id, vbry_north_id
      integer*4 tbry_north_id(NT)
      common /t_north_id/ tbry_north_id
      real zetabry_west(-1:Mm+2+padd_E),
     &    zetabry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(-1:Mm+2+padd_E),
     &    ubarbry_west_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_west(-1:Mm+2+padd_E),
     &    vbarbry_west_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real ubry_west(-1:Mm+2+padd_E,N),
     &    ubry_west_dt(-1:Mm+2+padd_E,N,2)
     &    ,vbry_west(-1:Mm+2+padd_E,N),
     &    vbry_west_dt(-1:Mm+2+padd_E,N,2)
      common /bry_u_west/ ubry_west, ubry_west_dt,
     &                    vbry_west, vbry_west_dt
      real tbry_west(-1:Mm+2+padd_E,N,NT),
     &    tbry_west_dt(-1:Mm+2+padd_E,N,2,NT)
      common /bry_t_west/ tbry_west, tbry_west_dt
      real zetabry_east(-1:Mm+2+padd_E),
     &    zetabry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real ubarbry_east(-1:Mm+2+padd_E),
     &    ubarbry_east_dt(-1:Mm+2+padd_E,2)
     &    ,vbarbry_east(-1:Mm+2+padd_E),
     &    vbarbry_east_dt(-1:Mm+2+padd_E,2)
      common /bry_ubar_east/ ubarbry_east, ubarbry_east_dt,
     &                       vbarbry_east, vbarbry_east_dt
      real ubry_east(-1:Mm+2+padd_E,N),
     &    ubry_east_dt(-1:Mm+2+padd_E,N,2)
     &    ,vbry_east(-1:Mm+2+padd_E,N),
     &    vbry_east_dt(-1:Mm+2+padd_E,N,2)
      common /bry_u_east/ ubry_east, ubry_east_dt,
     &                    vbry_east, vbry_east_dt
      real tbry_east(-1:Mm+2+padd_E,N,NT),
     &    tbry_east_dt(-1:Mm+2+padd_E,N,2,NT)
      common /bry_t_east/ tbry_east, tbry_east_dt
      real zetabry_south(-1:Lm+2+padd_X),
     &    zetabry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real ubarbry_south(-1:Lm+2+padd_X),
     &    ubarbry_south_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_south(-1:Lm+2+padd_X),
     &    vbarbry_south_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_south/ ubarbry_south, ubarbry_south_dt,
     &                        vbarbry_south, vbarbry_south_dt
      real ubry_south(-1:Lm+2+padd_X,N),
     &    ubry_south_dt(-1:Lm+2+padd_X,N,2)
     &    ,vbry_south(-1:Lm+2+padd_X,N),
     &    vbry_south_dt(-1:Lm+2+padd_X,N,2)
      common /bry_u_south/ ubry_south, ubry_south_dt,
     &                     vbry_south, vbry_south_dt
      real tbry_south(-1:Lm+2+padd_X,N,NT),
     &    tbry_south_dt(-1:Lm+2+padd_X,N,2,NT)
      common /bry_t_south/ tbry_south, tbry_south_dt
      real zetabry_north(-1:Lm+2+padd_X),
     &    zetabry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
      real ubarbry_north(-1:Lm+2+padd_X),
     &    ubarbry_north_dt(-1:Lm+2+padd_X,2)
     &    ,vbarbry_north(-1:Lm+2+padd_X),
     &    vbarbry_north_dt(-1:Lm+2+padd_X,2)
      common /bry_ubar_north/ ubarbry_north, ubarbry_north_dt,
     &                        vbarbry_north, vbarbry_north_dt
      real ubry_north(-1:Lm+2+padd_X,N),
     &    ubry_north_dt(-1:Lm+2+padd_X,N,2)
     &    ,vbry_north(-1:Lm+2+padd_X,N),
     &    vbry_north_dt(-1:Lm+2+padd_X,N,2)
      common /bry_u_north/ ubry_north, ubry_north_dt,
     &                     vbry_north, vbry_north_dt
      real tbry_north(-1:Lm+2+padd_X,N,NT),
     &    tbry_north_dt(-1:Lm+2+padd_X,N,2,NT)
      common /bry_t_north/ tbry_north, tbry_north_dt
      integer*4 Istr,Iend,Jstr,Jend,indx,itrc,i,j,k
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    eps, cff,
     &        cx,cy, dft,dfx,dfy, tau,tau_in,tau_out
      parameter (eps=1.D-20)
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
      tau_in=dt*tauT_in
      tau_out=dt*tauT_out
      if (.not.WEST_INTER) then
        do k=1,N
          do j=Jstr,Jend
              t(Istr-1,j,k,indx,itrc)=tbry_west(j,k,itrc)
     &                                       *rmask(Istr-1,j)
         enddo
        enddo
      endif
      if (.not.EAST_INTER) then
        do k=1,N
          do j=Jstr,Jend
            t(Iend+1,j,k,indx,itrc)=tbry_east(j,k,itrc)
     &                                   *rmask(Iend+1,j)
          enddo
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do k=1,N
          do i=Istr,Iend
            t(i,Jstr-1,k,indx,itrc)=tbry_south(i,k,itrc)
     &                                   *rmask(i,Jstr-1)
          enddo
        enddo
      endif
      if (.not.NORTH_INTER) then
        do k=1,N
          do i=Istr,Iend
            t(i,Jend+1,k,indx,itrc)=tbry_north(i,k,itrc)
     &                                   *rmask(i,Jend+1)
          enddo
        enddo
      endif
      if (.not.SOUTH_INTER .and. .not.WEST_INTER) then
        do k=1,N
          t(Istr-1,Jstr-1,k,indx,itrc)=0.5D0*
     &                          ( t(Istr,Jstr-1,k,indx,itrc)
     &                           +t(Istr-1,Jstr,k,indx,itrc))
     &                                 *rmask(Istr-1,Jstr-1)
        enddo
      endif
      if (.not.SOUTH_INTER .and. .not.EAST_INTER) then
        do k=1,N
          t(Iend+1,Jstr-1,k,indx,itrc)=0.5D0*
     &                          (t(Iend,Jstr-1,k,indx,itrc)
     &                          +t(Iend+1,Jstr,k,indx,itrc))
     &                                 *rmask(Iend+1,Jstr-1)
        enddo
      endif
      if (.not.NORTH_INTER .and. .not.WEST_INTER) then
        do k=1,N
          t(Istr-1,Jend+1,k,indx,itrc)=0.5D0*
     &                          ( t(Istr,Jend+1,k,indx,itrc)
     &                           +t(Istr-1,Jend,k,indx,itrc))
     &                                 *rmask(Istr-1,Jend+1)
        enddo
      endif
      if (.not.NORTH_INTER .and. .not.EAST_INTER) then
        do k=1,N
          t(Iend+1,Jend+1,k,indx,itrc)=0.5D0*
     &                          ( t(Iend,Jend+1,k,indx,itrc)
     &                           +t(Iend+1,Jend,k,indx,itrc))
     &                                 *rmask(Iend+1,Jend+1)
        enddo
      endif
      return
      end
