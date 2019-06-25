      subroutine v2dbc_tile(Istr,Iend,Jstr,Jend,grad)
      implicit none
      integer*4  LLm,Lm,MMm,Mm,N, LLm0,MMm0
      parameter (LLm0=512,  MMm0=512,  N=64)
      parameter (LLm=LLm0,  MMm=MMm0)
      integer*4 Lmmpi,Mmmpi,iminmpi,imaxmpi,jminmpi,jmaxmpi
      common /comm_setup_mpi1/ Lmmpi,Mmmpi
      common /comm_setup_mpi2/ iminmpi,imaxmpi,jminmpi,jmaxmpi
      integer*4 NSUB_X, NSUB_E, NPP
      integer*4 NP_XI, NP_ETA, NNODES
      parameter (NP_XI=8,  NP_ETA=4,  NNODES=NP_XI*NP_ETA)
      parameter (NPP=1)
      parameter (NSUB_X=1, NSUB_E=1)
      integer*4 NWEIGHT
      parameter (NWEIGHT=1000)
      integer*4 Msrc
      parameter (Msrc=6000)
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
      integer*4 Istr,Iend,Jstr,Jend, i,j
      real    grad(Istr-2:Iend+2,Jstr-2:Jend+2)
      real    eps,cff, cx,cy,
     &        dft,dfx,dfy, tau,tau_in,tau_out,hx,zx
      parameter (eps=1.D-20)
      real zetabry_west(-2:Mm+3+padd_E),
     &    zetabry_west_dt(-2:Mm+3+padd_E,2)
      common /bry_zeta_west/ zetabry_west, zetabry_west_dt
      real ubarbry_west(-2:Mm+3+padd_E),
     &    ubarbry_west_dt(-2:Mm+3+padd_E,2)
     &    ,vbarbry_west(-2:Mm+3+padd_E),
     &    vbarbry_west_dt(-2:Mm+3+padd_E,2)
      common /bry_ubar_west/ ubarbry_west, ubarbry_west_dt,
     &                       vbarbry_west, vbarbry_west_dt
      real ubry_west(-2:Mm+3+padd_E,N),
     &    ubry_west_dt(-2:Mm+3+padd_E,N,2)
     &    ,vbry_west(-2:Mm+3+padd_E,N),
     &    vbry_west_dt(-2:Mm+3+padd_E,N,2)
      common /bry_u_west/ ubry_west, ubry_west_dt,
     &                    vbry_west, vbry_west_dt
      real wbry_west(-2:Mm+3+padd_E,N),
     &    wbry_west_dt(-2:Mm+3+padd_E,N,2)
      common /bry_w_west/ wbry_west, wbry_west_dt
      real zetabry_east(-2:Mm+3+padd_E),
     &    zetabry_east_dt(-2:Mm+3+padd_E,2)
      common /bry_zeta_east/ zetabry_east, zetabry_east_dt
      real ubarbry_east(-2:Mm+3+padd_E),
     &    ubarbry_east_dt(-2:Mm+3+padd_E,2)
     &    ,vbarbry_east(-2:Mm+3+padd_E),
     &    vbarbry_east_dt(-2:Mm+3+padd_E,2)
      common /bry_ubar_east/ ubarbry_east, ubarbry_east_dt,
     &                       vbarbry_east, vbarbry_east_dt
      real ubry_east(-2:Mm+3+padd_E,N),
     &    ubry_east_dt(-2:Mm+3+padd_E,N,2)
     &    ,vbry_east(-2:Mm+3+padd_E,N),
     &    vbry_east_dt(-2:Mm+3+padd_E,N,2)
      common /bry_u_east/ ubry_east, ubry_east_dt,
     &                    vbry_east, vbry_east_dt
      real wbry_east(-2:Mm+3+padd_E,N),
     &    wbry_east_dt(-2:Mm+3+padd_E,N,2)
      common /bry_w_east/ wbry_east, wbry_east_dt
      real zetabry_south(-2:Lm+3+padd_X),
     &    zetabry_south_dt(-2:Lm+3+padd_X,2)
      common /bry_zeta_south/ zetabry_south, zetabry_south_dt
      real ubarbry_south(-2:Lm+3+padd_X),
     &    ubarbry_south_dt(-2:Lm+3+padd_X,2)
     &    ,vbarbry_south(-2:Lm+3+padd_X),
     &    vbarbry_south_dt(-2:Lm+3+padd_X,2)
      common /bry_ubar_south/ ubarbry_south, ubarbry_south_dt,
     &                        vbarbry_south, vbarbry_south_dt
      real ubry_south(-2:Lm+3+padd_X,N),
     &    ubry_south_dt(-2:Lm+3+padd_X,N,2)
     &    ,vbry_south(-2:Lm+3+padd_X,N),
     &    vbry_south_dt(-2:Lm+3+padd_X,N,2)
      common /bry_u_south/ ubry_south, ubry_south_dt,
     &                     vbry_south, vbry_south_dt
      real wbry_south(-2:Lm+3+padd_X,N),
     &    wbry_south_dt(-2:Lm+3+padd_X,N,2)
      common /bry_w_south/ wbry_south, wbry_south_dt
      real zetabry_north(-2:Lm+3+padd_X),
     &    zetabry_north_dt(-2:Lm+3+padd_X,2)
      common /bry_zeta_north/ zetabry_north, zetabry_north_dt
      real ubarbry_north(-2:Lm+3+padd_X),
     &    ubarbry_north_dt(-2:Lm+3+padd_X,2)
     &    ,vbarbry_north(-2:Lm+3+padd_X),
     &    vbarbry_north_dt(-2:Lm+3+padd_X,2)
      common /bry_ubar_north/ ubarbry_north, ubarbry_north_dt,
     &                        vbarbry_north, vbarbry_north_dt
      real ubry_north(-2:Lm+3+padd_X,N),
     &    ubry_north_dt(-2:Lm+3+padd_X,N,2)
     &    ,vbry_north(-2:Lm+3+padd_X,N),
     &    vbry_north_dt(-2:Lm+3+padd_X,N,2)
      common /bry_u_north/ ubry_north, ubry_north_dt,
     &                     vbry_north, vbry_north_dt
      real wbry_north(-2:Lm+3+padd_X,N),
     &    wbry_north_dt(-2:Lm+3+padd_X,N,2)
      common /bry_w_north/ wbry_north, wbry_north_dt
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
      tau_in=dtfast*tauM_in
      tau_out=dtfast*tauM_out
      grad = 1.D0
      if (.not.WEST_INTER .and. .not.SOUTH_INTER) then
        grad(Istr,Jstr) = 0.5D0
      endif
      if (.not.EAST_INTER .and. .not.SOUTH_INTER) then
        grad(Iend,Jstr) = 0.5D0
      endif
      if (.not.WEST_INTER .and. .not.NORTH_INTER) then
        grad(Istr,Jend+1) = 0.5D0
      endif
      if (.not.EAST_INTER .and. .not.NORTH_INTER) then
        grad(Iend,Jend+1) = 0.5D0
      endif
      if (.not.SOUTH_INTER) then
        do i=Istr,Iend
          cff=-sqrt(2.D0*g/(h(i,Jstr)+h(i,Jstr-1)+
     &                    zeta(i,Jstr,kstp)+zeta(i,Jstr-1,kstp)))
          vbar(i,Jstr,knew)=cff
     &               *(0.5D0*(zeta(i,Jstr-1,knew)+zeta(i,Jstr,knew))
     &                                           -zetabry_south(i)
     &             )*grad(i,Jstr)
     &                                           +vbarbry_south(i)
          vbar(i,Jstr,knew)=vbar(i,Jstr,knew)*vmask(i,Jstr)
        enddo
      endif
      if (.not.NORTH_INTER) then
        do i=Istr,Iend
        cff=sqrt(2.D0*g/(h(i,Jend)+h(i,Jend+1)
     &                +zeta(i,Jend,kstp)+zeta(i,Jend+1,kstp)))
          vbar(i,Jend+1,knew)=cff
     &               *(0.5D0*( zeta(i,Jend,knew)+zeta(i,Jend+1,knew))
     &                                            -zetabry_north(i)
     &               )*grad(i,Jend+1)
     &                                            +vbarbry_north(i)
          vbar(i,Jend+1,knew)=vbar(i,Jend+1,knew)*vmask(i,Jend+1)
        enddo
      endif
      if (.not.WEST_INTER) then
        do j=JstrV,Jend
          cff=sqrt(0.5D0*g*(h(Istr-1,j-1)+h(Istr-1,j)+
     &                    zeta(Istr-1,j-1,kstp)+zeta(Istr-1,j,kstp)))
          cx=dtfast*cff*0.5D0*(pm(Istr-1,j-1)+pm(Istr-1,j))
          vbar(Istr-1,j,knew)=( vbar(Istr-1,j,kstp)
     &                         +cx*vbar(Istr,j,knew) )/(1.D0+cx)
     &                        *vmask(Istr-1,j)
        enddo
      endif
      if (.not.EAST_INTER) then
        do j=JstrV,Jend
          cff=sqrt(0.5D0*g*(h(Iend+1,j-1)+h(Iend+1,j)+
     &                    zeta(Iend+1,j-1,kstp)+zeta(Iend+1,j,kstp)))
          cx=dtfast*cff*0.5D0*(pm(Iend+1,j-1)+pm(Iend+1,j))
          vbar(Iend+1,j,knew)=( vbar(Iend+1,j,kstp)
     &                         +cx*vbar(Iend,j,knew) )/(1.D0+cx)
     &                        *vmask(Iend+1,j)
        enddo
      endif
      if (.not.WEST_INTER .and. .not.SOUTH_INTER) then
        vbar(Istr-1,Jstr,knew)=0.5D0*( vbar(Istr-1,Jstr+1,knew)
     &                              +vbar(Istr  ,Jstr  ,knew))
     &                         *vmask(Istr-1,Jstr)
      endif
      if (.not.EAST_INTER .and. .not.SOUTH_INTER) then
        vbar(Iend+1,Jstr,knew)=0.5D0*( vbar(Iend+1,Jstr+1,knew)
     &                              +vbar(Iend  ,Jstr  ,knew))
     &                         *vmask(Iend+1,Jstr)
      endif
      if (.not.WEST_INTER .and. .not.NORTH_INTER) then
        vbar(Istr-1,Jend+1,knew)=0.5D0*( vbar(Istr-1,Jend,knew)
     &                                +vbar(Istr,Jend+1,knew))
     &                         *vmask(Istr-1,Jend+1)
      endif
      if (.not.EAST_INTER .and. .not.NORTH_INTER) then
        vbar(Iend+1,Jend+1,knew)=0.5D0*( vbar(Iend+1,Jend,knew)
     &                                +vbar(Iend,Jend+1,knew))
     &                         *vmask(Iend+1,Jend+1)
      endif
      return
      end
