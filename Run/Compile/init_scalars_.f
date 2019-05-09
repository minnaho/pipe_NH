      subroutine init_scalars (ierr)
      implicit none
      integer*4 ierr, i,j,itrc, lvar, lenstr
      character*20 nametrc, unitt
      character*60 vname1, vname3
      integer*4 omp_get_num_threads
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
      integer*4 filetype_his, filetype_avg
     &       ,filetype_dia, filetype_dia_avg
     &       ,filetype_diaM, filetype_diaM_avg
     &       ,filetype_diabio, filetype_diabio_avg
      parameter (filetype_his=1, filetype_avg=2,
     &           filetype_dia=3, filetype_dia_avg=4,
     &           filetype_diaM=5, filetype_diaM_avg=6,
     &           filetype_diabio=7,filetype_diabio_avg=8)
      integer*4 iloop, indextemp
      integer*4 indxTime, indxZ, indxUb, indxVb
      parameter (indxTime=1, indxZ=2, indxUb=3, indxVb=4)
      integer*4 indxU, indxV, indxT
      parameter (indxU=6, indxV=7, indxT=8)
      integer*4 indxS
      parameter (indxS=indxT+1)
      integer*4 indxBSD, indxBSS
      parameter (indxBSD=indxT+ntrc_salt+ntrc_pas+ntrc_bio+1,
     &           indxBSS=101)
      integer*4 indxO, indxW, indxR, indxVisc, indxDiff, indxAkv, 
     &                                                           indxAkt
      parameter (indxO=indxT+ntrc_salt+ntrc_pas+ntrc_bio+ntrc_sed
     &                      +ntrc_diats+ntrc_diauv+ntrc_diabio+1,
     &           indxW=indxO+1, indxR=indxO+2, indxVisc=indxO+3,
     &           indxDiff=indxO+4,indxAkv=indxO+5, indxAkt=indxO+6)
      integer*4 indxAks
      parameter (indxAks=indxAkt+4)
      integer*4 indxHbl
      parameter (indxHbl=indxAkt+5)
      integer*4 indxHbbl
      parameter (indxHbbl=indxAkt+6)
      integer*4 indxSSH
      parameter (indxSSH=indxAkt+12)
      integer*4 indxSUSTR, indxSVSTR
      parameter (indxSUSTR=indxSSH+1, indxSVSTR=indxSSH+2)
      integer*4 indxTime2
      parameter (indxTime2=indxSSH+3)
      integer*4 indxShflx, indxShflx_rsw
      parameter (indxShflx=indxSSH+4)
      integer*4 indxSwflx
      parameter (indxSwflx=indxShflx+1, indxShflx_rsw=indxShflx+2)
      integer*4 indxSST, indxdQdSST
      parameter (indxSST=indxShflx_rsw+1, indxdQdSST=indxShflx_rsw+2)
      integer*4 indxSSS
      parameter (indxSSS=indxSST+2)
      integer*4 indxWstr
      parameter (indxWstr=indxSUSTR+21)
      integer*4 indxUWstr
      parameter (indxUWstr=indxSUSTR+22)
      integer*4 indxVWstr
      parameter (indxVWstr=indxSUSTR+23)
      integer*4 indxBostr
      parameter (indxBostr=indxSUSTR+24)
      integer*4 indxWWA,indxWWD,indxWWP,indxWEB
      parameter (indxWWA=indxSUSTR+32, indxWWD=indxWWA+1,
     &           indxWWP=indxWWA+2
     &                             )
      integer*4 r2dvar, u2dvar, v2dvar, p2dvar, r3dvar,
     &                u3dvar, v3dvar, p3dvar, w3dvar, b3dvar
      parameter (r2dvar=0, u2dvar=1, v2dvar=2, p2dvar=3,
     & r3dvar=4, u3dvar=5, v3dvar=6, p3dvar=7, w3dvar=8,b3dvar=12)
      integer*4 xi_rho,xi_u, eta_rho,eta_v
      integer*4 ncidfrc, ncidbulk, ncidclm,  ntsms ,
     &        ntsrf,  ntssh,  ntsst, ntsss, ntuclm,
     &        ntbulk, ncidqbar, ntqbar, ntww
      integer*4 nttclm(NT), ntstf(NT), nttsrc(NT)
      integer*4 ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV
      integer*4 rstT(NT)
      integer*4  ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,   hisV,   hisR,    hisHbl, hisHbbl
     &      , hisO,   hisW,   hisVisc, hisDiff
     &      , hisAkv, hisAkt, hisAks
      integer*4 hisT(NT)
      integer*4 ncidavg, nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ, avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUwstr, avgVwstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgU,   avgV,   avgR,    avgHbl, avgHbbl
     &      , avgO,   avgW,   avgVisc, avgDiff
     &      , avgAkv, avgAkt, avgAks
      integer*4 avgT(NT)
      logical wrthis(500+NT)
     &      , wrtavg(500+NT)
      common/incscrum/
     &        ncidfrc, ncidbulk,ncidclm, ntsms, ntsrf, ntssh, ntsst
     &      , ntuclm, ntsss, ntbulk, ncidqbar, ntqbar, ntww
     &      , xi_rho,  xi_u
     &      , eta_rho, eta_v
     &                        ,  nttclm, ntstf, nttsrc
     &      , ncidrst, nrecrst,  nrpfrst
     &      , rstTime, rstTime2, rstTstep, rstZ,    rstUb,  rstVb
     &                         , rstU,    rstV,   rstT
     &      , ncidhis, nrechis,  nrpfhis
     &      , hisTime, hisTime2, hisTstep, hisZ,    hisUb,  hisVb
     &      , hisBostr, hisWstr, hisUWstr, hisVWstr
     &      , hisShflx, hisSwflx, hisShflx_rsw
     &      , hisU,    hisV,     hisT,    hisR
     &      , hisO,    hisW,     hisVisc, hisDiff
     &      , hisAkv,  hisAkt,   hisAks
     &      , hisHbl,  hisHbbl
     &      , ncidavg,  nrecavg,  nrpfavg
     &      , avgTime, avgTime2, avgTstep, avgZ,    avgUb,  avgVb
     &      , avgBostr, avgWstr, avgUWstr, avgVWstr
     &      , avgShflx, avgSwflx, avgShflx_rsw
     &      , avgU,    avgV,     avgT,     avgR
     &      , avgO,    avgW,     avgVisc,  avgDiff
     &      , avgAkv,  avgAkt,   avgAks
     &      , avgHbl,  avgHbbl
     &      , wrthis
     &      , wrtavg
      character*80 date_str, title, start_date
      character*80 ininame,  grdname,  hisname
     &         ,   rstname,  frcname,  bulkname,  usrname
     &         ,   qbarname, tsrcname
     &                                ,   avgname
     &                                ,   bry_file
      character*75  vname(20, 500)
      common /cncscrum/       date_str,   title,  start_date
     &         ,   ininame,  grdname, hisname
     &         ,   rstname,  frcname, bulkname,  usrname
     &         ,   qbarname, tsrcname
     &                                ,  avgname
     &                                ,   bry_file
     &                      ,  vname
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
C$OMP PARALLEL
C$OMP CRITICAL (isca_cr_rgn)
      numthreads=omp_get_num_threads()
C$OMP END CRITICAL (isca_cr_rgn)
C$OMP END PARALLEL
      if (mynode.eq.0) write(stdout,'(1x,A,3(1x,A,I3),A)') 'NUMBER',
     &    'OF THREADS:',numthreads,'BLOCKING:',NSUB_X,'x',NSUB_E,'.'
      if (numthreads.gt.NPP) then
        if (mynode.eq.0) write(stdout,'(/1x,A,I3/)')
     &    'ERROR: Requested number of threads exceeds NPP =', NPP
        ierr=ierr+1
      elseif (mod(NSUB_X*NSUB_E,numthreads).ne.0) then
        if (mynode.eq.0) write(stdout,
     &                '(/1x,A,1x,A,I3,4x,A,I3,4x,A,I4,A)') 'ERROR:',
     &                'wrong choice of numthreads =', numthreads,
     &                'NSUB_X =', NSUB_X, 'NSUB_E =', NSUB_E, '.'
        ierr=ierr+1
      endif
      time=0.D0
      tdays=0.D0
      PREDICTOR_2D_STEP=.FALSE.
      iic=0
      kstp=1
      krhs=1
      knew=1
      ntstart=1
      nstp=1
      nrhs=1
      nnew=1
      nfast=1
      synchro_flag=.true.
      first_time=0
      may_day_flag=0
      do j=0,NPP
        do i=0,31
          proc(i,j)=0
          CPU_time(i,j)=0.D0
        enddo
      enddo
      trd_count=0
      nrecrst=0
      nrechis=0
      nrecavg=0
      tile_count=0
      bc_count=0
      avgke=0.D0
      avgpe=0.D0
      avgkp=0.D0
      volume=0.D0
      hmin=+1.D+20
      hmax=-1.D+20
      grdmin=+1.D+20
      grdmax=-1.D+20
      Cu_min=+1.D+20
      Cu_max=-1.D+20
      lonmin=+1.D+20
      lonmax=-1.D+20
      latmin=+1.D+20
      latmax=-1.D+20
      bc_crss=0.D0
      rx0=-1.D+20
      rx1=-1.D+20
      ncidrst=-1
      ncidhis=-1
      ncidavg=-1
      ncidfrc=-1
      ncidbulk=-1
      ncidclm=-1
      ncidqbar=-1
       bry_id=-1
      call get_date (date_str)
      csmooth=dt/(86400.D0)
      ismooth=int(1.D0/csmooth)
      vname(:,:)='  '
      vname(1,indxTime)='ocean_time'
      vname(2,indxTime)='time since initialization'
      vname(3,indxTime)='second'
      vname(4,indxTime)='time, scalar, series'
      vname(5,indxTime)='time'
      vname(6,indxTime)='    '
      vname(7,indxTime)='T'
      vname(8,indxTime)='.NO.'
      vname(9,indxTime)=' '
      vname(1,indxTime2)='time'
      vname(2,indxTime2)='time since initialization'
      vname(3,indxTime2)='second'
      vname(4,indxTime2)='time, scalar, series'
      vname(5,indxTime2)='time'
      vname(6,indxTime2)='  '
      vname(7,indxTime2)='T '
      vname(8,indxTime2)='.NO.'
      vname(9,indxTime2)=' '
      vname(1,indxZ)='zeta                                        '
      vname(2,indxZ)='free-surface                                '
      vname(3,indxZ)='meter                                       '
      vname(4,indxZ)='free-surface, scalar, series                '
      vname(5,indxZ)='sea_surface_height                          '
      vname(6,indxZ)='lat_rho lon_rho                             '
      vname(7,indxZ)='                                            '
      vname(1,indxUb)='ubar                                       '
      vname(2,indxUb)='vertically integrated u-momentum component '
      vname(3,indxUb)='meter second-1                             '
      vname(4,indxUb)='ubar-velocity, scalar, series              '
      vname(5,indxUb)='barotropic_sea_water_x_'
     &               // 'velocity_at_u_location'
      vname(6,indxUb)='lat_u lon_u                                '
      vname(7,indxUb)='                                           '
      vname(1,indxVb)='vbar                                       '
      vname(2,indxVb)='vertically integrated v-momentum component '
      vname(3,indxVb)='meter second-1                             '
      vname(4,indxVb)='vbar-velocity, scalar, series              '
      vname(5,indxVb)='barotropic_sea_water_y_velocity_at_'
     &                // 'v_location'
      vname(6,indxVb)='lat_v lon_v                                '
      vname(7,indxVb)='                                           '
      vname(11,indxVb)=' '
      vname(1,indxBostr)='bostr                                   '
      vname(2,indxBostr)='Kinematic bottom stress                 '
      vname(3,indxBostr)='N/m2                                    '
      vname(4,indxBostr)='                                        '
      vname(5,indxBostr)='                                        '
      vname(6,indxBostr)='lat_rho lon_rho                         '
      vname(7,indxBostr)='                                        '
      vname(1,indxWstr)='wstr                                     '
      vname(2,indxWstr)='Kinematic wind stress                    '
      vname(3,indxWstr)='N/m2                                     '
      vname(4,indxWstr)='                                         '
      vname(5,indxWstr)='magnitude_of_surface_downward_stress     '
      vname(6,indxWstr)='lat_rho lon_rho                          '
      vname(7,indxWstr)='                                         '
      vname(1,indxUWstr)='sustr                                   '
      vname(2,indxUWstr)='Kinematic u wind stress component       '
      vname(3,indxUWstr)='N/m2                                    '
      vname(4,indxUWstr)='                                        '
      vname(5,indxUWstr)='surface_downward_eastward_stress        '
      vname(6,indxUWstr)='lat_u lon_u                             '
      vname(7,indxUWstr)='                                        '
      vname(8,indxUWstr)='                                        '
      vname(1,indxVWstr)='svstr                                   '
      vname(2,indxVWstr)='Kinematic v wind stress component       '
      vname(3,indxVWstr)='N/m2                                    '
      vname(4,indxVWstr)='                                        '
      vname(5,indxVWstr)='surface_downward_northward_stress       '
      vname(6,indxVWstr)='lat_v lon_v                             '
      vname(7,indxVWstr)='                                        '
      vname(1,indxU)='u                                           '
      vname(2,indxU)='u-momentum component                        '
      vname(3,indxU)='meter second-1                              '
      vname(4,indxU)='u-velocity, scalar, series                  '
      vname(5,indxU)='sea_water_x_velocity_at_u_location          '
      vname(6,indxU)='lat_u lon_u                                 '
      vname(7,indxU)='                                            '
      vname(1,indxV)='v                                           '
      vname(2,indxV)='v-momentum component                        '
      vname(3,indxV)='meter second-1                              '
      vname(4,indxV)='v-velocity, scalar, series                  '
      vname(5,indxV)='sea_water_y_velocity_at_v_location          '
      vname(6,indxV)='lat_v lon_v                                 '
      vname(7,indxV)='                                            '
      vname(1,indxT)='temp                                        '
      vname(2,indxT)='potential temperature                       '
      vname(3,indxT)='Celsius                                     '
      vname(4,indxT)='temperature, scalar, series                 '
      vname(5,indxT)='sea_water_potential_temperature             '
      vname(6,indxT)='lat_rho lon_rho                             '
      vname(7,indxT)='                                            '
      vname(1,indxS)='salt                                        '
      vname(2,indxS)='salinity                                    '
      vname(3,indxS)='PSU                                         '
      vname(4,indxS)='salinity, scalar, series                    '
      vname(5,indxS)='sea_water_salinity                          '
      vname(6,indxS)='lat_rho lon_rho                             '
      vname(7,indxS)='                                            '
      vname(1,indxShflx)='shflux                                  '
      vname(2,indxShflx)='surface net heat flux                   '
      vname(3,indxShflx)='Watts meter-2                           '
      vname(4,indxShflx)='surface heat flux, scalar, series       '
      vname(6,indxShflx)='lat_rho lon_rho                         '
      vname(7,indxShflx)='                                        '
      vname(1,indxSwflx)='swflux                                  '
      vname(2,indxSwflx)='surface freshwater flux (E-P)           '
      vname(3,indxSwflx)='centimeter day-1                        '
      vname(4,indxSwflx)='surface freshwater flux, scalar, series '
      vname(5,indxSwflx)='                                        '
      vname(6,indxSwflx)='lat_rho lon_rho                         '
      vname(7,indxSwflx)='                                        '
      vname(8,indxSwflx)='                                        '
      vname(1,indxShflx_rsw)='swrad                               '
      vname(2,indxShflx_rsw)='Short-wave surface radiation        '
      vname(3,indxShflx_rsw)='Watts meter-2                       '
      vname(4,indxShflx_rsw)='                                    '
      vname(5,indxShflx_rsw)='                                    '
      vname(6,indxShflx_rsw)='lat_rho lon_rho                     '
      vname(7,indxShflx_rsw)='                                    '
      vname(1,indxO)='omega                                       '
      vname(2,indxO)='S-coordinate vertical momentum component    '
      vname(3,indxO)='meter second-1                              '
      vname(4,indxO)='omega, scalar, series                       '
      vname(5,indxO)='                                            '
      vname(6,indxO)='lat_rho lon_rho                             '
      vname(7,indxO)='                                            '
      vname(1,indxW)='w                                           '
      vname(2,indxW)='vertical momentum component                 '
      vname(3,indxW)='meter second-1                              '
      vname(4,indxW)='w-velocity, scalar, series                  '
      vname(5,indxW)='upward_sea_water_velocity                   '
      vname(6,indxW)='lat_rho lon_rho                             '
      vname(1,indxR)='rho                                         '
      vname(2,indxR)='density anomaly                             '
      vname(3,indxR)='kilogram meter-3                            '
      vname(4,indxR)='density, scalar, series                     '
      vname(5,indxR)='sea_water_sigma_t                           '
      vname(7,indxR)='                                            '
      vname(1,indxDiff)='diff3d                                   '
      vname(2,indxDiff)='horizontal diffusivity coefficient       '
      vname(3,indxDiff)='meter4 second-1                          '
      vname(5,indxDiff)='ocean_tracer_xy_biharmonic_diffusivity   '
      vname(4,indxDiff)='diff3d, scalar, series                   '
      vname(6,indxDiff)='lat_rho lon_rho                          '
      vname(7,indxDiff)='                                         '
      vname(1,indxAkv)='AKv                                       '
      vname(2,indxAkv)='vertical viscosity coefficient            '
      vname(3,indxAkv)='meter2 second-1                           '
      vname(4,indxAkv)='AKv, scalar, series                       '
      vname(5,indxAkv)='ocean_vertical_momentum_diffusivity_'
     &                 / / 'at_w_location                         '
      vname(6,indxAkv)='lat_rho lon_rho                           '
      vname(7,indxAkv)='                                          '
      vname(1,indxAkt)='AKt                                       '
      vname(2,indxAkt)='temperature vertical diffusion coefficient'
      vname(3,indxAkt)='meter2 second-1                           '
      vname(4,indxAkt)='AKt, scalar, series                       '
      vname(5,indxAkt)='ocean_vertical_heat_diffusivity_'
     &                               / /  'at_w_location          '
      vname(6,indxAkt)='lat_rho lon_rho                           '
      vname(7,indxAkt)='                                          '
      vname(1,indxAks)='AKs                                       '
      vname(2,indxAks)='salinity vertical diffusion coefficient   '
      vname(3,indxAks)='meter2 second-1                           '
      vname(4,indxAks)='AKs, scalar, series                       '
      vname(5,indxAks)='ocean_vertical_salt_diffusivity_'
     &               / /  'at_w_location                          '
      vname(6,indxAks)='lat_rho lon_rho                           '
      vname(7,indxAks)='                                          '
      vname(1,indxHbl)='hbl                                       '
      vname(2,indxHbl)='depth of planetary boundary layer         '
      vname(3,indxHbl)='meter                                     '
      vname(4,indxHbl)='hbl, scalar, series                       '
      vname(5,indxHbl)='ocean_mixed_layer_thickness_defined_'
     &                / / 'by_mixing_scheme                       '
      vname(6,indxHbl)='lat_rho lon_rho                           '
      vname(7,indxHbl)='                                          '
      vname(1,indxHbbl)='hbbl                                     '
      vname(2,indxHbbl)='depth of bottom boundary layer           '
      vname(3,indxHbbl)='meter                                    '
      vname(4,indxHbbl)='hbbl, scalar, series                     '
      vname(5,indxHbbl)='                                         '
      vname(6,indxHbbl)='lat_rho lon_rho                          '
      vname(7,indxHbbl)='                                         '
      vname(1,indxSSH)='SSH                                       '
      vname(2,indxSSH)='sea surface height                        '
      vname(3,indxSSH)='meter                                     '
      vname(4,indxSSH)='SSH, scalar, series                       '
      vname(5,indxSSH)='sea_surface_height_above_sea_level        '
      vname(6,indxSSH)='lat_rho lon_rho                           '
      vname(7,indxSSH)='                                          '
      vname(1,indxSUSTR)='sustr                                   '
      vname(2,indxSUSTR)='surface u-momentum stress               '
      vname(3,indxSUSTR)='Newton meter-2                          '
      vname(4,indxSUSTR)='surface u-mom. stress, scalar, series   '
      vname(5,indxSUSTR)='surface_downward_x_stress               '
      vname(6,indxSUSTR)='lat_u lon_u                             '
      vname(7,indxSUSTR)='                                        '
      vname(1,indxSVSTR)='svstr                                   '
      vname(2,indxSVSTR)='surface v-momentum stress               '
      vname(3,indxSVSTR)='Newton meter-2                          '
      vname(4,indxSVSTR)='surface v-mom. stress, scalar, series   '
      vname(5,indxSVSTR)='surface_downward_y_stress               '
      vname(6,indxSVSTR)='lat_v lon_v                             '
      vname(7,indxSVSTR)='                                        '
      vname(1,indxSST)='SST                                       '
      vname(2,indxSST)='sea surface temperature                   '
      vname(3,indxSST)='Celsius                                   '
      vname(4,indxSST)='SST, scalar, series                       '
      vname(5,indxSST)='sea_surface_temperature                   '
      vname(6,indxSST)='lat_rho lon_rho                           '
      vname(7,indxSST)='                                          '
      vname(8,indxSST)=' '
      vname(1,indxdQdSST)='dQdSST                                 '
      vname(2,indxdQdSST)='surface heat flux sensitivity to SST   '
      vname(3,indxdQdSST)='Watts meter-2 Celsius-1                '
      vname(4,indxdQdSST)='dQdSST, scalar, series                 '
      vname(5,indxdQdSST)='heat_flux_correction                   '
      vname(6,indxdQdSST)='lat_rho lon_rho                        '
      vname(7,indxdQdSST)='                                       '
      vname(1,indxSSS)='SSS                                       '
      vname(2,indxSSS)='sea surface salinity                      '
      vname(3,indxSSS)='PSU                                       '
      vname(4,indxSSS)='SSS, scalar, series                       '
      vname(5,indxSSS)='sea_surface_salinity                      '
      vname(6,indxSSS)='lat_rho lon_rho                           '
      vname(7,indxSSS)='                                          '
      return
      end
