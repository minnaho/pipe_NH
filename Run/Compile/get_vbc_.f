      subroutine get_vbc
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
      call get_smflux
      call get_stflux (itemp)
      call get_sst
      call get_stflux (isalt)
      call get_sss
      call get_srflux
      return
      end
      subroutine set_vbc (tile)
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
      real A2d(N2d,NSA,0:NPP-1), A3d(N3d,5,0:NPP-1)
      integer*4 B2d(N2d,0:NPP-1)
      common /private_scratch/ A2d,A3d
      integer*4 tile, trd, omp_get_thread_num
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
      call set_vbc_tile (Istr,Iend,Jstr,Jend,
     &         A2d(1,1,trd),A2d(1,2,trd),A2d(1,3,trd))
      return
      end
      subroutine set_vbc_tile (Istr,Iend,Jstr,Jend, wrk,wrk1,wrk2)
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
!$AGRIF_DO_NOT_TREAT
      INTEGER*4 :: ocean_grid_comm
      common /cpl_comm/ ocean_grid_comm
!$AGRIF_END_DO_NOT_TREAT
      integer*4 i,j,is, Istr,Iend,Jstr,Jend
      real cff,cff1
      real wrk(Istr-2:Iend+2,Jstr-2:Jend+2)
      real wrk1(Istr-2:Iend+2,Jstr-2:Jend+2)
      real wrk2(Istr-2:Iend+2,Jstr-2:Jend+2)
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
      call set_smflux_tile (Istr,Iend,Jstr,Jend)
      call set_stflux_tile (Istr,Iend,Jstr,Jend, itemp)
      call set_srflux_tile (Istr,Iend,Jstr,Jend)
      call set_sst_tile (Istr,Iend,Jstr,Jend)
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,itemp)=stflx(i,j,itemp)+
     &                  dqdt(i,j)*(t(i,j,N,nrhs,itemp)-sst(i,j))
        enddo
      enddo
      call set_stflux_tile (Istr,Iend,Jstr,Jend, isalt)
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)*t(i,j,N,nrhs,isalt)
        enddo
      enddo
      call set_sss_tile (Istr,Iend,Jstr,Jend)
      do j=JstrR,JendR
        do i=IstrR,IendR
          stflx(i,j,isalt)=stflx(i,j,isalt)+
     &                  dqdt(i,j)*(t(i,j,N,nrhs,isalt)-sss(i,j))
        enddo
      enddo
      call ana_btflux_tile (Istr,Iend,Jstr,Jend, itemp)
      call ana_btflux_tile (Istr,Iend,Jstr,Jend, isalt)
      do j=JstrR,JendR
        do i=IstrR,IendR
          btflx(i,j,isalt)=btflx(i,j,isalt)*t(i,j,1,nrhs,isalt)
        enddo
      enddo
      if (Zob.ne.0.D0) then
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            cff=vonKar/LOG((z_r(i,j,1)-z_w(i,j,0))/Zob)
            wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25D0*(v(i  ,j,1,nrhs)+v(i  ,j+1,1,nrhs)+
     &                v(i-1,j,1,nrhs)+v(i-1,j+1,1,nrhs))
            bustr(i,j)=0.5D0*(wrk(i-1,j)+wrk(i,j))*u(i,j,1,nrhs)*
     &                 SQRT(u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff*cff)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25D0*(u(i,j  ,1,nrhs)+u(i+1,j,1,nrhs)+
     &                u(i,j-1,1,nrhs)+u(i+1,j-1,1,nrhs))
            bvstr(i,j)=0.5D0*(wrk(i,j-1)+wrk(i,j))*v(i,j,1,nrhs)*
     &                 SQRT(cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs))
          enddo
        enddo
      elseif (rdrg2.gt.0.D0) then
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25D0*(v(i,j,1,nrhs)+v(i,j+1,1,nrhs)+v(i-1,j,1,nrhs)+
     &                                            v(i-1,j+1,1,nrhs))
            bustr(i,j)=u(i,j,1,nrhs)*(rdrg2*sqrt(
     &                         u(i,j,1,nrhs)*u(i,j,1,nrhs)+cff*cff
     &                                                            ))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25D0*(u(i,j,1,nrhs)+u(i+1,j,1,nrhs)+u(i,j-1,1,nrhs)+
     &                                            u(i+1,j-1,1,nrhs))
            bvstr(i,j)=v(i,j,1,nrhs)*(rdrg2*sqrt(
     &                         cff*cff+v(i,j,1,nrhs)*v(i,j,1,nrhs)
     &                                                            ))
          enddo
        enddo
      else
        do j=Jstr,Jend
          do i=Istr,Iend
            bustr(i,j)=rdrg*u(i,j,1,nrhs)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            bvstr(i,j)=rdrg*v(i,j,1,nrhs)
          enddo
        enddo
      endif
      cff=0.75D0/dt
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1=cff*0.5D0*(Hz(i-1,j,1)+Hz(i,j,1))
          bustr(i,j)=SIGN(1.D0, bustr(i,j))*
     &               MIN(ABS(bustr(i,j)),
     &                   ABS(u(i,j,1,nrhs))*cff1)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1=cff*0.5D0*(Hz(i,j-1,1)+Hz(i,j,1))
          bvstr(i,j)=SIGN(1.D0, bvstr(i,j))*
     &               MIN(ABS(bvstr(i,j)),
     &                   ABS(v(i,j,1,nrhs))*cff1)
        enddo
      enddo
      IF (.not.EAST_INTER) THEN
        DO j=Jstr,Jend
          bustr(Iend+1,j)=bustr(Iend,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Iend+1,j)=bvstr(Iend,j)
        END DO
      END IF
      IF (.not.WEST_INTER) THEN
        DO j=Jstr,Jend
          bustr(IstrU-1,j)=bustr(IstrU,j)
        END DO
        DO j=JstrV,Jend
          bvstr(Istr-1,j)=bvstr(Istr,j)
        END DO
      END IF
      IF (.not.NORTH_INTER) THEN
        DO i=IstrU,Iend
          bustr(i,Jend+1)=bustr(i,Jend)
        END DO
        DO i=Istr,Iend
          bvstr(i,Jend+1)=bvstr(i,Jend)
        END DO
      END IF
      IF (.not.SOUTH_INTER) THEN
        DO i=IstrU,Iend
          bustr(i,Jstr-1)=bustr(i,Jstr)
        END DO
        DO i=Istr,Iend
          bvstr(i,JstrV-1)=bvstr(i,JstrV)
        END DO
      END IF
      IF (.not.SOUTH_INTER.and..not.WEST_INTER) THEN
        bustr(Istr,Jstr-1)=0.5D0*(bustr(Istr+1,Jstr-1)+bustr(Istr,Jstr))
        bvstr(Istr-1,Jstr)=0.5D0*(bvstr(Istr,Jstr)+bvstr(Istr-1,Jstr+1))
      END IF
      IF (.not.SOUTH_INTER.and..not.EAST_INTER) THEN
        bustr(Iend+1,Jstr-1)=0.5D0*(bustr(Iend+1,Jstr)+bustr(Iend,
     &                                                          Jstr-1))
        bvstr(Iend+1,Jstr)=0.5D0*(bvstr(Iend+1,Jstr+1)+bvstr(Iend,Jstr))
      END IF
      IF (.not.NORTH_INTER.and..not.WEST_INTER) THEN
        bustr(Istr,Jend+1)=0.5D0*(bustr(Istr,Jend)+bustr(Istr+1,Jend+1))
        bvstr(Istr-1,Jend+1)=0.5D0*(bvstr(Istr-1,Jend)+bvstr(Istr,
     &                                                          Jend+1))
      END IF
      IF (.not.NORTH_INTER.and..not.EAST_INTER) THEN
        bustr(Iend+1,Jend+1)=0.5D0*(bustr(Iend+1,Jend)+bustr(Iend,
     &                                                          Jend+1))
        bvstr(Iend+1,Jend+1)=0.5D0*(bvstr(Iend+1,Jend)+bvstr(Iend,
     &                                                          Jend+1))
      END IF
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,bustr(-1,-1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,bvstr(-1,-1))
      call exchange_u2d_tile (Istr,Iend,Jstr,Jend,sustr(-1,-1))
      call exchange_v2d_tile (Istr,Iend,Jstr,Jend,svstr(-1,-1))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                      stflx(-1,-1,itemp))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,srflx(-1,-1))
      call exchange_r2d_tile (Istr,Iend,Jstr,Jend,
     &                                      stflx(-1,-1,isalt))
      return
      end
