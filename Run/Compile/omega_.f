      subroutine omega (tile)
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
      call omega_tile (Istr,Iend,Jstr,Jend, A2d(1,1,trd)
     &                                     ,A2d(1,2,trd)
     &                                                  )
      return
      end
      subroutine omega_tile (Istr,Iend,Jstr,Jend, wrk
     &      ,Cu_adv
     &                      )
      implicit none
      integer*4 Istr,Iend,Jstr,Jend, i,j,k
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
      real wrk   (Istr-2:Iend+2    )
      real wrk2  (Istr-2:Iend+2    )
      real cu_adv(Istr-2:Iend+2,0:N)
      real gamma,Cdt,cff,my_cu_adv,c_min,c_max,cu_cut,fcu
      integer*4 count
      parameter(c_min=0.6D0,c_max=0.95D0,cu_cut=2.D0*c_max-c_min,
     &          fcu=4.D0*c_max*(c_max-c_min))
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
      real   ,dimension(0:N) :: Q_int
      integer*4 is,iii,jjj
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
      gamma = 0.08333333333333D0
      if(nrhs.eq.3) then
        cdt =     dt
      elseif(iic.eq.ntstart) then
        cdt = 0.5D0*dt
      else
        cdt = (1.D0-2*gamma)*dt
      endif
      count = 0
      do j=Jstr,Jend
        do i=Istr,Iend
          Wi(i,j,0)=0.D0
          We(i,j,0)=0.D0
        enddo
        do k=1,N,+1
          do i=Istr,Iend
            We(i,j,k)=We(i,j,k-1) -Huon(i+1,j,k)+Huon(i,j,k)
     &                            -Hvom(i,j+1,k)+Hvom(i,j,k)
            Cu_adv(i,k)=max(HUon(i+1,j,k),0.D0)-min(Huon(i,j,k),0.D0)
     &                 +max(Hvom(i,j+1,k),0.D0)-min(Hvom(i,j,k),0.D0)
          enddo
        enddo
      do is=1,Nsrc
        iii=Isrc_mpi(is,mynode)
        jjj=Jsrc_mpi(is,mynode)
        if (Dsrc(is).eq.2) then
          if (Istr.le.iii .and. iii.le.Iend+1 .and.
     &        Jstr.le.jjj .and. jjj.le.Jend+1.and.j.eq.jjj) then
            Q_int(1) = Qsrc(is,1)
            do k=2,N
              Q_int(k) = Q_int(k-1) + Qsrc(is,k)
            enddo
            do k=1,N
              We(iii,jjj,k)=We(iii,jjj,k)+Q_int(k)
            enddo
            We(iii,jjj,0)=We(iii,jjj,0)+Q_int(1)
          endif
        endif
      enddo
        do i=Istr,Iend
          wrk(i)=We(i,j,N)/(z_w(i,j,N)-z_w(i,j,0))
        enddo
        do k=N-1,1,-1
          do i=Istr,Iend
            We(i,j,k)=We(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          enddo
        enddo
        do i=Istr,Iend
          We(i,j,N)=0.D0
          Wi(i,j,N)=0.D0
        enddo
        do k=1,N
          do i=Istr,Iend
            cff = cdt*pm(i,j)*pn(i,j)/Hz(i,j,k)
            Cu_adv(i,k)=cff*(
     &                   Cu_adv(i,k)
     &                 + max(We(i,j,k),0.D0)-min(We(i,j,k-1),0.D0)
     &                       )
          enddo
        enddo
        do k=1,N-1
          do i=istr,iend
             my_cu_adv = max(Cu_adv(i,k+1),Cu_adv(i,k))
             if(my_cu_adv.lt.c_min) then
                 cff = 0.D0
             elseif(my_cu_adv.lt.cu_cut) then
                 cff = (my_cu_adv-c_min)**2
                 cff = cff / (fcu+cff)
                 count = count +1
             else
                 cff = (my_cu_adv-c_max)/my_cu_adv
             endif
             Wi(i,j,k) =     cff *We(i,j,k)
             We(i,j,k) = (1.D0-cff)*We(i,j,k)
          enddo
        enddo
      enddo
      if (.not.WEST_INTER) then
        do k=0,N
          do j=Jstr,Jend
            We(0,j,k)=We(1,j,k)
            Wi(0,j,k)=Wi(1,j,k)
          enddo
        enddo
      endif
      if (.not.EAST_INTER) then
        do k=0,N
          do j=Jstr,Jend
            We(Lmmpi+1,j,k)=We(Lmmpi,j,k)
            Wi(Lmmpi+1,j,k)=Wi(Lmmpi,j,k)
          enddo
        enddo
      endif
      if (.not.SOUTH_INTER) then
        do k=0,N
          do i=Istr,Iend
            We(i,0,k)=We(i,1,k)
            Wi(i,0,k)=Wi(i,1,k)
          enddo
        enddo
      endif
      if (.not.NORTH_INTER) then
        do k=0,N
          do i=Istr,Iend
            We(i,Mmmpi+1,k)=We(i,Mmmpi,k)
            Wi(i,Mmmpi+1,k)=Wi(i,Mmmpi,k)
          enddo
        enddo
      endif
      if (.not.WEST_INTER .and. .not.SOUTH_INTER) then
        do k=0,N
          We(0,0,k)=We(1,1,k)
          Wi(0,0,k)=Wi(1,1,k)
        enddo
      endif
      if (.not.WEST_INTER .and. .not.NORTH_INTER) then
        do k=0,N
          We(0,Mmmpi+1,k)=We(1,Mmmpi,k)
          Wi(0,Mmmpi+1,k)=Wi(1,Mmmpi,k)
        enddo
      endif
      if (.not.EAST_INTER .and. .not.SOUTH_INTER) then
        do k=0,N
          We(Lmmpi+1,0,k)=We(Lmmpi,1,k)
          Wi(Lmmpi+1,0,k)=Wi(Lmmpi,1,k)
        enddo
      endif
      if (.not.EAST_INTER .and. .not.NORTH_INTER) then
        do k=0,N
          We(Lmmpi+1,Mmmpi+1,k)=We(Lmmpi,Mmmpi,k)
          Wi(Lmmpi+1,Mmmpi+1,k)=Wi(Lmmpi,Mmmpi,k)
        enddo
      endif
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, We(-2,-2,0))
      call exchange_w3d_tile (Istr,Iend,Jstr,Jend, Wi(-2,-2,0))
      return
      end
