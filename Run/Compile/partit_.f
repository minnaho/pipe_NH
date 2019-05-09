      program partit
      implicit none
      integer*4 stdout, max_buff_size, maxdims, maxvars, maxnodes
      parameter (stdout=6,           max_buff_size=300*300*30,
     &           maxdims=40,         maxvars=64,     maxnodes=512)
      character*120 ncname0, ncname(0:maxnodes-1), string
      character*32 dimname(maxdims), varname(maxvars)
      logical XiPeriodic, EtaPeriodic, part_switch(maxvars),
     &                                           series(maxvars)
      real*8 buff(max_buff_size)
      integer*4 narg,  NP_XI, NNODES,  xi_rho,  id_xi_rho,  id_xi_psi,
     &        arg,  NP_ETA,   node,  xi_u,    id_xi_u,    id_xi_v,
     &        iargc, subXI,     ii,  eta_rho, id_eta_rho, id_eta_psi,
     &        ierr, subETA,     jj,  eta_v,   id_eta_v,   id_eta_u,
     &        ndims, nvars, ngatts,  tsize,   unlimdimid, varatts,
     &        i,j,k,  lstr,   lvar,  lenstr,  rec,  size, ncid0,
     &  ncid(0:maxnodes-1), dimid(maxdims),   dimsize(maxdims),
     &      varid(maxvars), vartype(maxvars), vardims(maxvars),
     &      start(maxdims), count(maxdims),   start1(maxdims),
     &      dimids(maxdims,maxvars),          ibuff(maxdims)
      common /partit_main/  buff,  ncid,  dimid,  dimsize, varid,
     &    vartype, vardims, start, count, start1, dimids,  ibuff
      integer*4 LLm, MMm
      integer*4 chunk_size_X, margin_X
      integer*4 chunk_size_E, margin_E
      integer*4 istrmpi, jstrmpi
      integer*4 iendmpi, jendmpi
      integer*4 Lmmpi, Mmmpi
      integer*4 nf_byte
      integer*4 nf_int1
      integer*4 nf_char
      integer*4 nf_short
      integer*4 nf_int2
      integer*4 nf_int
      integer*4 nf_float
      integer*4 nf_real
      integer*4 nf_double
      parameter (nf_byte = 1)
      parameter (nf_int1 = nf_byte)
      parameter (nf_char = 2)
      parameter (nf_short = 3)
      parameter (nf_int2 = nf_short)
      parameter (nf_int = 4)
      parameter (nf_float = 5)
      parameter (nf_real = nf_float)
      parameter (nf_double = 6)
      integer*4           nf_fill_byte
      integer*4           nf_fill_int1
      integer*4           nf_fill_char
      integer*4           nf_fill_short
      integer*4           nf_fill_int2
      integer*4           nf_fill_int
      real              nf_fill_float
      real              nf_fill_real
      doubleprecision   nf_fill_double
      parameter (nf_fill_byte = -127)
      parameter (nf_fill_int1 = nf_fill_byte)
      parameter (nf_fill_char = 0)
      parameter (nf_fill_short = -32767)
      parameter (nf_fill_int2 = nf_fill_short)
      parameter (nf_fill_int = -2147483647)
      parameter (nf_fill_float = 9.9692099683868690D+36)
      parameter (nf_fill_real = nf_fill_float)
      parameter (nf_fill_double = 9.9692099683868690D+36)
      integer*4 nf_nowrite
      integer*4 nf_write
      integer*4 nf_clobber
      integer*4 nf_noclobber
      integer*4 nf_fill
      integer*4 nf_nofill
      integer*4 nf_lock
      integer*4 nf_share
      integer*4 nf_64bit_offset
      integer*4 nf_sizehint_default
      integer*4 nf_align_chunk
      integer*4 nf_format_classic
      integer*4 nf_format_64bit
      integer*4 nf_diskless
      integer*4 nf_mmap
      parameter (nf_nowrite = 0)
      parameter (nf_write = 1)
      parameter (nf_clobber = 0)
      parameter (nf_noclobber = 4)
      parameter (nf_fill = 0)
      parameter (nf_nofill = 256)
      parameter (nf_lock = 1024)
      parameter (nf_share = 2048)
      parameter (nf_64bit_offset = 512)
      parameter (nf_sizehint_default = 0)
      parameter (nf_align_chunk = -1)
      parameter (nf_format_classic = 1)
      parameter (nf_format_64bit = 2)
      parameter (nf_diskless = 8)
      parameter (nf_mmap = 16)
      integer*4 nf_unlimited
      parameter (nf_unlimited = 0)
      integer*4 nf_global
      parameter (nf_global = 0)
      integer*4 nf_max_dims
      integer*4 nf_max_attrs
      integer*4 nf_max_vars
      integer*4 nf_max_name
      integer*4 nf_max_var_dims
      parameter (nf_max_dims = 1024)
      parameter (nf_max_attrs = 8192)
      parameter (nf_max_vars = 8192)
      parameter (nf_max_name = 256)
      parameter (nf_max_var_dims = nf_max_dims)
      integer*4 nf_noerr
      integer*4 nf_ebadid
      integer*4 nf_eexist
      integer*4 nf_einval
      integer*4 nf_eperm
      integer*4 nf_enotindefine
      integer*4 nf_eindefine
      integer*4 nf_einvalcoords
      integer*4 nf_emaxdims
      integer*4 nf_enameinuse
      integer*4 nf_enotatt
      integer*4 nf_emaxatts
      integer*4 nf_ebadtype
      integer*4 nf_ebaddim
      integer*4 nf_eunlimpos
      integer*4 nf_emaxvars
      integer*4 nf_enotvar
      integer*4 nf_eglobal
      integer*4 nf_enotnc
      integer*4 nf_ests
      integer*4 nf_emaxname
      integer*4 nf_eunlimit
      integer*4 nf_enorecvars
      integer*4 nf_echar
      integer*4 nf_eedge
      integer*4 nf_estride
      integer*4 nf_ebadname
      integer*4 nf_erange
      integer*4 nf_enomem
      integer*4 nf_evarsize
      integer*4 nf_edimsize
      integer*4 nf_etrunc
      parameter (nf_noerr = 0)
      parameter (nf_ebadid = -33)
      parameter (nf_eexist = -35)
      parameter (nf_einval = -36)
      parameter (nf_eperm = -37)
      parameter (nf_enotindefine = -38)
      parameter (nf_eindefine = -39)
      parameter (nf_einvalcoords = -40)
      parameter (nf_emaxdims = -41)
      parameter (nf_enameinuse = -42)
      parameter (nf_enotatt = -43)
      parameter (nf_emaxatts = -44)
      parameter (nf_ebadtype = -45)
      parameter (nf_ebaddim = -46)
      parameter (nf_eunlimpos = -47)
      parameter (nf_emaxvars = -48)
      parameter (nf_enotvar = -49)
      parameter (nf_eglobal = -50)
      parameter (nf_enotnc = -51)
      parameter (nf_ests = -52)
      parameter (nf_emaxname = -53)
      parameter (nf_eunlimit = -54)
      parameter (nf_enorecvars = -55)
      parameter (nf_echar = -56)
      parameter (nf_eedge = -57)
      parameter (nf_estride = -58)
      parameter (nf_ebadname = -59)
      parameter (nf_erange = -60)
      parameter (nf_enomem = -61)
      parameter (nf_evarsize = -62)
      parameter (nf_edimsize = -63)
      parameter (nf_etrunc = -64)
      integer*4  nf_fatal
      integer*4 nf_verbose
      parameter (nf_fatal = 1)
      parameter (nf_verbose = 2)
      character*80   nf_inq_libvers
      external       nf_inq_libvers
      character*80   nf_strerror
      external       nf_strerror
      logical        nf_issyserr
      external       nf_issyserr
      integer*4         nf_inq_base_pe
      external        nf_inq_base_pe
      integer*4         nf_set_base_pe
      external        nf_set_base_pe
      integer*4         nf_create
      external        nf_create
      integer*4         nf__create
      external        nf__create
      integer*4         nf__create_mp
      external        nf__create_mp
      integer*4         nf_open
      external        nf_open
      integer*4         nf__open
      external        nf__open
      integer*4         nf__open_mp
      external        nf__open_mp
      integer*4         nf_set_fill
      external        nf_set_fill
      integer*4         nf_set_default_format
      external        nf_set_default_format
      integer*4         nf_redef
      external        nf_redef
      integer*4         nf_enddef
      external        nf_enddef
      integer*4         nf__enddef
      external        nf__enddef
      integer*4         nf_sync
      external        nf_sync
      integer*4         nf_abort
      external        nf_abort
      integer*4         nf_close
      external        nf_close
      integer*4         nf_delete
      external        nf_delete
      integer*4         nf_inq
      external        nf_inq
      integer*4 nf_inq_path
      external nf_inq_path
      integer*4         nf_inq_ndims
      external        nf_inq_ndims
      integer*4         nf_inq_nvars
      external        nf_inq_nvars
      integer*4         nf_inq_natts
      external        nf_inq_natts
      integer*4         nf_inq_unlimdim
      external        nf_inq_unlimdim
      integer*4         nf_inq_format
      external        nf_inq_format
      integer*4         nf_def_dim
      external        nf_def_dim
      integer*4         nf_inq_dimid
      external        nf_inq_dimid
      integer*4         nf_inq_dim
      external        nf_inq_dim
      integer*4         nf_inq_dimname
      external        nf_inq_dimname
      integer*4         nf_inq_dimlen
      external        nf_inq_dimlen
      integer*4         nf_rename_dim
      external        nf_rename_dim
      integer*4         nf_inq_att
      external        nf_inq_att
      integer*4         nf_inq_attid
      external        nf_inq_attid
      integer*4         nf_inq_atttype
      external        nf_inq_atttype
      integer*4         nf_inq_attlen
      external        nf_inq_attlen
      integer*4         nf_inq_attname
      external        nf_inq_attname
      integer*4         nf_copy_att
      external        nf_copy_att
      integer*4         nf_rename_att
      external        nf_rename_att
      integer*4         nf_del_att
      external        nf_del_att
      integer*4         nf_put_att_text
      external        nf_put_att_text
      integer*4         nf_get_att_text
      external        nf_get_att_text
      integer*4         nf_put_att_int1
      external        nf_put_att_int1
      integer*4         nf_get_att_int1
      external        nf_get_att_int1
      integer*4         nf_put_att_int2
      external        nf_put_att_int2
      integer*4         nf_get_att_int2
      external        nf_get_att_int2
      integer*4         nf_put_att_int
      external        nf_put_att_int
      integer*4         nf_get_att_int
      external        nf_get_att_int
      integer*4         nf_put_att_real
      external        nf_put_att_real
      integer*4         nf_get_att_real
      external        nf_get_att_real
      integer*4         nf_put_att_double
      external        nf_put_att_double
      integer*4         nf_get_att_double
      external        nf_get_att_double
      integer*4         nf_def_var
      external        nf_def_var
      integer*4         nf_inq_var
      external        nf_inq_var
      integer*4         nf_inq_varid
      external        nf_inq_varid
      integer*4         nf_inq_varname
      external        nf_inq_varname
      integer*4         nf_inq_vartype
      external        nf_inq_vartype
      integer*4         nf_inq_varndims
      external        nf_inq_varndims
      integer*4         nf_inq_vardimid
      external        nf_inq_vardimid
      integer*4         nf_inq_varnatts
      external        nf_inq_varnatts
      integer*4         nf_rename_var
      external        nf_rename_var
      integer*4         nf_copy_var
      external        nf_copy_var
      integer*4         nf_put_var_text
      external        nf_put_var_text
      integer*4         nf_get_var_text
      external        nf_get_var_text
      integer*4         nf_put_var_int1
      external        nf_put_var_int1
      integer*4         nf_get_var_int1
      external        nf_get_var_int1
      integer*4         nf_put_var_int2
      external        nf_put_var_int2
      integer*4         nf_get_var_int2
      external        nf_get_var_int2
      integer*4         nf_put_var_int
      external        nf_put_var_int
      integer*4         nf_get_var_int
      external        nf_get_var_int
      integer*4         nf_put_var_real
      external        nf_put_var_real
      integer*4         nf_get_var_real
      external        nf_get_var_real
      integer*4         nf_put_var_double
      external        nf_put_var_double
      integer*4         nf_get_var_double
      external        nf_get_var_double
      integer*4         nf_put_var1_text
      external        nf_put_var1_text
      integer*4         nf_get_var1_text
      external        nf_get_var1_text
      integer*4         nf_put_var1_int1
      external        nf_put_var1_int1
      integer*4         nf_get_var1_int1
      external        nf_get_var1_int1
      integer*4         nf_put_var1_int2
      external        nf_put_var1_int2
      integer*4         nf_get_var1_int2
      external        nf_get_var1_int2
      integer*4         nf_put_var1_int
      external        nf_put_var1_int
      integer*4         nf_get_var1_int
      external        nf_get_var1_int
      integer*4         nf_put_var1_real
      external        nf_put_var1_real
      integer*4         nf_get_var1_real
      external        nf_get_var1_real
      integer*4         nf_put_var1_double
      external        nf_put_var1_double
      integer*4         nf_get_var1_double
      external        nf_get_var1_double
      integer*4         nf_put_vara_text
      external        nf_put_vara_text
      integer*4         nf_get_vara_text
      external        nf_get_vara_text
      integer*4         nf_put_vara_int1
      external        nf_put_vara_int1
      integer*4         nf_get_vara_int1
      external        nf_get_vara_int1
      integer*4         nf_put_vara_int2
      external        nf_put_vara_int2
      integer*4         nf_get_vara_int2
      external        nf_get_vara_int2
      integer*4         nf_put_vara_int
      external        nf_put_vara_int
      integer*4         nf_get_vara_int
      external        nf_get_vara_int
      integer*4         nf_put_vara_real
      external        nf_put_vara_real
      integer*4         nf_get_vara_real
      external        nf_get_vara_real
      integer*4         nf_put_vara_double
      external        nf_put_vara_double
      integer*4         nf_get_vara_double
      external        nf_get_vara_double
      integer*4         nf_put_vars_text
      external        nf_put_vars_text
      integer*4         nf_get_vars_text
      external        nf_get_vars_text
      integer*4         nf_put_vars_int1
      external        nf_put_vars_int1
      integer*4         nf_get_vars_int1
      external        nf_get_vars_int1
      integer*4         nf_put_vars_int2
      external        nf_put_vars_int2
      integer*4         nf_get_vars_int2
      external        nf_get_vars_int2
      integer*4         nf_put_vars_int
      external        nf_put_vars_int
      integer*4         nf_get_vars_int
      external        nf_get_vars_int
      integer*4         nf_put_vars_real
      external        nf_put_vars_real
      integer*4         nf_get_vars_real
      external        nf_get_vars_real
      integer*4         nf_put_vars_double
      external        nf_put_vars_double
      integer*4         nf_get_vars_double
      external        nf_get_vars_double
      integer*4         nf_put_varm_text
      external        nf_put_varm_text
      integer*4         nf_get_varm_text
      external        nf_get_varm_text
      integer*4         nf_put_varm_int1
      external        nf_put_varm_int1
      integer*4         nf_get_varm_int1
      external        nf_get_varm_int1
      integer*4         nf_put_varm_int2
      external        nf_put_varm_int2
      integer*4         nf_get_varm_int2
      external        nf_get_varm_int2
      integer*4         nf_put_varm_int
      external        nf_put_varm_int
      integer*4         nf_get_varm_int
      external        nf_get_varm_int
      integer*4         nf_put_varm_real
      external        nf_put_varm_real
      integer*4         nf_get_varm_real
      external        nf_get_varm_real
      integer*4         nf_put_varm_double
      external        nf_put_varm_double
      integer*4         nf_get_varm_double
      external        nf_get_varm_double
      integer*4 nf_ubyte
      integer*4 nf_ushort
      integer*4 nf_uint
      integer*4 nf_int64
      integer*4 nf_uint64
      integer*4 nf_string
      integer*4 nf_vlen
      integer*4 nf_opaque
      integer*4 nf_enum
      integer*4 nf_compound
      parameter (nf_ubyte = 7)
      parameter (nf_ushort = 8)
      parameter (nf_uint = 9)
      parameter (nf_int64 = 10)
      parameter (nf_uint64 = 11)
      parameter (nf_string = 12)
      parameter (nf_vlen = 13)
      parameter (nf_opaque = 14)
      parameter (nf_enum = 15)
      parameter (nf_compound = 16)
      integer*4           nf_fill_ubyte
      integer*4           nf_fill_ushort
      parameter (nf_fill_ubyte = 255)
      parameter (nf_fill_ushort = 65535)
      integer*4 nf_format_netcdf4
      parameter (nf_format_netcdf4 = 3)
      integer*4 nf_format_netcdf4_classic
      parameter (nf_format_netcdf4_classic = 4)
      integer*4 nf_netcdf4
      parameter (nf_netcdf4 = 4096)
      integer*4 nf_classic_model
      parameter (nf_classic_model = 256)
      integer*4 nf_chunk_seq
      parameter (nf_chunk_seq = 0)
      integer*4 nf_chunk_sub
      parameter (nf_chunk_sub = 1)
      integer*4 nf_chunk_sizes
      parameter (nf_chunk_sizes = 2)
      integer*4 nf_endian_native
      parameter (nf_endian_native = 0)
      integer*4 nf_endian_little
      parameter (nf_endian_little = 1)
      integer*4 nf_endian_big
      parameter (nf_endian_big = 2)
      integer*4 nf_chunked
      parameter (nf_chunked = 0)
      integer*4 nf_contiguous
      parameter (nf_contiguous = 1)
      integer*4 nf_nochecksum
      parameter (nf_nochecksum = 0)
      integer*4 nf_fletcher32
      parameter (nf_fletcher32 = 1)
      integer*4 nf_noshuffle
      parameter (nf_noshuffle = 0)
      integer*4 nf_shuffle
      parameter (nf_shuffle = 1)
      integer*4 nf_szip_ec_option_mask
      parameter (nf_szip_ec_option_mask = 4)
      integer*4 nf_szip_nn_option_mask
      parameter (nf_szip_nn_option_mask = 32)
      integer*4 nf_mpiio
      parameter (nf_mpiio = 8192)
      integer*4 nf_mpiposix
      parameter (nf_mpiposix = 16384)
      integer*4 nf_pnetcdf
      parameter (nf_pnetcdf = 32768)
      integer*4 nf_independent
      parameter (nf_independent = 0)
      integer*4 nf_collective
      parameter (nf_collective = 1)
      integer*4 nf_ehdferr
      parameter (nf_ehdferr = -101)
      integer*4 nf_ecantread
      parameter (nf_ecantread = -102)
      integer*4 nf_ecantwrite
      parameter (nf_ecantwrite = -103)
      integer*4 nf_ecantcreate
      parameter (nf_ecantcreate = -104)
      integer*4 nf_efilemeta
      parameter (nf_efilemeta = -105)
      integer*4 nf_edimmeta
      parameter (nf_edimmeta = -106)
      integer*4 nf_eattmeta
      parameter (nf_eattmeta = -107)
      integer*4 nf_evarmeta
      parameter (nf_evarmeta = -108)
      integer*4 nf_enocompound
      parameter (nf_enocompound = -109)
      integer*4 nf_eattexists
      parameter (nf_eattexists = -110)
      integer*4 nf_enotnc4
      parameter (nf_enotnc4 = -111)
      integer*4 nf_estrictnc3
      parameter (nf_estrictnc3 = -112)
      integer*4 nf_enotnc3
      parameter (nf_enotnc3 = -113)
      integer*4 nf_enopar
      parameter (nf_enopar = -114)
      integer*4 nf_eparinit
      parameter (nf_eparinit = -115)
      integer*4 nf_ebadgrpid
      parameter (nf_ebadgrpid = -116)
      integer*4 nf_ebadtypid
      parameter (nf_ebadtypid = -117)
      integer*4 nf_etypdefined
      parameter (nf_etypdefined = -118)
      integer*4 nf_ebadfield
      parameter (nf_ebadfield = -119)
      integer*4 nf_ebadclass
      parameter (nf_ebadclass = -120)
      integer*4 nf_emaptype
      parameter (nf_emaptype = -121)
      integer*4 nf_elatefill
      parameter (nf_elatefill = -122)
      integer*4 nf_elatedef
      parameter (nf_elatedef = -123)
      integer*4 nf_edimscale
      parameter (nf_edimscale = -124)
      integer*4 nf_enogrp
      parameter (nf_enogrp = -125)
      integer*4 nf_create_par
      external nf_create_par
      integer*4 nf_open_par
      external nf_open_par
      integer*4 nf_var_par_access
      external nf_var_par_access
      integer*4 nf_inq_ncid
      external nf_inq_ncid
      integer*4 nf_inq_grps
      external nf_inq_grps
      integer*4 nf_inq_grpname
      external nf_inq_grpname
      integer*4 nf_inq_grpname_full
      external nf_inq_grpname_full
      integer*4 nf_inq_grpname_len
      external nf_inq_grpname_len
      integer*4 nf_inq_grp_parent
      external nf_inq_grp_parent
      integer*4 nf_inq_grp_ncid
      external nf_inq_grp_ncid
      integer*4 nf_inq_grp_full_ncid
      external nf_inq_grp_full_ncid
      integer*4 nf_inq_varids
      external nf_inq_varids
      integer*4 nf_inq_dimids
      external nf_inq_dimids
      integer*4 nf_def_grp
      external nf_def_grp
      integer*4 nf_rename_grp
      external nf_rename_grp
      integer*4 nf_def_var_deflate
      external nf_def_var_deflate
      integer*4 nf_inq_var_deflate
      external nf_inq_var_deflate
      integer*4 nf_def_var_fletcher32
      external nf_def_var_fletcher32
      integer*4 nf_inq_var_fletcher32
      external nf_inq_var_fletcher32
      integer*4 nf_def_var_chunking
      external nf_def_var_chunking
      integer*4 nf_inq_var_chunking
      external nf_inq_var_chunking
      integer*4 nf_def_var_fill
      external nf_def_var_fill
      integer*4 nf_inq_var_fill
      external nf_inq_var_fill
      integer*4 nf_def_var_endian
      external nf_def_var_endian
      integer*4 nf_inq_var_endian
      external nf_inq_var_endian
      integer*4 nf_inq_typeids
      external nf_inq_typeids
      integer*4 nf_inq_typeid
      external nf_inq_typeid
      integer*4 nf_inq_type
      external nf_inq_type
      integer*4 nf_inq_user_type
      external nf_inq_user_type
      integer*4 nf_def_compound
      external nf_def_compound
      integer*4 nf_insert_compound
      external nf_insert_compound
      integer*4 nf_insert_array_compound
      external nf_insert_array_compound
      integer*4 nf_inq_compound
      external nf_inq_compound
      integer*4 nf_inq_compound_name
      external nf_inq_compound_name
      integer*4 nf_inq_compound_size
      external nf_inq_compound_size
      integer*4 nf_inq_compound_nfields
      external nf_inq_compound_nfields
      integer*4 nf_inq_compound_field
      external nf_inq_compound_field
      integer*4 nf_inq_compound_fieldname
      external nf_inq_compound_fieldname
      integer*4 nf_inq_compound_fieldindex
      external nf_inq_compound_fieldindex
      integer*4 nf_inq_compound_fieldoffset
      external nf_inq_compound_fieldoffset
      integer*4 nf_inq_compound_fieldtype
      external nf_inq_compound_fieldtype
      integer*4 nf_inq_compound_fieldndims
      external nf_inq_compound_fieldndims
      integer*4 nf_inq_compound_fielddim_sizes
      external nf_inq_compound_fielddim_sizes
      integer*4 nf_def_vlen
      external nf_def_vlen
      integer*4 nf_inq_vlen
      external nf_inq_vlen
      integer*4 nf_free_vlen
      external nf_free_vlen
      integer*4 nf_def_enum
      external nf_def_enum
      integer*4 nf_insert_enum
      external nf_insert_enum
      integer*4 nf_inq_enum
      external nf_inq_enum
      integer*4 nf_inq_enum_member
      external nf_inq_enum_member
      integer*4 nf_inq_enum_ident
      external nf_inq_enum_ident
      integer*4 nf_def_opaque
      external nf_def_opaque
      integer*4 nf_inq_opaque
      external nf_inq_opaque
      integer*4 nf_put_att
      external nf_put_att
      integer*4 nf_get_att
      external nf_get_att
      integer*4 nf_put_var
      external nf_put_var
      integer*4 nf_put_var1
      external nf_put_var1
      integer*4 nf_put_vara
      external nf_put_vara
      integer*4 nf_put_vars
      external nf_put_vars
      integer*4 nf_get_var
      external nf_get_var
      integer*4 nf_get_var1
      external nf_get_var1
      integer*4 nf_get_vara
      external nf_get_vara
      integer*4 nf_get_vars
      external nf_get_vars
      integer*4 nf_put_var1_int64
      external nf_put_var1_int64
      integer*4 nf_put_vara_int64
      external nf_put_vara_int64
      integer*4 nf_put_vars_int64
      external nf_put_vars_int64
      integer*4 nf_put_varm_int64
      external nf_put_varm_int64
      integer*4 nf_put_var_int64
      external nf_put_var_int64
      integer*4 nf_get_var1_int64
      external nf_get_var1_int64
      integer*4 nf_get_vara_int64
      external nf_get_vara_int64
      integer*4 nf_get_vars_int64
      external nf_get_vars_int64
      integer*4 nf_get_varm_int64
      external nf_get_varm_int64
      integer*4 nf_get_var_int64
      external nf_get_var_int64
      integer*4 nf_get_vlen_element
      external nf_get_vlen_element
      integer*4 nf_put_vlen_element
      external nf_put_vlen_element
      integer*4 nf_set_chunk_cache
      external nf_set_chunk_cache
      integer*4 nf_get_chunk_cache
      external nf_get_chunk_cache
      integer*4 nf_set_var_chunk_cache
      external nf_set_var_chunk_cache
      integer*4 nf_get_var_chunk_cache
      external nf_get_var_chunk_cache
      integer*4 nccre
      integer*4 ncopn
      integer*4 ncddef
      integer*4 ncdid
      integer*4 ncvdef
      integer*4 ncvid
      integer*4 nctlen
      integer*4 ncsfil
      external nccre
      external ncopn
      external ncddef
      external ncdid
      external ncvdef
      external ncvid
      external nctlen
      external ncsfil
      integer*4 ncrdwr
      integer*4 nccreat
      integer*4 ncexcl
      integer*4 ncindef
      integer*4 ncnsync
      integer*4 nchsync
      integer*4 ncndirty
      integer*4 nchdirty
      integer*4 nclink
      integer*4 ncnowrit
      integer*4 ncwrite
      integer*4 ncclob
      integer*4 ncnoclob
      integer*4 ncglobal
      integer*4 ncfill
      integer*4 ncnofill
      integer*4 maxncop
      integer*4 maxncdim
      integer*4 maxncatt
      integer*4 maxncvar
      integer*4 maxncnam
      integer*4 maxvdims
      integer*4 ncnoerr
      integer*4 ncebadid
      integer*4 ncenfile
      integer*4 nceexist
      integer*4 nceinval
      integer*4 nceperm
      integer*4 ncenotin
      integer*4 nceindef
      integer*4 ncecoord
      integer*4 ncemaxds
      integer*4 ncename
      integer*4 ncenoatt
      integer*4 ncemaxat
      integer*4 ncebadty
      integer*4 ncebadd
      integer*4 ncests
      integer*4 nceunlim
      integer*4 ncemaxvs
      integer*4 ncenotvr
      integer*4 nceglob
      integer*4 ncenotnc
      integer*4 ncfoobar
      integer*4 ncsyserr
      integer*4 ncfatal
      integer*4 ncverbos
      integer*4 ncentool
      integer*4 ncbyte
      integer*4 ncchar
      integer*4 ncshort
      integer*4 nclong
      integer*4 ncfloat
      integer*4 ncdouble
      parameter(ncbyte = 1)
      parameter(ncchar = 2)
      parameter(ncshort = 3)
      parameter(nclong = 4)
      parameter(ncfloat = 5)
      parameter(ncdouble = 6)
      parameter(ncrdwr = 1)
      parameter(nccreat = 2)
      parameter(ncexcl = 4)
      parameter(ncindef = 8)
      parameter(ncnsync = 16)
      parameter(nchsync = 32)
      parameter(ncndirty = 64)
      parameter(nchdirty = 128)
      parameter(ncfill = 0)
      parameter(ncnofill = 256)
      parameter(nclink = 32768)
      parameter(ncnowrit = 0)
      parameter(ncwrite = ncrdwr)
      parameter(ncclob = nf_clobber)
      parameter(ncnoclob = nf_noclobber)
      integer*4 ncunlim
      parameter(ncunlim = 0)
      parameter(ncglobal  = 0)
      parameter(maxncop = 64)
      parameter(maxncdim = 1024)
      parameter(maxncatt = 8192)
      parameter(maxncvar = 8192)
      parameter(maxncnam = 256)
      parameter(maxvdims = maxncdim)
      parameter(ncnoerr = nf_noerr)
      parameter(ncebadid = nf_ebadid)
      parameter(ncenfile = -31)
      parameter(nceexist = nf_eexist)
      parameter(nceinval = nf_einval)
      parameter(nceperm = nf_eperm)
      parameter(ncenotin = nf_enotindefine )
      parameter(nceindef = nf_eindefine)
      parameter(ncecoord = nf_einvalcoords)
      parameter(ncemaxds = nf_emaxdims)
      parameter(ncename = nf_enameinuse)
      parameter(ncenoatt = nf_enotatt)
      parameter(ncemaxat = nf_emaxatts)
      parameter(ncebadty = nf_ebadtype)
      parameter(ncebadd = nf_ebaddim)
      parameter(nceunlim = nf_eunlimpos)
      parameter(ncemaxvs = nf_emaxvars)
      parameter(ncenotvr = nf_enotvar)
      parameter(nceglob = nf_eglobal)
      parameter(ncenotnc = nf_enotnc)
      parameter(ncests = nf_ests)
      parameter (ncentool = nf_emaxname)
      parameter(ncfoobar = 32)
      parameter(ncsyserr = -31)
      parameter(ncfatal = 1)
      parameter(ncverbos = 2)
      integer*4 filbyte
      integer*4 filchar
      integer*4 filshort
      integer*4 fillong
      real filfloat
      doubleprecision fildoub
      parameter (filbyte = -127)
      parameter (filchar = 0)
      parameter (filshort = -32767)
      parameter (fillong = -2147483647)
      parameter (filfloat = 9.9692099683868690D+36)
      parameter (fildoub = 9.9692099683868690D+36)
      narg=iargc()
      if (narg .lt. 3) then
        write(stdout,'(/1x,A,1x,A/32x,A/)') 'Usage of partit',
     &    'should be:', 'partit NP_X NP_E ncname1 ... ncnameN'
        stop
      endif
      call getarg(1,string)
      lvar=lenstr(string)
      NP_XI=0
      do i=1,lvar
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_XI=10*NP_XI+j
        else
          write(stdout,'(/8x,2(A,1x),2A/)') 'ERROR: illegal first',
     &     'argument', string(1:lvar), ', must be an integer number.'
          stop
        endif
      enddo
      call getarg(2,string)
      lvar=lenstr(string)
      NP_ETA=0
      do i=1,lvar
        j=ichar(string(i:i))-48
        if (j.ge.0 .and. j.le.9) then
          NP_ETA=10*NP_ETA+j
        else
          write(stdout,'(/8x,2(A,1x),2A/)') 'ERROR: illegal second',
     &      'argument', string(1:lvar), ', must be an integer number.'
          stop
        endif
      enddo
      NNODES=NP_XI*NP_ETA
      if (NNODES.gt.maxnodes) then
        write(stdout,'(/8x,A,1x,A,I3,1x,A/15x,A/)') 'ERROR:',
     &      'requested number of nodes',NNODES,'exceeds limit.',
     &      'Increase parameter maxnodes in file "partit.F".'
        stop
      endif
      write(stdout,'(/4x,2(4x,A,I3)/)') 'NP_XI =',  NP_XI,
     &                                  'NP_ETA =', NP_ETA
      do arg=3,narg
        call getarg(arg,ncname0)
        lstr=lenstr(ncname0)
        ierr=nf_open (ncname0(1:lstr), nf_nowrite, ncid0)
        if (ierr .ne. nf_noerr) then
          write(stdout,'(/8x,A,1x,A,1x,A,A/)') 'ERROR: Cannot',
     &                  'open netCDF file', ncname0(1:lstr),'.'
          goto 97
        endif
        ierr=nf_inq_att (ncid0, nf_global, 'partition', i,j)
        if (ierr .eq. nf_noerr) then
          write(stdout,'(/8x,3A/17x,A/)') 'WARNING: netCDF file ''',
     &                 ncname0(1:lstr), ''' is already partitioned',
     &                 'file. It cannot be partitioned any further.'
          goto 97
        endif
        write(stdout,'(8x,3A)') 'Processing netCDF file ''',
     &                               ncname0(1:lstr), '''...'
        ierr=nf_inq (ncid0, ndims, nvars, ngatts, unlimdimid)
        if (ierr .ne. nf_noerr) then
           write(stdout,'(/8x,A,1x,A/15x,A,1x,A,A/)') 'ERROR:',
     &        'Cannot determine number of dimensions, variables',
     &        'and attributes in netCDF file',ncname0(1:lstr),'.'
          goto 97
        elseif (ndims .gt. maxdims) then
          write(stdout,'(/8x,A,I4,1x,4A/15x,A,1x,A/)')
     &        'ERROR: number of dimensions', ndims,  'in netCDF',
     &        'file ''', ncname0(1:lstr), '''', 'exceeds limit.',
     &        'Increase parameter maxdims in file "partit.F".'
          goto 97
         elseif (nvars .gt. maxvars) then
          write(stdout,'(/8x,A,I4,1x,4A/15x,A,1x,A/)')
     &        'ERROR: number of variables',  nvars,  'in netCDF',
     &        'file ''', ncname0(1:lstr), '''', 'exceeds limit.',
     &        'Increase parameter maxvars in file "partit.F".'
          goto 97
        endif
        tsize=1
        do i=1,ndims
          ierr=nf_inq_dimname (ncid0, i, dimname(i))
          if (ierr .ne. nf_noerr) then
             write(stdout,'(/8x,A,I3/15x,3A/)')
     &           'ERROR: Cannot determine name for dimension ID =',
     &            i,  'in netCDF file ''',  ncname0(1:lstr),  '''.'
             goto 97
          endif
          ierr=nf_inq_dimlen  (ncid0, i, dimsize(i))
          if (ierr .ne. nf_noerr) then
             lvar=lenstr(dimname(i))
             write(stdout,'(/8x,A,A,A/15x,3A/)')
     &            'ERROR: Cannot determine length of dimension ''',
     &             dimname(i)(1:lvar),  '''',  'in netCDF file ''',
     &                                       ncname0(1:lstr), '''.'
             goto 97
          endif
          if (i.eq. unlimdimid) then
            tsize=dimsize(i)
            dimsize(i)=nf_unlimited
          endif
        enddo
        id_xi_rho=0
        id_xi_u=0
        id_eta_rho=0
        id_eta_v=0
        id_xi_psi=0
        id_xi_v=0
        id_eta_psi=0
        id_eta_u=0
        do i=1,ndims
          lvar=lenstr(dimname(i))
          if (lvar.eq.6 .and. dimname(i)(1:lvar).eq.'xi_rho') then
            id_xi_rho=i
            xi_rho=dimsize(i)
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_u') then
            id_xi_u=i
            xi_u=dimsize(i)
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_rho') then
            id_eta_rho=i
            eta_rho=dimsize(i)
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_v') then
            id_eta_v=i
            eta_v=dimsize(i)
          elseif (lvar.eq.6 .and.dimname(i)(1:lvar).eq.'xi_psi') then
            id_xi_psi=i
          elseif (lvar.eq.4 .and. dimname(i)(1:lvar).eq.'xi_v') then
            id_xi_v=i
          elseif (lvar.eq.7.and.dimname(i)(1:lvar).eq.'eta_psi') then
            id_eta_psi=i
          elseif (lvar.eq.5 .and. dimname(i)(1:lvar).eq.'eta_u') then
            id_eta_u=i
          endif
        enddo
        if (id_xi_rho.eq.0  .or. id_xi_u.eq.0 .or.
     &      id_eta_rho.eq.0 .or. id_eta_v.eq.0) then
          write(stdout,'(/8x,2A/15x,3A/)') 'ERROR: not all ',
     &            'partitionable dimensions are found',
     &            'in netCDF file ''', ncname0(1:lstr), '''.'
          goto 97
        endif
        LLm = xi_rho - 2
        MMm = eta_rho - 2
        if (xi_rho .eq. xi_u) then
          XiPeriodic=.true.
          subXI=xi_rho/NP_XI
          if (subXI*NP_XI .ne. xi_rho) then
            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
     &         'ERROR: Cannot partition XI-direction:', 'xi_rho =',
     &          xi_rho,  'is',  'not divisible by NP_XI =',  NP_XI,
     &         'in XI-periodic case.'
            goto 97
          endif
        elseif (xi_rho .eq. xi_u+1) then
          XiPeriodic=.false.
          subXI=(xi_rho-2)/NP_XI
          chunk_size_X=(LLm+NP_XI-1)/NP_XI
          margin_X=(NP_XI*chunk_size_X-LLm)/2
        else
          write(stdout,'(/8x,2A/)') 'ERROR: inconsistent ',
     &                'dimensions ''xi_rho'' and ''xi_u''.'
          goto 97
        endif
        if (eta_rho .eq. eta_v) then
          EtaPeriodic=.true.
          subETA=eta_rho/NP_ETA
          if (subETA*NP_ETA .ne. eta_rho) then
            write(stdout,'(/8x,A,1x,A,I4,1x,A/15x,A,I3,1x,A/)')
     &         'ERROR: Cannot partition ETA-direction:', 'eta_rho =',
     &          eta_rho,  'is',  'not divisible by NP_ETA =', NP_ETA,
     &         'in ETA-periodic case.'
            goto 97
          endif
        elseif (eta_rho .eq. eta_v+1) then
          EtaPeriodic=.false.
          subETA=(eta_rho-2)/NP_ETA
          chunk_size_E=(MMm+NP_ETA-1)/NP_ETA
          margin_E=(NP_ETA*chunk_size_E-MMm)/2
        else
          write(stdout,'(/8x,A,1x,A/)') 'ERROR: inconsistent',
     &                 'dimensions ''eta_rho'' and ''eta_v''.'
          goto 97
        endif
        do node=0,NNODES-1
          lstr=lenstr(ncname0)
          ncname(node)=ncname0
          ierr=0
          call insert_node (ncname(node), lstr, node, NNODES, ierr)
          if (ierr. ne. 0) goto 97
          ierr=nf_create (ncname(node)(1:lstr),nf_64bit_offset,
     &                                                       ncid(node))
          if (ierr .eq. nf_noerr) then
            write(stdout,'(12x,3A)') 'Created partitioned file ''',
     &                                  ncname(node)(1:lstr), '''.'
          else
            write(stdout,'(/8x,A,1x,3A/)') 'ERROR: cannot create',
     &              'netCDF file ''', ncname(node)(1:lstr), '''.'
            goto 97
          endif
          jj=node/NP_XI
          ii=node-jj*NP_XI
      istrmpi=1+ii*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)
      jstrmpi=1+jj*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)
      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1
          do i=1,ndims
            size=dimsize(i)
            if (i .eq. id_xi_rho) then
              size=Lmmpi
              if (.not.XiPeriodic) then
                if (ii.eq.0      ) size=size+1
                if (ii.eq.NP_XI-1) size=size+1
              endif
            elseif (i .eq. id_xi_u) then
              size=Lmmpi
              if (.not.XiPeriodic) then
                if (ii.eq.NP_XI-1) size=size+1
              endif
            elseif (i .eq. id_eta_rho) then
              size=Mmmpi
              if (.not.EtaPeriodic) then
                if (jj.eq.0       ) size=size+1
                if (jj.eq.NP_ETA-1) size=size+1
              endif
            elseif (i .eq. id_eta_v) then
              size=Mmmpi
              if (.not.EtaPeriodic) then
                if (jj.eq.NP_ETA-1) size=size+1
              endif
            endif
            if (i.ne.id_xi_psi   .and.  i.ne.id_xi_v  .and.
     &          i.ne.id_eta_psi  .and.  i.ne.id_eta_u)  then
              lvar=lenstr(dimname(i))
              ierr=nf_def_dim (ncid(node), dimname(i)(1:lvar),
     &                                         size, dimid(i))
              if (ierr .ne. nf_noerr) then
                write(stdout,'(/8x,4A/15x,A,I4,A)') 'ERROR: ',
     &            'Cannot define dimension ''', dimname(i)(1:lvar),
     &            '''.',    'netCDF ettor status =',       ierr, '.'
              endif
            else
              dimid(i)=0
            endif
          enddo
          ibuff(1)=ii
          ibuff(2)=jj
          ibuff(3)=NP_XI
          ibuff(4)=NP_ETA
          ierr=nf_put_att_int (ncid(node), nf_global, 'partition',
     &                                           nf_int, 4, ibuff)
        enddo
        do i=1,ngatts
          ierr=nf_inq_attname (ncid0, nf_global, i, string)
          if (ierr. eq. nf_noerr) then
            lvar=lenstr(string)
            do node=0,NNODES-1
              ierr=nf_copy_att (ncid0, nf_global, string(1:lvar),
     &                                     ncid(node), nf_global)
              if (ierr. ne. nf_noerr) then
                lstr=lenstr(ncname(node))
                write(stdout,'(/8x,4A/15x,3A/)')  'ERROR: Cannot ',
     &            'copy global attribute ''', string(1:lvar), '''',
     &            'into netCDF file ''', ncname(node)(1:lstr),'''.'
                goto 97
              endif
            enddo
          else
            lstr=lenstr(ncname(0))
            write(stdout,'(/8x,2A,I3/15x,3A/)') 'ERROR: Cannot ',
     &         'determine mame of global attribute with ID =', i,
     &         'from netCDF file ''',    ncname0(1:lstr),   '''.'
            goto 97
          endif
        enddo
        do i=1,nvars
          ierr=nf_inq_var (ncid0,   i, varname(i),  vartype(i),
     &                      vardims(i), dimids(1,i),   varatts)
          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_psi) then
              dimids(j,i)=id_xi_u
            elseif (dimids(j,i).eq.id_xi_v) then
              dimids(j,i)=id_xi_rho
            elseif (dimids(j,i).eq.id_eta_psi) then
              dimids(j,i)=id_eta_v
            elseif (dimids(j,i).eq.id_eta_u) then
              dimids(j,i)=id_eta_rho
            endif
          enddo
          series(i)=.false.
          part_switch(i)=.false.
          do j=1,vardims(i)
            if (dimids(j,i).eq.id_xi_rho .or.
     &          dimids(j,i).eq.id_xi_u    .or.
     &          dimids(j,i).eq.id_eta_rho .or.
     &          dimids(j,i).eq.id_eta_v) then
              part_switch(i)=.true.
            elseif (dimids(j,i).eq.unlimdimid) then
              series(i)=.true.
            endif
          enddo
          do j=1,vardims(i)
            do k=1,ndims
              if (dimids(j,i).eq.k) ibuff(j)=dimid(k)
            enddo
          enddo
          lvar=lenstr(varname(i))
          do node=0,NNODES-1
            ierr=nf_def_var (ncid(node), varname(i)(1:lvar),
     &              vartype(i), vardims(i), ibuff, varid(i))
          enddo
          do j=1,varatts
            ierr=nf_inq_attname (ncid0, varid(i), j, string)
            lvar=lenstr(string)
            do node=0,NNODES-1
            ierr=nf_copy_att (ncid0, i, string(1:lvar),
     &                            ncid(node), varid(i))
            enddo
          enddo
        enddo
        do node=0,NNODES-1
          ierr=nf_enddef(ncid(node))
        enddo
        do rec=1,tsize
          if (tsize.gt.1) write(stdout,'(16x,A,I4,A)')
     &                 'Processing record', rec, '...'
          do i=1,nvars
            if (series(i) .or. rec.eq.1) then
              if (.not.part_switch(i) .and. .not.series(i)) then
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_var_text (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_var_int    (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_real) then
                  ierr=nf_get_var_real   (ncid0, i, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_var_double (ncid0, i, buff)
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A/)') 'ERROR: scalar variable',
     &              ' ''', varname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,NNODES-1
                    if (vartype(i) .eq. nf_char) then
                     ierr=nf_put_var_text (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_int) then
                     ierr=nf_put_var_int   (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_real) then
                     ierr=nf_put_var_real  (ncid(node),varid(i),buff)
                    elseif (vartype(i) .eq. nf_double) then
                     ierr=nf_put_var_double(ncid(node),varid(i),buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A/15x,3A,I4,A/)')
     &                    'ERROR: Cannot write scalar variable ''',
     &                     varname(i)(1:lvar), ''' into netCDF file',
     &                    '''',  ncname(node)(1:lstr),
     &                    '''.     netCDF error code =',  ierr,   '.'
                      goto 97
                    endif
                  enddo
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A/)') 'ERROR: Cannot read ',
     &             'scalar variable ''', varname(i)(1:lvar), '''.'
                  goto 97
                endif
              elseif (.not.part_switch(i)) then
                size=1
                do j=1,vardims(i)
                  if (dimids(j,i).eq.unlimdimid) then
                    start(j)=rec
                    count(j)=1
                  else
                    start(j)=1
                    count(j)=dimsize(dimids(j,i))
                  endif
                  size=size*count(j)
                enddo
                if (vartype(i) .eq. nf_char) then
                  size=size*1
                elseif (vartype(i) .eq. nf_int) then
                  size=size*4
                elseif (vartype(i) .eq. nf_real) then
                  size=size*4
                elseif (vartype(i) .eq. nf_double) then
                  size=size*8
                else
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,3A/)') 'ERROR: variable ''',
     &                 varname(i)(1:lvar), ''' has unknown type.'
                  goto 97
                endif
                if (size .gt. 8*max_buff_size) then
                  write(stdout,'(/8x,A,3(/15x,A,I10,1x,A)/)')
     &              'ERROR: unsufficient buffer size in "partit.F":',
     &              'requested:',         size,      'Bytes,',
     &              'available:',   8*max_buff_size, 'Bytes.',
     &              'Increase parameter max_buff_size and recompile.'
                  goto 97
                endif
                if (vartype(i) .eq. nf_char) then
                  ierr=nf_get_vara_text   (ncid0, i, start,
     &                                          count, buff)
                elseif (vartype(i) .eq. nf_int) then
                  ierr=nf_get_vara_int    (ncid0, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_real) then
                  ierr=nf_get_vara_real   (ncid0, i, start,
     &                                         count, buff)
                elseif (vartype(i) .eq. nf_double) then
                  ierr=nf_get_vara_double (ncid0, i, start,
     &                                         count, buff)
                endif
                if (ierr .eq. nf_noerr) then
                  do node=0,NNODES-1
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_real) then
                      ierr=nf_put_vara_real   (ncid(node), varid(i),
     &                                           start, count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (ncid(node), varid(i),
     &                                           start, count, buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot write variable ''',
     &                   varname(i)(1:lvar),''' for time record',rec,
     &                  'into netCDF file ''',  ncname(node)(1:lstr),
     &                  '''. netCDF error code =', ierr, '.'
                      goto 97
                    endif
                  enddo
                else
                  lstr=lenstr(ncname0)
                  lvar=lenstr(varname(i))
                  write(stdout,'(/8x,4A,I3/15x,3A,I4/)') 'ERROR: ',
     &              'Cannot read variable ''',  varname(i)(1:lvar),
     &              ''' for time record',rec,'from netCDF file ''',
     &              ncname0(1:lstr),'''. netCDF error code =',ierr
                  goto 97
                endif
              elseif (part_switch(i)) then
                do node=0,NNODES-1
                  jj=node/NP_XI
                  ii=node-jj*NP_XI
      istrmpi=1+ii*chunk_size_X-margin_X
      iendmpi=istrmpi+chunk_size_X-1
      istrmpi=max(istrmpi,1)
      iendmpi=min(iendmpi,LLm)
      jstrmpi=1+jj*chunk_size_E-margin_E
      jendmpi=jstrmpi+chunk_size_E-1
      jstrmpi=max(jstrmpi,1)
      jendmpi=min(jendmpi,MMm)
      Lmmpi=iendmpi-istrmpi+1
      Mmmpi=jendmpi-jstrmpi+1
                  size=1
                  do j=1,vardims(i)
                    if (dimids(j,i).eq.id_xi_rho) then
                      start(j)=istrmpi
                      count(j)=Lmmpi
                      if (.not.XiPeriodic) then
                        if (ii.gt.0      ) start(j)=start(j)+1
                        if (ii.eq.0      ) count(j)=count(j)+1
                        if (ii.eq.NP_XI-1) count(j)=count(j)+1
                      endif
                      start1(j)=1
                    elseif (dimids(j,i).eq.id_xi_u) then
                      start(j)=istrmpi
                      count(j)=Lmmpi
                      if (.not.XiPeriodic) then
                        if (ii.eq.NP_XI-1) count(j)=count(j)+1
                      endif
                      start1(j)=1
                    elseif (dimids(j,i).eq.id_eta_rho) then
                      start(j)=jstrmpi
                      count(j)=Mmmpi
                      if (.not.EtaPeriodic) then
                        if (jj.gt.0       ) start(j)=start(j)+1
                        if (jj.eq.0       ) count(j)=count(j)+1
                        if (jj.eq.NP_ETA-1) count(j)=count(j)+1
                      endif
                      start1(j)=1
                    elseif (dimids(j,i).eq.id_eta_v) then
                      start(j)=jstrmpi
                      count(j)=Mmmpi
                      if (.not.EtaPeriodic) then
                        if (jj.eq.NP_ETA-1) count(j)=count(j)+1
                      endif
                      start1(j)=1
                    elseif (dimids(j,i).eq.unlimdimid) then
                      start(j)=rec
                      count(j)=1
                      start1(j)=rec
                    else
                      start(j)=1
                      count(j)=dimsize(dimids(j,i))
                      start1(j)=1
                    endif
                    size=size*count(j)
                  enddo
                  if (vartype(i) .eq. nf_char) then
                    size=size*1
                  elseif (vartype(i) .eq. nf_int) then
                    size=size*4
                  elseif (vartype(i) .eq. nf_real) then
                    size=size*4
                  elseif (vartype(i) .eq. nf_double) then
                    size=size*8
                  else
                    lvar=lenstr(varname(i))
                    write(stdout,'(/8x,4A/)') 'ERROR: variable ''',
     &                   varname(i)(1:lvar), ''' has unknown type.'
                    goto 97
                  endif
                  if (size .gt. 8*max_buff_size) then
                    write(stdout,'(/8x,A,3(/15x,A,I10,1x,A)/)')
     &              'ERROR: unsufficient buffer size in "partit.F":',
     &              'requested:',         size,      'Bytes,',
     &              'available:',   8*max_buff_size, 'Bytes.',
     &              'Increase parameter max_buff_size and recompile.'
                    goto 97
                  endif
                  if (vartype(i) .eq. nf_char) then
                    ierr=nf_get_vara_text   (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_int) then
                    ierr=nf_get_vara_int    (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_real) then
                    ierr=nf_get_vara_real   (ncid0, i, start,
     &                                           count, buff)
                  elseif (vartype(i) .eq. nf_double) then
                    ierr=nf_get_vara_double (ncid0, i, start,
     &                                           count, buff)
                  endif
                  if (ierr .eq. nf_noerr) then
                    if (vartype(i) .eq. nf_char) then
                      ierr=nf_put_vara_text   (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_int) then
                      ierr=nf_put_vara_int    (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_real) then
                      ierr=nf_put_vara_real   (ncid(node), varid(i),
     &                                          start1, count, buff)
                    elseif (vartype(i) .eq. nf_double) then
                      ierr=nf_put_vara_double (ncid(node), varid(i),
     &                                          start1, count, buff)
                    endif
                    if (ierr .ne. nf_noerr) then
                      lvar=lenstr(varname(i))
                      lstr=lenstr(ncname(node))
                      write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot write partitioned array ''',
     &                   varname(i)(1:lvar),''' for time record',rec,
     &                  'into netCDF file ''',  ncname(node)(1:lstr),
     &                  '''. netCDF error code =', ierr, '.'
                      goto 97
                    endif
                  else
                    lstr=lenstr(ncname0)
                    lvar=lenstr(varname(i))
                    write(stdout,'(/8x,3A,I3/15x,3A,I4,A/)')
     &                  'ERROR: Cannot read partitioned array ''',
     &                   varname(i)(1:lvar),   ''' for time record',
     &                   rec, 'from netCDF file ''',ncname0(1:lstr),
     &                  '''. netCDF error code =',   ierr,    '.'
                    goto 97
                  endif
                enddo
              endif
            endif
          enddo
        enddo
  97    ierr=nf_close (ncid0)
        do node=0,NNODES-1
         ierr=nf_close (ncid(node))
        enddo
      enddo
      stop
      end
