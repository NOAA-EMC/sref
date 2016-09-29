!-------------------------------------------------------------------------------
program postgp
!$$$  Main program documentation block
!
! Main program: postgp       Transform sigma to pressure grib
!   Prgmmr: Iredell          Org: np23        Date: 1999-10-18
!
! Abstract: Program transforms sigma input to pressure grib output.
!   The output consists of data on a regular lat/lon grid.
!   Geopotential height, wind components, relative humidity,
!   temperature, vertical velocity and absolute vorticity
!   are output on mandatory pressures levels at 50 mb resolution.
!   Constant pressure thickness fields are also output.
!   Constant height above sea level fields are output too.
!   Also output are sundry fields consisting of
!   precipitable water, sea level pressure, two lifted indices,
!   surface orography, temperature, pressure, and pressure tendency,
!   tropopause temperature, pressure, wind components and shear,
!   and maxwind level temperature, pressure and wind components.
!   First nampgb namelist is read to determine output format.
!   Then a sigma (grid or spectral) file is read from unit 11 and
!   the program produces and writes a pressure grib1 file to unit 51.
!   Then a sigma file is read from unit 12 and
!   the program produces and writes a pressure grib1 file to unit 52.
!   The program continues until an empty input file is encountered.
!
! Program history log:
!   92-10-31  Iredell
!   98-02-01  Iredell  new levs 975,950,900,800,750,650,600,550,450,350,
!                      deleted all abs vor but 500 mb,
!                      smoothed abs vor to t72,
!                      included 1000 to 330 rh,
!                      made y2k-compliant,
!                      allow for tracers on the gaussian grid,
!                      reset fhour units for precip after day 10,
!                      corrected gds for gaussian grid option.
!   98-04-08  Iredell  extend relative humidity to 100 mb and
!                      restrict ozone to 100 mb and above.
! 1999-05-08  Iredell  SP version
! 1999-10-18  Iredell  Documented
! 2000-09-11  Iredell  Added freezing level
! 2000-09-27  Iredell  Smoothed prmsl to t80
! 2001-10-02  Iredell  Extrapolated freezing level underground
! 2001-11-05  Iredell  Included potential vorticity levels
! 2001-11-15  Iredell  Allowed potential algorithm computation of PMSL
! 2003-02-01  Iredell  Added level type 126 for isobaric levels above 1 mb
! 2003-02-01  Iredell  Added relative humidity to constant height level and
!                      added constant height level 4267.
! 2003-03-01  Iredell  Included potential temperature levels
! 2003-03-01  Iredell  Added flexibility in fields to output
!                      with an optional filter list 
! 2005-02-02  Iredell  Added flux fields, changed maxwind and pv output
!
! Namelists:
!   nampgb:      parameters determining output format
!     ddsig      character(255) input sigma filename
!                (default: 'postgp.inp.sig')
!     ddflx      character(255) input flux filename
!                (default: 'postgp.inp.flx')
!     ddpgb      character(255) output pressure grib filename
!                (default: 'postgp.out.pgb')
!     idrtc      integer compute grid flag (0 for latlon or 4 for gaussian)
!                (default: idrt)
!     ioc        integer number of compute longitude points (default: io)
!     joc        integer number of compute latitude points (default: jo)
!     idrt       integer output grid flag (0 for latlon or 4 for gaussian)
!                (default: 0 for latlon)
!     io         integer number of output longitude points
!                (default: usually 360)
!     jo         integer number of output latitude points
!                (default: usually 181)
!     lpo        integer pressure level set id
!                (10=mandatory up to 10, 1=every 50 up to 1, 2=every 25 up to 1)
!                (default: usually 1)
!     kpo        integer number of pressure levels (default: 26)
!     po         real (kpo) pressures in mb (default: every 50 mb
!                plus 975 mb plus mandatory levels up to 10 mb)
!     kpto       integer number of pressure thickness layers output (default: 6)
!     kpt        integer number of pressure thickness layers computed
!                (default: 6)
!     pt         real (kpt) pressure thickness in mb (default: 30.)
!     kzz        integer number of constant height levels (default: 8)
!     zz         real (kzz) constant heights in m
!                (default: 305.,457.,610.,914.,1829.,2743.,3658.,4572.)
!     kth        integer number of potential temperature levels (default: 0)
!     th         real (kth) potential temperatures in K (default: )
!     kpv        integer number of potential vorticity levels (default: 8)
!     pv         real (kpv) potential vorticities in PV units
!                (default: 0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0)
!     pvpt       real (kpv) top pressures for PV search (default: 8*50.)
!     pvsb       real (kpv) bottom sigmas for PV search (default: 8*0.8)
!     ncpus      integer number of parallel processes (default: 1)
!     mxbit      integer maximum number of bits to pack data (default: 16)
!     ids        integer (255,255) decimal scaling of packed data
!                (default: set by subprogram idsset)
!     pob        real (255) lowest level pressure in mb to output data
!                as a function of parameter indicator
!                (default: 500 for 5-wave height, 100 for ozone)
!     pot        real (255) highest level pressure in mb to output data
!                as a function of parameter indicator
!                (default: 500 for 5-wave height, 100 for rh or cloud,
!                          100 for omega, 0 otherwise)
!     omin       real (255) minimum value in output data
!                as a function of parameter indicator
!                (default: 0. for humidity and ozone fields, -huge otherwise)
!     omax       real (255) maximum value in output data
!                as a function of parameter indicator
!                (default: 100. for rh, +huge otherwise)
!     moo        integer smoothing truncation for pressure level fields
!                (set to 255 for no smoothing and no abs vorticity)
!                (default: 2*joc/3 if idrt=0, jcap if idrt=4)
!     mow        integer smoothing width for pressure level fields
!                (default: joc/6 if idrt=0, 0 if idrt=4)
!     mooa       integer smoothing truncation for pressure level abs vorticity
!                (default: moo)
!     mowa       integer smoothing truncation for pressure level abs vorticity
!                (default: mow)
!     mooslp     integer smoothing truncation for sea level pressure
!                (default: 80)
!     mowslp     integer smoothing truncation for sea level pressure
!                (default: 20)
!     icen       integer forecast center identifier (default: 7)
!     icen2      integer forecast sub-center identifier
!                (default: from sigma file)
!     igen       integer model generating code (default: from sigma file)
!     ienst      integer ensemble type (default: from sigma file)
!     iensi      integer ensemble identification
!                (default: from sigma file)
!     lmw        integer multiple write flag (default: 1 for yes)
!     nkplist    integer number of entries in the filter list (default: 1)
!     kplist     integer (4,nkplist) filter list of fields to output
!                in iptv,ipu,itl,il1*256+il2 order (-1 means wildcard)
!                (default: kplist(:,1)=-1 and kplist(:,2:)=0)
!
! Input files:
!   unit   11-?  sigma file(s)
!
! Output files:
!   unit   51-?  pressure grib1 file(s)
!
! Modules used:
!   postgp_module  Shared data for postgp
!   sigio_r_module Sigma file I/O
!   funcphys       Physical functions
!
! Subprograms called:
!   baopenr        byte-addressable open for read
!   baopenwt       byte-addressable open for write with truncate
!   errmsg         write error message to stderr
!   gfuncphys      initialize physics functions
!   idsset         set defaults for decimal scaling
!   instrument     collect wall-clock timings
!   makglgds       create global GDS
!   mpabort        abort a distributed program
!   mpfilli4       gather grid decomposition
!   mpfillr4       gather grid decomposition
!   mpstart        initialize a distributed program
!   mpstop         finalize a distributed program
!   mptgen         generate grid decomposition dimensions
!   mptran1l1      transpose grid decompositions
!   mptran1r4      transpose grid decompositions
!   posta          post subprogram to vertically interpolate
!   posto          post subprogram to pack and write output
!   postx          post subprogram to filter horizontally
!   postx1         post subprogram to filter horizontally
!   rtflx          read and transform flux file
!   rtsig          read and transform sigma file
!   sigio_rrhead   read sigma file header
!   sigio_rropen   open sigma file
!   w3tagb         document beginning of program
!   w3tage         document end of program
!
! Attributes:
!   Language: cray fortran
!
!$$$
  use postgp_module
  use sigio_module
  use sigio_r_module
  use funcphys
  use machine,only:kint_mpi
  implicit none
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Input namelist
  integer,parameter:: komax=200,nomax=2000
  character(255) ddsig,ddflx,ddpgb
  integer idrtc,ioc,joc,idrt,io,jo
  integer lpo,kpo,kpto,kpt,kzz,kth,kpv
  integer ncpus,mxbit,moo,mow,mooa,mowa,mooslp,mowslp
  integer icen,icen2,igen,ienst,iensi
  integer lmw
  real po(komax),pt(komax),zz(komax),th(komax),pv(komax),pvpt(komax),pvsb(komax)
  integer ids(255,255)
  real pob(255),pot(255)
  real omin(255),omax(255)
  integer nkplist,kplist(4,nomax)
  namelist/nampgb/ ddsig,ddflx,ddpgb,&
                   idrtc,ioc,joc,idrt,io,jo,lpo,kpo,po,kpto,kpt,pt,kzz,zz,&
                   kth,th,kpv,pv,pvpt,pvsb,&
                   ncpus,mxbit,ids,pob,pot,omin,omax,&
                   moo,mow,mooa,mowa,mooslp,mowslp,&
                   icen,icen2,igen,ienst,iensi,lmw,&
                   nkplist,kplist
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Other local data
  integer iret,iretflx
  integer(sigio_intkind):: irets
  integer igridc,igrido
  integer kgdsc(200),kgdso(200),kenso(5)
  integer ipopt(20)
  real popa(komax),ptpa(komax),pvptpa(komax)
  integer(sigio_intkind),parameter:: lusig=11
  integer,parameter:: luflx=31,lupgb=51
  type(sigio_head):: sighead
  integer jcap,levs,mo3,mct
  integer kpdsflxs(200,nflxs),kpdsflxv(200,nflxv)
  integer kpdsfoxs(200,nfoxs),kpdsfoxv(200,nfoxv)
  real hrflx,fthour
  integer ijo,ijoc,nfpos,nfpov,nfpts,nfptv,nfzzs,nfzzv,nfths,nfthv,nfpvs,nfpvv
  integer(kint_mpi):: comm,rank,size
  integer k1a,k2a,kma,kna,kxa
  integer ns1b,ns2b,nsxb,nsmb,nsnb
  integer nv1b,nv2b,nvxb,nvmb,nvnb
  integer ij1c,ij2c,ijxc,ijmc,ijnc
  integer k1d,k2d,kxd,kmd,knd
  integer nfpos1e,nfpos2e,nfposxe,nfposme,nfposne
  integer nfpov1e,nfpov2e,nfpovxe,nfpovme,nfpovne
  integer nfpts1e,nfpts2e,nfptsxe,nfptsme,nfptsne
  integer nfptv1e,nfptv2e,nfptvxe,nfptvme,nfptvne
  integer nfzzs1e,nfzzs2e,nfzzsxe,nfzzsme,nfzzsne
  integer nfzzv1e,nfzzv2e,nfzzvxe,nfzzvme,nfzzvne
  integer nfths1e,nfths2e,nfthsxe,nfthsme,nfthsne
  integer nfthv1e,nfthv2e,nfthvxe,nfthvme,nfthvne
  integer nfpvs1e,nfpvs2e,nfpvsxe,nfpvsme,nfpvsne
  integer nfpvv1e,nfpvv2e,nfpvvxe,nfpvvme,nfpvvne
  integer nsuns1e,nsuns2e,nsunsxe,nsunsme,nsunsne
  integer nsunv1e,nsunv2e,nsunvxe,nsunvme,nsunvne
  integer nfoxs1e,nfoxs2e,nfoxsxe,nfoxsme,nfoxsne
  integer nfoxv1e,nfoxv2e,nfoxvxe,nfoxvme,nfoxvne
  logical*1 ldum(1)
  real gdum(1)
  integer iyr,imo,idy,ihr
  integer iprev,nfw
  integer i,nk,k,n,ng
  integer kall
  real ttot,tmin,tmax
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Allocatable data
  real,allocatable:: ha(:),pa(:),pxa(:),pya(:)
  real,allocatable:: ta(:,:),txa(:,:),tya(:,:)
  real,allocatable:: ua(:,:),va(:,:),da(:,:),za(:,:),sha(:,:),o3a(:,:),cta(:,:)
  integer,allocatable:: kpdsflxsb(:,:),kpdsflxvb(:,:)
  logical*1,allocatable:: lflxsb(:,:),lflxvb(:,:)
  real,allocatable:: fflxsb(:,:),fflxub(:,:),fflxvb(:,:)
  real,allocatable:: hc(:),pc(:),pxc(:),pyc(:)
  real,allocatable:: tc(:,:),txc(:,:),tyc(:,:)
  real,allocatable:: uc(:,:),vc(:,:),dc(:,:),zc(:,:),shc(:,:),o3c(:,:),ctc(:,:)
  logical*1,allocatable:: lflxsc(:,:),lflxvc(:,:)
  real,allocatable:: fflxsc(:,:),fflxuc(:,:),fflxvc(:,:)
  real,allocatable:: fposc(:,:),fpouc(:,:),fpovc(:,:)
  real,allocatable:: fptsc(:,:),fptuc(:,:),fptvc(:,:)
  real,allocatable:: fzzsc(:,:),fzzuc(:,:),fzzvc(:,:)
  logical*1,allocatable:: lthsc(:,:),lthvc(:,:)
  real,allocatable:: fthsc(:,:),fthuc(:,:),fthvc(:,:)
  logical*1,allocatable:: lpvsc(:,:),lpvvc(:,:)
  real,allocatable:: fpvsc(:,:),fpvuc(:,:),fpvvc(:,:)
  real,allocatable:: fsunsc(:,:),fsunuc(:,:),fsunvc(:,:)
  logical*1,allocatable:: lfoxsc(:,:),lfoxvc(:,:)
  real,allocatable:: ffoxsc(:,:),ffoxuc(:,:),ffoxvc(:,:)
  real,allocatable:: fposd(:,:,:),fpoud(:,:,:),fpovd(:,:,:)
  integer,allocatable:: ibmse(:),iptve(:),ipue(:),ipve(:)
  integer,allocatable:: itle(:),il1e(:),il2e(:)
  integer,allocatable:: ifue(:),ip1e(:),ip2e(:),itre(:)
  character(16) cinstep(30)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Initial setup
  call mpstart(comm,rank,size,iret)
  if(rank.eq.0) then
    call w3tagb('postgp  ',0000,0000,0000,'np23   ')
  endif
  call instrument(30,kall,ttot,tmin,tmax)
  call gfuncphys
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Set defaults and read namelist
  ddsig='postgp.inp.sig'
  ddflx='postgp.inp.flx'
  ddpgb='postgp.out.pgb'
  idrtc=0
  ioc=0
  joc=0
  idrt=0
  io=0
  jo=0
  lpo=0
  kpo=0
  po=0.
  kpto=6
  kpt=6
  pt=(/30.,60.,90.,120.,150.,180.,(0.,k=6+1,komax)/)
  kzz=8
  zz=(/305.,457.,610.,914.,1829.,2743.,3658.,4572.,(0.,k=8+1,komax)/)
  kth=0
  kpv=8
  pv=(/0.5,-0.5,1.0,-1.0,1.5,-1.5,2.0,-2.0,(0.,k=kpv+1,komax)/)
  pvpt=(/(50.,k=1,kpv),(0.,k=kpv+1,komax)/)
  pvsb=(/(0.8,k=1,kpv),(0.,k=kpv+1,komax)/)
  ncpus=1
  mxbit=20
  call idsset(ids)
  pob=1000.
! pob(051)=0.
! pob(153)=0.
  pob(154)=100.
  pob(222)=500.
  pot=0.
  pot(039)=100.
  pot(051)=100.
  pot(052)=100.
  pot(153)=100.
  pot(222)=500.
  omin=-huge(omin)
  omin(51)=0.
  omin(52)=0.
  omin(153)=0.
  omin(154)=0.
  omax=+huge(omax)
  omax(52)=100.
  moo=0
  mow=0
  mooa=0
  mowa=0
  mooslp=80
  mowslp=20
  icen=7
  icen2=0
  igen=0
  ienst=0
  iensi=0
  lmw=1
  nkplist=1
  kplist(:,1)=-1
  kplist(:,2:)=0
  read(*,nampgb,iostat=iret)
! write(*,nampgb)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Open sigma file
  call sigio_rropen(lusig,ddsig,irets)
  if(irets.ne.0) then
    call errmsg('Error opening input sigma file '//ddsig)
    call mpabort(11)
  endif
  if(rank.eq.0)&
  print '("postgp: opened sigma file ",a)',ddsig
  call sigio_rrhead(lusig,sighead,irets)
  if(irets.ne.0) then
    call errmsg('Error reading input sigma file '//ddsig)
    call mpabort(11)
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Set more defaults
  jcap=sighead%jcap
  levs=sighead%levs
  if(icen2.eq.0) icen2=sighead%icen2
  if(igen.eq.0) igen=sighead%igen
  if(ienst.eq.0) ienst=sighead%iens(1)
  if(iensi.eq.0) iensi=sighead%iens(2)
  if(io.eq.0.or.jo.eq.0) then
    idrt=0
    if(jcap.eq.62) then
      io=144
      jo=73
    elseif(jcap.ge.382) then
      io=720
      jo=361
    else
      io=360
      jo=181
    endif
  endif
  if(ioc.eq.0.or.joc.eq.0) then
    idrtc=idrt
    ioc=io
    joc=jo
  endif
  if(moo.eq.0) then
    if(idrtc.eq.0) then
      moo=2*joc/3
      mow=joc/6
    else
      moo=jcap
      mow=0
    endif
  endif
  if(mooa.eq.0) then
    mooa=moo
    mowa=mow
  endif
  if(mooslp.eq.0) then
    mooslp=moo
    mowslp=mow
  endif
  if(mooslp.gt.moo) then
    mooslp=moo
    mowslp=mow
  endif
  if(kpo.eq.0) then
    if(lpo.eq.0) then
      if(levs.lt.64) then
        lpo=1
      else
        lpo=2
      endif
    endif
    if(lpo.eq.1) then
      kpo=31
      po=(/1000.,975.,950.,925.,900.,850.,800.,750.,700.,650.,600.,&
           550.,500.,450.,400.,350.,300.,250.,200.,150.,100.,&
           70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,&
           (0.,k=31+1,komax)/)
    elseif(lpo.eq.2) then
      kpo=47
      po=(/1000.,975.,950.,925.,900.,875.,850.,825.,800.,775.,750.,725.,700.,&
                 675.,650.,625.,600.,575.,550.,525.,500.,475.,450.,425.,400.,&
                 375.,350.,325.,300.,275.,250.,225.,200.,175.,150.,125.,100.,&
                 70.,50.,30.,20.,10.,7.,5.,3.,2.,1.,&
           (0.,k=47+1,komax)/)
    else
      kpo=15
      po=(/1000.,850.,700.,500.,400.,300.,250.,200.,150.,100.,&
           70.,50.,30.,20.,10.,&
           (0.,k=15+1,komax)/)
    endif
  endif
  nfpos=npos*kpo
  nfpov=npov*kpo
  nfpts=npts*kpt
  nfptv=nptv*kpt
  nfzzs=nzzs*kzz
  nfzzv=nzzv*kzz
  nfths=nths*kth
  nfthv=nthv*kth
  nfpvs=npvs*kpv
  nfpvv=npvv*kpv
  call makglgds(idrtc,ioc,joc,igridc,kgdsc)
  call makglgds(idrt,io,jo,igrido,kgdso)
  ijo=io*jo
  ijoc=ioc*joc
  popa=po*1.e+2
  ptpa=pt*1.e+2
  pvptpa=pvpt*1.e+2
  if(rank.eq.0) then
    print '("postgp: sigma file resolution T ",i4,"   L ",i4)',jcap,levs
    if(idrtc.eq.0) then
      print '("postgp: compute grid is equal-spaced ",i4," by ",i4)',ioc,joc
    else
      print '("postgp: compute grid is Gaussian ",i4," by ",i4)',ioc,joc
    endif
    if(idrt.eq.0) then
      print '("postgp: output grid is equal-spaced ",i4," by ",i4)',io,jo
    else
      print '("postgp: output grid is Gaussian ",i4," by ",i4)',io,jo
    endif
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Open flux file
  call baopenr(luflx,ddflx,iretflx)
  if(iretflx.ne.0) then
    if(rank.eq.0)&
    print '("postgp: error opening flux file ",a)',ddflx
  else
    if(rank.eq.0)&
    print '("postgp: opened flux file ",a)',ddflx
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Open postgp file
  if(lmw.ne.0.or.rank.eq.0) then
    call baopenwt(lupgb,ddpgb,iret)
    if(iret.ne.0) then
      call errmsg('Error opening output postgp file '//ddpgb)
      call mpabort(51)
    endif
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Determine grid decompositions
  call mptgen(rank,size,1,1,levs,k1a,k2a,kxa,kma,kna)
  call mptgen(rank,size,1,1,nflxs,ns1b,ns2b,nsxb,nsmb,nsnb)
  call mptgen(rank,size,1,1,nflxv,nv1b,nv2b,nvxb,nvmb,nvnb)
  call mptgen(rank,size,1,1,ijoc,ij1c,ij2c,ijxc,ijmc,ijnc)
  call mptgen(rank,size,1,1,kpo,k1d,k2d,kxd,kmd,knd)
  call mptgen(rank,size,1,1,nfpos,nfpos1e,nfpos2e,nfposxe,nfposme,nfposne)
  call mptgen(rank,size,1,1,nfpov,nfpov1e,nfpov2e,nfpovxe,nfpovme,nfpovne)
  call mptgen(rank,size,1,1,nfpts,nfpts1e,nfpts2e,nfptsxe,nfptsme,nfptsne)
  call mptgen(rank,size,1,1,nfptv,nfptv1e,nfptv2e,nfptvxe,nfptvme,nfptvne)
  call mptgen(rank,size,1,1,nfzzs,nfzzs1e,nfzzs2e,nfzzsxe,nfzzsme,nfzzsne)
  call mptgen(rank,size,1,1,nfzzv,nfzzv1e,nfzzv2e,nfzzvxe,nfzzvme,nfzzvne)
  call mptgen(rank,size,1,1,nfths,nfths1e,nfths2e,nfthsxe,nfthsme,nfthsne)
  call mptgen(rank,size,1,1,nfthv,nfthv1e,nfthv2e,nfthvxe,nfthvme,nfthvne)
  call mptgen(rank,size,1,1,nfpvs,nfpvs1e,nfpvs2e,nfpvsxe,nfpvsme,nfpvsne)
  call mptgen(rank,size,1,1,nfpvv,nfpvv1e,nfpvv2e,nfpvvxe,nfpvvme,nfpvvne)
  call mptgen(rank,size,1,1,nsuns,nsuns1e,nsuns2e,nsunsxe,nsunsme,nsunsne)
  call mptgen(rank,size,1,1,nsunv,nsunv1e,nsunv2e,nsunvxe,nsunvme,nsunvne)
  call mptgen(rank,size,1,1,nfoxs,nfoxs1e,nfoxs2e,nfoxsxe,nfoxsme,nfoxsne)
  call mptgen(rank,size,1,1,nfoxv,nfoxv1e,nfoxv2e,nfoxvxe,nfoxvme,nfoxvne)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Allocate step c input arrays
  allocate(hc(ij1c:ij2c),&
           pc(ij1c:ij2c),&
           pxc(ij1c:ij2c),&
           pyc(ij1c:ij2c),&
           tc(levs,ij1c:ij2c),&
           txc(levs,ij1c:ij2c),&
           tyc(levs,ij1c:ij2c),&
           uc(levs,ij1c:ij2c),&
           vc(levs,ij1c:ij2c),&
           dc(levs,ij1c:ij2c),&
           zc(levs,ij1c:ij2c),&
           shc(levs,ij1c:ij2c),&
           o3c(levs,ij1c:ij2c),&
           ctc(levs,ij1c:ij2c),&
           lflxsc(nflxs,ij1c:ij2c),&
           lflxvc(nflxv,ij1c:ij2c),&
           fflxsc(nflxs,ij1c:ij2c),&
           fflxuc(nflxv,ij1c:ij2c),&
           fflxvc(nflxv,ij1c:ij2c))
  call instrument(1,kall,ttot,tmin,tmax);cinstep(1)='setup'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Step A: read and transform sigma file
  if(rank.eq.0)&
  print *,'Begin Step A.'
  allocate(ha(ijoc),&
           pa(ijoc),&
           pxa(ijoc),&
           pya(ijoc),&
           ta(ijoc,k1a:k2a),&
           txa(ijoc,k1a:k2a),&
           tya(ijoc,k1a:k2a),&
           ua(ijoc,k1a:k2a),&
           va(ijoc,k1a:k2a),&
           da(ijoc,k1a:k2a),&
           za(ijoc,k1a:k2a),&
           sha(ijoc,k1a:k2a),&
           o3a(ijoc,k1a:k2a),&
           cta(ijoc,k1a:k2a))
  call rtsig(lusig,sighead,k1a,k2a,kgdsc,ijoc,1,ha,pa,pxa,pya,&
             ta,txa,tya,ua,va,da,za,sha,o3a,cta,mo3,mct,iret)
  if(iret.ne.0) then
    call errmsg('Error reading input sigma file '//ddsig)
    call mpabort(11)
  endif
  if(mo3.eq.0) o3a=0
  if(mct.eq.0) cta=0
  call instrument(2,kall,ttot,tmin,tmax);cinstep(2)='step a rtsig'
  hc=ha(ij1c:ij2c)
  deallocate(ha)
  pc=pa(ij1c:ij2c)
  deallocate(pa)
  pxc=pxa(ij1c:ij2c)
  deallocate(pxa)
  pyc=pya(ij1c:ij2c)
  deallocate(pya)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,ta,tc)
  deallocate(ta)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,txa,txc)
  deallocate(txa)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,tya,tyc)
  deallocate(tya)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,ua,uc)
  deallocate(ua)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,va,vc)
  deallocate(va)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,da,dc)
  deallocate(da)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,za,zc)
  deallocate(za)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,sha,shc)
  deallocate(sha)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,o3a,o3c)
  deallocate(o3a)
  call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kma,kxa,levs,levs,cta,ctc)
  deallocate(cta)
  call instrument(3,kall,ttot,tmin,tmax);cinstep(3)='step a tran'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Step B: read an transform flux file
  if(rank.eq.0)&
  print *,'Begin Step B.'
  if(iretflx.eq.0) then
    allocate(kpdsflxsb(200,ns1b:ns2b),&
             kpdsflxvb(200,nv1b:nv2b),&
             lflxsb(ijoc,ns1b:ns2b),&
             lflxvb(ijoc,nv1b:nv2b),&
             fflxsb(ijoc,ns1b:ns2b),&
             fflxub(ijoc,nv1b:nv2b),&
             fflxvb(ijoc,nv1b:nv2b))
    call rtflx(luflx,nsxb,nvxb,kgdsc,ijoc,&
               jpdsflxs(1:24,ns1b:ns2b),&
               jpdsflxu(1:24,nv1b:nv2b),jpdsflxv(1:24,nv1b:nv2b),&
               kpdsflxsb,kpdsflxvb,lflxsb,lflxvb,fflxsb,fflxub,fflxvb)
    call instrument(4,kall,ttot,tmin,tmax);cinstep(4)='step b rtflx'
    call mpfilli4(comm,size,200,200,200,nsmb,nsxb,nflxs,kpdsflxsb,kpdsflxs)
    deallocate(kpdsflxsb)
    call mpfilli4(comm,size,200,200,200,nvmb,nvxb,nflxv,kpdsflxvb,kpdsflxv)
    deallocate(kpdsflxvb)
    call mptran1l1(comm,size,ijmc,ijoc,ijxc,ijoc,nsmb,nsxb,nflxs,nflxs,&
                   lflxsb,lflxsc)
    deallocate(lflxsb)
    call mptran1l1(comm,size,ijmc,ijoc,ijxc,ijoc,nvmb,nvxb,nflxv,nflxv,&
                   lflxvb,lflxvc)
    deallocate(lflxvb)
    call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,nsmb,nsxb,nflxs,nflxs,&
                   fflxsb,fflxsc)
    deallocate(fflxsb)
    call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,nvmb,nvxb,nflxv,nflxv,&
                   fflxub,fflxuc)
    deallocate(fflxub)
    call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,nvmb,nvxb,nflxv,nflxv,&
                   fflxvb,fflxvc)
    deallocate(fflxvb)
    hrflx=0.
    if(kpdsflxs(16,iflxs_prate_sfc).eq.3)&
    hrflx=fthour(kpdsflxs(13,iflxs_prate_sfc),&
                 kpdsflxs(15,iflxs_prate_sfc)-kpdsflxs(14,iflxs_prate_sfc))
  else
    lflxsc=.false.
    lflxvc=.false.
    fflxsc=0.
    fflxuc=0.
    fflxvc=0.
    hrflx=0.
  endif
  call instrument(5,kall,ttot,tmin,tmax);cinstep(5)='step b tran'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Step C: vertically interpolate to output levels
  if(rank.eq.0)&
  print *,'Begin Step C.'
  allocate(fposc(nfpos,ij1c:ij2c),&
           fpouc(nfpov,ij1c:ij2c),&
           fpovc(nfpov,ij1c:ij2c),&
           fptsc(nfpts,ij1c:ij2c),&
           fptuc(nfptv,ij1c:ij2c),&
           fptvc(nfptv,ij1c:ij2c),&
           fzzsc(nfzzs,ij1c:ij2c),&
           fzzuc(nfzzv,ij1c:ij2c),&
           fzzvc(nfzzv,ij1c:ij2c),&
           lthsc(nfths,ij1c:ij2c),&
           lthvc(nfthv,ij1c:ij2c),&
           fthsc(nfths,ij1c:ij2c),&
           fthuc(nfthv,ij1c:ij2c),&
           fthvc(nfthv,ij1c:ij2c),&
           lpvsc(nfpvs,ij1c:ij2c),&
           lpvvc(nfpvv,ij1c:ij2c),&
           fpvsc(nfpvs,ij1c:ij2c),&
           fpvuc(nfpvv,ij1c:ij2c),&
           fpvvc(nfpvv,ij1c:ij2c),&
           fsunsc(nsuns,ij1c:ij2c),&
           fsunuc(nsunv,ij1c:ij2c),&
           fsunvc(nsunv,ij1c:ij2c),&
           lfoxsc(nfoxs,ij1c:ij2c),&
           lfoxvc(nfoxv,ij1c:ij2c),&
           ffoxsc(nfoxs,ij1c:ij2c),&
           ffoxuc(nfoxv,ij1c:ij2c),&
           ffoxvc(nfoxv,ij1c:ij2c))
  call posta(ijxc,levs,&
             sighead%idvc,sighead%idsl,sighead%nvcoord,sighead%vcoord,&
             hc,pc,pxc,pyc,tc,txc,tyc,uc,vc,dc,zc,shc,o3c,ctc,&
             hrflx,kpdsflxs,kpdsflxv,lflxsc,lflxvc,fflxsc,fflxuc,fflxvc,&
             kpo,popa,kpt,ptpa,kzz,zz,kth,th,kpv,pv,pvptpa,pvsb,&
             nfpos,nfpov,nfpts,nfptv,nfzzs,nfzzv,&
             nfths,nfthv,nfpvs,nfpvv,&
             kpdsfoxs,kpdsfoxv,lfoxsc,lfoxvc,&
             fposc,fpouc,fpovc,fptsc,fptuc,fptvc,fzzsc,fzzuc,fzzvc,&
             lthsc,lthvc,fthsc,fthuc,fthvc,lpvsc,lpvvc,fpvsc,fpvuc,fpvvc,&
             fsunsc,fsunuc,fsunvc,ffoxsc,ffoxuc,ffoxvc)
  deallocate(hc,pc,pxc,pyc,tc,txc,tyc,uc,vc,dc,zc,shc,o3c,ctc,&
             lflxsc,lflxvc,fflxsc,fflxuc,fflxvc)
  call instrument(6,kall,ttot,tmin,tmax);cinstep(6)='step c'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Step D: filter horizontally
  if(rank.eq.0)&
  print *,'Begin Step D.'
  if(moo.ne.255) then
    allocate(fposd(ijoc,k1d:k2d,npos),&
             fpoud(ijoc,k1d:k2d,npov),&
             fpovd(ijoc,k1d:k2d,npov))
    do n=1,npos
      call mptran1r4(comm,size,kmd,kpo,kxd,nfpos,ijmc,ijxc,ijoc,ijoc,&
                     fposc(1+(n-1)*kpo,ij1c),fposd(1,k1d,n))
    enddo
    do n=1,npov
      call mptran1r4(comm,size,kmd,kpo,kxd,nfpov,ijmc,ijxc,ijoc,ijoc,&
                     fpouc(1+(n-1)*kpo,ij1c),fpoud(1,k1d,n))
      call mptran1r4(comm,size,kmd,kpo,kxd,nfpov,ijmc,ijxc,ijoc,ijoc,&
                     fpovc(1+(n-1)*kpo,ij1c),fpovd(1,k1d,n))
    enddo
    call instrument(8,kall,ttot,tmin,tmax);cinstep(8)='step d tran'
    call postx(idrtc,ioc,joc,kxd,moo,mow,mooa,mowa,fposd,fpoud,fpovd)
    call instrument(7,kall,ttot,tmin,tmax);cinstep(7)='step d postx'
    do n=1,npos
      call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kmd,kxd,kpo,nfpos,&
                     fposd(1,k1d,n),fposc(1+(n-1)*kpo,ij1c))
    enddo
    do n=1,npov
      call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kmd,kxd,kpo,nfpov,&
                     fpoud(1,k1d,n),fpouc(1+(n-1)*kpo,ij1c))
      call mptran1r4(comm,size,ijmc,ijoc,ijxc,ijoc,kmd,kxd,kpo,nfpov,&
                     fpovd(1,k1d,n),fpovc(1+(n-1)*kpo,ij1c))
    enddo
    deallocate(fposd,fpoud,fpovd)
    call instrument(8,kall,ttot,tmin,tmax);cinstep(8)='step d tran'
    if(mooslp.ne.0) then
      allocate(fposd(ijoc,1,1))
      call mpfillr4(comm,size,1,1,1,ijmc,ijxc,ijoc,&
                    fsunsc(isuns_prmsl_msl,ij1c:ij2c),fposd)
      call instrument(8,kall,ttot,tmin,tmax);cinstep(8)='step d tran'
      call postx1(idrtc,ioc,joc,mooslp,mowslp,fposd)
      call instrument(7,kall,ttot,tmin,tmax);cinstep(7)='step d postx'
      fsunsc(isuns_prmsl_msl,ij1c:ij2c)=fposd(ij1c:ij2c,1,1)
      deallocate(fposd)
      call instrument(8,kall,ttot,tmin,tmax);cinstep(8)='step d tran'
    endif
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Step E: pack and write output file
  if(rank.eq.0)&
  print *,'Begin Step E.'
  iyr=sighead%idate(4)
  imo=sighead%idate(2)
  idy=sighead%idate(3)
  ihr=sighead%idate(1)
  iprev=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output pressure level scalar
  allocate(ibmse(nfpos1e:nfpos2e),&
           iptve(nfpos1e:nfpos2e),&
           ipue(nfpos1e:nfpos2e),&
           ipve(nfpos1e:nfpos2e),&
           itle(nfpos1e:nfpos2e),&
           il1e(nfpos1e:nfpos2e),&
           il2e(nfpos1e:nfpos2e),&
           ifue(nfpos1e:nfpos2e),&
           ip1e(nfpos1e:nfpos2e),&
           ip2e(nfpos1e:nfpos2e),&
           itre(nfpos1e:nfpos2e))
  do nk=nfpos1e,nfpos2e
    k=mod(nk-1,kpo)+1
    n=(nk-1)/kpo+1
    ibmse(nk)=0
    iptve(nk)=jpdspos(19,n)
    ipue(nk)=jpdspos(5,n)
    ipve(nk)=-1
    itle(nk)=100
    il1e(nk)=0
    il2e(nk)=nint(po(k))
    if(abs(po(k)-nint(po(k))).gt.0.005.and.po(k).lt.655.) then
      itle(nk)=126
      il1e(nk)=0
      il2e(nk)=nint(po(k)*1.e2)
    endif
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
    if(n.eq.ipos_o3mr_prs.and.mo3.eq.0) ipue(nk)=0
    if(n.eq.ipos_clwmr_prs.and.mct.eq.0) ipue(nk)=0
    if(po(k).gt.pob(ipue(nk)).or.po(k).lt.pot(ipue(nk))) ipue(nk)=0
  enddo
  call kpfilt(nkplist,kplist,nfposxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfpos,nfpos1e,nfpos2e,nfposme,nfposxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fposc,gdum,0,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fposc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," pressure level scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(9,kall,ttot,tmin,tmax);cinstep(9)='step e pos'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output pressure level vector
  allocate(ibmse(nfpov1e:nfpov2e),&
           iptve(nfpov1e:nfpov2e),&
           ipue(nfpov1e:nfpov2e),&
           ipve(nfpov1e:nfpov2e),&
           itle(nfpov1e:nfpov2e),&
           il1e(nfpov1e:nfpov2e),&
           il2e(nfpov1e:nfpov2e),&
           ifue(nfpov1e:nfpov2e),&
           ip1e(nfpov1e:nfpov2e),&
           ip2e(nfpov1e:nfpov2e),&
           itre(nfpov1e:nfpov2e))
  do nk=nfpov1e,nfpov2e
    k=mod(nk-1,kpo)+1
    n=(nk-1)/kpo+1
    ibmse(nk)=0
    iptve(nk)=jpdspov(19,n)
    ipue(nk)=jpdspou(5,n)
    ipve(nk)=jpdspov(5,n)
    itle(nk)=100
    il1e(nk)=0
    il2e(nk)=nint(po(k))
    if(abs(po(k)-nint(po(k))).gt.0.005.and.po(k).lt.655.) then
      itle(nk)=126
      il1e(nk)=0
      il2e(nk)=nint(po(k)*1.e2)
    endif
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
    if(po(k).gt.pob(ipue(nk)).or.po(k).lt.pot(ipue(nk))) ipue(nk)=0
    if(po(k).gt.pob(ipve(nk)).or.po(k).lt.pot(ipve(nk))) ipve(nk)=0
  enddo
  call kpfilt(nkplist,kplist,nfpovxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfpovxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfpov,nfpov1e,nfpov2e,nfpovme,nfpovxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fpouc,fpovc,0,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fpouc,fpovc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," pressure level vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(10,kall,ttot,tmin,tmax);cinstep(10)='step e pov'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output pressure layer scalar
  allocate(ibmse(nfpts1e:nfpts2e),&
           iptve(nfpts1e:nfpts2e),&
           ipue(nfpts1e:nfpts2e),&
           ipve(nfpts1e:nfpts2e),&
           itle(nfpts1e:nfpts2e),&
           il1e(nfpts1e:nfpts2e),&
           il2e(nfpts1e:nfpts2e),&
           ifue(nfpts1e:nfpts2e),&
           ip1e(nfpts1e:nfpts2e),&
           ip2e(nfpts1e:nfpts2e),&
           itre(nfpts1e:nfpts2e))
  do nk=nfpts1e,nfpts2e
    k=mod(nk-1,kpt)+1
    n=(nk-1)/kpt+1
    ibmse(nk)=0
    iptve(nk)=jpdspts(19,n)
    ipue(nk)=jpdspts(5,n)
    ipve(nk)=-1
    itle(nk)=116
    il1e(nk)=nint(pt(k))
    il2e(nk)=0
    if(k.gt.1) il2e(nk)=nint(pt(k-1))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
    if(k.gt.kpto) ipue(nk)=0
  enddo
  call kpfilt(nkplist,kplist,nfptsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfpts,nfpts1e,nfpts2e,nfptsme,nfptsxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fptsc,gdum,0,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fptsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," pressure layer scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(11,kall,ttot,tmin,tmax);cinstep(11)='step e pts'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output pressure layer vector
  allocate(ibmse(nfptv1e:nfptv2e),&
           iptve(nfptv1e:nfptv2e),&
           ipue(nfptv1e:nfptv2e),&
           ipve(nfptv1e:nfptv2e),&
           itle(nfptv1e:nfptv2e),&
           il1e(nfptv1e:nfptv2e),&
           il2e(nfptv1e:nfptv2e),&
           ifue(nfptv1e:nfptv2e),&
           ip1e(nfptv1e:nfptv2e),&
           ip2e(nfptv1e:nfptv2e),&
           itre(nfptv1e:nfptv2e))
  do nk=nfptv1e,nfptv2e
    k=mod(nk-1,kpt)+1
    n=(nk-1)/kpt+1
    ibmse(nk)=0
    iptve(nk)=jpdsptv(19,n)
    ipue(nk)=jpdsptu(5,n)
    ipve(nk)=jpdsptv(5,n)
    itle(nk)=116
    il1e(nk)=nint(pt(k))
    il2e(nk)=0
    if(k.gt.1) il2e(nk)=nint(pt(k-1))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
    if(k.gt.kpto) ipue(nk)=0
    if(k.gt.kpto) ipve(nk)=0
  enddo
  call kpfilt(nkplist,kplist,nfptvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfptvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfptv,nfptv1e,nfptv2e,nfptvme,nfptvxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fptuc,fptvc,0,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fptuc,fptvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," pressure layer vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(12,kall,ttot,tmin,tmax);cinstep(12)='step e ptv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output height above sea level scalar
  allocate(ibmse(nfzzs1e:nfzzs2e),&
           iptve(nfzzs1e:nfzzs2e),&
           ipue(nfzzs1e:nfzzs2e),&
           ipve(nfzzs1e:nfzzs2e),&
           itle(nfzzs1e:nfzzs2e),&
           il1e(nfzzs1e:nfzzs2e),&
           il2e(nfzzs1e:nfzzs2e),&
           ifue(nfzzs1e:nfzzs2e),&
           ip1e(nfzzs1e:nfzzs2e),&
           ip2e(nfzzs1e:nfzzs2e),&
           itre(nfzzs1e:nfzzs2e))
  do nk=nfzzs1e,nfzzs2e
    k=mod(nk-1,kzz)+1
    n=(nk-1)/kzz+1
    ibmse(nk)=0
    iptve(nk)=jpdszzs(19,n)
    ipue(nk)=jpdszzs(5,n)
    ipve(nk)=-1
    itle(nk)=103
    il1e(nk)=0
    il2e(nk)=nint(zz(k))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
     if(ipue(nk).eq.52) ipue(nk)=0  ! skip rh for now
  enddo
  call kpfilt(nkplist,kplist,nfzzsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfzzs,nfzzs1e,nfzzs2e,nfzzsme,nfzzsxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fzzsc,gdum,0,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fzzsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," height level scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(13,kall,ttot,tmin,tmax);cinstep(13)='step e zzs'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output height above sea level vector
  allocate(ibmse(nfzzv1e:nfzzv2e),&
           iptve(nfzzv1e:nfzzv2e),&
           ipue(nfzzv1e:nfzzv2e),&
           ipve(nfzzv1e:nfzzv2e),&
           itle(nfzzv1e:nfzzv2e),&
           il1e(nfzzv1e:nfzzv2e),&
           il2e(nfzzv1e:nfzzv2e),&
           ifue(nfzzv1e:nfzzv2e),&
           ip1e(nfzzv1e:nfzzv2e),&
           ip2e(nfzzv1e:nfzzv2e),&
           itre(nfzzv1e:nfzzv2e))
  do nk=nfzzv1e,nfzzv2e
    k=mod(nk-1,kzz)+1
    n=(nk-1)/kzz+1
    ibmse(nk)=0
    iptve(nk)=jpdszzv(19,n)
    ipue(nk)=jpdszzu(5,n)
    ipve(nk)=jpdszzv(5,n)
    itle(nk)=103
    il1e(nk)=0
    il2e(nk)=nint(zz(k))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
  enddo
  call kpfilt(nkplist,kplist,nfzzvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfzzvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfzzv,nfzzv1e,nfzzv2e,nfzzvme,nfzzvxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fzzuc,fzzvc,0,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fzzuc,fzzvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," height level vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(14,kall,ttot,tmin,tmax);cinstep(14)='step e zzv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output potential temperature level scalar
  allocate(ibmse(nfths1e:nfths2e),&
           iptve(nfths1e:nfths2e),&
           ipue(nfths1e:nfths2e),&
           ipve(nfths1e:nfths2e),&
           itle(nfths1e:nfths2e),&
           il1e(nfths1e:nfths2e),&
           il2e(nfths1e:nfths2e),&
           ifue(nfths1e:nfths2e),&
           ip1e(nfths1e:nfths2e),&
           ip2e(nfths1e:nfths2e),&
           itre(nfths1e:nfths2e))
  do nk=nfths1e,nfths2e
    k=mod(nk-1,kth)+1
    n=(nk-1)/kth+1
    ibmse(nk)=1
    iptve(nk)=jpdsths(19,n)
    ipue(nk)=jpdsths(5,n)
    ipve(nk)=-1
    itle(nk)=113
    il1e(nk)=0
    il2e(nk)=nint(th(k))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
  enddo
  call kpfilt(nkplist,kplist,nfthsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfths,nfths1e,nfths2e,nfthsme,nfthsxe,&
             ijmc,ijxc,ijoc,ijo,lthsc,fthsc,gdum,1,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(lthsc,fthsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," isentropic level scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(15,kall,ttot,tmin,tmax);cinstep(13)='step e ths'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output potential temperature level vector
  allocate(ibmse(nfthv1e:nfthv2e),&
           iptve(nfthv1e:nfthv2e),&
           ipue(nfthv1e:nfthv2e),&
           ipve(nfthv1e:nfthv2e),&
           itle(nfthv1e:nfthv2e),&
           il1e(nfthv1e:nfthv2e),&
           il2e(nfthv1e:nfthv2e),&
           ifue(nfthv1e:nfthv2e),&
           ip1e(nfthv1e:nfthv2e),&
           ip2e(nfthv1e:nfthv2e),&
           itre(nfthv1e:nfthv2e))
  do nk=nfthv1e,nfthv2e
    k=mod(nk-1,kth)+1
    n=(nk-1)/kth+1
    ibmse(nk)=1
    iptve(nk)=jpdsthv(19,n)
    ipue(nk)=jpdsthu(5,n)
    ipve(nk)=jpdsthv(5,n)
    itle(nk)=113
    il1e(nk)=0
    il2e(nk)=nint(th(k))
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
  enddo
  call kpfilt(nkplist,kplist,nfthvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfthvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfthv,nfthv1e,nfthv2e,nfthvme,nfthvxe,&
             ijmc,ijxc,ijoc,ijo,lthvc,fthuc,fthvc,1,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(lthvc,fthuc,fthvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," isentropic level vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(16,kall,ttot,tmin,tmax);cinstep(14)='step e thv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output potential vorticity level scalar
  allocate(ibmse(nfpvs1e:nfpvs2e),&
           iptve(nfpvs1e:nfpvs2e),&
           ipue(nfpvs1e:nfpvs2e),&
           ipve(nfpvs1e:nfpvs2e),&
           itle(nfpvs1e:nfpvs2e),&
           il1e(nfpvs1e:nfpvs2e),&
           il2e(nfpvs1e:nfpvs2e),&
           ifue(nfpvs1e:nfpvs2e),&
           ip1e(nfpvs1e:nfpvs2e),&
           ip2e(nfpvs1e:nfpvs2e),&
           itre(nfpvs1e:nfpvs2e))
  do nk=nfpvs1e,nfpvs2e
    k=mod(nk-1,kpv)+1
    n=(nk-1)/kpv+1
    ibmse(nk)=1
    iptve(nk)=jpdspvs(19,n)
    ipue(nk)=jpdspvs(5,n)
    ipve(nk)=-1
    itle(nk)=117
    il1e(nk)=0
    il2e(nk)=nint(pv(k)*1.e3)
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
  enddo
  call kpfilt(nkplist,kplist,nfpvsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfpvs,nfpvs1e,nfpvs2e,nfpvsme,nfpvsxe,&
             ijmc,ijxc,ijoc,ijo,lpvsc,fpvsc,gdum,1,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(lpvsc,fpvsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," PV level scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(15,kall,ttot,tmin,tmax);cinstep(13)='step e pvs'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output potential vorticity level vector
  allocate(ibmse(nfpvv1e:nfpvv2e),&
           iptve(nfpvv1e:nfpvv2e),&
           ipue(nfpvv1e:nfpvv2e),&
           ipve(nfpvv1e:nfpvv2e),&
           itle(nfpvv1e:nfpvv2e),&
           il1e(nfpvv1e:nfpvv2e),&
           il2e(nfpvv1e:nfpvv2e),&
           ifue(nfpvv1e:nfpvv2e),&
           ip1e(nfpvv1e:nfpvv2e),&
           ip2e(nfpvv1e:nfpvv2e),&
           itre(nfpvv1e:nfpvv2e))
  do nk=nfpvv1e,nfpvv2e
    k=mod(nk-1,kpv)+1
    n=(nk-1)/kpv+1
    ibmse(nk)=1
    iptve(nk)=jpdspvv(19,n)
    ipue(nk)=jpdspvu(5,n)
    ipve(nk)=jpdspvv(5,n)
    itle(nk)=117
    il1e(nk)=0
    il2e(nk)=nint(pv(k)*1.e3)
    ifue(nk)=1
    ip1e(nk)=nint(sighead%fhour)
    ip2e(nk)=0
    itre(nk)=10
  enddo
  call kpfilt(nkplist,kplist,nfpvvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfpvvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfpvv,nfpvv1e,nfpvv2e,nfpvvme,nfpvvxe,&
             ijmc,ijxc,ijoc,ijo,lpvvc,fpvuc,fpvvc,1,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(lpvvc,fpvuc,fpvvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," PV level vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(16,kall,ttot,tmin,tmax);cinstep(14)='step e pvv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output sundry scalar
  allocate(ibmse(nsuns1e:nsuns2e),&
           iptve(nsuns1e:nsuns2e),&
           ipue(nsuns1e:nsuns2e),&
           ipve(nsuns1e:nsuns2e),&
           itle(nsuns1e:nsuns2e),&
           il1e(nsuns1e:nsuns2e),&
           il2e(nsuns1e:nsuns2e),&
           ifue(nsuns1e:nsuns2e),&
           ip1e(nsuns1e:nsuns2e),&
           ip2e(nsuns1e:nsuns2e),&
           itre(nsuns1e:nsuns2e))
  do n=nsuns1e,nsuns2e
    ibmse(n)=0
    iptve(n)=jpdssuns(19,n)
    ipue(n)=jpdssuns(5,n)
    ipve(n)=-1
    itle(n)=jpdssuns(6,n)
    il1e(n)=jpdssuns(7,n)/256
    il2e(n)=mod(jpdssuns(7,n),256)
    ifue(n)=1
    ip1e(n)=nint(sighead%fhour)
    ip2e(n)=0
    itre(n)=10
    if(n.eq.isuns_tozne_clm.and.mo3.eq.0) ipue(n)=0
    if(n.eq.isuns_cwat_clm.and.mct.eq.0) ipue(n)=0
  enddo
  call kpfilt(nkplist,kplist,nsunsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nsuns,nsuns1e,nsuns2e,nsunsme,nsunsxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fsunsc,gdum,0,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fsunsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," sundry scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(17,kall,ttot,tmin,tmax);cinstep(15)='step e suns'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output sundry vector
  allocate(ibmse(nsunv1e:nsunv2e),&
           iptve(nsunv1e:nsunv2e),&
           ipue(nsunv1e:nsunv2e),&
           ipve(nsunv1e:nsunv2e),&
           itle(nsunv1e:nsunv2e),&
           il1e(nsunv1e:nsunv2e),&
           il2e(nsunv1e:nsunv2e),&
           ifue(nsunv1e:nsunv2e),&
           ip1e(nsunv1e:nsunv2e),&
           ip2e(nsunv1e:nsunv2e),&
           itre(nsunv1e:nsunv2e))
  do n=nsunv1e,nsunv2e
    ibmse(n)=0
    iptve(n)=jpdssunv(19,n)
    ipue(n)=jpdssunu(5,n)
    ipve(n)=jpdssunv(5,n)
    itle(n)=jpdssunv(6,n)
    il1e(n)=jpdssunv(7,n)/256
    il2e(n)=mod(jpdssunv(7,n),256)
    ifue(n)=1
    ip1e(n)=nint(sighead%fhour)
    ip2e(n)=0
    itre(n)=10
  enddo
  call kpfilt(nkplist,kplist,nsunvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nsunvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nsunv,nsunv1e,nsunv2e,nsunvme,nsunvxe,&
             ijmc,ijxc,ijoc,ijo,ldum,fsunuc,fsunvc,0,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(fsunuc,fsunvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," sundry vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(18,kall,ttot,tmin,tmax);cinstep(16)='step e sunv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output flux scalar
  allocate(ibmse(nfoxs1e:nfoxs2e),&
           iptve(nfoxs1e:nfoxs2e),&
           ipue(nfoxs1e:nfoxs2e),&
           ipve(nfoxs1e:nfoxs2e),&
           itle(nfoxs1e:nfoxs2e),&
           il1e(nfoxs1e:nfoxs2e),&
           il2e(nfoxs1e:nfoxs2e),&
           ifue(nfoxs1e:nfoxs2e),&
           ip1e(nfoxs1e:nfoxs2e),&
           ip2e(nfoxs1e:nfoxs2e),&
           itre(nfoxs1e:nfoxs2e))
  do n=nfoxs1e,nfoxs2e
    ibmse(n)=1
    iptve(n)=jpdsfoxs(19,n)
    ipue(n)=jpdsfoxs(5,n)
    ipve(n)=-1
    itle(n)=jpdsfoxs(6,n)
    il1e(n)=jpdsfoxs(7,n)/256
    il2e(n)=mod(jpdsfoxs(7,n),256)
    ifue(n)=kpdsfoxs(13,n)
    ip1e(n)=kpdsfoxs(14,n)
    ip2e(n)=kpdsfoxs(15,n)
    itre(n)=kpdsfoxs(16,n)
    if(any(kpdsfoxs(13:16,n).lt.0)) ipue(n)=0
  enddo
  call kpfilt(nkplist,kplist,nfoxsxe,iptve,ipue,itle,il1e,il2e)
  call posto(comm,size,nfoxs,nfoxs1e,nfoxs2e,nfoxsme,nfoxsxe,&
             ijmc,ijxc,ijoc,ijo,lfoxsc,ffoxsc,gdum,1,0,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(ffoxsc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," flux scalar",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(19,kall,ttot,tmin,tmax);cinstep(17)='step e foxs'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Output flux vector
  allocate(ibmse(nfoxv1e:nfoxv2e),&
           iptve(nfoxv1e:nfoxv2e),&
           ipue(nfoxv1e:nfoxv2e),&
           ipve(nfoxv1e:nfoxv2e),&
           itle(nfoxv1e:nfoxv2e),&
           il1e(nfoxv1e:nfoxv2e),&
           il2e(nfoxv1e:nfoxv2e),&
           ifue(nfoxv1e:nfoxv2e),&
           ip1e(nfoxv1e:nfoxv2e),&
           ip2e(nfoxv1e:nfoxv2e),&
           itre(nfoxv1e:nfoxv2e))
  do n=nfoxv1e,nfoxv2e
    ibmse(n)=1
    iptve(n)=jpdsfoxv(19,n)
    ipue(n)=jpdsfoxu(5,n)
    ipve(n)=jpdsfoxv(5,n)
    itle(n)=jpdsfoxv(6,n)
    il1e(n)=jpdsfoxv(7,n)/256
    il2e(n)=mod(jpdsfoxv(7,n),256)
    ifue(n)=kpdsfoxv(13,n)
    ip1e(n)=kpdsfoxv(14,n)
    ip2e(n)=kpdsfoxv(15,n)
    itre(n)=kpdsfoxv(16,n)
    if(any(kpdsfoxv(13:16,n).lt.0)) ipue(n)=0
    if(any(kpdsfoxv(13:16,n).lt.0)) ipve(n)=0
  enddo
  call kpfilt(nkplist,kplist,nfoxvxe,iptve,ipue,itle,il1e,il2e)
  call kpfilt(nkplist,kplist,nfoxvxe,iptve,ipve,itle,il1e,il2e)
  call posto(comm,size,nfoxv,nfoxv1e,nfoxv2e,nfoxvme,nfoxvxe,&
             ijmc,ijxc,ijoc,ijo,lfoxvc,ffoxuc,ffoxvc,1,1,lmw,&
             ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre,&
             idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
             icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
             ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
  deallocate(ffoxuc,ffoxvc)
  deallocate(ibmse,iptve,ipue,ipve,itle,il1e,il2e,ifue,ip1e,ip2e,itre)
  if(iret.ne.0) then
    call errmsg('Error packing or writing to output file '//ddpgb)
    call mpabort(51)
  endif
  if(rank.eq.0)&
  print '(i6," flux vector",t30,"fields written.",&
         &i15," total bytes written.")',nfw,iprev
  call instrument(20,kall,ttot,tmin,tmax);cinstep(18)='step e foxv'
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Finish up
  call mpstop(iret)
  call instrument(21,kall,ttot,tmin,tmax);cinstep(21)='mpstop'
  if(rank.eq.0) then
    do n=1,30
      call instrument(-n,kall,ttot,tmin,tmax)
      if(kall.gt.0)&
       print '("TIME SPENT IN STEP ",i2,x,a,g16.6)',n,cinstep(n),ttot
    enddo
    call w3tage('postgp  ')
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end program
!-------------------------------------------------------------------------------
subroutine rtsig(lusig,head,k1,k2,kgds,ijo,nct,&
                 h,p,px,py,t,tx,ty,u,v,d,z,sh,o3,ct,mo3,mct,iret)
!$$$  Subprogram documentation block
!
! Subprogram: rtsig      Read and transform sigma file
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram reads a sigma file and transforms
!   the fields to a designated global grid.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call rtsig(lusig,head,k1,k2,kgds,ijo,nct,&
!                    h,p,px,py,t,tx,ty,u,v,d,z,sh,o3,ct,mo3,mct,iret)
!   Input argument list:
!     lusig    integer(sigio_intkind) sigma file unit number
!     head     type(sigio_head) sigma file header
!     k1       integer first model level to return
!     k2       integer last model level to return
!     kgds     integer (200) GDS to which to transform
!     ijo      integer dimension of output fields
!     nct      integer number of output cloud types
!   Output argument list:
!     h        real (ijo) surface orography (m)
!     p        real (ijo) surface pressure (Pa)
!     px       real (ijo) log surface pressure x-gradient (1/m)
!     py       real (ijo) log surface pressure y-gradient (1/m)
!     t        real (ijo,k1:k2) temperature (K)
!     tx       real (ijo,k1:k2) virtual temperature x-gradient (K/m)
!     ty       real (ijo,k1:k2) virtual temperature y-gradient (K/m)
!     u        real (ijo,k1:k2) x-component wind (m/s)
!     v        real (ijo,k1:k2) y-component wind (m/s)
!     d        real (ijo,k1:k2) wind divergence (1/s)
!     z        real (ijo,k1:k2) absolute vorticity (1/s)
!     sh       real (ijo,k1:k2) specific humidity (kg/kg)
!     o3       real (ijo,k1:k2) specific ozone (kg/kg)
!     ct       real (ijo,k1:k2,nct) specific cloud water (kg/kg)
!     mo3      integer number of ozone fields returned
!     mct      integer number of cloud water fields returned
!     iret     integer return code
!
! Modules used:
!   sigio_r_module sigma file I/O
!
! Subprograms called:
!   sigio_aldatm   allocate sigma upper air data
!   sigio_aldats   allocate sigma surface data
!   sigio_axdatm   deallocate sigma upper air data
!   sigio_axdats   deallocate sigma surface data
!   sigio_rrdatm   read sigma upper air data
!   sigio_rrdats   read sigma surface data
!   sptez          scalar spectral transform
!   sptezd         gradient spectral transform
!   sptezm         multiple scalar spectral transform
!   sptezmv        multiple vector spectral transform
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use sigio_module
  use sigio_r_module
  use physcons
  implicit none
  integer(sigio_intkind),intent(in):: lusig
  type(sigio_head),intent(in):: head
  integer,intent(in):: k1,k2,kgds(200),ijo,nct
  real,dimension(ijo),intent(out):: h,p,px,py
  real,dimension(ijo,k1:k2),intent(out):: t,tx,ty,u,v,d,z,sh,o3
  real,dimension(ijo,k1:k2,nct),intent(out):: ct
  integer,intent(out):: mo3,mct,iret
  integer idrt,io,jo
  integer(sigio_intkind):: irets
  type(sigio_dats):: dats
  type(sigio_datm):: datm
  integer io3,ict,jct
  real pmean,tmean(k1:k2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Determine output grid
  idrt=kgds(1)
  if(kgds(1).eq.0.and.kgds(4).lt.90000) idrt=256
  io=kgds(2)
  jo=kgds(3)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Read and transform surface fields
  iret=1
  call sigio_aldats(head,dats,irets)
  call sigio_rrdats(lusig,head,dats,irets)
  if(irets.ne.0) return
  call sptez(0,head%jcap,idrt,io,jo,dats%hs,h,1)
  call sptez(0,head%jcap,idrt,io,jo,dats%ps,p,1)
  p=1.e3*exp(p)
  call sptezd(0,head%jcap,idrt,io,jo,dats%ps,pmean,px,py,1)
  call sigio_axdats(dats,irets)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Determine types of tracer fields
  mo3=0
  mct=0
  if(head%idvt.eq.0) then
    io3=2
    ict=3
  else
    io3=mod(head%idvt,10)+1
    ict=mod(head%idvt/10,10)+1
  endif
  if(io3.gt.1.and.io3.le.head%ntrac) mo3=1
  if(ict.gt.1.and.ict.le.head%ntrac) mct=min(nct,head%ntrac-ict+1,head%ncldt)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Read and transform fields on levels k1 through k2
  iret=2
  if(k2.ge.k1) then
    call sigio_aldatm(head,k1,k2,datm,irets)
    call sigio_rrdatm(lusig,head,datm,irets)
    if(irets.ne.0) return
    call sptezm(0,head%jcap,idrt,io,jo,k2-k1+1,datm%t,t,1)
    call sptezmd(0,head%jcap,idrt,io,jo,k2-k1+1,datm%t,tmean,tx,ty,1)
    call sptezmv(0,head%jcap,idrt,io,jo,k2-k1+1,datm%d,datm%z,u,v,1)
    call sptezm(0,head%jcap,idrt,io,jo,k2-k1+1,datm%d,d,1)
    datm%z(3,:)=datm%z(3,:)+2*con_omega/sqrt(1.5)
    call sptezm(0,head%jcap,idrt,io,jo,k2-k1+1,datm%z,z,1)
    call sptezm(0,head%jcap,idrt,io,jo,k2-k1+1,datm%q,sh,1)
    if(mo3.gt.0) then
      call sptezm(0,head%jcap,idrt,io,jo,mo3*(k2-k1+1),datm%q(1,k1,io3),o3,1)
    endif
    if(mct.gt.0) then
      call sptezm(0,head%jcap,idrt,io,jo,mct*(k2-k1+1),datm%q(1,k1,ict),ct,1)
    endif
    t=t/(1+con_fvirt*sh)
    call sigio_axdatm(datm,irets)
  endif
  iret=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine rtflx(luflx,ns,nv,kgdsc,ijoc,jpdss,jpdsu,jpdsv,&
                 kpdss,kpdsv,ls,lv,fs,fu,fv)
!$$$  Subprogram documentation block
!
! Subprogram: rtflx      Read and transform flux file
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram reads a flux file and transforms
!   the fields to a designated global grid.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call rtflx(luflx,ns,nv,kgdsc,ijoc,jpdss,jpdsu,jpdsv,&
!                    kpdss,kpdsv,ls,lv,fs,fu,fv)
!   Input argument list:
!     luflx    integer flux file unit number
!     ns       integer number of scalar fields to return
!     nv       integer number of vector fields to return
!     kgdsc    integer (200) GDS to which to transform
!     ijoc     integer dimension of output fields
!     jpdss    integer (1:24,ns) PDS type and level for scalar fields
!     jpdsu    integer (1:24,nv) PDS type and level for vector x-component fields
!     jpdsv    integer (1:24,nv) PDS type and level for vector y-component fields
!   Output argument list:
!     kpdss    integer (200,ns) PDS from scalar fields
!     kpdsv    integer (200,ns) PDS from vector y-component fields
!     ls       logical*1 (ijoc,ns) bitmaps for scalar fields
!     lv       logical*1 (ijoc,nv) bitmaps for vector fields
!     fs       real (ijoc,ns) scalar fields
!     fu       real (ijoc,ns) vector x-component fields
!     fv       real (ijoc,ns) vector y-component fields
!
! Subprograms called:
!   getgb          read and unpack GRIB message
!   getgbh         read and unpack GRIB header
!   ipolates       general scalar interpolation
!   ipolatev       general vector interpolation
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
  integer,intent(in):: luflx,ns,nv,kgdsc(200),ijoc
  integer,intent(in):: jpdss(1:24,ns),jpdsu(1:24,nv),jpdsv(1:24,nv)
  integer,intent(out):: kpdss(200,ns),kpdsv(200,nv)
  logical*1,intent(out):: ls(ijoc,ns),lv(ijoc,nv)
  real,intent(out):: fs(ijoc,ns),fu(ijoc,nv),fv(ijoc,nv)
  integer jpds1(200),jgds(200),kpds(200),kgds(200)
  integer kfa,kg,kf,k,iret,ko,n
  integer ibis(ns),ibiv(nv),ibos(ns),ibov(nv)
  logical*1,allocatable:: lrs(:,:),lrv(:,:)
  real,allocatable:: frs(:,:),fru(:,:),frv(:,:)
  integer ipopt(20)
  real rlat(ijoc),rlon(ijoc),crot(ijoc),srot(ijoc)
   integer kall,ttot,tmin,tmax
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(ns.le.0.and.nv.le.0) return
  jpds1=-1
  jgds=-1
  call getgbh(luflx,0,0,jpds1,jgds,kg,kfa,k,kpds,kgds,iret)
  if(iret.ne.0) then
    kpdss=-9999
    kpdsv=-9999
    ls=.false.
    lv=.false.
    fs=0.
    fv=0.
    return
  endif
  jpds1(3)=kpds(3)
  jgds=kgds
  allocate(lrs(kfa,ns),lrv(kfa,nv),frs(kfa,ns),fru(kfa,nv),frv(kfa,nv))
  do n=1,ns
    jpds1(19)=jpdss(19,n)
    jpds1(5:7)=jpdss(5:7,n)
    call getgb(luflx,0,kfa,0,jpds1,jgds,&
               kf,k,kpdss(1,n),kgds,lrs(1,n),frs(1,n),iret)
    if(iret.ne.0) then
      kpdss(:,n)=-9999
      lrs(:,n)=.false.
      frs(:,n)=0.
      ibis(n)=1
    else
      ibis(n)=mod(kpdss(4,n)/64,2)
    endif
  enddo
  do n=1,nv
    jpds1(19)=jpdsu(19,n)
    jpds1(5:7)=jpdsu(5:7,n)
    call getgb(luflx,0,kfa,0,jpds1,jgds,&
               kf,k,kpdsv(1,n),kgds,lrv(1,n),fru(1,n),iret)
    if(iret.eq.0) then
      jpds1(19)=jpdsv(19,n)
      jpds1(5:7)=jpdsv(5:7,n)
      call getgb(luflx,0,kfa,0,jpds1,jgds,&
                 kf,k,kpdsv(1,n),kgds,lrv(1,n),frv(1,n),iret)
    endif
    if(iret.ne.0) then
      kpdsv(:,n)=-9999
      lrv(:,n)=.false.
      frv(:,n)=0.
      ibiv(n)=1
    else
      ibiv(n)=mod(kpdsv(4,n)/64,2)
    endif
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ipopt=-1
  ipopt(1)=1
  do n=1,ns
    if(ibis(n).eq.0) then
      ibis(n)=1
      lrs(:,n)=frs(:,n).ne.0.
    endif
  enddo
  call ipolates(1,ipopt,kgds,kgdsc,kf,ijoc,ns,ibis,lrs,frs,&
                ko,rlat,rlon,ibos,ls,fs,iret)
  do n=1,ns
    ls(:ko,n)=ls(:ko,n).or.mod(kpdss(4,n)/64,2).eq.0
  enddo
  call ipolatev(1,ipopt,kgds,kgdsc,kfa,ijoc,nv,ibiv,lrv,fru,frv,&
                ko,rlat,rlon,crot,srot,ibov,lv,fu,fv,iret)
  deallocate(lrs,lrv,frs,fru,frv)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
real function fthour(iftu,ift)
!$$$  Subprogram documentation block
!
! Subprogram: fthour     Convert GRIB forecast time units to hours
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This function returns a forecast hour as encoded in the GRIB PDS.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  ...=fthour(iftu,ift)
!   Input argument list:
!     iftu     integer GRIB PDS forecast time unit
!     ift      integer GRIB PDS forecast time
!   Output argument list:
!     fthour   integer forecast hour
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
  integer,intent(in):: iftu,ift
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(iftu.eq.0) then              ! convert from minutes
    fthour=ift/60.
  elseif(iftu.eq.1) then          ! convert from hours
    fthour=ift
  elseif(iftu.eq.2) then          ! convert from days
    fthour=ift*24.
  elseif(iftu.eq.10) then         ! convert from 3-hours
    fthour=ift*3.
  elseif(iftu.eq.11) then         ! convert from 6-hours
    fthour=ift*6.
  elseif(iftu.eq.12) then         ! convert from 12-hours
    fthour=ift*12.
  elseif(iftu.eq.254) then        ! convert from seconds
    fthour=ift/3600.
  else                            ! unknown or imprecise forecast time unit
    fthour=-1.
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end function
!-------------------------------------------------------------------------------
subroutine makglgds(idrt,imax,jmax,igrid,kgds)
!$$$  Subprogram documentation block
!
! Subprogram: makglgds   Create global GDS
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram creates a global Grid Description Section.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call makglgds(idrt,imax,jmax,igrid,kgds)
!   Input argument list:
!     idrt     integer data representation type (0 for latlon, 4 for Gaussian)
!     imax     integer zonal dimension
!     jmax     integer meridional dimension
!   Output argument list:
!     igrid    integer NCEP grid identifier
!     kgds     integer (200) unpacked GDS
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
  integer,intent(in):: idrt,imax,jmax
  integer,intent(out):: igrid,kgds(200)
  real slat(jmax),wlat(jmax)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  igrid=255
  if(idrt.eq.0.and.imax.eq.144.and.jmax.eq.73) igrid=2
  if(idrt.eq.0.and.imax.eq.360.and.jmax.eq.181) igrid=3
  if(idrt.eq.0.and.imax.eq.720.and.jmax.eq.361) igrid=4
  if(idrt.eq.4.and.imax.eq.192.and.jmax.eq.94) igrid=98
  if(idrt.eq.4.and.imax.eq.384.and.jmax.eq.192) igrid=126
  if(idrt.eq.4.and.imax.eq.512.and.jmax.eq.256) igrid=170
  if(idrt.eq.4.and.imax.eq.768.and.jmax.eq.384) igrid=127
  kgds(1)=modulo(idrt,256)
  kgds(2)=imax
  kgds(3)=jmax
  select case(idrt)
  case(0)
    kgds(4)=90000
  case(4)
    call splat(4,jmax,slat,wlat)
    kgds(4)=nint(180000./acos(-1.)*asin(slat(1)))
  case(256)
    kgds(4)=90000-nint(0.5*180000./jmax)
  end select
  kgds(5)=0
  kgds(6)=128
  kgds(7)=-kgds(4)
  kgds(8)=-nint(360000./imax)
  kgds(9)=-kgds(8)
  select case(idrt)
  case(0)
    kgds(10)=nint(180000./(jmax-1))
  case(4)
    kgds(10)=jmax/2
  case(256)
    kgds(10)=nint(180000./jmax)
  end select
  kgds(11:19)=0
  kgds(20)=255
  kgds(21:200)=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine makglpds(iptv,icen,igen,igrid,ibms,ipu,itl,il1,il2,&
                    iyr,imo,idy,ihr,iftu,ip1,ip2,itr,&
                    ina,inm,icen2,ids,kpds)
!$$$  Subprogram documentation block
!
! Subprogram: makglpds   Create global PDS
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram creates a global Product Description Section.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call makglpds(iptv,icen,igen,igrid,ibms,ipu,itl,il1,il2,&
!                       iyr,imo,idy,ihr,iftu,ip1,ip2,itr,&
!                       ina,inm,icen2,ids,kpds)
!   Input argument list:
!     iptv     integer parameter table version
!     icen     integer center
!     igen     integer model generating code
!     igrid    integer center grid identifier
!     ibms     integer bitmap section flag
!     ipu      integer parameter identifier
!     itl      integer type of level
!     il1      integer level 1
!     il2      integer level 2
!     iyr      integer year
!     imo      integer month
!     idy      integer day
!     ihr      integer hour
!     iftu     integer forecast time unit
!     ip1      integer time 1
!     ip2      integer time 2
!     itr      integer time representation
!     ina      integer number in average
!     inm      integer number missing
!     icen2    integer subcenter
!     ids      integer decimal scaling
!   Output argument list:
!     kpds     integer (200) unpacked PDS
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
  integer,intent(in):: iptv,icen,igen,igrid,ibms,ipu,itl,il1,il2,&
                       iyr,imo,idy,ihr,iftu,ip1,ip2,itr,&
                       ina,inm,icen2,ids
  integer,intent(out):: kpds(200)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  kpds(01)=icen
  kpds(02)=igen
  kpds(03)=igrid
  kpds(04)=128+64*ibms
  kpds(05)=ipu
  kpds(06)=itl
  kpds(07)=256*il1+il2
  kpds(08)=mod(iyr-1,100)+1
  kpds(09)=imo
  kpds(10)=idy
  kpds(11)=ihr
  kpds(12)=0
  kpds(13)=iftu
  kpds(14)=ip1
  kpds(15)=ip2
  kpds(16)=itr
  kpds(17)=ina
  kpds(18)=1
  kpds(19)=iptv
  kpds(20)=inm
  kpds(21)=(iyr-1)/100+1
  kpds(22)=ids
  kpds(23)=icen2
  kpds(24)=0
  kpds(25:)=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine idsset(ids)
!$$$  subprogram documentation block
!
! subprogram: idsset     sets default decimal scalings
!   prgmmr: iredell      org: w/nmc23     date: 92-10-31
!
! abstract: sets decimal scalings defaults for various parameters.
!   a decimal scaling of -3 means data is packed in kilo-si units.
!
! program history log:
!   92-10-31  iredell
!
! usage:call idsset(ids)
!   output arguments:
! ids          integer (255,255) decimal scalings
!              (unknown decimal scalings will not be set)
!
! attributes:
!   language: cray fortran
!
!$$$
  implicit none
  integer,intent(out):: ids(255,255)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  standard fields
  ids(1:8,001)=0      ! pressure (pa)
  ids(1:8,002)=0      ! sea-level pressure (pa)
  ids(1:8,003)=3      ! pressure tendency (pa/s)
  ids(1:8,004)=1      ! potential vorticity (Km2/kg/s)
                      !
  ids(1:8,006)=0      ! geopotential (m2/s2)
  ids(1:8,007)=1      ! geopotential height (m)
  ids(1:8,008)=1      ! geometric height (m)
  ids(1:8,009)=1      ! standard deviation of height (m)
  ids(1:8,010)=1      ! total ozone (dobson)
  ids(1:8,011)=1      ! temperature (k)
  ids(1:8,012)=1      ! virtual temperature (k)
  ids(1:8,013)=1      ! potential temperature (k)
  ids(1:8,014)=1      ! pseudo-adiabatic potential temperature (k)
  ids(1:8,015)=1      ! maximum temperature (k)
  ids(1:8,016)=1      ! minimum temperature (k)
  ids(1:8,017)=1      ! dewpoint temperature (k)
  ids(1:8,018)=1      ! dewpoint depression (k)
  ids(1:8,019)=4      ! temperature lapse rate (k/m)
  ids(1:8,020)=0      ! visibility (m)
                      ! radar spectra 1 ()
                      ! radar spectra 2 ()
                      ! radar spectra 3 ()
                      !
  ids(1:8,025)=1      ! temperature anomaly (k)
  ids(1:8,026)=0      ! pressure anomaly (pa)
  ids(1:8,027)=1      ! geopotential height anomaly (m)
                      ! wave spectra 1 ()
                      ! wave spectra 2 ()
                      ! wave spectra 3 ()
  ids(1:8,031)=0      ! wind direction (degrees)
  ids(1:8,032)=1      ! wind speed (m/s)
  ids(1:8,033)=1      ! zonal wind (m/s)
  ids(1:8,034)=1      ! meridional wind (m/s)
  ids(1:8,035)=-4     ! streamfunction (m2/s)
  ids(1:8,036)=-4     ! velocity potential (m2/s)
  ids(1:8,037)=0      ! montgomery stream function (m2/s2)
  ids(1:8,038)=8      ! sigma vertical velocity (1/s)
  ids(1:8,039)=3      ! pressure vertical velocity (pa/s)
  ids(1:8,040)=4      ! geometric vertical velocity (m/s)
  ids(1:8,041)=6      ! absolute vorticity (1/s)
  ids(1:8,042)=6      ! absolute divergence (1/s)
  ids(1:8,043)=6      ! relative vorticity (1/s)
  ids(1:8,044)=6      ! relative divergence (1/s)
  ids(1:8,045)=4      ! vertical u shear (1/s)
  ids(1:8,046)=4      ! vertical v shear (1/s)
  ids(1:8,047)=0      ! direction of current (degrees)
                      ! speed of current (m/s)
                      ! u of current (m/s)
                      ! v of current (m/s)
  ids(1:8,051)=5      ! specific humidity (kg/kg)
  ids(1:8,052)=0      ! relative humidity (percent)
  ids(1:8,053)=5      ! humidity mixing ratio (kg/kg)
  ids(1:8,054)=1      ! precipitable water (kg/m2)
  ids(1:8,055)=0      ! vapor pressure (pa)
  ids(1:8,056)=0      ! saturation deficit (pa)
  ids(1:8,057)=1      ! evaporation (kg/m2)
  ids(1:8,058)=1      ! cloud ice (kg/m2)
  ids(1:8,059)=6      ! precipitation rate (kg/m2/s)
  ids(1:8,060)=0      ! thunderstorm probability (percent)
  ids(1:8,061)=1      ! total precipitation (kg/m2)
  ids(1:8,062)=1      ! large-scale precipitation (kg/m2)
  ids(1:8,063)=1      ! convective precipitation (kg/m2)
  ids(1:8,064)=6      ! water equivalent snowfall rate (kg/m2/s)
  ids(1:8,065)=0      ! water equivalent of snow depth (kg/m2)
  ids(1:8,066)=2      ! snow depth (m)
                      ! mixed-layer depth (m)
                      ! transient thermocline depth (m)
                      ! main thermocline depth (m)
                      ! main thermocline anomaly (m)
  ids(1:8,071)=0      ! total cloud cover (percent)
  ids(1:8,072)=0      ! convective cloud cover (percent)
  ids(1:8,073)=0      ! low cloud cover (percent)
  ids(1:8,074)=0      ! middle cloud cover (percent)
  ids(1:8,075)=0      ! high cloud cover (percent)
  ids(1:8,076)=2      ! cloud water (kg/m2)
                      !
  ids(1:8,078)=1      ! convective snow (kg/m2)
  ids(1:8,079)=1      ! large scale snow (kg/m2)
  ids(1:8,080)=1      ! water temperature (k)
  ids(1:8,081)=0      ! sea-land mask ()
                      ! deviation of sea level from mean (m)
  ids(1:8,083)=5      ! roughness (m)
  ids(1:8,084)=1      ! albedo (percent)
  ids(1:8,085)=1      ! soil temperature (k)
  ids(1:8,086)=0      ! soil wetness (kg/m2)
  ids(1:8,087)=0      ! vegetation (percent)
                      ! salinity (kg/kg)
  ids(1:8,089)=4      ! density (kg/m3)
  ids(1:8,090)=1      ! runoff (kg/m2)
  ids(1:8,091)=2      ! ice concentration ()
  ids(1:8,092)=2      ! ice thickness (m)
  ids(1:8,093)=0      ! direction of ice drift (degrees)
                      ! speed of ice drift (m/s)
                      ! u of ice drift (m/s)
                      ! v of ice drift (m/s)
                      ! ice growth (m)
                      ! ice divergence (1/s)
  ids(1:8,099)=1      ! snow melt (kg/m2)
                      ! sig height of waves and swell (m)
  ids(1:8,101)=0      ! direction of wind waves (degrees)
                      ! sig height of wind waves (m)
                      ! mean period of wind waves (s)
  ids(1:8,104)=0      ! direction of swell waves (degrees)
                      ! sig height of swell waves (m)
                      ! mean period of swell waves (s)
  ids(1:8,107)=0      ! primary wave direction (degrees)
                      ! primary wave mean period (s)
  ids(1:8,109)=0      ! secondary wave direction (degrees)
                      ! secondary wave mean period (s)
  ids(1:8,111)=0      ! net solar radiative flux at surface (w/m2)
  ids(1:8,112)=0      ! net longwave radiative flux at surface (w/m2)
  ids(1:8,113)=0      ! net solar radiative flux at top (w/m2)
  ids(1:8,114)=0      ! net longwave radiative flux at top (w/m2)
  ids(1:8,115)=0      ! net longwave radiative flux (w/m2)
  ids(1:8,116)=0      ! net solar radiative flux (w/m2)
  ids(1:8,117)=0      ! total radiative flux (w/m2)
                      !
                      !
                      !
  ids(1:8,121)=0      ! latent heat flux (w/m2)
  ids(1:8,122)=0      ! sensible heat flux (w/m2)
  ids(1:8,123)=0      ! boundary layer dissipation (w/m2)
  ids(1:8,124)=3      ! u wind stress (n/m2)
  ids(1:8,125)=3      ! v wind stress (n/m2)
                      ! wind mixing energy (j)
                      ! image data ()
  ids(1:8,128)=0      ! mean sea-level pressure (stdatm) (pa)
  ids(1:8,129)=0      ! mean sea-level pressure (maps) (pa)
  ids(1:8,130)=0      ! mean sea-level pressure (eta) (pa)
  ids(1:8,131)=1      ! surface lifted index (k)
  ids(1:8,132)=1      ! best lifted index (k)
  ids(1:8,133)=1      ! k index (k)
  ids(1:8,134)=1      ! sweat index (k)
  ids(1:8,135)=10     ! horizontal moisture divergence (kg/kg/s)
  ids(1:8,136)=4      ! speed shear (1/s)
  ids(1:8,137)=3      ! 3-hr pressure tendency (pa/s)
  ids(1:8,138)=6      ! brunt-vaisala frequency squared (1/s2)
  ids(1:8,139)=11     ! potential vorticity (mass-weighted) (1/s/m)
  ids(1:8,140)=0      ! rain mask ()
  ids(1:8,141)=0      ! freezing rain mask ()
  ids(1:8,142)=0      ! ice pellets mask ()
  ids(1:8,143)=0      ! snow mask ()
  ids(1:8,144)=3      ! volumetric soil moisture content (fraction)
  ids(1:8,145)=0      ! potential evaporation rate (w/m2)
  ids(1:8,146)=0      ! cloud workfunction (j/kg)
  ids(1:8,147)=3      ! u gravity wave stress (n/m2)
  ids(1:8,148)=3      ! v gravity wave stress (n/m2)
  ids(1:8,149)=10     ! potential vorticity (m2/s/kg)
                      ! covariance between v and u (m2/s2)
                      ! covariance between u and t (k*m/s)
                      ! covariance between v and t (k*m/s)
  ids(1:8,153)=6      ! cloud water mixing ratio (kg/kg)
  ids(1:8,154)=9      ! ozone mixing ratio (kg/kg)
  ids(1:8,155)=0      ! ground heat flux (w/m2)
  ids(1:8,156)=0      ! convective inhibition (j/kg)
  ids(1:8,157)=0      ! convective ape (j/kg)
  ids(1:8,158)=0      ! turbulent ke (j/kg)
  ids(1:8,159)=0      ! condensation pressure of lifted parcel (pa)
  ids(1:8,160)=0      ! clear sky upward solar flux (w/m2)
  ids(1:8,161)=0      ! clear sky downward solar flux (w/m2)
  ids(1:8,162)=0      ! clear sky upward longwave flux (w/m2)
  ids(1:8,163)=0      ! clear sky downward longwave flux (w/m2)
  ids(1:8,164)=0      ! cloud forcing net solar flux (w/m2)
  ids(1:8,165)=0      ! cloud forcing net longwave flux (w/m2)
  ids(1:8,166)=0      ! visible beam downward solar flux (w/m2)
  ids(1:8,167)=0      ! visible diffuse downward solar flux (w/m2)
  ids(1:8,168)=0      ! near ir beam downward solar flux (w/m2)
  ids(1:8,169)=0      ! near ir diffuse downward solar flux (w/m2)
                      !
                      !
  ids(1:8,172)=3      ! momentum flux (n/m2)
  ids(1:8,173)=0      ! mass point model surface ()
  ids(1:8,174)=0      ! velocity point model surface ()
  ids(1:8,175)=0      ! sigma layer number ()
  ids(1:8,176)=2      ! latitude (degrees)
  ids(1:8,177)=2      ! east longitude (degrees)
                      !
                      !
                      !
  ids(1:8,181)=9      ! x-gradient log pressure (1/m)
  ids(1:8,182)=9      ! y-gradient log pressure (1/m)
  ids(1:8,183)=5      ! x-gradient height (m/m)
  ids(1:8,184)=5      ! y-gradient height (m/m)
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
                      !
  ids(1:8,201)=0      ! ice-free water surcace (percent)
                      !
                      !
  ids(1:8,204)=0      ! downward solar radiative flux (w/m2)
  ids(1:8,205)=0      ! downward longwave radiative flux (w/m2)
                      !
  ids(1:8,207)=0      ! moisture availability (percent)
                      ! exchange coefficient (kg/m2/s)
  ids(1:8,209)=0      ! number of mixed layer next to sfc ()
                      !
  ids(1:8,211)=0      ! upward solar radiative flux (w/m2)
  ids(1:8,212)=0      ! upward longwave radiative flux (w/m2)
  ids(1:8,213)=0      ! non-convective cloud cover (percent)
  ids(1:8,214)=6      ! convective precipitation rate (kg/m2/s)
  ids(1:8,215)=7      ! total diabatic heating rate (k/s)
  ids(1:8,216)=7      ! total radiative heating rate (k/s)
  ids(1:8,217)=7      ! total diabatic nonradiative heating rate (k/s)
  ids(1:8,218)=2      ! precipitation index (fraction)
  ids(1:8,219)=1      ! std dev of ir t over 1x1 deg area (k)
  ids(1:8,220)=4      ! natural log of surface pressure over 1 kpa ()
  ids(1:8,221)=1      ! planetary boundary layer height (m)
  ids(1:8,222)=1      ! 5-wave geopotential height (m)
  ids(1:8,223)=1      ! plant canopy surface water (kg/m2)
                      !
                      !
                      ! blackadars mixing length (m)
                      ! asymptotic mixing length (m)
  ids(1:8,228)=1      ! potential evaporation (kg/m2)
  ids(1:8,229)=0      ! snow phase-change heat flux (w/m2)
                      !
  ids(1:8,231)=3      ! convective cloud mass flux (pa/s)
  ids(1:8,232)=0      ! downward total radiation flux (w/m2)
  ids(1:8,233)=0      ! upward total radiation flux (w/m2)
  ids(1:8,234)=1      ! baseflow-groundwater runoff (kg/m2)
  ids(1:8,235)=1      ! storm surface runoff (kg/m2)
                      !
  ids(1:8,237)=6      ! total ozone (kg/m2)
  ids(1:8,238)=0      ! snow cover (percent)
  ids(1:8,239)=1      ! snow temperature (k)
                      !
  ids(1:8,241)=7      ! large scale condensation heating rate (k/s)
  ids(1:8,242)=7      ! deep convective heating rate (k/s)
  ids(1:8,243)=10     ! deep convective moistening rate (kg/kg/s)
  ids(1:8,244)=7      ! shallow convective heating rate (k/s)
  ids(1:8,245)=10     ! shallow convective moistening rate (kg/kg/s)
  ids(1:8,246)=7      ! vertical diffusion heating rate (kg/kg/s)
  ids(1:8,247)=7      ! vertical diffusion zonal acceleration (m/s/s)
  ids(1:8,248)=7      ! vertical diffusion merid acceleration (m/s/s)
  ids(1:8,249)=10     ! vertical diffusion moistening rate (kg/kg/s)
  ids(1:8,250)=7      ! solar radiative heating rate (k/s)
  ids(1:8,251)=7      ! longwave radiative heating rate (k/s)
                      ! drag coefficient ()
                      ! friction velocity (m/s)
                      ! richardson number ()
                      !
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  other fields
  ids(129,200)=2      ! UV-B downward solar flux (W/m2)
  ids(129,201)=2      ! clear sky UV-B downward solar flux (W/m2)
  ids(130,066)=2      ! snow depth (m)
  ids(130,160)=3      ! liquid volumetric soil moisture (non-frozen) ()
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine posta(im,levs,idvc,idsl,nvcoord,vcoord,&
                 hs,ps,psx,psy,t,tx,ty,u,v,d,z,sh,o3,ct,&
                 hrflx,kpdsflxs,kpdsflxv,lflxs,lflxv,fflxs,fflxu,fflxv,&
                 kpo,po,kpt,pt,kzz,zz,kth,th,kpv,pv,pvpt,pvsb,&
                 nfpos,nfpov,nfpts,nfptv,nfzzs,nfzzv,&
                 nfths,nfthv,nfpvs,nfpvv,&
                 kpdsfoxs,kpdsfoxv,lfoxs,lfoxv,&
                 fpos,fpou,fpov,fpts,fptu,fptv,fzzs,fzzu,fzzv,&
                 lths,lthv,fths,fthu,fthv,lpvs,lpvv,fpvs,fpvu,fpvv,&
                 fsuns,fsunu,fsunv,ffoxs,ffoxu,ffoxv)
!$$$  Subprogram documentation block
!
! Subprogram: posta      Post subprogram to vertically interpolate
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes profiles of posted products
!           from profiles of model variables, from soup to nuts.
!           Relative humidity and height are calculated on model surfaces.
!           Then fields are vertically interpolated to specified
!           pressure levels, specified pressure layers, and specified
!           height levels.  Sundry single level fields are also computed.
!           Additional flux fields are calculated if necessary.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call posta(im,levs,idvc,idsl,nvcoord,vcoord,&
!                    hs,ps,psx,psy,t,tx,ty,u,v,d,z,sh,o3,ct,&
!                    hrflx,kpdsflxs,kpdsflxv,lflxs,lflxv,fflxs,fflxu,fflxv,&
!                    kpo,po,kpt,pt,kzz,zz,kth,th,kpv,pv,pvpt,pvsb,&
!                    nfpos,nfpov,nfpts,nfptv,nfzzs,nfzzv,&
!                    nfths,nfthv,nfpvs,nfpvv,&
!                    kpdsfoxs,kpdsfoxv,lfoxs,lfoxv,&
!                    fpos,fpou,fpov,fpts,fptu,fptv,fzzs,fzzu,fzzv,&
!                    lths,lthv,fths,fthu,fthv,lpvs,lpvv,fpvs,fpvu,fpvv,&
!                    fsuns,fsunu,fsunv,ffoxs,ffoxu,ffoxv)
!   Input argument list:
!     im       integer number of profiles
!     levs     integer number of model levels
!     idvc     integer vertical coordinate id (1 for sigma and 2 for hybrid)
!     idsl     integer type of sigma structure (1 for phillips or 2 for mean)
!     nvcoord  integer number of vertical coordinates
!     vcoord   real (km+1,nvcoord) vertical coordinates
!     hs       real (im) surface orography (m)
!     ps       real (im) surface pressure (Pa)
!     psx      real (im) log surface pressure x-gradient (1/m)
!     psy      real (im) log surface pressure y-gradient (1/m)
!     t        real (levs,im) temperature (K)
!     tx       real (levs,im) temperature x-gradient (K/m)
!     ty       real (levs,im) temperature y-gradient (K/m)
!     u        real (levs,im) x-component wind (m/s)
!     v        real (levs,im) y-component wind (m/s)
!     d        real (levs,im) wind divergence (1/s)
!     z        real (levs,im) absolute vorticity (1/s)
!     sh       real (levs,im) specific humidity (kg/kg)
!     o3       real (levs,im) specific ozone (kg/kg)
!     ct       real (levs,im) specific cloud water (kg/kg)
!     hrflx    real hours over which fluxes are averaged
!     kpdsflxs integer (200,nflxs) input scalar flux PDS
!     kpdsflxv integer (200,nflxv) input vector flux PDS
!     lflxs    logical*1 (nflxs,im) scalar flux bitmaps
!     lflxv    logical*1 (nflxv,im) vector flux bitmaps
!     fflxs    real (nflxs,im) scalar flux fields
!     fflxu    real (nflxv,im) vector flux x-component fields
!     fflxv    real (nflxv,im) vector flux y-component fields
!     kpo      integer number of pressure levels
!     po       real (kpo) pressure levels (Pa)
!     kpt      integer number of pressure layers
!     pt       real (kpt) pressure layer edges above the surface (Pa)
!     kzz      integer number of height levels
!     zz       real (kpt) heights above sea level (m)
!     kth      integer number of potential temperature levels
!     th       real (kpv) potential temperatures (K)
!     kpv      integer number of potential vorticity levels
!     pv       real (kpv) potential vorticities (10**-6*K*m**2/kg/s)
!     pvpt     real (kpv) top pressures for PV search (Pa)
!     pvsb     real (kpv) bottom sigmas for PV search ()
!     nfpos    integer first dimension of output pressure level scalar data
!     nfpov    integer first dimension of output pressure level vector data
!     nfpts    integer first dimension of output pressure layer scalar data
!     nfptv    integer first dimension of output pressure layer vector data
!     nfzzs    integer first dimension of output height level scalar data
!     nfzzv    integer first dimension of output height level vector data
!     nfths    integer first dimension of output potential temperature level scalar data
!     nfthv    integer first dimension of output potential temperature level vector data
!     nfpvs    integer first dimension of output potential vorticity level scalar data
!     nfpvv    integer first dimension of output potential vorticity level vector data
!   Output argument list:
!     kpdsfoxs integer (200,nfoxs) output scalar flux PDS
!     kpdsfoxv integer (200,nfoxv) output vector flux PDS
!     lfoxs    logical*1 (nfoxs,im) scalar flux bitmaps
!     lfoxv    logical*1 (nfoxv,im) vector flux bitmaps
!     fpos     real (nfpos,im) pressure level scalar fields
!     fpou     real (nfpos,im) pressure level vector x-component fields
!     fpov     real (nfpos,im) pressure level vector y-component fields
!     fpts     real (nfpts,im) pressure layer scalar fields
!     fptu     real (nfpts,im) pressure layer vector x-component fields
!     fptv     real (nfpts,im) pressure layer vector y-component fields
!     fzzs     real (nfzzs,im) height level scalar fields
!     fzzu     real (nfzzs,im) height level vector x-component fields
!     fzzv     real (nfzzs,im) height level vector y-component fields
!     lths     logical*1 (nfpvs,im) potential temperature level scalar bitmaps
!     lthv     logical*1 (nfpvv,im) potential temperature level vector bitmaps
!     fths     real (nfpvs,im) potential temperature level scalar fields
!     fthu     real (nfpvv,im) potential temperature level vector x-component fields
!     fthv     real (nfpvv,im) potential temperature level vector y-component fields
!     lpvs     logical*1 (nfpvs,im) potential vorticity level scalar bitmaps
!     lpvv     logical*1 (nfpvv,im) potential vorticity level vector bitmaps
!     fpvs     real (nfpvs,im) potential vorticity level scalar fields
!     fpvu     real (nfpvv,im) potential vorticity level vector x-component fields
!     fpvv     real (nfpvv,im) potential vorticity level vector y-component fields
!     fsuns    real (nsuns,im) scalar sundry fields
!     fsunu    real (nsunv,im) vector sundry x-component fields
!     fsunv    real (nsunv,im) vector sundry y-component fields
!     ffoxs    real (nfoxs,im) scalar flux fields
!     ffoxu    real (nfoxv,im) vector flux x-component fields
!     ffoxv    real (nfoxv,im) vector flux y-component fields
!
! Modules used:
!   postgp_module  Shared data for postgp
!
! Subprograms called:
!   flxcnv1        Copy flux metadata
!   modstuff       Compute model coordinate dependent functions
!   getrh          Compute saturation humidity and relative humidity
!   hydro          Compute geopotential heights
!   p2po           Interpolate to pressure level
!   p2pt           Interpolate to pressure layer
!   p2zz           Interpolate to height level
!   sundry         Compute sundry single-level posted fields
!   calwxt1        Compute precipitation type
!   flxcnv         Copy flux data
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use postgp_module
  implicit none
  integer,intent(in):: im,levs,idvc,idsl,nvcoord,kpo,kpt,kzz,kth,kpv
  integer,intent(in):: nfpos,nfpov,nfpts,nfptv,nfzzs,nfzzv
  integer,intent(in):: nfths,nfthv,nfpvs,nfpvv
  real,intent(in):: hrflx
  real,intent(in):: vcoord(levs+1,nvcoord)
  real,intent(in):: hs(im),ps(im),psx(im),psy(im)
  real,intent(in):: t(levs,im),tx(levs,im),ty(levs,im)
  real,intent(in):: u(levs,im),v(levs,im),d(levs,im),z(levs,im)
  real,intent(in):: sh(levs,im),o3(levs,im),ct(levs,im)
  integer,intent(in):: kpdsflxs(200,nflxs),kpdsflxv(200,nflxv)
  logical*1,intent(in):: lflxs(nflxs,im),lflxv(nflxv,im)
  real,intent(in):: fflxs(nflxs,im),fflxu(nflxv,im),fflxv(nflxv,im)
  real,intent(in):: po(kpo),pt(kpt),zz(kzz),th(kth),pv(kpv),pvpt(kpv),pvsb(kpv)
  integer,intent(out):: kpdsfoxs(200,nfoxs),kpdsfoxv(200,nfoxv)
  logical*1,intent(out):: lfoxs(nfoxs,im),lfoxv(nfoxv,im)
  real,intent(out):: fpos(nfpos,im),fpou(nfpov,im),fpov(nfpov,im)
  real,intent(out):: fpts(nfpts,im),fptu(nfptv,im),fptv(nfptv,im)
  real,intent(out):: fzzs(nfzzs,im),fzzu(nfzzv,im),fzzv(nfzzv,im)
  logical*1,intent(out):: lths(nfths,im),lthv(nfthv,im)
  real,intent(out):: fths(nfths,im),fthu(nfthv,im),fthv(nfthv,im)
  logical*1,intent(out):: lpvs(nfpvs,im),lpvv(nfpvv,im)
  real,intent(out):: fpvs(nfpvs,im),fpvu(nfpvv,im),fpvv(nfpvv,im)
  real,intent(out):: fsuns(nsuns,im),fsunu(nsunv,im),fsunv(nsunv,im)
  real,intent(out):: ffoxs(nfoxs,im),ffoxu(nfoxv,im),ffoxv(nfoxv,im)
  real,dimension(kpo):: apo
  real,dimension(levs):: pd,pl,apl,om,px,py,shs,rh,tv,hl
  real,dimension(levs):: hm,se,bvf2,pvn,theta,sigma,pvu
  real,dimension(levs+1):: pi
  real os,aps
  integer ipoh,ipou,ipov,ipot,ipow,ipor,ipoa,ipos,ipoo,ipoc
  integer iptu,iptv,iptt,iptr,iptq
  integer izzu,izzv,izzt,izzr
  integer ithu,ithv,ithh,itht,ithz
  integer ipvu,ipvv,ipvh,ipvt,ipvp,ipvs
  integer isun,ifox
  integer i,k,iwx
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute indices
  ipoh=1+kpo*(ipos_hgt_prs-1)
  ipou=1+kpo*(ipou_ugrd_prs-1)
  ipov=1+kpo*(ipov_vgrd_prs-1)
  ipot=1+kpo*(ipos_tmp_prs-1)
  ipow=1+kpo*(ipos_vvel_prs-1)
  ipor=1+kpo*(ipos_rh_prs-1)
  ipoa=1+kpo*(ipos_absv_prs-1)
  ipos=1+kpo*(ipos_spfh_prs-1)
  ipoo=1+kpo*(ipos_o3mr_prs-1)
  ipoc=1+kpo*(ipos_clwmr_prs-1)
  iptu=1+kpt*(iptu_ugrd_plg-1)
  iptv=1+kpt*(iptv_vgrd_plg-1)
  iptt=1+kpt*(ipts_tmp_plg-1)
  iptr=1+kpt*(ipts_rh_plg-1)
  iptq=1+kpt*(ipts_spfh_plg-1)
  izzu=1+kzz*(izzu_ugrd_hml-1)
  izzv=1+kzz*(izzv_vgrd_hml-1)
  izzt=1+kzz*(izzs_tmp_hml-1)
  izzr=1+kzz*(izzs_rh_hml-1)
  ithu=1+kth*(ithu_ugrd_thel-1)
  ithv=1+kth*(ithv_vgrd_thel-1)
  ithh=1+kth*(iths_mntsf_thel-1)
  itht=1+kth*(iths_tmp_thel-1)
  ithz=1+kth*(iths_pvort_thel-1)
  ipvu=1+kpv*(ipvu_ugrd_pvl-1)
  ipvv=1+kpv*(ipvv_vgrd_pvl-1)
  ipvh=1+kpv*(ipvs_hgt_pvl-1)
  ipvt=1+kpv*(ipvs_tmp_pvl-1)
  ipvp=1+kpv*(ipvs_pres_pvl-1)
  ipvs=1+kpv*(ipvs_vwsh_pvl-1)
  isun=1
  ifox=1
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute log pressures and copy flux PDS data
  apo=log(po)
  call flxcnv1(kpdsflxs,kpdsflxv,kpdsfoxs,kpdsfoxv)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  For every profile, compute model coordinate information,
!  compute relative humidity and heights on model levels,
!  interpolate to pressure level, pressure thickness and height level,
!  compute sundry fields, compute precipitation type, and copy flux fields.
  do i=1,im
    call modstuff(levs,idvc,idsl,nvcoord,vcoord,&
                  ps(i),psx(i),psy(i),d(1,i),u(1,i),v(1,i),&
                  pd,pi,pl,aps,apl,os,om,px,py)
    call getrh(levs,pl,sh(1,i),t(1,i),shs,rh)
    call hydro(levs,hs(i),aps,apl,t(1,i),sh(1,i),tv,hl)
    call pvetc(levs,pl,px,py,tv,tx(1,i),ty(1,i),hl,u(1,i),v(1,i),z(1,i),&
               hm,se,bvf2,pvn,theta,sigma,pvu)
    call p2po(levs,apl,u(1,i),v(1,i),om,hl,t(1,i),rh,sh(1,i),o3(1,i),ct(1,i),&
              kpo,apo,&
              fpou(ipou,i),fpov(ipov,i),fpos(ipow,i),fpos(ipoh,i),&
              fpos(ipot,i),fpos(ipor,i),fpos(ipos,i),fpos(ipoo,i),fpos(ipoc,i))
    call p2pt(levs,ps(i),pd,u(1,i),v(1,i),t(1,i),sh(1,i),shs,kpt,pt,&
              fptu(iptu,i),fptv(iptv,i),fpts(iptt,i),fpts(iptq,i),fpts(iptr,i))
    call p2zz(levs,hl,u(1,i),v(1,i),t(1,i),rh,kzz,zz,&
              fzzu(izzu,i),fzzv(izzv,i),fzzs(izzt,i),fzzs(izzr,i))
    call p2th(levs,theta,u(1,i),v(1,i),hm,t(1,i),pvu,kth,th,&
              lths(ithh,i),fthu(ithu,i),fthv(ithv,i),&
              fths(ithh,i),fths(itht,i),fths(ithz,i))
    lthv(ithu:ithu+kth-1,i)=lths(ithh:ithh+kth-1,i)
    lths(itht:itht+kth-1,i)=lths(ithh:ithh+kth-1,i)
    lths(ithz:ithz+kth-1,i)=lths(ithh:ithh+kth-1,i)
    call p2pv(levs,pvu,hl,t(1,i),pl,u(1,i),v(1,i),kpv,pv,pvpt,pvsb*ps(i),&
              lpvs(ipvh,i),fpvu(ipvu,i),fpvv(ipvv,i),&
              fpvs(ipvh,i),fpvs(ipvt,i),fpvs(ipvp,i),fpvs(ipvs,i))
    lpvv(ipvu:ipvu+kpv-1,i)=lpvs(ipvh:ipvh+kpv-1,i)
    lpvs(ipvt:ipvt+kpv-1,i)=lpvs(ipvh:ipvh+kpv-1,i)
    lpvs(ipvp:ipvp+kpv-1,i)=lpvs(ipvh:ipvh+kpv-1,i)
    lpvs(ipvs:ipvs+kpv-1,i)=lpvs(ipvh:ipvh+kpv-1,i)
    call sundry(hs(i),ps(i),levs,pd,pl,u(1,i),v(1,i),t(1,i),&
                hl,om,rh,sh(1,i),shs,o3(1,i),ct(1,i),&
                kpo,po,fpos(ipoh,i),kpt,pt,fpts(iptt,i),fpts(iptq,i),&
                fsuns(1,i),fsunu(1,i),fsunv(1,i))
    call calwxt1(levs,fflxs(iflxs_prate_sfc,i),&
                 t(levs:1:-1,i),sh(levs:1:-1,i),pl(levs:1:-1),pi(levs+1:1:-1),&
                 iwx)
    call flxcnv(hrflx,lflxs(1,i),lflxv(1,i),fflxs(1,i),fflxu(1,i),fflxv(1,i),&
                iwx,kpdsfoxs,kpdsfoxv,&
                lfoxs(1,i),lfoxv(1,i),ffoxs(1,i),ffoxu(1,i),ffoxv(1,i))
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
contains
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine modstuff(km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
                      pd,pi,pm,aps,apm,os,om,px,py)
!$$$  Subprogram documentation block
!
! Subprogram: modstuff   Compute model coordinate dependent functions
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes fields which depend on the model coordinate
!           such as pressure thickness and vertical velocity.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call modstuff(km,idvc,idsl,nvcoord,vcoord,ps,psx,psy,d,u,v,&
!                       pd,pi,pm,aps,apm,os,om,px,py)
!   Input argument list:
!     km       integer number of levels
!     idvc     integer vertical coordinate id (1 for sigma and 2 for hybrid)
!     idsl     integer type of sigma structure (1 for phillips or 2 for mean)
!     nvcoord  integer number of vertical coordinates
!     vcoord   real (km+1,nvcoord) vertical coordinates
!     ps       real surface pressure (Pa)
!     psx      real log surface pressure x-gradient (1/m)
!     psy      real log surface pressure y-gradient (1/m)
!     d        real (km) wind divergence (1/s)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!   Output argument list:
!     pd       real (km) pressure thickness (Pa)
!     pi       real (km+1) interface pressure (Pa)
!     pm       real (km+1) mid-layer pressure (Pa)
!     aps      real log surface pressure ()
!     apm      real (km+1) log mid-layer pressure ()
!     os       real (km) surface pressure tendency (Pa/s)
!     om       real (km) vertical velocity (Pa/s)
!     px       real (km) mid-layer pressure x-gradient (Pa/m)
!     py       real (km) mid-layer pressure y-gradient (Pa/m)
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use sigio_module
    implicit none
    integer,intent(in):: km,idvc,idsl,nvcoord
    real,intent(in):: vcoord(km+1,nvcoord)
    real,intent(in):: ps,psx,psy
    real,intent(in):: u(km),v(km),d(km)
    real,intent(out):: pd(km),pi(km+1),pm(km)
    real,intent(out):: aps,apm(km),os,om(km),px(km),py(km)
    real dpmdps(km),dpddps(km),dpidps(km+1),vgradp
    integer k,iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call sigio_modpr(1,1,km,nvcoord,idvc,idsl,vcoord,iret,&
                     ps=(/ps/),&
                     pm=pm,pd=pd,dpmdps=dpmdps,dpddps=dpddps)
    pi(1)=ps
    dpidps(1)=1.
    do k=1,km
      pi(k+1)=pi(k)-pd(k)
      dpidps(k+1)=dpidps(k)-dpddps(k)
    enddo
    aps=log(ps)
    apm=log(pm)
    os=0
    do k=km,1,-1
      vgradp=u(k)*psx+v(k)*psy
      os=os-vgradp*ps*(dpmdps(k)-dpidps(k+1))-d(k)*(pm(k)-pi(k+1))
      om(k)=vgradp*ps*dpmdps(k)+os
      os=os-vgradp*ps*(dpidps(k)-dpmdps(k))-d(k)*(pi(k)-pm(k))
    enddo
    px=ps*dpmdps*psx
    py=ps*dpmdps*psy
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine modpr1(km,idvc,idsl,si,ak,bk,ps,pi,pm)
!$$$  subprogram documentation block
!
! subprogram:    modpr1      compute model pressures
!   prgmmr: iredell          org: w/nmc23     date: 92-10-31
!
! abstract: compute model pressures.
!
! program history log:
! 2001-07-25  mark iredell
!
! usage:    call modpr1(km,idvc,idsl,si,ak,bk,ps,pi,pm)
!   input argument list:
!     km           integer number of levels
!     idvc         integer vertical coordinate id
!                  (1 for sigma and 2 for hybrid)
!     idsl         integer type of sigma structure
!                  (1 for phillips or 2 for mean)
!     si           real (km+1) sigma interface values (idvc=1)
!     ak           real (km+1) hybrid interface a (idvc=2)
!     bk           real (km+1) hybrid interface b (idvc=2)
!     ps           real surface pressure (Pa)
!   output argument list:
!     pi           real (km+1) interface pressure (Pa)
!     pm           real (km) mid-layer pressure (Pa)
!
! Subprograms called:
!   fpkap          compute pressure to the kappa
!   frkap          compute pressure to the 1/kappa
!
! attributes:
!   language: fortran
!
!$$$
    use funcphys
    use physcons
    implicit none
    integer,intent(in):: km,idvc,idsl
    real,intent(in):: si(km+1),ak(km+1),bk(km+1),ps
    real,intent(out):: pi(km+1),pm(km)
    integer k
    real(krealfp) pd,pu,pkd,pku
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(idvc.eq.2) then
      do k=1,km+1
        pi(k)=ak(k)+bk(k)*ps
      enddo
    else
      do k=1,km+1
        pi(k)=si(k)*ps
      enddo
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if(idsl.eq.2) then
      do k=1,km
        pm(k)=(pi(k)+pi(k+1))/2
      enddo
    else
      pd=pi(1)
      pkd=fpkap(pd)
      do k=1,km
        if(pi(k+1).gt.0) then
          pu=pi(k+1)
          pku=fpkap(pu)
          pm(k)=frkap((pd*pkd-pu*pku)/((con_rocp+1)*(pd-pu)))
        else
          pu=0
          pku=0
          pm(k)=frkap(pkd/(con_rocp+1))
        endif
        pd=pu
        pkd=pku
      enddo
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine getrh(km,p,sh,t,shs,rh)
!$$$  Subprogram documentation block
!
! Subprogram: getrh      Compute saturation humidity and relative humidity
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes the saturation specific humidity and the
!           relative humidity.  The relative humidity is constrained to be
!           between 0 and 100.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call getrh(km,p,sh,t,shs,rh)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     sh       real (km) specific humidity (kg/kg)
!     t        real (km) temperature (K)
!   Output argument list:
!     shs      real (km) saturation specific humidity (kg/kg)
!     rh       real (km) relative humidity (percent)
!
! Modules used:
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   fpvs           compute saturation vapor pressure
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use funcphys
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in):: p(km),sh(km),t(km)
    real,intent(out):: shs(km),rh(km)
    real(krealfp) pr,tr,es
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      pr=p(k)
      tr=t(k)
      es=fpvs(tr)
      es=min(es,pr)
      shs(k)=con_eps*es/(pr+con_epsm1*es)
      rh(k)=1.e2*min(max(sh(k)/shs(k),0.),1.)
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine hydro(km,hs,aps,ap,t,sh,tv,h)
!$$$  Subprogram documentation block
!
! Subprogram: hydro      Compute geopotential heights
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes geopotential heights by integrating
!           the hydrostatic equation.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call hydro(km,hs,aps,ap,t,sh,tv,h)
!   Input argument list:
!     km       integer number of levels
!     hs       real surface height (m)
!     aps      real log surface pressure ()
!     ap       real (km) log pressure ()
!     t        real (km) temperature (K)
!     sh       real (km) specific humidity (kg/kg)
!   Output argument list:
!     tv       real (km) virtual temperature (K)
!     h        real (km) height (m)
!
! Files included:
!   physcons.h     Physical constants
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in):: hs,aps,ap(km),t(km),sh(km)
    real,intent(out):: tv(km),h(km)
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    tv(1)=t(1)*(1+con_fvirt*sh(1))
    h(1)=hs-con_rog*tv(1)*(ap(1)-aps)
    do k=2,km
      tv(k)=t(k)*(1+con_fvirt*sh(k))
      h(k)=h(k-1)-con_rog*0.5*(tv(k-1)+tv(k))*(ap(k)-ap(k-1))
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine pvetc(km,p,px,py,t,tx,ty,h,u,v,av,hm,s,bvf2,pvn,theta,sigma,pvu)
!$$$  Subprogram documentation block
!
! Subprogram: pvetc      Compute potential vorticity, etc
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes
!             the Montgomery streamfunction
!               hm=cp*t+g*z
!             the specific entropy
!               s=cp*log(t/t0)-r*log(p/p0)
!             the Brunt-Vaisala frequency squared
!               bvf2=g/cp*ds/dz
!             the potential vorticity defined as
!               pvn=(av*ds/dz-dv/dz*ds/dx+du/dz*ds/dy)/rho/cp
!             the potential temperature
!               theta=t0*exp(s/cp)
!             the static stability
!               sigma=t/g*bvf2
!             and the potential vorticity in PV units
!               pvu=10**-6*theta*pvn
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call pvetc(km,p,px,py,t,tx,ty,h,u,v,av,s,bvf2,pvn,theta,sigma,pvu)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     px       real (km) pressure x-gradient (Pa/m)
!     py       real (km) pressure y-gradient (Pa/m)
!     t        real (km) (virtual) temperature (K)
!     tx       real (km) (virtual) temperature x-gradient (K/m)
!     ty       real (km) (virtual) temperature y-gradient (K/m)
!     h        real (km) height (m)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     av       real (km) absolute vorticity (1/s)
!   Output argument list:
!     hm       real (km) Montgomery streamfunction (m**2/s**2)
!     s        real (km) specific entropy (J/K/kg)
!     bvf2     real (km) Brunt-Vaisala frequency squared (1/s**2)
!     pvn      real (km) potential vorticity (m**2/kg/s)
!     theta    real (km) (virtual) potential temperature (K)
!     sigma    real (km) static stability (K/m)
!     pvu      real (km) potential vorticity (10**-6*K*m**2/kg/s)
!
! Modules used:
!   physcons       Physical constants
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,px,py,t,tx,ty,h,u,v,av
    real,intent(out),dimension(km):: hm,s,bvf2,pvn,theta,sigma,pvu
!   real,parameter:: hhmin=500.,t0=2.e2,p0=1.e5
    real,parameter:: hhmin=5.,t0=2.e2,p0=1.e5
    integer k,kd,ku,k2(2)
    real cprho,sx,sy,sz,uz,vz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,km
      hm(k)=con_cp*t(k)+con_g*h(k)
      s(k)=con_cp*log(t(k)/t0)-con_rd*log(p(k)/p0)
    enddo
    do k=1,km
      call rsearch1(km,h,2,(/h(k)-hhmin,h(k)+hhmin/),k2)
      kd=max(k2(1),1)
      ku=min(k2(2)+1,km)
      cprho=p(k)/(con_rocp*t(k))
      sx=con_cp*tx(k)/t(k)-con_rd*px(k)/p(k)
      sy=con_cp*ty(k)/t(k)-con_rd*py(k)/p(k)
      sz=(s(ku)-s(kd))/(h(ku)-h(kd))
      uz=(u(ku)-u(kd))/(h(ku)-h(kd))
      vz=(v(ku)-v(kd))/(h(ku)-h(kd))
      bvf2(k)=con_g/con_cp*sz
      pvn(k)=(av(k)*sz-vz*sx+uz*sy)/cprho
      theta(k)=t0*exp(s(k)/con_cp)
      sigma(k)=t(k)/con_g*bvf2(k)
      pvu(k)=1.e6*theta(k)*pvn(k)
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2po(km,ap,u,v,om,h,t,rh,sh,o3,ct,&
                  kpo,apo,up,vp,omp,hp,tp,rhp,shp,o3p,ctp)
!$$$  Subprogram documentation block
!
! Subprogram: p2po       Interpolate to pressure level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given pressure levels.
!   The interpolation is linear in log pressure within the domain.
!   The height is then integrated hydrostatically by assuming that
!   the virtual temperature is linear in log pressure.
!   Outside the domain, fields are held constant with the following exceptions.
!   Above the top, the height is integrated hydrostatically.
!   Below the bottom, the height and temperature are found by the Shuell method,
!   which limits the errors in extrapolation.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2po(km,ap,u,v,om,h,t,rh,sh,o3,ct,&
!                   kpo,apo,up,vp,omp,hp,tp,rhp,shp,o3p,ctp)
!   Input argument list:
!     km       integer number of levels
!     ap       real (km) log pressure ()
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     om       real (km) vertical velocity (Pa/s)
!     h        real (km) height (m)
!     t        real (km) temperature (K)
!     rh       real (km) relative humidity (percent)
!     sh       real (km) specific humidity (kg/kg)
!     o3       real (km) specific ozone (kg/kg)
!     ct       real (km) specific cloud water (kg/kg)
!     kpo      integer number of pressure levels
!     apo      real (kpo) log pressure levels ()
!   Output argument list:
!     up       real (kpo) x-component wind (m/s)
!     vp       real (kpo) y-component wind (m/s)
!     omp      real (kpo) vertical velocity (Pa/s)
!     hp       real (kpo) height (m)
!     tp       real (kpo) temperature (K)
!     rhp      real (kpo) relative humidity (percent)
!     shp      real (kpo) specific humidity (kg/kg)
!     o3p      real (kpo) specific ozone (kg/kg)
!     ctp      real (kpo) specific cloud water (kg/kg)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km,kpo
    real,intent(in),dimension(km):: ap,u,v,om,h,t,rh,sh,o3,ct
    real,intent(in):: apo(kpo)
    real,intent(out),dimension(kpo):: up,vp,omp,hp,tp,rhp,shp,o3p,ctp
    real,parameter:: gammam=-6.5e-3,zshul=75.,tvshul=290.66
    real w,tvd,tvu,gammas,part
    integer loc(kpo),l
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call rsearch1(km,ap(1),kpo,apo(1),loc(1))
    do k=1,kpo
      l=loc(k)
      if(l.eq.0) then
! below bottom
        up(k)=u(1)
        vp(k)=v(1)
        omp(k)=om(1)
        rhp(k)=rh(1)
        shp(k)=sh(1)
        o3p(k)=o3(1)
        ctp(k)=ct(1)
        tvu=t(1)*(1.+con_fvirt*sh(1))
        if(h(1).gt.zshul) then
          tvd=tvu-gammam*h(1)
          if(tvd.gt.tvshul) then
            if(tvu.gt.tvshul) then
              tvd=tvshul-5.e-3*(tvu-tvshul)**2
            else
              tvd=tvshul
            endif
          endif
          gammas=(tvu-tvd)/h(1)
        else
          gammas=0.
        endif
        part=con_rog*(apo(k)-ap(1))
        hp(k)=h(1)-tvu*part/(1.+0.5*gammas*part)
!       tp(k)=t(1)+gammas*(hp(k)-h(1))
        tp(k)=t(1)+gammam*(hp(k)-h(1))
! within domain
      elseif(l.lt.km) then
        w=(apo(k)-ap(l))/(ap(l+1)-ap(l))
        up(k)=u(l)+w*(u(l+1)-u(l))
        vp(k)=v(l)+w*(v(l+1)-v(l))
        omp(k)=om(l)+w*(om(l+1)-om(l))
        tp(k)=t(l)+w*(t(l+1)-t(l))
        rhp(k)=rh(l)+w*(rh(l+1)-rh(l))
        shp(k)=sh(l)+w*(sh(l+1)-sh(l))
        o3p(k)=o3(l)+w*(o3(l+1)-o3(l))
        ctp(k)=ct(l)+w*(ct(l+1)-ct(l))
        tvd=t(l)*(1+con_fvirt*sh(l))
        tvu=tp(k)*(1+con_fvirt*shp(k))
        hp(k)=h(l)-con_rog*0.5*(tvd+tvu)*(apo(k)-ap(l))
! above top
      else
        up(k)=u(km)
        vp(k)=v(km)
        omp(k)=om(km)
        tp(k)=t(km)
        rhp(k)=rh(km)
        shp(k)=sh(km)
        o3p(k)=o3(km)
        ctp(k)=ct(km)
        tvd=t(km)*(1+con_fvirt*sh(km))
        hp(k)=h(km)-con_rog*tvd*(apo(k)-ap(km))
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2pt(km,ps,pd,u,v,t,sh,shs,kpt,pt,upt,vpt,tpt,shpt,rhpt)
!$$$  Subprogram documentation block
!
! Subprogram: p2pt       Interpolate to pressure layer
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given pressure layers.
!   The pressure layers begin abutting the surface and lay contiguously
!   up into the boundary layer.  The interpolation is actually an integrated
!   mean value in pressure through each output layer.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2pt(km,ps,pd,u,v,t,sh,shs,kpt,pt,upt,vpt,tpt,shpt,rhpt)
!   Input argument list:
!     km       integer number of levels
!     ps       real surface pressure (Pa)
!     pd       real (km) pressure thickness (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     sh       real (km) specific humidity (kg/kg)
!     shs      real (km) saturation specific humidity (kg/kg)
!     kpt      integer number of pressure layers
!     pt       real (kpt) pressure layer edges above surface (Pa)
!   Output argument list:
!     upt      real (kpt) x-component wind (m/s)
!     vpt      real (kpt) y-component wind (m/s)
!     omp      real (kpt) vertical velocity (Pa/s)
!     tpt      real (kpt) temperature (K)
!     shpt     real (kpt) specific humidity (kg/kg)
!     rhpt     real (kpt) relative humidity (percent)
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    implicit none
    integer,intent(in):: km,kpt
    real,intent(in):: ps
    real,intent(in),dimension(km):: pd,u,v,t,sh,shs
    real,intent(in):: pt(kpt)
    real,intent(out),dimension(kpt):: upt,vpt,tpt,shpt,rhpt
    real pi(km+1),pta(0:kpt),ptd
    integer loc(kpt),l1,l2,l
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    pi(1)=ps
    do k=2,km+1
      pi(k)=pi(k-1)-pd(k-1)
    enddo
    pta(0)=ps
    pta(1:kpt)=ps-pt
    call rsearch1(km,pi(1),kpt,pta(1),loc(1))
    l1=1
    do k=1,kpt
      l2=loc(k)
      upt(k)=0.
      vpt(k)=0.
      tpt(k)=0.
      shpt(k)=0.
      rhpt(k)=0.
      do l=l1,l2
        ptd=max(min(pi(l),pta(k-1))-max(pi(l+1),pta(k)),0.)
        upt(k)=upt(k)+ptd*u(l)
        vpt(k)=vpt(k)+ptd*v(l)
        tpt(k)=tpt(k)+ptd*t(l)
        shpt(k)=shpt(k)+ptd*sh(l)
        rhpt(k)=rhpt(k)+ptd*shs(l)
      enddo
      l1=l2
      rhpt(k)=1.e2*min(max(shpt(k)/rhpt(k),0.),1.)
      upt(k)=upt(k)/(pta(k-1)-pta(k))
      vpt(k)=vpt(k)/(pta(k-1)-pta(k))
      tpt(k)=tpt(k)/(pta(k-1)-pta(k))
      shpt(k)=shpt(k)/(pta(k-1)-pta(k))
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2zz(km,h,u,v,t,r,kzz,zz,uzz,vzz,tzz,rzz)
!$$$  Subprogram documentation block
!
! Subprogram: p2zz       Interpolate to height level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given height levels.
!   The output levels are constant height above sea level.  The interpolation
!   is linear in height.  Outside the domain fields are held constant.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2001-05-04  Mark Iredell  Added relative humidity
!
! Usage:  call p2zz(km,h,u,v,t,r,kzz,zz,uzz,vzz,tzz,rzz)
!   Input argument list:
!     km       integer number of levels
!     h        real (km) height (m)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     r        real (km) relative humidity (percent)
!     kzz      integer number of height levels
!     zz       real (kzz) height levels (m)
!   Output argument list:
!     uzz      real (kzz) x-component wind (m/s)
!     vzz      real (kzz) y-component wind (m/s)
!     tzz      real (kzz) temperature (K)
!     rzz      real (kzz) relative humidity (percent)
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    implicit none
    integer,intent(in):: km,kzz
    real,intent(in),dimension(km):: h,u,v,t,r
    real,intent(in):: zz(kzz)
    real,intent(out),dimension(kzz):: uzz,vzz,tzz,rzz
    real w
    integer loc(kzz),l
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call rsearch1(km,h(1),kzz,zz(1),loc(1))
    do k=1,kzz
      l=loc(k)
      if(l.eq.0) then
! below bottom
        uzz(k)=u(1)
        vzz(k)=v(1)
        tzz(k)=t(1)
        rzz(k)=r(1)
! within domain
      elseif(l.lt.km) then
        w=(zz(k)-h(l))/(h(l+1)-h(l))
        uzz(k)=u(l)+w*(u(l+1)-u(l))
        vzz(k)=v(l)+w*(v(l+1)-v(l))
        tzz(k)=t(l)+w*(t(l+1)-t(l))
        rzz(k)=r(l)+w*(r(l+1)-r(l))
! above top
      else
        uzz(k)=u(km)
        vzz(k)=v(km)
        tzz(k)=t(km)
        rzz(k)=r(km)
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2th(km,theta,u,v,h,t,pvu,kth,th,lth,uth,vth,hth,tth,zth)
!$$$  Subprogram documentation block
!
! Subprogram: p2th       Interpolate to isentropic level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given isentropic levels.
!   The interpolation is linear in entropy.
!   Outside the domain the bitmap is set to false.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2th(km,theta,u,v,h,t,puv,kth,th,uth,vth,tth)
!   Input argument list:
!     km       integer number of levels
!     theta    real (km) potential temperature (K)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     h        real (km) height (m)
!     t        real (km) temperature (K)
!     pvu      real (km) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!     kth      integer number of isentropic levels
!     th       real (kth) isentropic levels (K)
!   Output argument list:
!     lpv      logical*1 (kth) bitmap
!     uth      real (kth) x-component wind (m/s)
!     vth      real (kth) y-component wind (m/s)
!     hth      real (kth) height (m)
!     tth      real (kth) temperature (K)
!     zth      real (kth) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    implicit none
    integer,intent(in):: km,kth
    real,intent(in),dimension(km):: theta,u,v,h,t,pvu
    real,intent(in):: th(kth)
    logical*1,intent(out),dimension(kth):: lth
    real,intent(out),dimension(kth):: uth,vth,hth,tth,zth
    real w
    integer loc(kth),l
    integer k
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    call rsearch1(km,theta(1),kth,th(1),loc(1))
    do k=1,kth
      l=loc(k)
      lth(k)=l.gt.0.and.l.lt.km
      if(lth(k)) then
        w=log(th(k)/theta(l))/log(theta(l+1)/theta(l))
        uth(k)=u(l)+w*(u(l+1)-u(l))
        vth(k)=v(l)+w*(v(l+1)-v(l))
        hth(k)=h(l)+w*(h(l+1)-h(l))
        tth(k)=t(l)+w*(t(l+1)-t(l))
        zth(k)=pvu(l)+w*(pvu(l+1)-pvu(l))
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine p2pv(km,pvu,h,t,p,u,v,kpv,pv,pvpt,pvpb,&
                  lpv,upv,vpv,hpv,tpv,ppv,spv)
!$$$  Subprogram documentation block
!
! Subprogram: p2pv       Interpolate to potential vorticity level
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram interpolates fields to given potential vorticity
!   levels within given pressure limits.
!   The output level is the first  encountered from the top pressure limit.
!   If the given potential vorticity level is not found, the outputs are zero
!   and the bitmap is false. The interpolation is linear in potential vorticity.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call p2pv(km,pvu,h,t,p,u,v,kpv,pv,pvpt,pvpb,&
!                   lpv,upv,vpv,hpv,tpv,ppv,spv)
!   Input argument list:
!     km       integer number of levels
!     pvu      real (km) potential vorticity in PV units (10**-6*K*m**2/kg/s)
!     h        real (km) height (m)
!     t        real (km) temperature (K)
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     kpv      integer number of potential vorticity levels
!     pv       real (kpv) potential vorticity levels (10**-6*K*m**2/kg/s)
!     pvpt     real (kpv) top pressures for PV search (Pa)
!     pvpb     real (kpv) bottom pressures for PV search (Pa)
!   Output argument list:
!     lpv      logical*1 (kpv) bitmap
!     upv      real (kpv) x-component wind (m/s)
!     vpv      real (kpv) y-component wind (m/s)
!     hpv      real (kpv) temperature (K)
!     tpv      real (kpv) temperature (K)
!     ppv      real (kpv) pressure (Pa)
!     spv      real (kpv) wind speed shear (1/s)
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km,kpv
    real,intent(in),dimension(km):: pvu,h,t,p,u,v
    real,intent(in):: pv(kpv),pvpt(kpv),pvpb(kpv)
    logical*1,intent(out),dimension(kpv):: lpv
    real,intent(out),dimension(kpv):: upv,vpv,hpv,tpv,ppv,spv
    real,parameter:: pd=2500.
    real w,spdu,spdd
    integer k,l1,l2,lu,ld,l
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do k=1,kpv
      call rsearch1(km,p,1,pvpb(k),l1)
      call rsearch1(km,p,1,pvpt(k),l2)
      l1=l1+1
      l=0
      if(pv(k).ge.0.) then
        do lu=l2-1,l1,-1
          if(pv(k).lt.pvu(lu+1).and.pv(k).ge.pvu(lu)) then
            call rsearch1(km,p,1,p(lu)+pd,ld)
            if(all(pv(k).ge.pvu(ld:lu-1))) then
              l=lu
              exit
            endif
          endif
        enddo
      else
        do lu=l2-1,l1,-1
          if(pv(k).gt.pvu(lu+1).and.pv(k).le.pvu(lu)) then
            call rsearch1(km,p,1,p(lu)+pd,ld)
            if(all(pv(k).le.pvu(ld:lu-1))) then
              l=lu
              exit
            endif
          endif
        enddo
      endif
      lpv(k)=l.gt.0
      if(lpv(k)) then
        w=(pv(k)-pvu(l))/(pvu(l+1)-pvu(l))
        upv(k)=u(l)+w*(u(l+1)-u(l))
        vpv(k)=v(l)+w*(v(l+1)-v(l))
        hpv(k)=h(l)+w*(h(l+1)-h(l))
        tpv(k)=t(l)+w*(t(l+1)-t(l))
        ppv(k)=p(l)*exp((h(l)-hpv(k))*(1-0.5*(tpv(k)/t(l)-1))/(con_rog*t(l)))
        spdu=sqrt(u(l+1)**2+v(l+1)**2)
        spdd=sqrt(u(l)**2+v(l)**2)
        spv(k)=(spdu-spdd)/(h(l+1)-h(l))
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine sundry(hs,ps,km,pd,p,u,v,t,h,om,rh,sh,shs,o3,ct,&
                    kpo,po,hpo,kpt,pt,tpt,shpt,suns,sunu,sunv)
!$$$  Subprogram documentation block
!
! Subprogram: sundry     Compute sundry single-level posted fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes sundry single-level posted fields.
!   The fields returned are surface orography and pressure; tropopause wind,
!   temperature, height, and vertical shear; both surface and best lifted index,
!   convective available potential energy and convective inhibition; maximum
!   wind level wind, temperature, and height; sea level pressure; column
!   precipitable water and average relative humidities; bottom sigma fields;
!   and total ozone and total cloud water.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call sundry(hs,ps,km,pd,p,u,v,t,h,om,rh,sh,shs,o3,ct,&
!                     kpo,po,hpo,kpt,pt,tpt,shpt,suns,sunu,sunv)
!   Input argument list:
!     hs       real surface height (m)
!     ps       real surface pressure (Pa)
!     km       integer number of levels
!     pd       real (km) pressure thickness (Pa)
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!     om       real (km) vertical velocity (Pa/s)
!     rh       real (km) relative humidity (percent)
!     sh       real (km) specific humidity (kg/kg)
!     shs      real (km) saturation specific humidity (kg/kg)
!     o3       real (km) specific ozone (kg/kg)
!     ct       real (km) specific cloud water (kg/kg)
!     kpo      integer number of pressure levels
!     po       real (kpo) pressure levels (Pa)
!     hpo      real (kpo) height (m)
!     kpt      integer number of pressure layers
!     pt       real (kpt) pressure layer edges above surface (Pa)
!     tpt      real (kpt) temperature (K)
!     shpt     real (kpt) specific humidity (kg/kg)
!   Output argument list:
!     suns     real (nsuns) sundry scalar fields
!     sunu     real (nsunv) sundry vector x-component fields
!     sunv     real (nsunv) sundry vector y-component fields
!
! Modules used:
!   postgp_module  Shared data for postgp
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   tpause         compute tropopause level fields
!   liftix         compute lifting index, cape and cin
!   mxwind         compute maximum wind level fields
!   freeze         compute freezing level fields
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use postgp_module
    use funcphys
    use physcons
    implicit none
    integer,intent(in):: km,kpo,kpt
    real,intent(in):: hs,ps
    real,intent(in),dimension(km):: pd,p,u,v,t,h,om,rh,sh,shs,o3,ct
    real,intent(in),dimension(kpo):: po,hpo
    real,intent(in),dimension(kpt):: pt,tpt,shpt
    real,intent(out):: suns(nsuns),sunu(nsunv),sunv(nsunv)
    real,parameter:: pslp(2)=(/1000.e+2,500.e+2/),&
                     pm1=1.e5,tm1=287.45,hm1=113.,hm2=5572.,&
                     fslp=con_g*(hm2-hm1)/(con_rd*tm1)
    real,parameter:: strh1=0.44,strh2=0.72,strh3=0.44,strh4=0.33,&
                     sbrh1=1.00,sbrh2=0.94,sbrh3=0.72,sbrh4=1.00
    real,parameter:: sl1=0.9950
    integer k,kslp(2)
    real sumtn,sumtd,sum1n,sum1d,sum2n,sum2d,sum3n,sum3d,sum4n,sum4d
    real pid,piu,dp1,dp2,dp3,dp4
    real sumo3,sumct,p1,p1k,f2
    real hfac
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  surface orography and surface pressure
    suns(isuns_hgt_sfc)=hs
    suns(isuns_pres_sfc)=ps
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  sundry tropopause fields
    call tpause(km,p,u,v,t,h,&
                suns(isuns_pres_trp),sunu(isunu_ugrd_trp),sunv(isunv_vgrd_trp),&
                suns(isuns_tmp_trp),suns(isuns_hgt_trp),suns(isuns_vwsh_trp))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  lifted index, cape and cin
    call liftix(ps,kpt,pt,tpt,shpt,km,p,t,sh,h,&
                suns(isuns_lftx_sfc),suns(isuns_cape_sfc),suns(isuns_cin_sfc),&
                suns(isuns_blftx_sfc),&
                suns(isuns_cape_plg_180_0),suns(isuns_cin_plg_180_0))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  sundry maxwind fields
    call mxwind(km,p,u,v,t,h,&
                suns(isuns_pres_mwl),sunu(isunu_ugrd_mwl),sunv(isunv_vgrd_mwl),&
                suns(isuns_tmp_mwl),suns(isuns_hgt_mwl))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  sundry freezing fields
    call freeze(km,hs,p,t,h,rh,&
                suns(isuns_hgt_zdeg),suns(isuns_rh_zdeg),&
                suns(isuns_hgt_htfl),suns(isuns_rh_htfl))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  sea level pressure
    call rsearch1(kpo,po(1),2,pslp(1),kslp(1))
    if(kslp(1).gt.0.and.po(kslp(1)).eq.pslp(1).and.&
       kslp(2).gt.0.and.po(kslp(2)).eq.pslp(2)) then
      hfac=hpo(kslp(1))/(hpo(kslp(2))-hpo(kslp(1)))
      suns(isuns_prmsl_msl)=pm1*exp(fslp*hfac)
    else
      suns(isuns_prmsl_msl)=0
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  column precipitable water and average relative humidities
    sumtn=0
    sumtd=0
    sum1n=0
    sum1d=0
    sum2n=0
    sum2d=0
    sum3n=0
    sum3d=0
    sum4n=0
    sum4d=0
    pid=ps
    do k=1,km
      sumtn=sumtn+sh(k)*pd(k)
      sumtd=sumtd+shs(k)*pd(k)
      piu=pid-pd(k)
      dp1=max(min(pid,sbrh1*ps)-max(piu,strh1*ps),0.)
      sum1n=sum1n+sh(k)*dp1
      sum1d=sum1d+shs(k)*dp1
      dp2=max(min(pid,sbrh2*ps)-max(piu,strh2*ps),0.)
      sum2n=sum2n+sh(k)*dp2
      sum2d=sum2d+shs(k)*dp2
      dp3=max(min(pid,sbrh3*ps)-max(piu,strh3*ps),0.)
      sum3n=sum3n+sh(k)*dp3
      sum3d=sum3d+shs(k)*dp3
      dp4=max(min(pid,sbrh4*ps)-max(piu,strh4*ps),0.)
      sum4n=sum4n+sh(k)*dp4
      sum4d=sum4d+shs(k)*dp4
      pid=piu
    enddo
    suns(isuns_pwat_clm)=max(sumtn,0.)/con_g
    suns(isuns_rh_clm)=1.e2*min(max(sumtn/sumtd,0.),1.)
    suns(isuns_rh_slr_044_100)=1.e2*min(max(sum1n/sum1d,0.),1.)
    suns(isuns_rh_slr_072_094)=1.e2*min(max(sum2n/sum2d,0.),1.)
    suns(isuns_rh_slr_044_072)=1.e2*min(max(sum3n/sum3d,0.),1.)
    suns(isuns_rh_slr_033_100)=1.e2*min(max(sum4n/sum4d,0.),1.)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  bottom sigma fields interpolated from first two model layers
    p1=sl1*ps
    f2=log(p(1)/p1)/log(p(1)/p(2))
    p1k=fpkap(real(p1,krealfp))
    suns(isuns_tmp_sig_9950)=t(1)+f2*(t(2)-t(1))
    suns(isuns_pot_sig_9950)=suns(isuns_tmp_sig_9950)/p1k
    suns(isuns_vvel_sig_9950)=om(1)+f2*(om(2)-om(1))
    suns(isuns_rh_sig_9950)=rh(1)+f2*(rh(2)-rh(1))
    sunu(isunu_ugrd_sig_9950)=u(1)+f2*(u(2)-u(1))
    sunv(isunv_vgrd_sig_9950)=v(1)+f2*(v(2)-v(1))
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  total ozone
    sumo3=0
    do k=1,km
      sumo3=sumo3+o3(k)*pd(k)
    enddo
!  convert ozone from kg/m2 to dobson units, which give the depth of the
!  ozone layer in 1e-5 m if brought to natural temperature and pressure.
    suns(isuns_tozne_clm)=max(sumo3,0.)/(con_g*2.14e-5)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  total cloud water
    sumct=0
    do k=1,km
      sumct=sumct+ct(k)*pd(k)
    enddo
    suns(isuns_cwat_clm)=max(sumct,0.)/con_g
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!$$$  Subprogram documentation block
!
! Subprogram: tpause     Compute tropopause level fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram finds the tropopause level and computes fields 
!   at the tropopause level.  The tropopause is defined as the lowest level
!   above 500 mb which has a temperature lapse rate of less than 2 K/km.
!   The lapse rate must average less than 2 K/km over a 2 km depth.
!   If no such level is found below 50 mb, the tropopause is set to 50 mb.
!   The tropopause fields are interpolated linearly in lapse rate.
!   The tropopause pressure is found hydrostatically.
!   The tropopause wind shear is computed as the partial derivative
!   of wind speed with respect to height at the tropopause level.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call tpause(km,p,u,v,t,h,ptp,utp,vtp,ttp,htp,shrtp)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!   Output argument list:
!     ptp      real tropopause pressure (Pa)
!     utp      real tropopause x-component wind (m/s)
!     vtp      real tropopause y-component wind (m/s)
!     ttp      real tropopause temperature (K)
!     htp      real tropopause height (m)
!     shrtp    real tropopause wind shear (1/s)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: ptp,utp,vtp,ttp,htp,shrtp
    real,parameter:: ptplim(2)=(/500.e+2,50.e+2/),gamtp=2.e-3,hd=2.e+3
    real gamu,gamd,td,gami,wtp,spdu,spdd
    integer klim(2),k,kd,ktp
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find tropopause level
    call rsearch1(km-2,p(2),2,ptplim(1),klim(1))
    klim(1)=klim(1)+2
    klim(2)=klim(2)+1
    gamd=1.e+9
    ktp=klim(2)
    wtp=0
    do k=klim(1),klim(2)
      gamu=(t(k-1)-t(k+1))/(h(k+1)-h(k-1))
      if(gamu.le.gamtp) then
        call rsearch1(km-k-1,h(k+1),1,h(k)+hd,kd)
        td=t(k+kd)+(h(k)+hd-h(k+kd))/(h(k+kd+1)-h(k+kd))*(t(k+kd+1)-t(k+kd))
        gami=(t(k)-td)/hd
        if(gami.le.gamtp) then
          ktp=k
          wtp=(gamtp-gamu)/(max(gamd,gamtp+0.1e-3)-gamu)
          exit
        endif
      endif
      gamd=gamu
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute tropopause level fields
    utp=u(ktp)-wtp*(u(ktp)-u(ktp-1))
    vtp=v(ktp)-wtp*(v(ktp)-v(ktp-1))
    ttp=t(ktp)-wtp*(t(ktp)-t(ktp-1))
    htp=h(ktp)-wtp*(h(ktp)-h(ktp-1))
    ptp=p(ktp)*exp((h(ktp)-htp)*(1-0.5*(ttp/t(ktp)-1))/(con_rog*t(ktp)))
    spdu=sqrt(u(ktp)**2+v(ktp)**2)
    spdd=sqrt(u(ktp-1)**2+v(ktp-1)**2)
    shrtp=(spdu-spdd)/(h(ktp)-h(ktp-1))
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine mxwind(km,p,u,v,t,h,pmw,umw,vmw,tmw,hmw)
!$$$  Subprogram documentation block
!
! Subprogram: mxwind     Compute maximum wind level fields
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram finds the maximum wind level and computes fields 
!   at the maximum wind level.  The maximum wind level is searched for
!   between 500 mb and 100 mb.  The height and wind speed at the maximum wind
!   speed level is calculated by assuming the wind speed varies quadratically
!   in height in the neighborhood of the maximum wind level.  The other fields
!   are interpolated linearly in height to the maximum wind level.
!   The maximum wind level pressure is found hydrostatically.
!
! Program history log:
!   1999-10-18  Mark Iredell
!   2005-02-02  Mark Iredell  changed upper limit to 100 mb
!
! Usage:  call mxwind(km,p,u,v,t,h,pmw,umw,vmw,tmw,hmw)
!   Input argument list:
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     u        real (km) x-component wind (m/s)
!     v        real (km) y-component wind (m/s)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!   Output argument list:
!     pmw      real maximum wind level pressure (Pa)
!     umw      real maximum wind level x-component wind (m/s)
!     vmw      real maximum wind level y-component wind (m/s)
!     tmw      real maximum wind level temperature (K)
!     hmw      real maximum wind level height (m)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in),dimension(km):: p,u,v,t,h
    real,intent(out):: pmw,umw,vmw,tmw,hmw
    real,parameter:: pmwlim(2)=(/500.e+2,100.e+2/)
    integer klim(2),k,kmw
    real spd(km),spdmw,wmw,dhd,dhu,shrd,shru,dhmw,ub,vb,spdb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find maximum wind level
    call rsearch1(km,p(1),2,pmwlim(1),klim(1))
    klim(1)=klim(1)+1
    spd(klim(1):klim(2))=sqrt(u(klim(1):klim(2))**2+v(klim(1):klim(2))**2)
    spdmw=spd(klim(1))
    kmw=klim(1)
    do k=klim(1)+1,klim(2)
      if(spd(k).gt.spdmw) then
        spdmw=spd(k)
        kmw=k
      endif
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find speed and height at the maximum wind level
    if(kmw.eq.klim(1).or.kmw.eq.klim(2)) then
      hmw=h(kmw)
      spdmw=spd(kmw)
      wmw=0.
    else
      dhd=h(kmw)-h(kmw-1)
      dhu=h(kmw+1)-h(kmw)
      shrd=(spd(kmw)-spd(kmw-1))/(h(kmw)-h(kmw-1))
      shru=(spd(kmw)-spd(kmw+1))/(h(kmw+1)-h(kmw))
      dhmw=(shrd*dhu-shru*dhd)/(2*(shrd+shru))
      hmw=h(kmw)+dhmw
      spdmw=spd(kmw)+dhmw**2*(shrd+shru)/(dhd+dhu)
      if(dhmw.gt.0) kmw=kmw+1
      wmw=(h(kmw)-hmw)/(h(kmw)-h(kmw-1))
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute maximum wind level fields
    ub=u(kmw)-wmw*(u(kmw)-u(kmw-1))
    vb=v(kmw)-wmw*(v(kmw)-v(kmw-1))
    spdb=max(sqrt(ub**2+vb**2),1.e-6)
    umw=ub*spdmw/spdb
    vmw=vb*spdmw/spdb
    tmw=t(kmw)-wmw*(t(kmw)-t(kmw-1))
    pmw=p(kmw)*exp((h(kmw)-hmw)*(1-0.5*(tmw/t(kmw)-1))/(con_rog*t(kmw)))
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine freeze(km,hs,p,t,h,rh,hfz1,rhfz1,hfz2,rhfz2)
!$$$  Subprogram documentation block
!
! Subprogram: freeze     Compute freezing level fields
!   Prgmmr: Iredell      Org: np23        Date: 2000-09-11
!
! Abstract: This subprogram finds the lowest and highest freezing levels
!   and interpolates the height and relative humidity to those levels.
!   The freezing levels are the respectively the lowest and highest levels
!   below 50 mb where temperature crosses the freezing point.
!   Fields are interpolated linearly in temperature.
!   If temperature is all below freezing, then the relative humidities
!   are set to the bottom relative humidity and the heights are set
!   to where the freezing level would be if the temperature is assumed
!   to be increasing at a lapse rate of 6.5 K/km below the bottom level,
!   with a minimum limit of 1 km below sea level.
!
! Program history log:
!   2000-09-11  Mark Iredell
!   2001-10-02  Mark Iredell  Below ground algorithm changed
!
! Usage:  call freeze(km,hs,p,t,h,rh,hfz1,rhfz1,hfz2,rhfz2)
!   Input argument list:
!     km       integer number of levels
!     hs       real surface height (m)
!     p        real (km) pressure (Pa)
!     t        real (km) temperature (K)
!     h        real (km) height (m)
!     rh       real (km) relative humidity (percent)
!   Output argument list:
!     hfz1     real lowest freezing level height (m)
!     rhfz1    real lowest freezing level relative humidity (percent)
!     hfz2     real highest freezing level height (m)
!     rhfz2    real highest freezing level relative humidity (percent)
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use physcons
    implicit none
    integer,intent(in):: km
    real,intent(in):: hs
    real,intent(in),dimension(km):: p,t,h,rh
    real,intent(out):: hfz1,rhfz1,hfz2,rhfz2
    real tc(km)
    real,parameter:: pfzlim=50.e+2,gammam=-6.5e-3,zfzlim=-1.e+3
    integer klim,k,kfz1,kfz2
    real wfz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find the lowest and highest freezing levels
    call rsearch1(km,p(1),1,pfzlim,klim)
    kfz1=klim+1
    kfz2=1
    do k=1,klim
      tc(k)=t(k)-con_t0c
      if(kfz1.gt.klim.and.tc(k).le.0) kfz1=k
      if(tc(k).gt.0) kfz2=k+1
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  interpolate fields to the lowest freezing level
    if(kfz1.eq.1) then
      hfz1=max(h(1)-tc(1)/gammam,zfzlim)
      rhfz1=rh(1)
    elseif(kfz1.eq.klim+1) then
      hfz1=h(klim)
      rhfz1=rh(klim)
    else
      wfz=-tc(kfz1)/(tc(kfz1-1)-tc(kfz1))
      hfz1=h(kfz1)-wfz*(h(kfz1)-h(kfz1-1))
      rhfz1=rh(kfz1)-wfz*(rh(kfz1)-rh(kfz1-1))
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  interpolate fields to the highest freezing level
    if(kfz2.eq.1) then
      hfz2=max(h(1)-tc(1)/gammam,zfzlim)
      rhfz2=rh(1)
    elseif(kfz2.eq.klim+1) then
      hfz2=h(klim)
      rhfz2=rh(klim)
    else
      wfz=-tc(kfz2)/(tc(kfz2-1)-tc(kfz2))
      hfz2=h(kfz2)-wfz*(h(kfz2)-h(kfz2-1))
      rhfz2=rh(kfz2)-wfz*(rh(kfz2)-rh(kfz2-1))
    endif
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine liftix(ps,kpt,pt,tpt,shpt,km,p,t,sh,h,&
                    sli,scape,scin,bli,bcape,bcin)
!$$$  Subprogram documentation block
!
! Subprogram: liftix     Compute lifting index, cape and cin
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes both surface and best lifting index,
!   convective available potential energy and convective inhibition.
!   The surface quantities result from lifting a parcel from a pressure
!   layer at the surface.  The so-called best quantities result from lifting
!   the pressure layer parcel with the warmest equivalent potential temperature.
!   The lifted index is the parcel temperature minus the environment temperature
!   when lifted to 500 mb.  The convective available potential energy (CAPE)
!   and the convective inhibition (CIN) are found by integrating buoyant energy
!   between the lifting condensation level and the highest buoyant level.
!   The CAPE is the integral of the positively buoyant region and
!   the CIN is the integral of the negatively buoyant region.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call liftix(ps,kpt,pt,tpt,shpt,km,p,t,sh,h,&
!                     sli,scape,scin,bli,bcape,bcin)
!   Input argument list:
!     ps       real surface pressure (Pa)
!     kpt      integer number of pressure layers
!     pt       real (kpt) pressure layer edges above surface (Pa)
!     tpt      real (kpt) temperature (K)
!     shpt     real (kpt) specific humidity (kg/kg)
!     km       integer number of levels
!     p        real (km) pressure (Pa)
!     t        real (km) temperature (K)
!     sh       real (km) specific humidity (kg/kg)
!     h        real (km) height (m)
!   Output argument list:
!     sli      real surface lifted index (K)
!     scape    real surface cape (J/kg)
!     scin     real surface cin (J/kg)
!     bli      real best lifted index (K)
!     bcape    real best cape (J/kg)
!     bcin     real best cin (J/kg)
!
! Modules used:
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   rsearch1       search for a surrounding real interval
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use funcphys
    use physcons
    implicit none
    integer,intent(in):: kpt,km
    real,intent(in):: ps
    real,intent(in),dimension(kpt):: pt,tpt,shpt
    real,intent(in),dimension(km):: p,t,sh,h
    real,intent(out):: sli,scape,scin,bli,bcape,bcin
    real,parameter:: plift=500.e+2,ptop=40.e+2
    integer k
    real(krealfp) pm,pv,tr,tdpd,tlcl,pklcl,thelcl,pkmas,themas,pkmab,themab
    real(krealfp) tlift,pliftr,pliftk,pks,pkb,tma,shma,pr,pk
    real pta(0:kpt),tvenv,gdz,tvcld
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  select the warmest equivalent potential temperature
    pta(0)=ps
    pta(1:kpt)=ps-pt
    do k=1,kpt
      pm=0.5*(pta(k-1)+pta(k))
      pv=pm*shpt(k)/(con_eps-con_epsm1*shpt(k))
      tr=tpt(k)
      tdpd=max(tr-ftdp(pv),0._krealfp)
      tlcl=ftlcl(tr,tdpd)
      pklcl=fpkap(pm)*tlcl/tr
      thelcl=fthe(tlcl,pklcl)
      if(k.eq.1) then
        pkmas=pklcl
        themas=thelcl
        pkmab=pklcl
        themab=thelcl
      elseif(thelcl.gt.themab) then
        pkmab=pklcl
        themab=thelcl
      endif
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  lift the parcel to 500 mb along a dry adiabat below the lcl
!  or along a moist adiabat above the lcl.
!  the lifted index is the environment minus parcel temperature.
    call rsearch1(km,p(1),1,plift,k)
    if(k.gt.0.and.k.lt.km) then
      tlift=t(k)+log(p(k)/plift)/log(p(k)/p(k+1))*(t(k+1)-t(k))
      pliftr=plift
      pliftk=fpkap(pliftr)
      pks=min(pliftk,pkmas)
      call stma(themas,pks,tma,shma)
      sli=tlift-tma*pliftk/pks
      pkb=min(pliftk,pkmab)
      call stma(themab,pkb,tma,shma)
      bli=tlift-tma*pliftk/pkb
    else
      sli=0.
      bli=0.
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  lift the parcel to its highest buoyant layer (below 40 mb).
!  separately integrate between the lcl and this top layer
!  positively buoyant area (cape) and negatively buoyant area (cin).
    scape=0.
    scin=0.
    bcape=0.
    bcin=0.
    do k=km-1,2,-1
      if(p(k).gt.ptop) then
        pr=p(k)
        pk=fpkap(pr)
        tvenv=t(k)*(1.+con_fvirt*sh(k))
        gdz=con_g*0.5*(h(k+1)-h(k-1))
        if(pk.le.pkmas) then
          call stma(themas,pk,tma,shma)
          tvcld=tma*(1.+con_fvirt*shma)
          if(tvcld.gt.tvenv) then
            scape=scape+gdz*log(tvcld/tvenv)
          elseif(scape.gt.0.) then
            scin=scin+gdz*log(tvcld/tvenv)
          endif
        endif
        if(pk.le.pkmab) then
          call stma(themab,pk,tma,shma)
          tvcld=tma*(1.+con_fvirt*shma)
          if(tvcld.gt.tvenv) then
            bcape=bcape+gdz*log(tvcld/tvenv)
          elseif(bcape.gt.0.) then
            bcin=bcin+gdz*log(tvcld/tvenv)
          endif
        endif
      endif
    enddo
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine flxcnv1(kpdsflxs,kpdsflxv,kpdsfoxs,kpdsfoxv)
!$$$  Subprogram documentation block
!
! Subprogram: flxcnv1    Copy flux metadata
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram copies flux metadata from the input flux arrays
!   to the output flux arrays.  The metadata copied consist of words 13-16
!   of the integer PDS arrays, which refer to forecast hour information.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call flxcnv1(kpdsflxs,kpdsflxv,kpdsfoxs,kpdsfoxv)
!   Input argument list:
!     kpdsflxs integer (200,nflxs) input flux scalar metadata
!     kpdsflxv integer (200,nflxv) input flux vector metadata
!   Output argument list:
!     kpdsfoxs integer (200,nfoxs) output flux scalar metadata
!     kpdsfoxv integer (200,nfoxv) output flux vector metadata
!
! Modules used:
!   postgp_module  Shared data for postgp
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use postgp_module
    implicit none
    integer,intent(in):: kpdsflxs(200,nflxs),kpdsflxv(200,nflxv)
    integer,intent(out):: kpdsfoxs(200,nfoxs),kpdsfoxv(200,nfoxv)
    integer i
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    do i=13,16
      kpdsfoxs(i,ifoxs_shtfl_sfc         )=kpdsflxs(i,iflxs_shtfl_sfc         )
      kpdsfoxs(i,ifoxs_lhtfl_sfc         )=kpdsflxs(i,iflxs_lhtfl_sfc         )
      kpdsfoxs(i,ifoxs_tmp_sfc           )=kpdsflxs(i,iflxs_tmp_sfc           )
      kpdsfoxs(i,ifoxs_soilw_dlr_0_10    )=kpdsflxs(i,iflxs_soilw_dlr_0_10    )
      kpdsfoxs(i,ifoxs_soilw_dlr_10_40   )=kpdsflxs(i,iflxs_soilw_dlr_10_40   )
      kpdsfoxs(i,ifoxs_soilw_dlr_40_100  )=kpdsflxs(i,iflxs_soilw_dlr_40_100  )
      kpdsfoxs(i,ifoxs_soilw_dlr_100_200 )=kpdsflxs(i,iflxs_soilw_dlr_100_200 )
      kpdsfoxs(i,ifoxs_soilw_dlr_10_200  )=kpdsflxs(i,iflxs_soilw_dlr_10_200  )
      kpdsfoxs(i,ifoxs_tmp_dlr_0_10      )=kpdsflxs(i,iflxs_tmp_dlr_0_10      )
      kpdsfoxs(i,ifoxs_tmp_dlr_10_40     )=kpdsflxs(i,iflxs_tmp_dlr_10_40     )
      kpdsfoxs(i,ifoxs_tmp_dlr_40_100    )=kpdsflxs(i,iflxs_tmp_dlr_40_100    )
      kpdsfoxs(i,ifoxs_tmp_dlr_100_200   )=kpdsflxs(i,iflxs_tmp_dlr_100_200   )
      kpdsfoxs(i,ifoxs_tmp_dlr_10_200    )=kpdsflxs(i,iflxs_tmp_dlr_10_200    )
      kpdsfoxs(i,ifoxs_soill_dlr_0_10    )=kpdsflxs(i,iflxs_soill_dlr_0_10    )
      kpdsfoxs(i,ifoxs_soill_dlr_10_40   )=kpdsflxs(i,iflxs_soill_dlr_10_40   )
      kpdsfoxs(i,ifoxs_soill_dlr_40_100  )=kpdsflxs(i,iflxs_soill_dlr_40_100  )
      kpdsfoxs(i,ifoxs_soill_dlr_100_200 )=kpdsflxs(i,iflxs_soill_dlr_100_200 )
      kpdsfoxs(i,ifoxs_cnwat_sfc         )=kpdsflxs(i,iflxs_cnwat_sfc         )
      kpdsfoxs(i,ifoxs_snod_sfc          )=kpdsflxs(i,iflxs_snod_sfc          )
      kpdsfoxs(i,ifoxs_weasd_sfc         )=kpdsflxs(i,iflxs_weasd_sfc         )
      kpdsfoxs(i,ifoxs_dlwrf_sfc         )=kpdsflxs(i,iflxs_dlwrf_sfc         )
      kpdsfoxs(i,ifoxs_ulwrf_sfc         )=kpdsflxs(i,iflxs_ulwrf_sfc         )
      kpdsfoxs(i,ifoxs_ulwrf_toa         )=kpdsflxs(i,iflxs_ulwrf_toa         )
      kpdsfoxs(i,ifoxs_uswrf_toa         )=kpdsflxs(i,iflxs_uswrf_toa         )
      kpdsfoxs(i,ifoxs_uswrf_sfc         )=kpdsflxs(i,iflxs_uswrf_sfc         )
      kpdsfoxs(i,ifoxs_dswrf_sfc         )=kpdsflxs(i,iflxs_dswrf_sfc         )
      kpdsfoxs(i,ifoxs_duvb_sfc          )=kpdsflxs(i,iflxs_duvb_sfc          )
      kpdsfoxs(i,ifoxs_cduvb_sfc         )=kpdsflxs(i,iflxs_cduvb_sfc         )
      kpdsfoxs(i,ifoxs_tcdc_hcl          )=kpdsflxs(i,iflxs_tcdc_hcl          )
      kpdsfoxs(i,ifoxs_pres_hct          )=kpdsflxs(i,iflxs_pres_hct          )
      kpdsfoxs(i,ifoxs_pres_hcb          )=kpdsflxs(i,iflxs_pres_hcb          )
      kpdsfoxs(i,ifoxs_tmp_hct           )=kpdsflxs(i,iflxs_tmp_hct           )
      kpdsfoxs(i,ifoxs_tcdc_mcl          )=kpdsflxs(i,iflxs_tcdc_mcl          )
      kpdsfoxs(i,ifoxs_pres_mct          )=kpdsflxs(i,iflxs_pres_mct          )
      kpdsfoxs(i,ifoxs_pres_mcb          )=kpdsflxs(i,iflxs_pres_mcb          )
      kpdsfoxs(i,ifoxs_tmp_mct           )=kpdsflxs(i,iflxs_tmp_mct           )
      kpdsfoxs(i,ifoxs_tcdc_lcl          )=kpdsflxs(i,iflxs_tcdc_lcl          )
      kpdsfoxs(i,ifoxs_pres_lct          )=kpdsflxs(i,iflxs_pres_lct          )
      kpdsfoxs(i,ifoxs_pres_lcb          )=kpdsflxs(i,iflxs_pres_lcb          )
      kpdsfoxs(i,ifoxs_tmp_lct           )=kpdsflxs(i,iflxs_tmp_lct           )
      kpdsfoxs(i,ifoxs_prate_sfc         )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_cprat_sfc         )=kpdsflxs(i,iflxs_cprat_sfc         )
      kpdsfoxs(i,ifoxs_gflux_sfc         )=kpdsflxs(i,iflxs_gflux_sfc         )
      kpdsfoxs(i,ifoxs_land_sfc          )=kpdsflxs(i,iflxs_land_sfc          )
      kpdsfoxs(i,ifoxs_icec_sfc          )=kpdsflxs(i,iflxs_icec_sfc          )
      kpdsfoxs(i,ifoxs_icetk_sfc         )=kpdsflxs(i,iflxs_icetk_sfc         )
      kpdsfoxs(i,ifoxs_tmp_hag_2         )=kpdsflxs(i,iflxs_tmp_hag_2         )
      kpdsfoxs(i,ifoxs_spfh_hag_2        )=kpdsflxs(i,iflxs_spfh_hag_2        )
      kpdsfoxs(i,ifoxs_tmax_hag_2        )=kpdsflxs(i,iflxs_tmax_hag_2        )
      kpdsfoxs(i,ifoxs_tmin_hag_2        )=kpdsflxs(i,iflxs_tmin_hag_2        )
      kpdsfoxs(i,ifoxs_watr_sfc          )=kpdsflxs(i,iflxs_watr_sfc          )
      kpdsfoxs(i,ifoxs_pevpr_sfc         )=kpdsflxs(i,iflxs_pevpr_sfc         )
      kpdsfoxs(i,ifoxs_cwork_clm         )=kpdsflxs(i,iflxs_cwork_clm         )
      kpdsfoxs(i,ifoxs_hpbl_sfc          )=kpdsflxs(i,iflxs_hpbl_sfc          )
      kpdsfoxs(i,ifoxs_albdo_sfc         )=kpdsflxs(i,iflxs_albdo_sfc         )
      kpdsfoxs(i,ifoxs_tcdc_clm          )=kpdsflxs(i,iflxs_tcdc_clm          )
      kpdsfoxs(i,ifoxs_tcdc_cvc          )=kpdsflxs(i,iflxs_tcdc_cvc          )
      kpdsfoxs(i,ifoxs_pres_cvt          )=kpdsflxs(i,iflxs_pres_cvt          )
      kpdsfoxs(i,ifoxs_pres_cvb          )=kpdsflxs(i,iflxs_pres_cvb          )
      kpdsfoxs(i,ifoxs_tcdc_blc          )=kpdsflxs(i,iflxs_tcdc_blc          )
      kpdsfoxs(i,ifoxs_apcp_sfc          )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_acpcp_sfc         )=kpdsflxs(i,iflxs_cprat_sfc         )
      kpdsfoxs(i,ifoxs_crain_sfc         )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_cfrzr_sfc         )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_cicep_sfc         )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_csnow_sfc         )=kpdsflxs(i,iflxs_prate_sfc         )
      kpdsfoxs(i,ifoxs_rh_hag_2          )=kpdsflxs(i,iflxs_tmp_hag_2         )
      kpdsfoxv(i,ifoxu_uflx_sfc          )=kpdsflxv(i,iflxu_uflx_sfc          )
      kpdsfoxv(i,ifoxu_ugrd_hag_10       )=kpdsflxv(i,iflxu_ugrd_hag_10       )
      kpdsfoxv(i,ifoxu_ugwd_sfc          )=kpdsflxv(i,iflxu_ugwd_sfc          )
    enddo
    kpdsfoxs(16,ifoxs_apcp_sfc          )=4
    kpdsfoxs(16,ifoxs_acpcp_sfc         )=4
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine flxcnv(hrflx,lflxs,lflxv,fflxs,fflxu,fflxv,iwx,&
                    kpdsfoxs,kpdsfoxv,lfoxs,lfoxv,ffoxs,ffoxu,ffoxv)
!$$$  Subprogram documentation block
!
! Subprogram: flxcnv     Copy flux data
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram copies flux data from the input flux arrays
!   to the output flux arrays.  The accumulated precipitation fields are
!   computed from the average precipitation rate fields.  The categorical
!   precipitation types are decoded.  The 2-meter relative humidity is computed.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call flxcnv(hrflx,lflxs,lflxv,fflxs,fflxu,fflxv,iwx,&
!                     kpdsfoxs,kpdsfoxv,lfoxs,lfoxv,ffoxs,ffoxu,ffoxv)
!   Input argument list:
!     hrflx    real number of hours flux over which fields are averaged
!     lflxs    logical*1 (nflxs) input scalar bitmap
!     lflxv    logical*1 (nflxv) input vector bitmap
!     fflxs    real (nflxs) input scalar data
!     fflxu    real (nflxv) input vector x-component data
!     fflxv    real (nflxv) input vector y-component data
!     iwx      integer encoded categorical precipitation type
!     kpdsfoxs integer (200,nfoxs) output flux scalar metadata
!     kpdsfoxv integer (200,nfoxv) output flux vector metadata
!   Output argument list:
!     lfoxs    logical*1 (nfoxs) output scalar bitmap
!     lfoxv    logical*1 (nfoxv) output vector bitmap
!     ffoxs    real (nfoxs) output scalar data
!     ffoxu    real (nfoxv) output vector x-component data
!     ffoxv    real (nfoxv) output vector y-component data
!
! Modules used:
!   postgp_module  Shared data for postgp
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use postgp_module
    use funcphys
    use physcons
    implicit none
    real,intent(in):: hrflx
    logical*1,intent(in):: lflxs(nflxs),lflxv(nflxv)
    real,intent(in):: fflxs(nflxs),fflxu(nflxv),fflxv(nflxv)
    integer,intent(in):: iwx
    integer,intent(in):: kpdsfoxs(200,nfoxs),kpdsfoxv(200,nfoxv)
    logical*1,intent(out):: lfoxs(nfoxs),lfoxv(nfoxv)
    real,intent(out):: ffoxs(nfoxs),ffoxu(nfoxv),ffoxv(nfoxv)
    integer n
    real t2,q2,ps,p2,pvs2,qs2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Copy bitmaps.
    lfoxs(ifoxs_shtfl_sfc         )=lflxs(iflxs_shtfl_sfc         )
    lfoxs(ifoxs_lhtfl_sfc         )=lflxs(iflxs_lhtfl_sfc         )
    lfoxs(ifoxs_tmp_sfc           )=lflxs(iflxs_tmp_sfc           )
    lfoxs(ifoxs_soilw_dlr_0_10    )=lflxs(iflxs_soilw_dlr_0_10    )
    lfoxs(ifoxs_soilw_dlr_10_40   )=lflxs(iflxs_soilw_dlr_10_40   )
    lfoxs(ifoxs_soilw_dlr_40_100  )=lflxs(iflxs_soilw_dlr_40_100  )
    lfoxs(ifoxs_soilw_dlr_100_200 )=lflxs(iflxs_soilw_dlr_100_200 )
    lfoxs(ifoxs_soilw_dlr_10_200  )=lflxs(iflxs_soilw_dlr_10_200  )
    lfoxs(ifoxs_tmp_dlr_0_10      )=lflxs(iflxs_tmp_dlr_0_10      )
    lfoxs(ifoxs_tmp_dlr_10_40     )=lflxs(iflxs_tmp_dlr_10_40     )
    lfoxs(ifoxs_tmp_dlr_40_100    )=lflxs(iflxs_tmp_dlr_40_100    )
    lfoxs(ifoxs_tmp_dlr_100_200   )=lflxs(iflxs_tmp_dlr_100_200   )
    lfoxs(ifoxs_tmp_dlr_10_200    )=lflxs(iflxs_tmp_dlr_10_200    )
    lfoxs(ifoxs_soill_dlr_0_10    )=lflxs(iflxs_soill_dlr_0_10    )
    lfoxs(ifoxs_soill_dlr_10_40   )=lflxs(iflxs_soill_dlr_10_40   )
    lfoxs(ifoxs_soill_dlr_40_100  )=lflxs(iflxs_soill_dlr_40_100  )
    lfoxs(ifoxs_soill_dlr_100_200 )=lflxs(iflxs_soill_dlr_100_200 )
    lfoxs(ifoxs_cnwat_sfc         )=lflxs(iflxs_cnwat_sfc         )
    lfoxs(ifoxs_snod_sfc          )=lflxs(iflxs_snod_sfc          )
    lfoxs(ifoxs_weasd_sfc         )=lflxs(iflxs_weasd_sfc         )
    lfoxs(ifoxs_dlwrf_sfc         )=lflxs(iflxs_dlwrf_sfc         )
    lfoxs(ifoxs_ulwrf_sfc         )=lflxs(iflxs_ulwrf_sfc         )
    lfoxs(ifoxs_ulwrf_toa         )=lflxs(iflxs_ulwrf_toa         )
    lfoxs(ifoxs_uswrf_toa         )=lflxs(iflxs_uswrf_toa         )
    lfoxs(ifoxs_uswrf_sfc         )=lflxs(iflxs_uswrf_sfc         )
    lfoxs(ifoxs_dswrf_sfc         )=lflxs(iflxs_dswrf_sfc         )
    lfoxs(ifoxs_duvb_sfc          )=lflxs(iflxs_duvb_sfc          )
    lfoxs(ifoxs_cduvb_sfc         )=lflxs(iflxs_cduvb_sfc         )
    lfoxs(ifoxs_tcdc_hcl          )=lflxs(iflxs_tcdc_hcl          )
    lfoxs(ifoxs_pres_hct          )=lflxs(iflxs_pres_hct          )
    lfoxs(ifoxs_pres_hcb          )=lflxs(iflxs_pres_hcb          )
    lfoxs(ifoxs_tmp_hct           )=lflxs(iflxs_tmp_hct           )
    lfoxs(ifoxs_tcdc_mcl          )=lflxs(iflxs_tcdc_mcl          )
    lfoxs(ifoxs_pres_mct          )=lflxs(iflxs_pres_mct          )
    lfoxs(ifoxs_pres_mcb          )=lflxs(iflxs_pres_mcb          )
    lfoxs(ifoxs_tmp_mct           )=lflxs(iflxs_tmp_mct           )
    lfoxs(ifoxs_tcdc_lcl          )=lflxs(iflxs_tcdc_lcl          )
    lfoxs(ifoxs_pres_lct          )=lflxs(iflxs_pres_lct          )
    lfoxs(ifoxs_pres_lcb          )=lflxs(iflxs_pres_lcb          )
    lfoxs(ifoxs_tmp_lct           )=lflxs(iflxs_tmp_lct           )
    lfoxs(ifoxs_prate_sfc         )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_cprat_sfc         )=lflxs(iflxs_cprat_sfc         )
    lfoxs(ifoxs_gflux_sfc         )=lflxs(iflxs_gflux_sfc         )
    lfoxs(ifoxs_land_sfc          )=lflxs(iflxs_land_sfc          )
    lfoxs(ifoxs_icec_sfc          )=lflxs(iflxs_icec_sfc          )
    lfoxs(ifoxs_icetk_sfc         )=lflxs(iflxs_icetk_sfc         )
    lfoxs(ifoxs_tmp_hag_2         )=lflxs(iflxs_tmp_hag_2         )
    lfoxs(ifoxs_spfh_hag_2        )=lflxs(iflxs_spfh_hag_2        )
    lfoxs(ifoxs_tmax_hag_2        )=lflxs(iflxs_tmax_hag_2        )
    lfoxs(ifoxs_tmin_hag_2        )=lflxs(iflxs_tmin_hag_2        )
    lfoxs(ifoxs_watr_sfc          )=lflxs(iflxs_watr_sfc          )
    lfoxs(ifoxs_pevpr_sfc         )=lflxs(iflxs_pevpr_sfc         )
    lfoxs(ifoxs_cwork_clm         )=lflxs(iflxs_cwork_clm         )
    lfoxs(ifoxs_hpbl_sfc          )=lflxs(iflxs_hpbl_sfc          )
    lfoxs(ifoxs_albdo_sfc         )=lflxs(iflxs_albdo_sfc         )
    lfoxs(ifoxs_tcdc_clm          )=lflxs(iflxs_tcdc_clm          )
    lfoxs(ifoxs_tcdc_cvc          )=lflxs(iflxs_tcdc_cvc          )
    lfoxs(ifoxs_pres_cvt          )=lflxs(iflxs_pres_cvt          )
    lfoxs(ifoxs_pres_cvb          )=lflxs(iflxs_pres_cvb          )
    lfoxs(ifoxs_tcdc_blc          )=lflxs(iflxs_tcdc_blc          )
    lfoxs(ifoxs_apcp_sfc          )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_acpcp_sfc         )=lflxs(iflxs_cprat_sfc         )
    lfoxs(ifoxs_crain_sfc         )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_cfrzr_sfc         )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_cicep_sfc         )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_csnow_sfc         )=lflxs(iflxs_prate_sfc         )
    lfoxs(ifoxs_rh_hag_2          )=lflxs(iflxs_tmp_hag_2         ).and.&
                                    lflxs(iflxs_spfh_hag_2        ).and.&
                                    lflxs(iflxs_pres_sfc          )
    lfoxv(ifoxu_uflx_sfc          )=lflxv(iflxu_uflx_sfc          )
    lfoxv(ifoxu_ugrd_hag_10       )=lflxv(iflxu_ugrd_hag_10       )
    lfoxv(ifoxu_ugwd_sfc          )=lflxv(iflxu_ugwd_sfc          )
    do n=1,nfoxs
      if(kpdsfoxs(16,n).ge.2.and.kpdsfoxs(16,n).le.5.and.hrflx.le.0) then
        lfoxs(n)=.false.
      endif
    enddo
    do n=1,nfoxv
      if(kpdsfoxv(16,n).ge.2.and.kpdsfoxv(16,n).le.5.and.hrflx.le.0) then
        lfoxv(n)=.false.
      endif
    enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Copy scalar data.
    ffoxs(ifoxs_shtfl_sfc         )=fflxs(iflxs_shtfl_sfc         )
    ffoxs(ifoxs_lhtfl_sfc         )=fflxs(iflxs_lhtfl_sfc         )
    ffoxs(ifoxs_tmp_sfc           )=fflxs(iflxs_tmp_sfc           )
    ffoxs(ifoxs_soilw_dlr_0_10    )=fflxs(iflxs_soilw_dlr_0_10    )
    ffoxs(ifoxs_soilw_dlr_10_40   )=fflxs(iflxs_soilw_dlr_10_40   )
    ffoxs(ifoxs_soilw_dlr_40_100  )=fflxs(iflxs_soilw_dlr_40_100  )
    ffoxs(ifoxs_soilw_dlr_100_200 )=fflxs(iflxs_soilw_dlr_100_200 )
    ffoxs(ifoxs_soilw_dlr_10_200  )=fflxs(iflxs_soilw_dlr_10_200  )
    ffoxs(ifoxs_tmp_dlr_0_10      )=fflxs(iflxs_tmp_dlr_0_10      )
    ffoxs(ifoxs_tmp_dlr_10_40     )=fflxs(iflxs_tmp_dlr_10_40     )
    ffoxs(ifoxs_tmp_dlr_40_100    )=fflxs(iflxs_tmp_dlr_40_100    )
    ffoxs(ifoxs_tmp_dlr_100_200   )=fflxs(iflxs_tmp_dlr_100_200   )
    ffoxs(ifoxs_tmp_dlr_10_200    )=fflxs(iflxs_tmp_dlr_10_200    )
    ffoxs(ifoxs_soill_dlr_0_10    )=fflxs(iflxs_soill_dlr_0_10    )
    ffoxs(ifoxs_soill_dlr_10_40   )=fflxs(iflxs_soill_dlr_10_40   )
    ffoxs(ifoxs_soill_dlr_40_100  )=fflxs(iflxs_soill_dlr_40_100  )
    ffoxs(ifoxs_soill_dlr_100_200 )=fflxs(iflxs_soill_dlr_100_200 )
    ffoxs(ifoxs_cnwat_sfc         )=fflxs(iflxs_cnwat_sfc         )
    ffoxs(ifoxs_snod_sfc          )=fflxs(iflxs_snod_sfc          )
    ffoxs(ifoxs_weasd_sfc         )=fflxs(iflxs_weasd_sfc         )
    ffoxs(ifoxs_dlwrf_sfc         )=fflxs(iflxs_dlwrf_sfc         )
    ffoxs(ifoxs_ulwrf_sfc         )=fflxs(iflxs_ulwrf_sfc         )
    ffoxs(ifoxs_ulwrf_toa         )=fflxs(iflxs_ulwrf_toa         )
    ffoxs(ifoxs_uswrf_toa         )=fflxs(iflxs_uswrf_toa         )
    ffoxs(ifoxs_uswrf_sfc         )=fflxs(iflxs_uswrf_sfc         )
    ffoxs(ifoxs_dswrf_sfc         )=fflxs(iflxs_dswrf_sfc         )
    ffoxs(ifoxs_duvb_sfc          )=fflxs(iflxs_duvb_sfc          )
    ffoxs(ifoxs_cduvb_sfc         )=fflxs(iflxs_cduvb_sfc         )
    ffoxs(ifoxs_tcdc_hcl          )=fflxs(iflxs_tcdc_hcl          )
    ffoxs(ifoxs_pres_hct          )=fflxs(iflxs_pres_hct          )
    ffoxs(ifoxs_pres_hcb          )=fflxs(iflxs_pres_hcb          )
    ffoxs(ifoxs_tmp_hct           )=fflxs(iflxs_tmp_hct           )
    ffoxs(ifoxs_tcdc_mcl          )=fflxs(iflxs_tcdc_mcl          )
    ffoxs(ifoxs_pres_mct          )=fflxs(iflxs_pres_mct          )
    ffoxs(ifoxs_pres_mcb          )=fflxs(iflxs_pres_mcb          )
    ffoxs(ifoxs_tmp_mct           )=fflxs(iflxs_tmp_mct           )
    ffoxs(ifoxs_tcdc_lcl          )=fflxs(iflxs_tcdc_lcl          )
    ffoxs(ifoxs_pres_lct          )=fflxs(iflxs_pres_lct          )
    ffoxs(ifoxs_pres_lcb          )=fflxs(iflxs_pres_lcb          )
    ffoxs(ifoxs_tmp_lct           )=fflxs(iflxs_tmp_lct           )
    ffoxs(ifoxs_prate_sfc         )=fflxs(iflxs_prate_sfc         )
    ffoxs(ifoxs_cprat_sfc         )=fflxs(iflxs_cprat_sfc         )
    ffoxs(ifoxs_gflux_sfc         )=fflxs(iflxs_gflux_sfc         )
    ffoxs(ifoxs_land_sfc          )=fflxs(iflxs_land_sfc          )
    ffoxs(ifoxs_icec_sfc          )=fflxs(iflxs_icec_sfc          )
    ffoxs(ifoxs_icetk_sfc         )=fflxs(iflxs_icetk_sfc         )
    ffoxs(ifoxs_tmp_hag_2         )=fflxs(iflxs_tmp_hag_2         )
    ffoxs(ifoxs_spfh_hag_2        )=fflxs(iflxs_spfh_hag_2        )
    ffoxs(ifoxs_tmax_hag_2        )=fflxs(iflxs_tmax_hag_2        )
    ffoxs(ifoxs_tmin_hag_2        )=fflxs(iflxs_tmin_hag_2        )
    ffoxs(ifoxs_watr_sfc          )=fflxs(iflxs_watr_sfc          )
    ffoxs(ifoxs_pevpr_sfc         )=fflxs(iflxs_pevpr_sfc         )
    ffoxs(ifoxs_cwork_clm         )=fflxs(iflxs_cwork_clm         )
    ffoxs(ifoxs_hpbl_sfc          )=fflxs(iflxs_hpbl_sfc          )
    ffoxs(ifoxs_albdo_sfc         )=fflxs(iflxs_albdo_sfc         )
    ffoxs(ifoxs_tcdc_clm          )=fflxs(iflxs_tcdc_clm          )
    ffoxs(ifoxs_tcdc_cvc          )=fflxs(iflxs_tcdc_cvc          )
    ffoxs(ifoxs_pres_cvt          )=fflxs(iflxs_pres_cvt          )
    ffoxs(ifoxs_pres_cvb          )=fflxs(iflxs_pres_cvb          )
    ffoxs(ifoxs_tcdc_blc          )=fflxs(iflxs_tcdc_blc          )
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute accumulated precipitation and categorical precipitation types.
    ffoxs(ifoxs_apcp_sfc          )=3600*hrflx*fflxs(iflxs_prate_sfc         )
    ffoxs(ifoxs_acpcp_sfc         )=3600*hrflx*fflxs(iflxs_cprat_sfc         )
    ffoxs(ifoxs_crain_sfc         )=mod(iwx/8,2)
    ffoxs(ifoxs_cfrzr_sfc         )=mod(iwx/4,2)
    ffoxs(ifoxs_cicep_sfc         )=mod(iwx/2,2)
    ffoxs(ifoxs_csnow_sfc         )=mod(iwx/1,2)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute 2-meter relative humidity.
    if(lfoxs(ifoxs_rh_hag_2          )) then
      t2=fflxs(iflxs_tmp_hag_2         )
      q2=fflxs(iflxs_spfh_hag_2        )
      ps=fflxs(iflxs_pres_sfc          )
      p2=ps*(1.-2./con_rog/t2)
      pvs2=fpvs(real(t2,krealfp))
      qs2=con_eps*pvs2/(p2+con_epsm1*pvs2)
      ffoxs(ifoxs_rh_hag_2          )=1.e2*min(max(q2/qs2,0.),1.)
    endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Copy vector data.
    ffoxu(ifoxu_uflx_sfc          )=fflxu(iflxu_uflx_sfc          )
    ffoxv(ifoxv_vflx_sfc          )=fflxv(iflxv_vflx_sfc          )
    ffoxu(ifoxu_ugrd_hag_10       )=fflxu(iflxu_ugrd_hag_10       )
    ffoxv(ifoxv_vgrd_hag_10       )=fflxv(iflxv_vgrd_hag_10       )
    ffoxu(ifoxu_ugwd_sfc          )=fflxu(iflxu_ugwd_sfc          )
    ffoxv(ifoxv_vgwd_sfc          )=fflxv(iflxv_vgwd_sfc          )
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  subroutine calwxt1(lm,prec,t,q,p,pint,iwx)
!$$$  Subprogram documentation block
!
! Subprogram: calwxt1    Compute precipitation type
!   Prgmmr: Baldwin      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram computes precipitation type using
!           a decision tree approach that uses variables
!           such as integrated wet bulb temperature below freezing
!           and lowest layer temperature.
!           For more details, see Baldwin and Contorno preprint
!           from 13th Weather Analysis and Forecasting Conference
!           (or Baldwin et al. 10th NWP Conference preprint).
!
! Program history log:
!   1993-11-11  Baldwin
!   1994-09-30  Baldwin       set up new decision tree
!   1995-03-27  Iredell       modularized and vectorized
!   1999-10-18  Mark Iredell  fortran 90
!
! Usage:  call calwxt1(lm,prec,t,q,p,pint,iwx)
!   Input argument list:
!     lm           integer number of levels in this profile
!     prec         real precipitation (mm?)
!     t            real (lm) temperature (from top) (k)
!     q            real (lm) specific humidity (from top) (kg/kg)
!     p            real (lm) pressure (from top) (pa)
!     pint         real (lm+1) interface pressure (from top) (pa)
!   output arguments:
!     iwx          integer instantaneous weather type
!                  (the one's digit is for snow;
!                   the two's digit is for ice pellets;
!                   the four's digit is for freezing rain;
!                   the eight's digit is for rain.)
!
! Modules used:
!   funcphys       Physical functions
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   (ftdp)       compute dewpoint temperature
!   (ftcl)       compute lcl temperature
!   (fthe)       compute equivalent potential temperature
!   (ftma)       compute moist adiabat temperature
!
! Remarks: Weather type is only computed where there is precipitation.
!          All profiles must start at the top and go down.
!          For efficiency, inline all function calls.
!
! Attributes:
!   Language: Fortran 90
!
!$$$
    use funcphys
    use physcons
    implicit none
    integer,intent(in):: lm
    real,intent(in):: prec,t(lm),q(lm),p(lm),pint(lm+1)
    integer,intent(out):: iwx
    integer lice,l,l150,lfrz,lwrm,lone
    real tcold,twarm,tdchk
    real psm150,tv,dz150,surfw,surfc
    real dz(lm)
    real(krealfp) pv,pr,tdpd,tr,pk,tlcl,thelcl,twet(lm),qwet,areap4,areas8
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Initialize the weather type.
    iwx=0
!  If there is precipitation, find the coldest and warmest temperatures
!  in the saturated layer between 70 mb above ground and 500 mb.
!  Also find highest saturated layer in that range.
    if(prec.gt.0.) then
      tcold=t(lm)
      twarm=t(lm)
      lice=lm
      tdchk=2.
!  Try to find first 2 then 4 then 6 degree dewpoint depressions.
      do while(tcold.eq.t(lm).and.tdchk.le.6)
        do l=1,lm
          if(p(l).ge.500.e2.and.p(l).le.pint(lm+1)-70.e2) then
            pv=p(l)*q(l)/(con_eps-con_epsm1*q(l))
            tdpd=t(l)-ftdp(pv)
            if(tdpd.lt.tdchk) then
              tcold=min(tcold,t(l))
              twarm=max(twarm,t(l))
              lice=min(lice,l)
            endif
          endif
        enddo
        tdchk=tdchk+2
      enddo
!  Decide weather type if there is no cold saturated air.
      if(tcold.gt.con_t0c-4) then
!  Decide freezing rain if the surface is cold.
        if(t(lm).lt.con_t0c) then
          iwx=iwx+4
!  Otherwise decide rain.
        else
          iwx=iwx+8
        endif
!  If there is cold saturated air, then find lowest 150 mb level
!  and the lowest freezing level and warm level.
      else
        psm150=pint(lm)-150.e2
        l150=1
        lfrz=1
        lwrm=1
        do l=1,lm
          if(pint(l+1).le.psm150) l150=l+1
          if(t(l).lt.con_t0c) lfrz=l+1
          if(t(l).ge.twarm) lwrm=l+1
        enddo
!  Compute layer thickness and wet bulb temperature where needed.
        lone=min(lice,l150,lfrz,lwrm)
        do l=lone,lm
          if(pint(l).gt.0) then
            tv=t(l)*(1+con_fvirt*q(l))
            dz(l)=con_rog*tv*log(pint(l+1)/pint(l))
          else
            dz(l)=0
          endif
          pv=p(l)*q(l)/(con_eps-con_epsm1*q(l))
          tdpd=t(l)-ftdp(pv)
          if(tdpd.gt.0.) then
            pr=p(l)
            tr=t(l)
            pk=fpkap(pr)
            tlcl=ftlcl(tr,tdpd)
            thelcl=fthe(tlcl,pk*tlcl/tr)
            call stma(thelcl,pk,twet(l),qwet)
          else
            twet(l)=t(l)
          endif
        enddo
!  Integrate area of twet above -4c below highest saturated layer.
        areap4=0
        do l=lice,lm
          if(twet(l).gt.con_t0c-4) then
            areap4=areap4+(twet(l)-(con_t0c-4))*dz(l)
          endif
        enddo
!  Decide snow if there is scarce warm air.
        if(areap4.lt.3000.) then
          iwx=iwx+1
        else
!  Otherwise integrate net area of twet w.r.t. 0c in lowest 150 mb.
          l=l150
          tv=t(l)*(1+con_fvirt*q(l))
          dz150=con_rog*tv*log(pint(l+1)/psm150)
          areas8=(twet(l)-con_t0c)*dz150
          do l=l150+1,lm
            areas8=areas8+(twet(l)-con_t0c)*dz(l)
          enddo
!  Integrate area of twet above 0c below the freezing level.
          surfw=0
          do l=lfrz,lm
            if(twet(l).gt.con_t0c) then
              surfw=surfw+(twet(l)-con_t0c)*dz(l)
            endif
          enddo
!  Integrate area of twet below 0c below the warmest saturated level.
          surfc=0
          do l=lwrm,lm
            if(twet(l).lt.con_t0c) then
              surfc=surfc+(twet(l)-con_t0c)*dz(l)
            endif
          enddo
!  Decide ice pellets if there is yet plenty of cold air.
          if(surfc.lt.-3000..or.(areas8.lt.-3000..and.surfw.lt.50.)) then
            iwx=iwx+2
          else
!  Otherwise decide freezing rain if the surface is cold.
            if(t(lm).lt.con_t0c) then
              iwx=iwx+4
!  Otherwise decide rain.
            else
              iwx=iwx+8
            endif
          endif
        endif
      endif
    endif
  end subroutine
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine rsearch1(km1,z1,km2,z2,l2)
!$$$  subprogram documentation block
!
! subprogram:    rsearch1    search for a surrounding real interval
!   prgmmr: iredell    org: w/nmc23     date: 98-05-01
!
! abstract: this subprogram searches a monotonic sequences of real numbers
!   for intervals that surround a given search set of real numbers.
!   the sequences may be monotonic in either direction; the real numbers
!   may be single or double precision.
!
! program history log:
! 1999-01-05  mark iredell
! 2013-01-09  mark iredell  vanilla
!
! usage:    call rsearch1(km1,z1,km2,z2,l2)
!   input argument list:
!     km1    integer number of points in the sequence
!     z1     real (km1) sequence values to search
!            (z1 must be monotonic in either direction)
!     km2    integer number of points to search for
!     z2     real (km2) set of values to search for
!            (z2 need not be monotonic)
!
!   output argument list:
!     l2     integer (km2) interval locations from 0 to km1
!            (z2 will be between z1(l2) and z1(l2+1))
!
! remarks:
!   returned values of 0 or km1 indicate that the given search value
!   is outside the range of the sequence.
!
!   if a search value is identical to one of the sequence values
!   then the location returned points to the identical value.
!   if the sequence is not strictly monotonic and a search value is
!   identical to more than one of the sequence values, then the
!   location returned may point to any of the identical values.
!
!   if l2(k)=0, then z2(k) is less than the start point z1(1)
!   for ascending sequences (or greater than for descending sequences).
!   if l2(k)=km1, then z2(k) is greater than or equal to the end point
!   z1(km1) for ascending sequences (or less than or equal to for
!   descending sequences).  otherwise z2(k) is between the values
!   z1(l2(k)) and z1(l2(k+1)) and may equal the former.
!
! attributes:
!   language: fortran
!
!$$$
  implicit none
  integer,intent(in):: km1,km2
  real,intent(in):: z1(km1),z2(km2)
  integer,intent(out):: l2(km2)
  integer k1,k2
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  find the surrounding input interval for each output point.
  if(z1(1)<z1(km1)) then
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically ascending.
    do k2=1,km2
      l2(k2)=km1
      do k1=1,km1
        if(z2(k2)<z1(k1)) then
          l2(k2)=k1-1
          exit
        endif
      enddo
    enddo
  else
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input coordinate is monotonically descending.
    do k2=1,km2
      l2(k2)=km1
      do k1=1,km1
        if(z2(k2)>z1(k1)) then
          l2(k2)=k1-1
          exit
        endif
      enddo
    enddo
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine postx(id,im,jm,km,mt,mw,mta,mwa,fpos,fpou,fpov)
!$$$  Subprogram documentation block
!
! Subprogram: postx      Post subprogram to filter horizontally
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram horizontally filters pressure level posted fields.
!   The filtering is done in spectral space.  The filtering need not be a
!   straight spectral truncation, however, but can gradually filter over
!   a finite width in spectral space, which should ameliorate Gibbsing.
!   Also calculated are absolute vorticity and five-wave geopotential height.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call postx(id,im,jm,km,mt,mw,mta,mwa,fpos,fpou,fpov)
!   Input argument list:
!     id       integer data representation type
!     im       integer number of longitudes
!     jm       integer number of latitudes
!     mt       integer smoothing truncation
!     mw       integer smoothing width
!     mta      integer smoothing truncation for absolute vorticity
!     mwa      integer smoothing width for absolute vorticity
!     fpos     real (im,jm,km,npos) scalar fields
!     fpou     real (im,jm,km,npov) vector x-component fields
!     fpov     real (im,jm,km,npov) vector y-component fields
!   Output argument list:
!     fpos     real (im,jm,km,npos) scalar fields
!     fpou     real (im,jm,km,npov) vector x-component fields
!     fpov     real (im,jm,km,npov) vector y-component fields
!
! Modules used:
!   postgp_module  Shared data for postgp
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   sptezm         multiple scalar spectral transform
!   sptezmv        multiple vector spectral transform
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use postgp_module
  use physcons
  implicit none
  integer,intent(in):: id,im,jm,km,mt,mw,mta,mwa
  real,intent(inout):: fpos(im,jm,km,npos)
  real,intent(inout):: fpou(im,jm,km,npov),fpov(im,jm,km,npov)
  real s((mt+mw/2+1)*(mt+mw/2+2),km,7)
  integer mx,k,i,l,n
  real fac
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Transform to spectral
  mx=mt+mw/2
  call sptezm(0,mx,id,im,jm,km,s(1,1,1),fpos(1,1,1,ipos_hgt_prs),-1)
  call sptezmv(0,mx,id,im,jm,km,s(1,1,2),s(1,1,3),&
               fpou(1,1,1,ipou_ugrd_prs),fpov(1,1,1,ipov_vgrd_prs),-1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,4),fpos(1,1,1,ipos_tmp_prs),-1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,5),fpos(1,1,1,ipos_vvel_prs),-1)
!omp$ parallel do private(i,fac)
  do k=1,km
    i=0
    do l=0,mx
      do n=l,mx
!  Compute five-wave geopotential height
!       if(n.le.5) then              ! triangular 5
!       if(l.le.5.and.n-l.le.5) then ! rhomboidal 5
        if(l.le.5) then              ! zonal-wave 5
          s(i+1,k,6)=s(i+1,k,1)
          s(i+2,k,6)=s(i+2,k,1)
        else
          s(i+1,k,6)=0.
          s(i+2,k,6)=0.
        endif
!  Compute relative vorticity
        if(n.gt.mta-mwa/2) then
          fac=0.
          if(mwa.gt.0) fac=min(max(0.5+(mta-n)/float(mwa),0.),1.)
          s(i+1,k,7)=fac*s(i+1,k,3)
          s(i+2,k,7)=fac*s(i+2,k,3)
        else
          s(i+1,k,7)=s(i+1,k,3)
          s(i+2,k,7)=s(i+2,k,3)
        endif
!  Filter other fields
        if(n.gt.mt-mw/2) then
          fac=0.
          if(mw.gt.0) fac=min(max(0.5+(mt-n)/float(mw),0.),1.)
          s(i+1,k,1)=fac*s(i+1,k,1)
          s(i+2,k,1)=fac*s(i+2,k,1)
          s(i+1,k,2)=fac*s(i+1,k,2)
          s(i+2,k,2)=fac*s(i+2,k,2)
          s(i+1,k,3)=fac*s(i+1,k,3)
          s(i+2,k,3)=fac*s(i+2,k,3)
          s(i+1,k,4)=fac*s(i+1,k,4)
          s(i+2,k,4)=fac*s(i+2,k,4)
          s(i+1,k,5)=fac*s(i+1,k,5)
          s(i+2,k,5)=fac*s(i+2,k,5)
        endif
        i=i+2
      enddo
    enddo
!  Add planetary vorticity to get absolute vorticity
    s(3,k,7)=s(3,k,7)+2*con_omega/sqrt(1.5)
  enddo
!  Transform back to grid
  call sptezm(0,mx,id,im,jm,km,s(1,1,1),fpos(1,1,1,ipos_hgt_prs),+1)
  call sptezmv(0,mx,id,im,jm,km,s(1,1,2),s(1,1,3),&
               fpou(1,1,1,ipou_ugrd_prs),fpov(1,1,1,ipov_vgrd_prs),+1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,4),fpos(1,1,1,ipos_tmp_prs),+1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,5),fpos(1,1,1,ipos_vvel_prs),+1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,6),fpos(1,1,1,ipos_fwavh_prs),+1)
  call sptezm(0,mx,id,im,jm,km,s(1,1,7),fpos(1,1,1,ipos_absv_prs),+1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine postx1(id,im,jm,mt,mw,f)
!$$$  Subprogram documentation block
!
! Subprogram: postx1     Post subprogram to filter horizontally
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram horizontally filters a scalar field.
!   The filtering is done in spectral space.  The filtering need not be a
!   straight spectral truncation, however, but can gradually filter over
!   a finite width in spectral space, which should ameliorate Gibbsing.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call postx1(id,im,jm,mt,mw,f)
!   Input argument list:
!     id       integer data representation type
!     im       integer number of longitudes
!     jm       integer number of latitudes
!     mt       integer smoothing truncation
!     mw       integer smoothing width
!     fpos     real (im,jm) scalar field
!   Output argument list:
!     fpos     real (im,jm) scalar field
!
! Modules used:
!   postgp_module  Shared data for postgp
!
! Files included:
!   physcons.h     Physical constants
!
! Subprograms called:
!   sptez          scalar spectral transform
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use postgp_module
  implicit none
  integer,intent(in):: id,im,jm,mt,mw
  real,intent(inout):: f(im,jm)
  real s((mt+mw/2+1)*(mt+mw/2+2))
  integer mx,i,l,n
  real fac
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Transform to spectral
  mx=mt+mw/2
  call sptez(0,mx,id,im,jm,s,f,-1)
!  Filter
  i=0
  do l=0,mx
    do n=l,mx
      if(n.gt.mt-mw/2) then
        fac=0.
        if(mw.gt.0) fac=min(max(0.5+(mt-n)/float(mw),0.),1.)
        s(i+1)=fac*s(i+1)
        s(i+2)=fac*s(i+2)
      endif
      i=i+2
    enddo
  enddo
!  Transform back to grid
  call sptez(0,mx,id,im,jm,s,f,+1)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine kpfilt(nkplist,kplist,nfx,iptv,ipu,itl,il1,il2)
!$$$  Subprogram documentation block
!
! Subprogram: kpfilt     Filter output before packing and writing
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram makes the determination whether or not fields are
!   packed and written based on an input filter list.  If the parameter,
!   level type and level values are in the filter list, then the fields
!   will be written.  The filter list can have -1 as a wild card.
!
! Program history log:
!   2003-02-24  Mark Iredell
!
! Usage:  call kpfilt(nkplist,kplist,nfx,ipu,itl,il1,il2)
!   Input argument list:
!     nkplist  integer number of entries in the filter list
!     kplist   integer (4,nkplist) filter list iptv,ipu,itl,il1*256+il2 
!     nfx      integer local size of field decomposition
!     iptv     integer (nfx) parameter table version
!     ipu      integer (nfx) parameter identifier
!     itl      integer (nfx) type of level
!     il1      integer (nfx) level 1
!     il2      integer (nfx) level 2
!   Output argument list:
!     ipu      integer (nfx) parameter identifier
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  integer,intent(in):: nkplist,nfx
  integer,dimension(4,nkplist),intent(in):: kplist
  integer,dimension(nfx),intent(in):: iptv,itl,il1,il2
  integer,dimension(nfx),intent(inout):: ipu
  integer n,k,kpn(4)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Determine if each parameter is in the list
  field: do n=1,nfx
    kpn(1)=iptv(n)
    kpn(2)=ipu(n)
    kpn(3)=itl(n)
    kpn(4)=abs(il1(n)*256+il2(n))
    if(il2(n).lt.0) kpn(4)=kpn(4)+32768
    do k=1,nkplist
      if(all(kplist(:,k).eq.-1.or.kplist(:,k).eq.kpn)) cycle field
    enddo
    ipu(n)=0  ! not found
  enddo field
end subroutine
!-------------------------------------------------------------------------------
subroutine posto(comm,size,nf,nf1,nf2,nfm,nfx,ijmc,ijxc,ijoc,ijo,&
                 lc,fc,gc,lb,lv,lmw,&
                 ibms,iptv,ipu,ipv,itl,il1,il2,ifu,ip1,ip2,itr,&
                 idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
                 icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
                 ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
!$$$  Subprogram documentation block
!
! Subprogram: posto      Post subprogram to pack and write output
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram packs and writes the posted output data.
!   It also may first interpolate the data from a different compute grid.
!   The data is passed in decomposed over the compute grid, but is internally
!   transposed to the full horizontal dimension and partial field dimension
!   before any processing is done.  The packed data are written by all tasks
!   to a random access file, with coordination to avoid conflict.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call posto(comm,size,nf,nf1,nf2,nfm,nfx,ijmc,ijxc,ijoc,ijo,&
!                    lc,fc,gc,lb,lv,lmw,&
!                    ibms,iptv,ipu,ipv,itl,il1,il2,ifu,ip1,ip2,itr,&
!                    idrtc,ioc,joc,kgdsc,idrt,io,jo,kgdso,&
!                    icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi,&
!                    ids,omin,omax,mxbit,lupgb,iprev,nfw,iret)
!   Input argument list:
!     comm     integer message passing communicator
!     size     integer message passing size
!     nf       integer total number of fields
!     nf1      integer local first field index
!     nf2      integer local last field index
!     nfm      integer maximum size of field decomposition
!     nfx      integer local size of field decomposition
!     ijmc     integer maximum size of compute grid decomposition
!     ijxc     integer local size of compute grid decomposition
!     ijoc     integer total number of compute grid points
!     ijo      integer total number of output grid points
!     lc       logical*1 (nf,ijxc) bitmap on local compute grid
!     fc       real (nf,ijxc) data field 1 on local compute grid
!     gc       real (nf,ijxc) data field 2 on local compute grid
!     lb       integer bitmap switch
!     lv       integer vector switch
!     lmw      integer multiple write switch
!     ibms     integer (nfx) bitmap flag
!     iptv     integer (nfx) parameter table version
!     ipu      integer (nfx) parameter identifier for field 1
!     ipv      integer (nfx) parameter identifier for field 2
!     itl      integer (nfx) type of level
!     il1      integer (nfx) level 1
!     il2      integer (nfx) level 2
!     ifu      integer (nfx) forecast time unit
!     ip1      integer (nfx) time 1
!     ip2      integer (nfx) time 2
!     itr      integer (nfx) time representation
!     idrtc    integer compute grid data representation type
!     ioc      integer compute grid number of longitudes
!     joc      integer compute grid number of latitudes
!     kgdsc    integer (200) compute grid GDS
!     idrt     integer output grid data representation type
!     io       integer output grid number of longitudes
!     jo       integer output grid number of latitudes
!     kgdso    integer (200) output grid GDS
!     icen     integer center
!     igen     integer model generating code
!     igrido   integer center grid identifier
!     iyr      integer year
!     imo      integer month
!     idy      integer day
!     ihr      integer hour
!     icen2    integer subcenter
!     ienst    integer ensemble type
!     iensi    integer ensemble identification
!     ids      integer (255,255) decimal scaling
!     omin     real (255) minimum value in output data
!     omax     real (255) maximum value in output data
!     mxbit    integer maximum number of bits to pack data
!     lupgb    integer unit number of random access output file
!     iprev    integer number of bytes to written to output file
!   Output argument list:
!     iprev    integer number of bytes to written to output file
!     nfw      integer number of fields written
!     iret     integer return code
!
! Subprograms called:
!   mptran1l1      transpose grid decompositions
!   mptran1r4      transpose grid decompositions
!   ipolates       interpolate scalar fields
!   ipolatev       interpolate vector fields
!   makglgrb       make GRIB messages
!   mpfilli4       gather grid decomposition
!   bawrite        write random access record
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  integer(kint_mpi),intent(in):: comm,size
  integer,intent(in):: nf,nf1,nf2,nfm,nfx,ijmc,ijxc,ijoc,ijo
  logical*1,intent(in):: lc(nf,ijxc)
  real,intent(in):: fc(nf,ijxc),gc(nf,ijxc)
  integer,intent(in):: lb,lv,lmw
  integer,dimension(nfx),intent(in):: ibms,iptv,ipu,ipv,itl,il1,il2
  integer,dimension(nfx),intent(in):: ifu,ip1,ip2,itr
  integer,intent(in):: idrtc,ioc,joc,kgdsc(200),idrt,io,jo,kgdso(200)
  integer,intent(in):: icen,igen,igrido,iyr,imo,idy,ihr,icen2,ienst,iensi
  integer,intent(in):: ids(255,255)
  real,intent(in):: omin(255),omax(255)
  integer,intent(in):: mxbit,lupgb
  integer,intent(inout):: iprev
  integer,intent(out):: nfw,iret
  logical*1,allocatable:: lt(:,:)
  real,allocatable:: ft(:,:),gt(:,:)
  integer ibo(nf1:nf2)
  logical*1,allocatable:: lo(:,:)
  real,allocatable:: fo(:,:),go(:,:)
  integer nxo,nco(nf1:nf2),ngo(nf)
  character,allocatable:: co(:)
  logical lbm,luv
  integer ipopt(20),ko
  real rlat(ijo),rlon(ijo),crot(ijo),srot(ijo)
  integer iskip,iread,nread
  integer,allocatable:: igread(:)
  character,allocatable:: cgo(:)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Transpose distributed data to full horizontal fields
  nfw=0
  nxo=(1000+ijo*(mxbit+1)/8)*nfx
  lbm=lb.eq.1
  luv=lv.eq.1
  allocate(lt(ijoc,nf1:nf2))
  if(lbm) then
    call mptran1l1(comm,size,nfm,nf,nfx,nf,ijmc,ijxc,ijoc,ijoc,lc,lt)
  endif
  allocate(ft(ijoc,nf1:nf2))
  call mptran1r4(comm,size,nfm,nf,nfx,nf,ijmc,ijxc,ijoc,ijoc,fc,ft)
  if(luv) then
    allocate(gt(ijoc,nf1:nf2))
    call mptran1r4(comm,size,nfm,nf,nfx,nf,ijmc,ijxc,ijoc,ijoc,gc,gt)
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Interpolate to output grid
  ipopt(1)=1
  allocate(fo(ijo,nf1:nf2))
  allocate(lo(ijo,nf1:nf2))
  if(luv) allocate(go(ijo,nf1:nf2))
  if(idrt.eq.idrtc.and.io.eq.ioc.and.jo.eq.joc) then
    ibo=ibms
    if(lbm) then
      lo=lt
    else
      lo=.true.
    endif
    fo=ft
    if(luv) then
      go=gt
    endif
  elseif(luv) then
    call ipolatev(1,ipopt,kgdsc,kgdso,ijoc,ijo,nfx,ibms,lt,ft,gt,&
                  ko,rlat,rlon,crot,srot,ibo,lo,fo,go,iret)
    if(iret.ne.0) return
  else
    call ipolates(1,ipopt,kgdsc,kgdso,ijoc,ijo,nfx,ibms,lt,ft,&
                  ko,rlat,rlon,ibo,lo,fo,iret)
    if(iret.ne.0) return
  endif
  deallocate(lt)
  deallocate(ft)
  if(luv) deallocate(gt)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Make and write GRIB messages for field 1
  allocate(co(nxo))
  call makglgrb(icen,igen,igrido,iyr,imo,idy,ihr,0,0,icen2,&
                nfx,ibo,iptv,ipu,itl,il1,il2,ifu,ip1,ip2,itr,&
                ids,ienst,iensi,kgdso,omin,omax,&
                mxbit,ijo,ijo,nxo,lo,fo,nco,co,iret)
  if(.not.luv) deallocate(lo)
  deallocate(fo)
  if(iret.ne.0) return
  call mpfilli4(comm,size,1,1,1,nfm,nfx,nf,nco,ngo)
  iskip=iprev+sum(ngo(1:nf1-1))
  iread=sum(ngo(nf1:nf2))
  if(lmw.eq.0) then
    allocate(igread(size))
    call mpfilli4(comm,size,1,1,1,1,1,size,iread,igread)
    allocate(cgo(sum(igread)))
    call mpfillc1(comm,size,1,1,1,maxval(igread),iread,sum(igread),co,cgo)
    if(nf1.eq.1) then
      iread=sum(igread)
      call bawrite(lupgb,iskip,iread,nread,cgo)
    else
      iread=0
      nread=0
    endif
    deallocate(igread,cgo)
  else
    call bawrite(lupgb,iskip,iread,nread,co)
  endif
  if(nread.ne.iread) then
    iret=23
    return
  endif
  iprev=iprev+sum(ngo)
  nfw=count(ngo.gt.0)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Make and write GRIB messages for field 2
  if(luv) then
    call makglgrb(icen,igen,igrido,iyr,imo,idy,ihr,0,0,icen2,&
                  nfx,ibo,iptv,ipv,itl,il1,il2,ifu,ip1,ip2,itr,&
                  ids,ienst,iensi,kgdso,omin,omax,&
                  mxbit,ijo,ijo,nxo,lo,go,nco,co,iret)
    deallocate(lo)
    deallocate(go)
    if(iret.ne.0) return
    call mpfilli4(comm,size,1,1,1,nfm,nfx,nf,nco,ngo)
    iskip=iprev+sum(ngo(1:nf1-1))
    iread=sum(ngo(nf1:nf2))
    if(lmw.eq.0) then
      allocate(igread(size))
      call mpfilli4(comm,size,1,1,1,1,1,size,iread,igread)
      allocate(cgo(sum(igread)))
      call mpfillc1(comm,size,1,1,1,maxval(igread),iread,sum(igread),co,cgo)
      if(nf1.eq.1) then
        iread=sum(igread)
        call bawrite(lupgb,iskip,iread,nread,co)
      else
        iread=0
        nread=0
      endif
      deallocate(igread,cgo)
    else
      call bawrite(lupgb,iskip,iread,nread,co)
    endif
    call bawrite(lupgb,iskip,iread,nread,co)
    if(nread.ne.iread) then
      iret=23
      return
    endif
    iprev=iprev+sum(ngo)
    nfw=nfw+count(ngo.gt.0)
  endif
  deallocate(co)
end subroutine
!-------------------------------------------------------------------------------
subroutine makglgrb(icen,igen,igrid,iyr,imo,idy,ihr,ina,inm,icen2,&
                    na,ibms,iptv,ipu,itl,il1,il2,iftu,ip1,ip2,itr,&
                    ids,ienst,iensi,kgds,famin,famax,&
                    mxbit,ija,ida,idc,la,fa,ncg,cg,iret)
!$$$  Subprogram documentation block
!
! Subprogram: makglgrb   Make GRIB messages
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram creates several GRIB messages.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call makglgrb(icen,igen,igrid,iyr,imo,idy,ihr,ina,inm,icen2,&
!                       na,ibms,iptv,ipu,itl,il1,il2,iftu,ip1,ip2,itr,&
!                       ids,ienst,iensi,kgds,famin,famax,&
!                       mxbit,ija,ida,idc,la,fa,ncg,cg,iret)
!   Input argument list:
!     icen     integer center
!     igen     integer model generating code
!     igrid    integer center grid identifier
!     iyr      integer year
!     imo      integer month
!     idy      integer day
!     ihr      integer hour
!     ina      integer number in average
!     inm      integer number missing
!     icen2    integer subcenter
!     na       integer number of GRIB messages
!     ibms     integer (na) bitmap flag
!     iptv     integer (na) parameter table version
!     ipu      integer (na) parameter identifier for field 1
!              (if any are zero, that message is skipped)
!     itl      integer (na) type of level
!     il1      integer (na) level 1
!     il2      integer (na) level 2
!     iftu     integer (na) forecast time unit
!     ip1      integer (na) time 1
!     ip2      integer (na) time 2
!     itr      integer (na) time representation
!     ids      integer (255,255) decimal scaling
!     ienst    integer ensemble type
!     iensi    integer ensemble identification
!     kgds     integer (200) grid GDS
!     famin    real (255) minimum value in output data
!     famax    real (255) maximum value in output data
!     mxbit    integer maximum number of bits to pack data
!     ija      integer total number of grid points
!     ida      integer first dimension of input data
!     idc      integer size of output buffer in bytes
!     la       logical*1 (ida,na) bitmap to pack
!     fa       real (ida,na) data to pack
!   Output argument list:
!     ncg      integer (na) length of GRIB messages
!     cg       character (idc) output buffer
!     iret     integer return code
!
! Subprograms called:
!   makglpds       create global PDS
!   pakgb          make GRIB message
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  implicit none
  integer,intent(in):: icen,igen,igrid,iyr,imo,idy,ihr,ina,inm,icen2,na
  integer,intent(in):: ibms(na),iptv(na),ipu(na),itl(na),il1(na),il2(na)
  integer,intent(in):: iftu(na),ip1(na),ip2(na),itr(na)
  integer,intent(in):: ids(255,255),ienst,iensi,kgds(200)
  integer,intent(in):: mxbit,ija,ida,idc
  logical*1,intent(in):: la(ida,na)
  real,intent(in):: famin(255),famax(255)
  real,intent(in):: fa(ida,na)
  integer,intent(out):: ncg(na)
  character,intent(out):: cg(idc)
  integer ng,iskip,n,ip,ib,kbm,nbit,iret
  integer kpds(200),kens(5)
  real f(ija)
  real fmin,fmax
  character c1(1000+ija*(mxbit+1)/8)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  For each GRIB message, make PDS array, then make full GRIB message
  ng=1000+ija*(mxbit+1)/8
  iskip=0
  do n=1,na
    ncg(n)=0
    ip=ipu(n)
    if(ip.le.0) cycle
    ib=ibms(n)
    if(ib.eq.1) then
      if(.not.any(la(1:ija,n))) cycle
      if(all(la(1:ija,n))) ib=0
    endif
    call makglpds(iptv(n),icen,igen,igrid,ib,ip,itl(n),il1(n),il2(n),&
                  iyr,imo,idy,ihr,iftu(n),ip1(n),ip2(n),itr(n),&
                  ina,inm,icen2,ids(iptv(n),ip),kpds)
    kens(1)=1
    kens(2)=ienst
    kens(3)=iensi
    kens(4)=1
    kens(5)=255
    f=min(max(fa(1:ija,n),famin(ip)),famax(ip))
    call pakgb(ija,ng,0,mxbit,0,kpds,kgds,kens,la(1,n),f,&
               kbm,fmin,fmax,nbit,ncg(n),c1,iret)
    if(iret.ne.0) ncg(n)=0
    cg(iskip+1:iskip+ncg(n))=c1(1:ncg(n))
    iskip=iskip+ncg(n)
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine pakgb(kf,kg,minbit,maxbit,ibs,kpds,kgds,kens,lb,f,&
                 kbm,fmin,fmax,nbit,lg,g,iret)
!$$$  Subprogram documentation block
!
! Subprogram: pakgb      Make GRIB message
!   Prgmmr: Iredell      Org: np23        Date: 1999-10-18
!
! Abstract: This subprogram creates one GRIB message.
!
! Program history log:
!   1999-10-18  Mark Iredell
!
! Usage:  call pakgb(kf,kg,minbit,maxbit,ibs,kpds,kgds,kens,lb,f,&
!                    kbm,fmin,fmax,nbit,lg,g,iret)
!   Input argument list:
!     kf       integer length of data arrays
!     kg       integer length of output buffer
!     minbit   integer minimum number of bits to pack data
!     maxbit   integer maximum number of bits to pack data
!     ibs      integer binary scaling
!     kpds     integer (200) PDS parameters
!     kgds     integer (200) GDS parameters
!     kens     integer (5) ensemble parameters
!     lb       logical*1 (kf) bitmap to pack
!     fa       real (kf) data to pack
!   Output argument list:
!     kbm      integer number of packed valid data
!     fmin     real data packed data minimum
!     fmax     real data packed data maximum
!     nbit     integer number of bits used
!     lg       integer length of GRIB message
!     g        character (kg) GRIB message
!     iret     integer return code
!
! Subprograms called:
!   r63w72         map w3fi63 parameters onto w3fi72 parameters
!   getbit         get number of bits and round data
!   w3fi68         convert 25 word array to grib pds
!   pdsens         packs grib pds extension 41- for ensemble
!   w3fi72         pack grib
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  integer kpds(200),kgds(200),kens(5)
  logical*1 lb(kf)
  real f(kf)
  character g(kg)
  integer ibm(kf),ipds(200),igds(200),ibds(200)
  real fr(kf)
  character pds(400)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  get w3fi72 parameters
  call r63w72(kpds,kgds,ipds,igds)
  ibds=0
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  count valid data
  kbm=kf
  if(ipds(7).ne.0) then
    kbm=0
    do i=1,kf
      if(lb(i)) then
        ibm(i)=1
        kbm=kbm+1
      else
        ibm(i)=0
      endif
    enddo
    if(kbm.eq.kf) ipds(7)=0
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  get number of bits and round data
  if(kbm.eq.0) then
    do i=1,kf
      fr(i)=0.
    enddo
    fmin=0.
    fmax=0.
    nbit=0
  else
    call getbit(ipds(7),ibs,ipds(25),kf,ibm,f,fr,fmin,fmax,nbit)
    if(nbit.lt.minbit) then
      ibs2=ibs+(minbit-nbit)
      call getbit(ipds(7),ibs2,ipds(25),kf,ibm,f,fr,fmin,fmax,nbit)
    endif
    if(maxbit.gt.0) nbit=min(nbit,maxbit)
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  create product definition section
  call w3fi68(ipds,pds)
  if(ipds(24).eq.2) call pdsens(kens,kprob,xprob,kclust,kmembr,45,pds)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  pack grib data
  call w3fi72(0,fr,0,nbit,1,ipds,pds,1,255,igds,0,0,ibm,kf,ibds,kfo,g,lg,iret)
  if(iret.eq.0.and.lg.gt.kg) iret=98
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
