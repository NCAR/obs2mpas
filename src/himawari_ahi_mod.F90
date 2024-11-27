!----------------------------------------------------------------------
! Code adapted from 
!    https://github.com/NCAR/obs2ioda/blob/main/obs2ioda-v2/src/hsd.f90
!----------------------------------------------------------------------
!
module mod_himawari_ahi
!
! Purpose: Get Himawari Standard Data (HSD) lat, lon, bt, etc and pass it back to main
!          Convert HSD files to ioda-v1 format.
!          Currently only processes bands 7-16.
!
! input files:
!    (2) namelist.obs2model
!        &data_nml
!          list_files = 'flist.txt'     ! list of files
!          data_dir = '/data/himawari', ! path of the HSD files
!          data_id = 'HS',              ! prefix of the downloaded HSD files
!          sat_id = 'H08',              ! satellite ID
!          n_subsample = 1,             ! value use for thinning if write_iodav1 = .true.
!          write_iodav1 = .false.,      ! option to write out an iodav1 file (no superobbing)
!        /

  use netcdf_mod, only: open_netcdf_for_write, close_netcdf, &
                        def_netcdf_dims, def_netcdf_var, &
                        def_netcdf_end, put_netcdf_var

  !use utils_mod, only: get_julian_time

  use control_para !BJJ
  use utils_mod
   
  implicit none

  ! prefix for Cloud Property from ftp.ptree.jaxa.jp
  character(len=15), parameter :: CLP_id       = 'L2CLP'

  ! prefix for Cloud Mask (Binary Cloud Mask) from AWS NOAA
  ! prefix for Himawari-8 data on AWS (https://noaa-himawari8.s3.amazonaws.com/index.html)
  character(len=14), parameter :: BCM_id_old   = 'CLOUD_MASK'
  character(len=15), parameter :: HT_id_old    = 'CLOUD_HEIGHT'
  character(len=15), parameter :: Phase_id_old = 'CLOUD_PHASE'
  ! prefix for Himawari-9 data on AWS (https://noaa-himawari9.s3.amazonaws.com/index.html)
  character(len=14), parameter :: BCM_id_new   = 'CMSK'
  character(len=15), parameter :: HT_id_new    = 'CHGT'
  character(len=15), parameter :: Phase_id_new = 'CPHS'

  integer(i_kind) :: mmday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
  integer(i_kind), parameter :: npixel = 5500
  integer(i_kind), parameter :: nline  = 5500
  integer(i_kind), parameter :: nband = 10  ! number of infrared bands
  integer(i_kind)            :: band_start = 7
  integer(i_kind)            :: band_end   = 16

  real(r_kind)               :: brit(npixel, nline)
  real(r_double)             :: rlat, rlon, lon_diff, tmp1, theta1, theta2
  integer(i_kind)            :: ntotal, npix, nlin
  real(r_kind), allocatable  :: gsolzen(:,:)  ! satellite zenith angle (nx,ny)
  real(r_kind), allocatable  :: gsatzen(:,:)  ! solar zenith angle (nx,ny)
  
  integer(i_kind)    :: nfile
  integer(i_kind)    :: ifile
  character(len=512) :: ffname
  logical            :: isfile

  character(len=256) :: hs_list_file  ! the text file that contains a list of files to process
  character(len=256) :: data_dir
  character(len=18)  :: data_id
  character(len=3)   :: sat_id
  integer(i_kind)    :: n_subsample
  logical            :: write_iodav1

  character(len=256), allocatable :: hs_fnames(:)
  character(len=2)   :: finfo
  character(len=3)   :: fsat_id
  character(len=4)   :: region
  character(len=3)   :: resolution
  character(len=22)  :: file_time
  integer(i_kind)    :: iband, jday

  character(len=256) :: out_fname

  character(len=22), allocatable  :: scan_time(:) ! 2017-10-01T18:02:19.6Z
  integer(i_kind),   allocatable  :: fband_id(:)
  integer(i_kind),   allocatable  :: ftime_id(:)
  integer(i_kind),   allocatable  :: julianday(:)
  logical,           allocatable  :: valid(:), is_CLP(:), is_BCM(:), is_HT(:), is_Phase(:)

  character(len=22), allocatable  :: time_start(:) ! (ntime) 2017-10-01T18:02:19.6Z

  integer(i_kind) :: ncid, nf_status
  integer(i_kind) :: nx, ny
  integer(i_kind) :: it, ib, ii, i, j
  integer(i_kind) :: ntime
  integer(i_kind) :: t_index
  integer(i_kind) :: band_id

  real(r_kind)                    :: sdtb ! to be done
  logical                         :: found_time

  integer(i_kind), allocatable    :: cm_2d(:,:)   ! cloud_mask(nx,ny)

  type rad_type
     real(r_kind),    allocatable :: rad(:,:,:)  ! radiance(nband,nx,ny)
     real(r_kind),    allocatable :: bt(:,:,:)   ! brightness temperature(nband,nx,ny)
     integer(i_kind), allocatable :: qf(:,:,:)   ! quality flag(nband,nx,ny)
     real(r_kind),    allocatable :: sd(:)       ! std_dev(nband)
     integer(i_kind), allocatable :: cm(:,:)     ! cloud mask(nx,ny)
  end type rad_type
  type(rad_type),     allocatable :: rdata(:)    ! (ntime)

  type basic_info
    integer(i_byte)    :: headerNum      ! header block number = 1
    integer(i_short)   :: blockLen       ! block length = 282 bytes
    integer(i_short)   :: numHeader      ! total number of header blocks = 11
    integer(i_byte)    :: byteOrder      ! 0: little_endian, 1: big endian
    character(len=16)  :: satName        ! 'Himawari-8'
    character(len=16)  :: procCenter     ! 'MSC', 'OSK'
    character(len=4)   :: obsArea        ! 'FLDK'
    character(len=2)   :: dummy2
    integer(i_short)   :: hhnn           ! observation timeline
    real(r_double)     :: obsStartTime   ! Modified Julian Date
    real(r_double)     :: obsEndTime     ! Modified Julian Date
    real(r_double)     :: fileCreateTime ! Modified Julian Date
    integer(i_long)    :: totalHeaderLen
    integer(i_long)    :: dataLen
    integer(i_byte)    :: qcflag1
    integer(i_byte)    :: qcflag2
    integer(i_byte)    :: qcflag3
    integer(i_byte)    :: qcflag4
    character(len=32)  :: version
    character(len=128) :: fileName
    character(len=40)  :: dummy40
  end type basic_info

  type data_info
    integer(i_byte)    :: headerNum      ! header block number = 2
    integer(i_short)   :: blockLen       ! block length = 50 bytes
    integer(i_short)   :: bitPix         ! number of bits per pixel = 16
    integer(i_short)   :: nPix           ! number of columns (pixels east-west)
    integer(i_short)   :: nLin           ! number of lines (pixels north-south)
    integer(i_byte)    :: compression    ! 0: no compression (default), 1: gzip, 2: bzip2
    character(len=40)  :: dummy40
  end type data_info

  type proj_info
    integer(i_byte)    :: headerNum      ! header block number = 3
    integer(i_short)   :: blockLen       ! block length = 127 bytes
    real(r_double)     :: subLon         ! 140.7 degree
    integer(i_long)    :: cfac           ! column scaling factor
    integer(i_long)    :: lfac           ! line scaling factor
    real(r_single)     :: coff           ! column offset
    real(r_single)     :: loff           ! line offset
    real(r_double)     :: satDis         ! distance from earth's center to virtual satellite = 42164 km
    real(r_double)     :: eqtrRadius     ! earth's equatorial radius = 6378.1370 km
    real(r_double)     :: polrRadius     ! earth's polar radius = 6356.7523 km
    real(r_double)     :: projParam1     ! 0.00669438444
    real(r_double)     :: projParam2     ! 0.993305616
    real(r_double)     :: projParam3     ! 1.006739501
    real(r_double)     :: projParamSd    ! 1737122264
    integer(i_short)   :: resampleKind
    integer(i_short)   :: resampleSize
    character(len=40)  :: dummy40
  end type proj_info

  type navi_info
    integer(i_byte)    :: headerNum      ! header block number = 4
    integer(i_short)   :: blockLen       ! block length = 139 bytes
    real(r_double)     :: navTime        ! navigation information time in MJD
    real(r_double)     :: sspLon
    real(r_double)     :: sspLat
    real(r_double)     :: satDis
    real(r_double)     :: nadirLon
    real(r_double)     :: nadirLat
    real(r_double)     :: sunPos_x
    real(r_double)     :: sunPos_y
    real(r_double)     :: sunPos_z
    real(r_double)     :: moonPos_x
    real(r_double)     :: moonPos_y
    real(r_double)     :: moonPos_z
    character(len=40)  :: dummy40
  end type navi_info

  type calib_info
    integer(i_byte)    :: headerNum      ! header block number = 5
    integer(i_short)   :: blockLen       ! block length = 147 bytes
    integer(i_short)   :: bandNo
    real(r_double)     :: waveLen
    integer(i_short)   :: bitPix
    integer(i_short)   :: errorCount
    integer(i_short)   :: outCount
  ! count-radiance conversion equation
    real(r_double)     :: gain_cnt2rad
    real(r_double)     :: cnst_cnt2rad
  ! correction coefficient of sensor Planck function for converting radiance to brightness temperature
    real(r_double)     :: rad2btp_c0
    real(r_double)     :: rad2btp_c1
    real(r_double)     :: rad2btp_c2
  ! for converting brightness temperature to radiance
    real(r_double)     :: btp2rad_c0
    real(r_double)     :: btp2rad_c1
    real(r_double)     :: btp2rad_c2
    real(r_double)     :: lightSpeed
    real(r_double)     :: planckConst
    real(r_double)     :: bolzConst
    character(len=40)  :: dummy40
  end type calib_info

  type interCalib_info
    integer(i_byte)    :: headerNum      ! header block number = 6
    integer(i_short)   :: blockLen       ! block length = 259 bytes
    character(len=256) :: dummy256
  end type interCalib_info

  type segm_info
    integer(i_byte)    :: headerNum      ! header block number = 7
    integer(i_short)   :: blockLen       ! block length = 47 bytes
    integer(i_byte)    :: totalSegNum    ! total number of segments. 1: no divisions
    integer(i_byte)    :: segSeqNo       ! segment sequence number
    integer(i_short)   :: startLineNo    ! first line number of image segment
    character(len=40)  :: dummy40
  end type segm_info

  type naviCorr_info
    integer(i_byte)    :: headerNum      ! header block number = 8
    integer(i_short)   :: blockLen       ! block length
    real(r_single)     :: RoCenterColumn
    real(r_single)     :: RoCenterLine
    real(r_double)     :: RoCorrection
    integer(i_short)   :: correctNum
    integer(i_short), allocatable :: lineNo(:)      !(correctNum)
    real(r_single),   allocatable :: columnShift(:) !(correctNum)
    real(r_single),   allocatable :: lineShift(:)   !(correctNum)
    character(len=40)  :: dummy40
  end type naviCorr_info

  type obsTime_info
    integer(i_byte)    :: headerNum      ! header block number = 9
    integer(i_short)   :: blockLen       ! block length
    integer(i_short)   :: obsNum         ! number of observation times
    integer(i_short), allocatable :: lineNo(:) !(obsNum)
    real(r_double),   allocatable :: obsMJD(:) !(obsNum) observation time in MJD
    character(len=40)  :: dummy40
  end type obsTime_info

  type error_info
    integer(i_byte)    :: headerNum      ! header block number = 10
    integer(i_long)    :: blockLen       ! block length = 47
    integer(i_short)   :: errorNum       ! number of error information data = 0
  !  integer(i_short), allocatable :: lineNo(:)    !(errorNum)
  !  integer(i_short), allocatable :: errPixNum(:) !(errorNum)
    character(len=40)  :: dummy40
  end type error_info

  type dummy_info
    integer(i_byte)    :: headerNum      ! header block number = 11
    integer(i_short)   :: blockLen       ! block length = 259 bytes
    character(len=256) :: dummy256
  end type dummy_info

  type HSD_header
    type(basic_info)      :: basic
    type(data_info)       :: data
    type(proj_info)       :: proj
    type(navi_info)       :: navi
    type(calib_info)      :: calib
    type(interCalib_info) :: interCalib
    type(segm_info)       :: segm
    type(naviCorr_info)   :: navicorr
    type(obsTime_info)    :: obstime
    type(error_info)      :: error
    type(dummy_info)      :: dummy
  end type HSD_header

  type(HSD_header) :: header

  contains

subroutine Himawari_ReBroadcast_converter(glon_out, glat_out, F_out, varname_out, got_latlon_out)

   implicit none
   real(r_kind),      allocatable, intent(out) :: glon_out(:,:)
   real(r_kind),      allocatable, intent(out) :: glat_out(:,:)
   real(r_kind),      allocatable, intent(out) :: F_out(:,:,:)         ! (nx,ny, nfield), nfield=nfile, one file for each field
   character(len=64), allocatable, intent(out) :: varname_out(:)       ! (nfield)
   logical,           allocatable, intent(out) :: got_latlon_out(:,:)  ! (nx,ny)

   ! get namelist variables
   call get_namelist_vars(nfile, hs_fnames, hs_list_file, data_dir, data_id, sat_id, n_subsample, write_iodav1)

   allocate (ftime_id(nfile))
   allocate (scan_time(nfile))
   allocate (fband_id(nfile))
   allocate (julianday(nfile))
   allocate (valid(nfile))
   valid(:) = .false.

   allocate (is_CLP(nfile))
   allocate (is_BCM(nfile))
   allocate (is_Phase(nfile))
   allocate (is_HT(nfile))
   is_CLP(:)   = .false.
   is_BCM(:)   = .false.
   is_Phase(:) = .false.
   is_HT(:)    = .false.

   allocate (glat_out(npixel, nline))
   allocate (glon_out(npixel, nline))
   allocate (gsolzen(npixel, nline))
   allocate (gsatzen(npixel, nline))
   allocate (got_latlon_out(npixel, nline))
   allocate (F_out(npixel, nline, nfile))
   allocate (varname_out(nfile))
   glat_out(:,:)  = missing_r
   glon_out(:,:)  = missing_r
   gsolzen(:,:)   = missing_r
   gsatzen(:,:)   = missing_r
   F_out(:,:,:)   = missing_r
   varname_out(:) = ''

   ! parse the file list
   t_index = 0
   file_loop1: do ifile = 1, nfile
     ffname = trim(data_dir)//'/'//trim(hs_fnames(ifile))
     inquire(file=trim(ffname), exist=isfile)
     if ( .not. isfile ) then
        write(0,*) 'File not found: '//trim(ffname)
     else
        write(0,*) 'File found: '//trim(ffname)         
     end if

     ! retrieve some basic info from the filename itself
     call decode_himawari_name(trim(hs_fnames(ifile)), finfo, fband_id(ifile), fsat_id, scan_time(ifile), region, resolution, julianday(ifile), is_CLP(ifile), is_BCM(ifile), is_Phase(ifile), is_HT(ifile))

     if ( fsat_id /= sat_id ) then
        cycle file_loop1
     end if

     if ( .not. ( is_CLP(ifile) .or. is_BCM(ifile) .or. is_Phase(ifile) &
                .or. is_HT(ifile) ) ) then
        ! id of the file name must match specified data_id in namelist
        if ( finfo /= data_id ) then
           write(0,*) 'Satellite ID from namelist /= Satellite ID in flist.txt'
           cycle file_loop1
        else
           ! only process band 7-16
           if ( fband_id(ifile) < band_start .or. fband_id(ifile) > band_end ) then
              write(0,*) 'Infrared band ', fband_id(ifile), ' NOT supported'
              cycle file_loop1
           end if
        end if
     end if

     valid(ifile) = .true.

     ! group files of the same scan time
     if ( t_index == 0 ) then
        t_index = t_index + 1
        ftime_id(ifile) = t_index
     else
        found_time = .false.
        find_time_loop: do ii = ifile-1, 1, -1
           if ( valid(ii) ) then
              if ( scan_time(ifile) == scan_time(ii) ) then
                 ftime_id(ifile) = ftime_id(ii)
                 found_time = .true.
                 exit find_time_loop
              end if
           end if
        end do find_time_loop
        if ( .not. found_time ) then
           t_index = t_index + 1
           ftime_id(ifile) = t_index
        end if
     end if
     ntime = t_index

   end do file_loop1

   if ( ntime <= 0 ) then
      write(0,*) 'ntime = ', ntime
      write(0,*) 'No valid files found from hs_list_file '//trim(hs_list_file)
      stop
   end if

   allocate (time_start(ntime))
   allocate (rdata(ntime))

   file_loop2: do ifile = 1, nfile

     if ( valid(ifile) ) then
     
       ffname = trim(data_dir)//'/'//trim(hs_fnames(ifile))

       if ( .not. ( is_CLP(ifile) .or. is_BCM(ifile) .or. is_Phase(ifile) &
                  .or. is_HT(ifile) ) ) then
         call read_HSD(ffname, fband_id(ifile), fsat_id, julianday(ifile), glon_out, glat_out, brit, gsolzen, gsatzen, got_latlon_out)
         write(15,*) ifile, scan_time(ifile), julianday(ifile)
         F_out(:,:,ifile) = brit(:,:)

         it = ftime_id(ifile)
         ib = fband_id(ifile)

         if ( .not. allocated(rdata(it)%rad) ) allocate (rdata(it)%rad(nband, npixel, nline))
         if ( .not. allocated(rdata(it)%bt) )  allocate (rdata(it)%bt(nband, npixel, nline))
         if ( .not. allocated(rdata(it)%qf) )  allocate (rdata(it)%qf(nband, npixel, nline))
         if ( .not. allocated(rdata(it)%sd) )  allocate (rdata(it)%sd(nband))

         rdata(it)%bt(ib-band_start+1,:,:)  = brit(:,:)
         rdata(it)%rad(ib-band_start+1,:,:) = missing_r !rad_2d(i,j)
         rdata(it)%qf(ib-band_start+1,:,:)  = missing_r !qf_2d(i,j)
         rdata(it)%sd(ib-band_start+1)      = missing_r !sdtb

         write(varname_out(ifile),"(A,I2.2)") 'BT_'//fsat_id//'C', ib

       else if ( is_CLP(ifile) ) then
         nf_status = nf_OPEN(trim(ffname), nf_NOWRITE, ncid)
         if ( nf_status == 0 ) then
            write(0,*) 'Reading '//trim(ffname)
         else
            write(0,*) 'ERROR reading '//trim(ffname)
            cycle file_loop2
         end if
         call read_GRB_dims(ncid, 'latitude', 'longitude', nx, ny)
         allocate (cm_2d(nx, ny))
         call read_CLP(ncid, nx, ny, cm_2d)
         if ( .not. allocated(rdata(it)%cm) )  allocate (rdata(it)%cm(nx,ny))
         rdata(it)%cm(:,:) = cm_2d(:,:)

         !BJJ copy to output array: use ifile index
         F_out(:,:,ifile) = cm_2d(:,:)

         varname_out(ifile) = 'BCM_'//fsat_id

       else if ( is_BCM(ifile) ) then
         nf_status = nf_OPEN(trim(ffname), nf_NOWRITE, ncid)
         if ( nf_status == 0 ) then
            write(0,*) 'Reading '//trim(ffname)
         else
            write(0,*) 'ERROR reading '//trim(ffname)
            cycle file_loop2
         end if
         call read_GRB_dims(ncid, 'Rows', 'Columns', nx, ny)
         allocate (cm_2d(nx, ny))
         call read_L2_BCM(ncid, nx, ny, cm_2d, time_start(it))
         if ( is_empty_string(time_start(it)) ) then
            continue
         else if ( time_start(it) /= scan_time(ifile) ) then
            write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
            cycle file_loop2
         end if
         if ( .not. allocated(rdata(it)%cm) )  allocate (rdata(it)%cm(nx,ny))
         rdata(it)%cm(:,:) = cm_2d(:,:)

         !BJJ copy to output array: use ifile index
         F_out(:,:,ifile) = cm_2d(:,:)

         varname_out(ifile) = 'BCM_'//fsat_id

       else if ( is_HT(ifile) .or. is_Phase(ifile) ) then
       else
         write(0,*) 'ERROR: something is wrong. Check the files'
         stop
       end if

       nf_status = nf_CLOSE(ncid)

     end if

   end do file_loop2

   if ( allocated(cm_2d) )  deallocate(cm_2d)

   ! write IODAv1 file
   if ( write_iodav1 ) then
      do it = 1, ntime
         out_fname = trim(data_id)//'_'//fsat_id//'_'//scan_time(it)//'.nc4'
         write(0,*) 'Writing ', trim(out_fname)
         if ( allocated(rdata(it)%cm) ) then
            call output_iodav1(trim(out_fname), scan_time(it),                       &
                               npixel, nline, nband, n_subsample,                    &
                               got_latlon_out, glat_out, glon_out, gsatzen, gsolzen, &
                               rdata(it)%bt, rdata(it)%qf, rdata(it)%sd, rdata(it)%cm)
         else
            call output_iodav1(trim(out_fname), scan_time(it),                       &
                               npixel, nline, nband, n_subsample,                    &
                               got_latlon_out, glat_out, glon_out, gsatzen, gsolzen, &
                               rdata(it)%bt, rdata(it)%qf, rdata(it)%sd)
         end if
      end do
   end if

   do it = 1, ntime
      if ( allocated(rdata(it)%rad) ) deallocate (rdata(it)%rad)
      if ( allocated(rdata(it)%bt)  ) deallocate (rdata(it)%bt)
      if ( allocated(rdata(it)%qf)  ) deallocate (rdata(it)%qf)
      if ( allocated(rdata(it)%cm)  ) deallocate (rdata(it)%cm)
   end do
   if ( allocated(gsatzen) ) deallocate(gsatzen)
   if ( allocated(gsolzen) ) deallocate(gsolzen)
   deallocate(rdata)
   deallocate(hs_fnames)
   deallocate(ftime_id)
   deallocate(scan_time)
   deallocate(fband_id)
   deallocate(julianday)
   deallocate(valid)

   deallocate(is_CLP)
   deallocate(is_BCM)
   deallocate(is_Phase)
   deallocate(is_HT)

end subroutine Himawari_ReBroadcast_converter

subroutine read_CLP(ncid, nx, ny, cm)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   integer(i_kind),   intent(inout) :: cm(nx,ny)
   integer(i_short),  allocatable   :: itmp_short_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_kind)                  :: imiss = -999
   integer(i_short)                 :: ifill

   continue

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_short_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'CLTYPE', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   cm(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         if ( (itmp_short_2d(i,j) /= ifill) .or. (itmp_short_2d(i,j) /= 10) ) then
            if (itmp_short_2d(i,j) == 0 ) then ! clear pixel
               cm(i,j) = 0
            else
               cm(i,j) = 1
            end if
         end if
      end do
   end do
   deallocate(itmp_short_2d)

   return
end subroutine read_CLP

subroutine read_L2_BCM(ncid, nx, ny, cm, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   integer(i_kind),   intent(inout) :: cm(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: qf(nx,ny)
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'CloudMaskQualFlag', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   qf(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         qf(i,j) = itmp_byte_2d(i,j)
      end do
   end do
   deallocate(itmp_byte_2d)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'CloudMaskBinary', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   cm(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         if ( qf(i,j) == 0 ) then ! good quality
            cm(i,j) = itmp_byte_2d(i,j)
         end if
      end do
   end do
   deallocate(itmp_byte_2d)

   return
end subroutine read_L2_BCM

subroutine read_HSD(ffname, iband, satid, jday, longitude, latitude, brit, solzen, satzen, valid)

  implicit none

  character(len=256),intent(in)  :: ffname
  integer(i_kind),   intent(in)  :: iband
  character(len=3),  intent(in)  :: satid
  integer(i_kind),   intent(in)  :: jday
  real(r_single),   intent(out)  :: longitude(npixel, nline)
  real(r_single),   intent(out)  :: latitude(npixel, nline)
  real(r_single),   intent(out)  :: brit(npixel, nline)
  real(r_single),   intent(out)  :: solzen(npixel, nline)
  real(r_single),   intent(out)  :: satzen(npixel, nline)
  logical,          intent(out)  :: valid(npixel, nline)
   
  integer(i_short), allocatable  :: idata(:)
  integer(i_kind) :: numCorrect
  integer(i_kind) :: numObs
  integer(i_kind) :: ipixel, iline
  integer(i_kind) :: startLine, endLine
  integer(i_kind) :: radcount, i, ii, jj, ij, iv
  integer(i_kind) :: ierr
  integer(i_kind) :: nlocs, nvars, iloc
  integer(i_kind) :: ihh, imm, idd, flength, rvalue, offset
  integer(i_kind) :: iunit = 21
  real(r_double)  :: lon, lat
  real(r_double)  :: radiance, tbb
  real(r_double)  :: lon_sat, h_sat, r_eq
  integer(i_kind) :: nodivisionsegm = 1
  ! end of declaration
  continue

  longitude(:,:) = missing_r
  latitude(:,:)  = missing_r
  brit(:,:)      = missing_r
  solzen(:,:)    = missing_r
  satzen(:,:)    = missing_r
  valid(:,:)     = .false.

  open(iunit, file=trim(ffname), form='unformatted', action='read', access='stream', status='old', convert='little_endian')
  print*,'Reading from ', trim(ffname)

  read(iunit) header%basic%headerNum, &
              header%basic%blockLen, &
              header%basic%numHeader, &
              header%basic%byteOrder, &
              header%basic%satName, &
              header%basic%procCenter, &
              header%basic%obsArea, &
              header%basic%dummy2, &
              header%basic%hhnn, &
              header%basic%obsStartTime, &
              header%basic%obsEndTime, &
              header%basic%fileCreateTime, &
              header%basic%totalHeaderLen, &
              header%basic%dataLen, &
              header%basic%qcflag1, &
              header%basic%qcflag2, &
              header%basic%qcflag3, &
              header%basic%qcflag4, &
              header%basic%version, &
              header%basic%fileName, &
              header%basic%dummy40
  read(iunit) header%data%headerNum, &
              header%data%blockLen,&
              header%data%bitPix, &
              header%data%nPix, &
              header%data%nLin, &
              header%data%compression, &
              header%data%dummy40
  read(iunit) header%proj%headerNum, &
              header%proj%blockLen, &
              header%proj%subLon, &
              header%proj%cfac, &
              header%proj%lfac, &
              header%proj%coff, &
              header%proj%loff, &
              header%proj%satDis, &
              header%proj%eqtrRadius, &
              header%proj%polrRadius, &
              header%proj%projParam1, &
              header%proj%projParam2, &
              header%proj%projParam3, &
              header%proj%projParamSd, &
              header%proj%resampleKind, &
              header%proj%resampleSize, &
              header%proj%dummy40
  read(iunit) header%navi%headerNum, &
              header%navi%blockLen, &
              header%navi%navTime, &
              header%navi%sspLon, &
              header%navi%sspLat, &
              header%navi%satDis, &
              header%navi%nadirLon, &
              header%navi%nadirLat, &
              header%navi%sunPos_x, &
              header%navi%sunPos_y, &
              header%navi%sunPos_z, &
              header%navi%moonPos_x, &
              header%navi%moonPos_y, &
              header%navi%moonPos_z, &
              header%navi%dummy40
  read(iunit) header%calib%headerNum, &
              header%calib%blockLen, &
              header%calib%bandNo, &
              header%calib%waveLen, &
              header%calib%bitPix, &
              header%calib%errorCount, &
              header%calib%outCount, &
              header%calib%gain_cnt2rad, &
              header%calib%cnst_cnt2rad, &
              header%calib%rad2btp_c0, &
              header%calib%rad2btp_c1, &
              header%calib%rad2btp_c2, &
              header%calib%btp2rad_c0, &
              header%calib%btp2rad_c1, &
              header%calib%btp2rad_c2, &
              header%calib%lightSpeed, &
              header%calib%planckConst, &
              header%calib%bolzConst, &
              header%calib%dummy40
  read(iunit) header%interCalib%headerNum, &
              header%interCalib%blockLen, &
              header%interCalib%dummy256
  read(iunit) header%segm%headerNum, &
              header%segm%blockLen, &
              header%segm%totalSegNum, &
              header%segm%segSeqNo, &
              header%segm%startLineNo, &
              header%segm%dummy40
  read(iunit) header%navicorr%headerNum, &
              header%navicorr%blockLen, &
              header%navicorr%RoCenterColumn, &
              header%navicorr%RoCenterLine, &
              header%navicorr%RoCorrection, &
              header%navicorr%correctNum
  if ( header%navicorr%correctNum > 0 ) then
     numCorrect = header%navicorr%correctNum
     allocate(header%navicorr%lineNo(numCorrect))
     allocate(header%navicorr%columnShift(numCorrect))
     allocate(header%navicorr%lineShift(numCorrect))
  end if
  read(iunit) header%navicorr%lineNo(numCorrect), &
              header%navicorr%columnShift(numCorrect), &
              header%navicorr%lineShift(numCorrect), &
              header%navicorr%dummy40
  rewind(iunit)
  read(iunit) header%obstime%headerNum, &
              header%obstime%blockLen, &
              header%obstime%obsNum
  if ( header%obsTime%obsNum > 0 ) then
     numObs = header%obsTime%obsNum
     allocate(header%obsTime%lineNo(numObs))
     allocate(header%obsTime%obsMJD(numObs))
  end if
  read(iunit) header%obstime%lineNo, &
              header%obstime%obsMJD, &
              header%obstime%dummy40
  read(iunit) header%error%headerNum, &
              header%error%blockLen, &
              header%error%errorNum, &
              header%error%dummy40
  read(iunit) header%dummy%headerNum, &
              header%dummy%blockLen, &
              header%dummy%dummy256
  npix = header%data%nPix
  nlin = header%data%nLin
  ntotal = npix * nlin
  allocate(idata(ntotal))

  inquire(file=trim(ffname), size=flength)
  ! Offset relative to the beginning of the file
  offset = flength - (npix * nlin *2)
  ! Reposition the file to the offset value for reading
  call fseek(iunit, offset, 0, rvalue)
  read(iunit)idata(:)

  startLine = header%segm%startLineNo
  endLine = startLine + header%data%nLin - 1
  lon_sat = header%proj%subLon * deg2rad
  h_sat =  header%proj%satDis * 1000.0
  r_eq = header%proj%eqtrRadius * 1000.0
  do jj = 1, header%data%nLin
     do ii = 1, header%data%nPix

        tbb = missing_r
        lon = missing_r
        lat = missing_r

        ij = ii + (jj-1) * header%data%nPix
        iline = header%segm%startLineNo + jj - 1
        ipixel = ii

        if ( abs(latitude(ipixel, iline)) > 90.0 .or. &
             abs(longitude(ipixel, iline)) > 360.0 ) then
           call pixlin_to_lonlat(ipixel, iline, lon, lat, ierr)
           if ( ierr == 0 ) then
              valid(ipixel, iline) = .true.
              latitude(ipixel, iline) = lat
              longitude(ipixel, iline) = lon
              call calc_solar_zenith_angle_h( &
                 latitude(ipixel, iline), &
                 longitude(ipixel, iline), &
                 ihh, imm, jday, &
                 solzen(ipixel, iline))
              ! calculate geostationary satellite zenith angle
              rlat = lat * deg2rad ! in radian
              rlon = lon * deg2rad ! in radian
              lon_diff = abs(rlon-lon_sat)
              tmp1 = sqrt((2.0*r_eq*sin(lon_diff/2.)-r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2 &
                     +(2.0*r_eq*sin(rlat/2.))**2-(r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2)
              theta1 = 2.0*asin(tmp1/r_eq/2.)
              theta2 = atan(r_eq*sin(theta1)/((h_sat-r_eq)+r_eq*(1.0-sin(theta1))))
              satzen(ipixel, iline) = (theta1+theta2) * rad2deg
              if ( satzen(ipixel, iline) > 65.0 ) then
                 valid(ipixel, iline) = .false.
              end if
           end if
        end if

        radcount = idata(ij)
        if ( radcount /= header%calib%outCount .and. &
             radcount /= header%calib%errorCount .and. &
             radcount > 0 ) then
           ! convert count value to radiance
           radiance = radcount * header%calib%gain_cnt2rad + &
                      header%calib%cnst_cnt2rad
           radiance = radiance * 1000000.0  ! [ W/(m^2 sr micro m)] => [ W/(m^2 sr m)]
           ! convert radiance to physical value
           call hisd_radiance_to_tbb(radiance, tbb)
           ! visible or near infrared band
           brit(ipixel, iline) = tbb
        end if

     end do ! pixel
  end do ! line

  deallocate(idata)
  if ( header%obsTime%obsNum > 0 ) then
     deallocate(header%obsTime%lineNo)
     deallocate(header%obsTime%obsMJD)
  end if
  if ( header%navicorr%correctNum > 0 ) then
     deallocate(header%navicorr%lineNo)
     deallocate(header%navicorr%columnShift)
     deallocate(header%navicorr%lineShift)
  end if
  close(iunit)

  ! additional info for writing ioda at MPAS mesh !BJJ
  write(15,*) lon_sat, r_eq, h_sat

end subroutine read_HSD

subroutine pixlin_to_lonlat(pix, lin, lon, lat, ierr)

 implicit none

 integer(i_kind), intent(in)  :: pix, lin
 real(r_double),  intent(out) :: lon, lat
 integer(i_kind), intent(out) :: ierr

 real(r_double) :: SCLUNIT = 2.0**(-16)
 real(r_double) :: x, y
 real(r_double) :: Sd, Sn, S1, S2, S3, Sxy
 real(r_double) :: c, l

 ! initialize
 lon = missing_r
 lat = missing_r
 ierr = 0

 ! pix, lin -> c, l
 c = float(pix)
 l = float(lin)

 ! intermediate coordinates (x,y)
 ! Global Specification 4.4.4 Scaling Function
 ! https://www.cgms-info.org/wp-content/uploads/2021/10/cgms-lrit-hrit-global-specification-(v2-8-of-30-oct-2013).pdf
 !    c = COFF + nint(x * 2^-16 * CFAC)
 !    l = LOFF + nint(y * 2^-16 * LFAC)
 ! The intermediate coordinates (x,y) are as follows :
 !    x = (c -COFF) / (2^-16 * CFAC)
 !    y = (l -LOFF) / (2^-16 * LFAC)
 !    SCLUNIT = 2^-16
 x = deg2rad * ( c - header%proj%coff) / ( SCLUNIT * header%proj%cfac)
 y = deg2rad * ( l - header%proj%loff) / ( SCLUNIT * header%proj%lfac)

 ! longtitude,latitude
 ! Global Specification 4.4.3.2
 ! The invers projection function is as follows :
 !   lon = arctan(S2/S1) + sub_lon
 !   lat = arctan( (Req^2/Rpol^2) * S3 / Sxy )
 !
 ! Thererin the variables S1,S2,S3,Sxy are as follows :
 !    S1  = Rs - Sn * cos(x) * cos(y)
 !    S2  = Sn * sin(x) * cos(y)
 !    S3  =-Sn * sin(y)
 !    Sxy = sqrt(S1^2 + S2^2)
 !    Sn  =(Rs * cos(x) * cos(y) - Sd ) /
 !         (cos(y) * cos(y) + (Req^2/Rpol^2) * sin(y) * sin(y))
 !    Sd  =sqrt( (Rs * cos(x) * cos(y))^2
 !               - ( cos(y) * cos(y) + (Req^2/Rpol^2) * sin(y) * sin(y) )
 !               * (Rs^2 - Req^2)
 ! The variables Rs,Rpol,Req,(Req^2/Rpol^2),(Rs^2 - Req^2) are as follows :
 !    Rs  : distance from Earth center to satellite= head->proj->satDis
 !    Rpol: polar radius of the Earth              = head->proj->polrRadius
 !    Req : equator raidus of the Earth            = head->proj->eqtrRadius
 !    (Req^2/Rpol^2)                               = head->proj->projParam3
 !    (Rs^2 - Req^2)                               = head->proj->projParamSd
 Sd = (header%proj%satDis * cos(x) * cos(y)) * &
      (header%proj%satDis * cos(x) * cos(y)) - &
      (cos(y) * cos(y) + header%proj%projParam3 * sin(y) * sin(y)) * &
       header%proj%projParamSd
 if ( Sd < 0 ) then
    !write(*,*) 'Error in Sd'
    ierr = -1
    return
 else
    Sd = sqrt(Sd)
 end if
 Sn = (header%proj%satDis * cos(x) * cos(y) -Sd) / &
      (cos(y) * cos(y) + header%proj%projParam3 * sin(y) * sin(y))
 S1 = header%proj%satDis - (Sn * cos(x) * cos(y))
 S2 = Sn * sin(x) * cos(y)
 S3 =-Sn * sin(y)
 Sxy=sqrt( S1 * S1 + S2 * S2)

 lon = rad2deg * atan2(S2,S1) + header%proj%subLon
 lat = rad2deg * atan(header%proj%projParam3 * S3 / Sxy)

 ! check longtitude
 if ( lon >  180.0 ) lon = lon - 360.0
 if ( lon < -180.0 ) lon = lon + 360.0

 return
end subroutine pixlin_to_lonlat

subroutine hisd_radiance_to_tbb (radiance, tbb)

 implicit none

 real(r_double), intent(in)  :: radiance
 real(r_double), intent(out) :: tbb

 real(r_double) :: lambda
 real(r_double) :: planck_c1
 real(r_double) :: planck_c2

 real(r_double) :: effective_temperature

 ! central wave length
 lambda = header%calib%waveLen / 1000000.0 ! [micro m] => [m]

 ! radiance = radiance * 1000000.0  ! [ W/(m^2 sr micro m)] => [ W/(m^2 sr m)]

 ! planck_c1 = (2 * h * c^2 / lambda^5)
 planck_c1 = 2.0 * header%calib%planckConst *  &
             header%calib%lightSpeed ** 2 / &
             lambda ** 5

 ! planck_c2 = (h * c / k / lambda )
 planck_c2 = header%calib%planckConst * header%calib%lightSpeed / &
             header%calib%bolzConst / lambda

 if ( radiance > 0 ) then
    effective_temperature = planck_c2 / &
                            log( (planck_c1 / radiance ) + 1.0 )
    tbb = header%calib%rad2btp_c0 + &
          header%calib%rad2btp_c1 * effective_temperature + &
          header%calib%rad2btp_c2 * effective_temperature ** 2
 else
    tbb = missing_r
 end if
 return
end subroutine hisd_radiance_to_tbb

subroutine calc_solar_zenith_angle_h(xlat, xlon, gmt, minute, julian, solzen)

! the calulcation is adapted from subroutines radconst and calc_coszen in
! WRF phys/module_radiation_driver.F

 implicit none

 real(r_single),  intent(in)    :: xlat, xlon
 integer(i_kind), intent(in)    :: gmt, minute, julian
 real(r_single),  intent(inout) :: solzen

 real(r_single) :: obliq = 23.5
 real(r_single) :: deg_per_day = 360.0/365.0
 real(r_single) :: slon   ! longitude of the sun
 real(r_single) :: declin ! declination of the sun
 real(r_single) :: hrang, da, eot, xt, tloctm, rlat

 ! initialize to missing values
 solzen = missing_r

 ! calculate longitude of the sun from vernal equinox
 if ( julian >= 80 ) slon = (julian - 80 ) * deg_per_day
 if ( julian <  80 ) slon = (julian + 285) * deg_per_day

 declin = asin(sin(obliq*deg2rad)*sin(slon*deg2rad)) ! in radian

 da = 6.2831853071795862*(julian-1)/365.
 eot = (0.000075+0.001868*cos(da)-0.032077*sin(da) &
        -0.014615*cos(2.0*da)-0.04089*sin(2.0*da))*(229.18)
 xt = gmt + (minute + eot)/60.0

 if ( abs(xlon) > 360.0 .or. abs(xlat) > 90.0 ) return
 tloctm = xt + xlon/15.0
 hrang = 15.0*(tloctm-12.0) * deg2rad
 rlat = xlat * deg2rad
 solzen = acos( sin(rlat)*sin(declin) + &
                cos(rlat)*cos(declin)*cos(hrang) )
 solzen = solzen * rad2deg

 return
end subroutine calc_solar_zenith_angle_h

subroutine set_ahi_obserr(name_inst, nchan, obserrors)
   implicit none

   character(len=*), intent(in)  :: name_inst  ! instrument name
   integer(i_kind),  intent(in)  :: nchan      ! channel number
   real(r_kind),     intent(out) :: obserrors(nchan)
   obserrors(:) = missing_r
   if ( name_inst(1:3) == 'ahi' ) then
      select case ( trim(name_inst) )
         case ( 'ahi_himawari8' )
            obserrors = (/ 2.2, 3.0, 2.5, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2, 2.2 /)
         case default
            return
      end select
   else
      return
   end if
end subroutine set_ahi_obserr

subroutine decode_himawari_name(fname, finfo, iband, satid, file_time, region, resolution, jday, is_CLP, is_BCM, is_Phase, is_HT)

   implicit none

   character(len=256),intent(in)  :: fname
   character(len=2),  intent(out) :: finfo
   character(len=3),  intent(out) :: satid
   character(len=22), intent(out) :: file_time
   integer(i_kind),   intent(out) :: iband      ! 7-16
   character(len=4),  intent(out) :: region     ! 'FLDK'
   character(len=3),  intent(out) :: resolution ! 'R20, R05, etc'
   integer(i_kind),   intent(out) :: jday
   logical,           intent(out) :: is_CLP
   logical,           intent(out) :: is_BCM
   logical,           intent(out) :: is_HT
   logical,           intent(out) :: is_Phase

   integer(i_kind)   :: isegm
   character(len=4)  :: syear
   character(len=2)  :: smonth, sday, shour, sminute, ssec1
   character(len=1)  :: ssec2
   character(len=4)  :: version
   character(len=5)  :: pixelnumber, linenumber
   integer(i_kind)   :: year, month, day, hour, minute, sec, sec1, sec2

   !CLP_id:     NC_H08_20180512_1800_L2CLP010_FLDK.02401_02401.nc
   !BCM_id_old: Himawari8_AHI_FLDK_2019345_0000_00_CLOUD_MASK_EN.nc
   !            Himawari8_AHI_2KM_FLDK_2021023_0000_20_CLOUD_MASK_EN.nc
   !BCM_id_new: AHI-CMSK_v1r0_h08_s202103230300208_e202103230300542_c202103231944583.nc
   !            AHI-CMSK_v1r1_h09_s202308161430206_e202308161439400_c202308161452133.nc

   isegm = -99
   if ( fname(22:26) == CLP_id ) then
      is_CLP = .true.
      iband  = -99
   else if ( fname(36:45) == BCM_id_old .or. fname( 40:49) == BCM_id_old ) then
      is_BCM = .true.
      iband = -99
   else if ( fname(36:45) == HT_id_old .or. fname( 40:49) == HT_id_old ) then
      is_HT = .true.
      iband = -99
   else if ( fname(36:45) == Phase_id_old .or. fname( 40:49) == Phase_id_old ) then
      is_Phase = .true.
      iband = -99
   else if ( fname(5:8) == BCM_id_new ) then
      is_BCM = .true.
      iband = -99
   else if ( fname(5:8) == HT_id_new ) then
      is_HT = .true.
      iband = -99
   else if ( fname(5:8) == Phase_id_new ) then
      is_Phase = .true.
      iband = -99
   else
      is_CLP   = .false.
      is_Phase = .false.
      is_HT    = .false.
      is_BCM   = .false.
   end if

   if ( .not. ( is_CLP .or. is_BCM .or. is_Phase .or. is_HT ) ) then
      !HS_H08_20210123_0000_B07_FLDK_R20_S0110.DAT
      read(fname(1:2),   '(a2)')   finfo
      read(fname(4:6),   '(a3)')   satid
      read(fname(8:11),  '(i4)')   year
      read(fname(8:11),  '(a4)')   syear
      read(fname(12:13), '(i2)')   month
      read(fname(12:13), '(a2)')   smonth
      read(fname(14:15), '(i2)')   day
      read(fname(14:15), '(a2)')   sday
      read(fname(17:18), '(i2)')   hour
      read(fname(17:18), '(a2)')   shour
      read(fname(19:20), '(i2)')   minute
      read(fname(19:20), '(a2)')   sminute
      read(fname(23:24), '(i2)')   iband
      read(fname(26:29), '(a4)')   region
      read(fname(31:33), '(a3)')   resolution
      read(fname(36:37), '(i2)')   isegm
      ! 2017-10-01T18:02:00.0Z
      !print*, syear, smonth, sday, shour, sminute
      file_time = syear//'-'//smonth//'-'//sday//'T'//shour//':'//sminute//':00.0Z'

      jday = 0
      do i = 1, month - 1
        jday = jday + mmday(i)
      end do
      jday = jday + day

   else if ( is_CLP ) then
      !NC_H08_20180512_1800_L2CLP010_FLDK.02401_02401.nc
      finfo = 'HS'
      resolution = 'R50'
      read(fname(4:6),   '(a3)')   satid
      read(fname(8:11),  '(i4)')   year
      read(fname(8:11),  '(a4)')   syear
      read(fname(12:13), '(i2)')   month
      read(fname(12:13), '(a2)')   smonth
      read(fname(14:15), '(i2)')   day
      read(fname(14:15), '(a2)')   sday
      read(fname(17:18), '(i2)')   hour
      read(fname(17:18), '(a2)')   shour
      read(fname(19:20), '(i2)')   minute
      read(fname(19:20), '(a2)')   sminute
      read(fname(27:29), '(a3)')   version
      read(fname(31:34), '(a4)')   region
      read(fname(36:40), '(a5)')   pixelnumber
      read(fname(41:45), '(a5)')   linenumber
      ! 2017-10-01T18:02:00.0Z
      !print*, syear, smonth, sday, shour, sminute
      file_time = syear//'-'//smonth//'-'//sday//'T'//shour//':'//sminute//':00.0Z'

      jday = 0
      do i = 1, month - 1
        jday = jday + mmday(i)
      end do
      jday = jday + day

   else if ( (is_BCM .or. is_Phase .or. is_HT) .and. (fname(1:9) == 'Himawari8') ) then
      finfo = 'HS'
      if ( fname(15:17) == '2KM' ) then
         ! from 1610Z 2 July 2020 to 1450Z 22 March 2021
         !Himawari8_AHI_2KM_FLDK_2021023_0000_20_CLOUD_MASK_EN.nc
         resolution = 'R20'
         satid = 'H08'
         read(fname(19:22), '(a4)')   region
         read(fname(24:27), '(i4)')   year
         read(fname(24:27), '(a4)')   syear
         read(fname(28:30), '(i3)')   jday
         read(fname(32:33), '(i2)')   hour
         read(fname(32:33), '(a2)')   shour
         read(fname(34:35), '(i2)')   minute
         read(fname(34:35), '(a2)')   sminute
      else
         ! from 11 Dec 2019 to 1600Z 2 July 2020
         !Himawari8_AHI_FLDK_2019345_0000_00_CLOUD_MASK_EN.nc
         resolution = 'R50'
         satid = 'H08'
         read(fname(15:18), '(a4)')   region
         read(fname(20:23), '(i4)')   year
         read(fname(20:23), '(a4)')   syear
         read(fname(24:26), '(i3)')   jday
         read(fname(28:29), '(i2)')   hour
         read(fname(28:29), '(a2)')   shour
         read(fname(30:31), '(i2)')   minute
         read(fname(30:31), '(a2)')   sminute
      end if

      ! get month and day from julian day
      call get_date(year, jday, month, day)
      write(smonth, '(i2.2)') month
      write(sday, '(i2)') day
      ! 2017-10-01T18:02:00.0Z
      file_time = syear//'-'//smonth//'-'//sday//'T'//shour//':'//sminute//':00.0Z'

   else if ( (is_BCM .or. is_Phase .or. is_HT) .and. (fname(1:3) == 'AHI') ) then
      !AHI-CMSK_v1r0_h08_s202103230300208_e202103230300542_c202103231944583.nc
      !AHI-CMSK_v1r1_h09_s202308161430206_e202308161439400_c202308161452133.nc
      finfo = 'HS'
      resolution = 'R20'
      region = 'FLDK'
      read(fname(10:13), '(a4)')   version
      read(fname(15:17), '(a3)')   satid
      read(fname(20:23), '(i4)')   year
      read(fname(20:23), '(a4)')   syear
      read(fname(24:25), '(i2)')   month
      read(fname(24:25), '(a2)')   smonth
      read(fname(26:27), '(i2)')   day
      read(fname(26:27), '(a2)')   sday
      read(fname(28:29), '(i2)')   hour
      read(fname(28:29), '(a2)')   shour
      read(fname(30:31), '(i2)')   minute
      read(fname(30:31), '(a2)')   sminute
      read(fname(32:33), '(i2)')   sec1   ! integer part of second
      read(fname(32:33), '(a2)')   ssec1
      read(fname(34:34), '(a1)')   sec2   ! decimal part of second
      read(fname(34:34), '(a1)')   ssec2
      ! 2017-10-01T18:02:00.0Z
      !print*, syear, smonth, sday, shour, sminute
      file_time = syear//'-'//smonth//'-'//sday//'T'//shour//':'//sminute//':'//ssec1//'Z'

      jday = 0
      do i = 1, month - 1
        jday = jday + mmday(i)
      end do
      jday = jday + day

   else
      write(0,*) 'Error decode_himawari_name'
      stop
   end if

   return
end subroutine decode_himawari_name

end module mod_himawari_ahi
