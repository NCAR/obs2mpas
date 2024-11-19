!----------------------------------------------------------------------
! Code adapted from 
!    https://github.com/jamiebresch/obs2ioda/blob/main/goes_abi/src/goes_abi_converter.f90
!
! Yonggang G. Yu
! 27-Oct-2022
!----------------------------------------------------------------------
!
module  mod_goes_abi
!
! Purpose: Convert GOES ReBroadcast netCDF files to ioda-v1 format.
!          Currently only processes bands 7-16.
!
! input files:
!    (1) flist.txt: contains a list of nc files (exclude path) to be processed
!                     GoesReBroadcast file
!                     (optional) Clear Sky Mask output of cspp-geo-aitf package
!    (2) namelist.obs2model
!        &data_nml
!          list_files = 'flist.txt'        ! list of files
!          data_dir = '/data/goes',        ! path of the GRB nc files
!          data_id = 'OR_ABI-L1b-RadC-M3', ! prefix of the downloaded GRB nc files
!          sat_id = 'G16',                 ! satellite ID
!          n_subsample = 1,                ! value use for thinning if write_iodav1 = .true.
!          write_iodav1 = .false.,         ! option to write out an iodav1 file (no superobbing)
!        /

   use control_para !BJJ
   use utils_mod

   implicit none

   ! prefix of Clear Sky Mask (Binary Cloud Mask) output of cspp-geo-aitf package
   character(len=14), parameter :: BCM_id   = 'OR_ABI-L2-ACMF'
   character(len=15), parameter :: TEMP_id  = 'OR_ABI-L2-ACHTF'
   character(len=15), parameter :: Phase_id = 'OR_ABI-L2-ACTPF'
   character(len=15), parameter :: HT_id    = 'OR_ABI-L2-ACHAF'
   character(len=14), parameter :: PRES_id  = 'OR_ABI-L2-CTPF'

   integer(i_kind), parameter :: nband      = 10  ! IR bands 7-16
   integer(i_kind) :: band_start = 7
   integer(i_kind) :: band_end   = 16

   logical, allocatable :: got_latlon(:,:)
   real(r_kind), allocatable :: glat(:,:)    ! grid latitude (nx,ny)
   real(r_kind), allocatable :: glon(:,:)    ! grid longitude (nx,ny)
   real(r_kind), allocatable :: gzen(:,:)    ! satellite zenith angle (nx,ny)
   real(r_kind), allocatable :: solzen(:,:)  ! solar zenith angle (nx,ny)

   real(r_kind),    allocatable :: rad_2d(:,:)  ! radiance(nx,ny)
   real(r_kind),    allocatable :: bt_2d(:,:)   ! brightness temperature(nx,ny)
   integer(i_kind), allocatable :: qf_2d(:,:)   ! quality flag(nx,ny)
            ! qf (DQF, Data Quality Flag)
            ! 0:good, 1:conditionally_usable, 2:out_of_range, 3:no_value
   integer(i_kind), allocatable :: cm_2d(:,:)   ! cloud_mask(nx,ny)
   real(r_kind),    allocatable :: ctt_2d(:,:)  ! cloud top temperature(nx,ny) !BJJ
   real(r_kind),    allocatable :: ctph_2d(:,:) ! cloud top phase(nx,ny) !BJJ
   real(r_kind),    allocatable :: cth_2d(:,:)  ! cloud top height(nx,ny) !BJJ
   real(r_kind),    allocatable :: ctp_2d(:,:)  ! cloud top pressure(nx,ny) !BJJ

   type rad_type
      real(r_kind),    allocatable :: rad(:,:,:)  ! radiance(nband,nx,ny)
      real(r_kind),    allocatable :: bt(:,:,:)   ! brightness temperature(nband,nx,ny)
      integer(i_kind), allocatable :: qf(:,:,:)   ! quality flag(nband,nx,ny)
      real(r_kind),    allocatable :: sd(:)       ! std_dev(nband)
      integer(i_kind), allocatable :: cm(:,:)     ! cloud mask(nx,ny)
   end type rad_type
   type(rad_type), allocatable  :: rdata(:)  ! (ntime)

   character(len=22), allocatable :: time_start(:)  ! (ntime) 2017-10-01T18:02:19.6Z

   integer(i_kind) :: ncid, nf_status
   integer(i_kind) :: nx, ny
   integer(i_kind) :: it, ib, ii, i, j
   integer(i_kind) :: ntime
   integer(i_kind) :: t_index
   integer(i_kind) :: band_id

   real(r_kind)                    :: sdtb ! to be done
   integer(i_kind)                 :: ifile, nlen
   logical                         :: isfile
   logical                         :: found_time
   logical                         :: got_grid_info
   logical, allocatable            :: valid(:), is_BCM(:), is_TEMP(:), is_Phase(:), is_HT(:), is_PRES(:)
   character(len=256)              :: fname
   character(len=256)              :: out_fname
   character(len=18)               :: finfo
   character(len=2)                :: mode_id, scan_mode
   character(len=3)                :: fsat_id
   character(len=22), allocatable  :: scan_time(:) ! 2017-10-01T18:02:19.6Z
   integer(i_kind),   allocatable  :: fband_id(:)
   integer(i_kind),   allocatable  :: ftime_id(:)
   integer(i_kind),   allocatable  :: julianday(:)

   contains


subroutine Goes_ReBroadcast_converter(glon_out, glat_out, F_out, varname_out, got_latlon_out)

   implicit none
   real(r_kind),      allocatable, intent(out) :: glon_out(:,:)
   real(r_kind),      allocatable, intent(out) :: glat_out(:,:)
   real(r_kind),      allocatable, intent(out) :: F_out(:,:,:)         ! (nx,ny, nfield), nfield=nfile, one file for each field
   character(len=64), allocatable, intent(out) :: varname_out(:)       ! (nfield)
   logical,           allocatable, intent(out) :: got_latlon_out(:,:)  ! (nx,ny)
   ! loc
   integer :: ix

   integer(i_kind)                 :: nfile
   character(len=256), allocatable :: nc_fnames(:)

   character(len=256)              :: nc_list_file  ! the text file that contains a list of netcdf files to process
   character(len=256)              :: data_dir
   character(len=18)               :: data_id
   character(len=3)                :: sat_id
   integer(i_kind)                 :: n_subsample
   logical                         :: write_iodav1

   ! get namelist variables
   call get_namelist_vars(nfile, nc_fnames, nc_list_file, data_dir, data_id, sat_id, n_subsample, write_iodav1)

   allocate (ftime_id(nfile))
   allocate (scan_time(nfile))
   allocate (julianday(nfile))
   allocate (fband_id(nfile))
   allocate (valid(nfile))
   allocate (is_BCM(nfile))
   allocate (is_TEMP(nfile)) !BJJ
   allocate (is_Phase(nfile)) !BJJ
   allocate (is_HT(nfile)) !BJJ
   allocate (is_PRES(nfile)) !BJJ
   valid(:) = .false.
   is_BCM(:) = .false.
   is_TEMP(:) = .false.
   is_Phase(:) = .false.
   is_HT(:) = .false.
   is_PRES(:) = .false.

   nlen = len_trim(data_id)
   mode_id = data_id(nlen-1:nlen)

   ! parse the file list
   t_index = 0
   file_loop1: do ifile = 1, nfile
      fname = trim(data_dir)//'/'//trim(nc_fnames(ifile))
      inquire(file=trim(fname), exist=isfile)
      if ( .not. isfile ) then
         write(0,*) 'File not found: '//trim(fname)
         cycle file_loop1
      else
         write(0,*) 'File found: '//trim(fname)         
      end if

      ! retrieve some basic info from the netcdf filename itself
      call decode_nc_fname(trim(nc_fnames(ifile)),finfo, scan_mode, &
         is_BCM(ifile), is_TEMP(ifile), is_Phase(ifile), is_HT(ifile), is_PRES(ifile), &
         fband_id(ifile), fsat_id, scan_time(ifile), julianday(ifile))

      ! all files must be the same mode
      if ( scan_mode /= mode_id ) then
         cycle file_loop1
      end if

      if ( fsat_id /= sat_id ) then
         cycle file_loop1
      end if

      if ( .not. ( is_BCM(ifile) .or. is_TEMP(ifile) .or. is_Phase(ifile) & 
                 .or. is_HT(ifile) .or. is_PRES(ifile) ) ) then
         ! id of the file name must match specified data_id
         if ( finfo /= data_id ) then
            cycle file_loop1
         else
            ! only process band 7-16
            if ( fband_id(ifile) < band_start .or. fband_id(ifile) > band_end ) then
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
      write(0,*) 'No valid files found from nc_list_file '//trim(nc_list_file)
      stop
   end if

   allocate (time_start(ntime))
   allocate (rdata(ntime))

   got_grid_info = .false.
   file_loop2: do ifile = 1, nfile

      if ( valid(ifile) ) then

         fname = trim(data_dir)//'/'//trim(nc_fnames(ifile))
         nf_status = nf_OPEN(trim(fname), nf_NOWRITE, ncid)
         if ( nf_status == 0 ) then
            write(0,*) 'Reading '//trim(fname)
         else
            write(0,*) 'ERROR reading '//trim(fname)
            cycle file_loop2
         end if

         if ( .not. got_grid_info ) then
            call read_GRB_dims(ncid, nx, ny)
            allocate (glat(nx, ny))
            allocate (glon(nx, ny))
            allocate (gzen(nx, ny))
            allocate (solzen(nx, ny))
            allocate (got_latlon(nx, ny))
            glat(:,:) = missing_r
            glon(:,:) = missing_r
            gzen(:,:) = missing_r
            solzen(:,:) = missing_r
            !BJJ allocate output arrays
            allocate (glat_out(nx, ny))
            allocate (glon_out(nx, ny))
            allocate (F_out(nx, ny, nfile))
            allocate (varname_out(nfile))
            allocate (got_latlon_out(nx, ny))
            glat_out(:,:) = missing_r
            glon_out(:,:) = missing_r
            F_out(:,:,:) = missing_r
            varname_out = ''
            write(0,*) 'Calculating lat/lon from fixed grid x/y...'
            call read_GRB_grid(ncid, nx, ny, glat, glon, gzen, got_latlon)
            call calc_solar_zenith_angle(nx, ny, glat, glon, scan_time(ifile), julianday(ifile), solzen, got_latlon)
            ! additional info for writing ioda at MPAS mesh !BJJ
            write(15,*) ifile, scan_time(ifile), julianday(ifile)
            !BJJ copy to output array
            glat_out(:,:)=glat(:,:)
            glon_out(:,:)=glon(:,:)
            got_latlon_out(:,:)=got_latlon(:,:)
            got_grid_info = .true.
            allocate (rad_2d(nx, ny))
            allocate (bt_2d(nx, ny))
            allocate (qf_2d(nx, ny))
            allocate (cm_2d(nx, ny))   
         end if

         it = ftime_id(ifile)
         ib = fband_id(ifile)

         if ( .not. ( is_BCM(ifile) .or. is_TEMP(ifile) .or. is_Phase(ifile) .or. is_HT(ifile) .or. is_PRES(ifile) ) ) then

            call read_GRB(ncid, nx, ny, rad_2d, bt_2d, qf_2d, sdtb, band_id, time_start(it))

            if ( band_id /= ib ) then
               write(0,*) 'ERROR: band_id from the file name and the file content do not match.'
               cycle file_loop2
            end if

            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if

            if ( .not. allocated(rdata(it)%rad) ) allocate (rdata(it)%rad(nband,nx,ny))
            if ( .not. allocated(rdata(it)%bt) )  allocate (rdata(it)%bt(nband,nx,ny))
            if ( .not. allocated(rdata(it)%qf) )  allocate (rdata(it)%qf(nband,nx,ny))
            if ( .not. allocated(rdata(it)%sd) )  allocate (rdata(it)%sd(nband))

            do j = 1, ny
               do i = 1, nx
                  ! convert band id 7-16 to array index 1-10
                  rdata(it)%rad(ib-band_start+1,i,j) = rad_2d(i,j)
                  rdata(it)%bt(ib-band_start+1,i,j)  = bt_2d(i,j)
                  rdata(it)%qf(ib-band_start+1,i,j)  = qf_2d(i,j)
                  rdata(it)%sd(ib-band_start+1)      = sdtb
               end do
            end do

            !BJJ copy to output array: use ifile index
            !F_out(1:nx,1:ny,ifile)=rad_2d(1:nx,1:ny)  !Radiance
            !!varname_out(ifile)='Rad_'//fsat_id//fband_id(ifile)
            !write(varname_out(ifile),"(A,I2.2)") 'Rad_'//fsat_id//'C', fband_id(ifile)
            F_out(1:nx,1:ny,ifile)=bt_2d(1:nx,1:ny)  !Brightness Temperature
            write(varname_out(ifile),"(A,I2.2)") 'BT_'//fsat_id//'C', fband_id(ifile)

         elseif ( is_BCM(ifile) ) then
            call read_L2_BCM(ncid, nx, ny, cm_2d, time_start(it))
            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if
            if ( .not. allocated(rdata(it)%cm) )  allocate (rdata(it)%cm(nx,ny))
            rdata(it)%cm(:,:) = cm_2d(:,:)

            !BJJ copy to output array: use ifile index
            F_out(1:nx,1:ny,ifile)=cm_2d(1:nx,1:ny)
            varname_out(ifile)='BCM_'//fsat_id

         elseif ( is_TEMP(ifile) ) then
            allocate (ctt_2d(nx, ny))   
            call read_L2_TEMP(ncid, nx, ny, ctt_2d, time_start(it))
            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if
            !BJJ copy to output array: use ifile index
            F_out(1:nx,1:ny,ifile)=ctt_2d(1:nx,1:ny)
            varname_out(ifile)='TEMP_'//fsat_id
            deallocate (ctt_2d)
            
         elseif ( is_Phase(ifile) ) then
            allocate (ctph_2d(nx, ny))
            call read_L2_Phase(ncid, nx, ny, ctph_2d, time_start(it))
            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if
            !BJJ copy to output array: use ifile index
            F_out(1:nx,1:ny,ifile)=ctph_2d(1:nx,1:ny)
            varname_out(ifile)='Phase_'//fsat_id
            deallocate (ctph_2d)

         elseif ( is_HT(ifile) ) then
            allocate (cth_2d(nx, ny))
            call read_L2_HT(ncid, nx, ny, cth_2d, time_start(it))
            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if
            !BJJ copy to output array: use ifile index
            F_out(1:nx,1:ny,ifile)=cth_2d(1:nx,1:ny)
            varname_out(ifile)='HT_'//fsat_id
            deallocate (cth_2d)

         elseif ( is_PRES(ifile) ) then
            allocate (ctp_2d(nx, ny))   
            call read_L2_PRES(ncid, nx, ny, ctp_2d, time_start(it))
            if ( time_start(it) /= scan_time(ifile) ) then
               write(0,*) 'ERROR: scan start time from the file name and the file content do not match.'
               cycle file_loop2
            end if
            !BJJ copy to output array: use ifile index
            F_out(1:nx,1:ny,ifile)=ctp_2d(1:nx,1:ny)
            varname_out(ifile)='PRES_'//fsat_id
            deallocate (ctp_2d)
            
         else
            write(0,*) 'ERROR: something is wrong. check the files'
            stop
         end if         

         nf_status = nf_CLOSE(ncid)

      end if

   end do file_loop2

   if ( allocated(rad_2d) ) deallocate(rad_2d)
   if ( allocated(bt_2d) )  deallocate(bt_2d)
   if ( allocated(qf_2d) )  deallocate(qf_2d)
   if ( allocated(cm_2d) )  deallocate(cm_2d)

   if ( write_iodav1 ) then
      do it = 1, ntime
         out_fname = trim(data_id)//'_'//sat_id//'_'//time_start(it)//'.nc4'
         write(0,*) 'Writing ', trim(out_fname)
         if ( allocated(rdata(it)%cm) ) then
            call output_iodav1(trim(out_fname), time_start(it),      &
                               nx, ny, nband, n_subsample,           &
                               got_latlon, glat, glon, gzen, solzen, &
                               rdata(it)%bt, rdata(it)%qf, rdata(it)%sd, rdata(it)%cm)
         else
            call output_iodav1(trim(out_fname), time_start(it),      &
                               nx, ny, nband, n_subsample,           &
                               got_latlon, glat, glon, gzen, solzen, &
                               rdata(it)%bt, rdata(it)%qf, rdata(it)%sd)
         end if
      end do
   end if

   
   if ( allocated(glat) )   deallocate(glat)
   if ( allocated(glon) )   deallocate(glon)
   if ( allocated(gzen) )   deallocate(gzen)
   if ( allocated(solzen) ) deallocate(solzen)

   do it = 1, ntime
      if ( allocated(rdata(it)%rad) ) deallocate (rdata(it)%rad)
      if ( allocated(rdata(it)%bt)  ) deallocate (rdata(it)%bt)
      if ( allocated(rdata(it)%qf)  ) deallocate (rdata(it)%qf)
      if ( allocated(rdata(it)%cm)  ) deallocate (rdata(it)%cm)
   end do
   deallocate(rdata)
   deallocate(time_start)

   deallocate(nc_fnames)
   deallocate(ftime_id)
   deallocate(scan_time)
   deallocate(julianday)
   deallocate(fband_id)
   deallocate(valid)
   deallocate(is_BCM)
   deallocate(is_TEMP) !BJJ
   deallocate(is_Phase) !BJJ
   deallocate(is_HT) !BJJ
   deallocate(is_PRES) !BJJ

end subroutine Goes_ReBroadcast_converter

 
subroutine read_GRB_dims(ncid, nx, ny)
   implicit none
   integer(i_kind), intent(in)  :: ncid
   integer(i_kind), intent(out) :: nx, ny
   integer(i_kind)              :: dimid
   integer(i_kind)              :: nf_status(4)
   continue
   nf_status(1) = nf_INQ_DIMID(ncid, 'x', dimid)
   nf_status(2) = nf_INQ_DIMLEN(ncid, dimid, nx)
   nf_status(3) = nf_INQ_DIMID(ncid, 'y', dimid)
   nf_status(4) = nf_INQ_DIMLEN(ncid, dimid, ny)
   if ( any(nf_status /= 0) ) then
      write(0,*) 'Error reading dimensions'
      stop
   end if
   return
end subroutine read_GRB_dims

!NC_BYTE 8-bit signed integer
!NC_SHORT 16-bit signed integer
!NC_INT (or NC_LONG) 32-bit signed integer
!NC_FLOAT 32-bit floating point
!NC_DOUBLE 64-bit floating point

subroutine read_GRB_grid(ncid, nx, ny, glat, glon, gzen, got_latlon)
   implicit none
   integer(i_kind), intent(in)    :: ncid
   integer(i_kind), intent(in)    :: nx, ny
   real(r_kind),    intent(inout) :: glat(nx,ny)
   real(r_kind),    intent(inout) :: glon(nx,ny)
   real(r_kind),    intent(inout) :: gzen(nx,ny)
   logical,         intent(inout) :: got_latlon(nx,ny)
   integer(i_kind)                :: varid, i, j
   integer(i_kind)                :: nf_status
   integer(i_kind)                :: istart(1), icount(1)
   integer(i_short), allocatable  :: itmp_short_1d(:)
   real(r_kind),     allocatable  :: x(:)
   real(r_kind),     allocatable  :: y(:)
   real(r_single) :: scalef, offset
   real(r_double) :: dtmp
   real(r_double) :: r_eq    ! GRS80 semi-major axis of earth
   real(r_double) :: r_pol   ! GRS80 semi-minor axis of earth = (1-f)*r_eq
   real(r_double) :: lon_sat ! satellite longitude, longitude_of_projection_origin
   real(r_double) :: h_sat   ! satellite height
   real(r_double) :: a, b, c, rs, sx, sy, sz
   real(r_kind)   :: rlat, rlon, lon_diff, tmp1, theta1, theta2
   continue

!int goes_imager_projection ;
!  goes_imager_projection:long_name = "GOES-R ABI fixed grid projection" ;
!  goes_imager_projection:grid_mapping_name = "geostationary" ;
!  goes_imager_projection:perspective_point_height = 35786023. ;
!  goes_imager_projection:semi_major_axis = 6378137. ;
!  goes_imager_projection:semi_minor_axis = 6356752.31414 ;
!  goes_imager_projection:inverse_flattening = 298.2572221 ;
!  goes_imager_projection:latitude_of_projection_origin = 0. ;
!  goes_imager_projection:longitude_of_projection_origin = -89.5 ;
!  goes_imager_projection:sweep_angle_axis = "x" ;

   nf_status = nf_INQ_VARID(ncid, 'goes_imager_projection', varid)
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'semi_major_axis',  dtmp)
   r_eq = dtmp
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'semi_minor_axis',  dtmp)
   r_pol = dtmp
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'perspective_point_height',  dtmp)
   h_sat = dtmp + r_eq  ! perspective_point_height + semi_major_axis
   nf_status = nf_GET_ATT_DOUBLE(ncid, varid, 'longitude_of_projection_origin',  dtmp)
   lon_sat = dtmp * deg2rad

!short x(x) ;
!  x:scale_factor = 5.6e-05f ;
!  x:add_offset = -0.075012f ;
!  x:units = "rad" ;
!  x:axis = "X" ;
!  x:long_name = "GOES fixed grid projection x-coordinate" ;

   istart(1) = 1
   icount(1) = nx
   allocate(itmp_short_1d(nx))
   nf_status = nf_INQ_VARID(ncid, 'x', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:1), icount(1:1), itmp_short_1d(:))
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   allocate(x(nx))
   do i = 1, nx
      x(i) = offset + itmp_short_1d(i) * scalef
   end do
   deallocate(itmp_short_1d)

!short y(y) ;
!  y:scale_factor = -5.6e-05f ;
!  y:add_offset = 0.126532f ;
!  y:units = "rad" ;
!  y:axis = "Y" ;
!  y:long_name = "GOES fixed grid projection y-coordinate" ;
!  y:standard_name = "projection_y_coordinate" ;

   istart(1) = 1
   icount(1) = ny
   allocate(itmp_short_1d(ny))
   nf_status = nf_INQ_VARID(ncid, 'y', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:1), icount(1:1), itmp_short_1d(:))
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   allocate(y(ny))
   do i = 1, ny
      y(i) = offset + itmp_short_1d(i) * scalef
   end do
   deallocate(itmp_short_1d)
   ! Product Definition and User's Guide (PUG) Volume 3, pp. 19-21
   ! from fixed grid x/y to geodetic lat/lon
   got_latlon(1:nx,1:ny) = .true.
   do j = 1, ny
      do i = 1, nx
         a = sin(x(i))*sin(x(i)) + cos(x(i))*cos(x(i)) * &
             (cos(y(j))*cos(y(j))+(r_eq/r_pol)*(r_eq/r_pol)*sin(y(j))*sin(y(j)))
         b = -2.0*h_sat*cos(x(i))*cos(y(j))
         c = h_sat*h_sat - r_eq*r_eq
         if ( (b*b-4.0*a*c) < 0.0 ) then
            got_latlon(i,j) = .false.
            cycle
         end if
         rs = (-1.0*b - sqrt(b*b-4.0*a*c)) / (2.0*a)
         sx = rs * cos(x(i)) * cos(y(j))
         sy = -1.0 * rs * sin(x(i))
         sz = rs * cos(x(i)) * sin(y(j))
         !glat(i,j) = (atan((r_eq/r_pol)*(r_eq/r_pol)*(sz/sqrt((h_sat-sx)*(h_sat-sx)+sy*sy)))) * rad2deg
         !glon(i,j) = (lon_sat - atan(sy/(h_sat-sx))) * rad2deg
         glat(i,j) = atan((r_eq/r_pol)*(r_eq/r_pol)*(sz/sqrt((h_sat-sx)*(h_sat-sx)+sy*sy)))
         glon(i,j) = lon_sat - atan(sy/(h_sat-sx))
      end do
   end do

   deallocate(x)
   deallocate(y)

   ! calculate geostationary satellite zenith angle
   do j = 1, ny
      do i = 1, nx
         if ( .not. got_latlon(i,j) ) cycle
         ! glat, glon, gzen are in [radian] in this routine.
         call calc_geostationary_satellite_zenith_angle( \
              glat(i,j), glon(i,j), lon_sat, r_eq, h_sat, gzen(i,j) )
         glat(i,j) = glat(i,j) * rad2deg
         glon(i,j) = glon(i,j) * rad2deg
         gzen(i,j) = gzen(i,j) * rad2deg
      end do
   end do

   ! additional info for writing ioda at MPAS mesh !BJJ
   write(15,*) lon_sat, r_eq, h_sat

   return
end subroutine read_GRB_grid

subroutine read_GRB(ncid, nx, ny, rad, bt, qf, sd, band_id, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   integer(i_kind),   intent(out)   :: band_id
   real(r_kind),      intent(out)   :: sd
   real(r_kind),      intent(inout) :: rad(nx,ny)
   real(r_kind),      intent(inout) :: bt(nx,ny)
   integer(i_kind),   intent(inout) :: qf(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_1d(:)
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_short), allocatable    :: itmp_short_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_short)                 :: ifill
   real(r_single)                   :: rfill
   real(r_single)                   :: rtmp
   real(r_single)                   :: planck_fk1, planck_fk2
   real(r_single)                   :: planck_bc1, planck_bc2
   real(r_single)                   :: scalef, offset
   real(r_kind)                     :: rmiss = -999.0
   integer(i_kind)                  :: imiss = -999
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = 1
   allocate(itmp_byte_1d(1))
   nf_status = nf_INQ_VARID(ncid, 'band_id', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:1), icount(1:1), itmp_byte_1d(:))
   band_id = itmp_byte_1d(1)
   deallocate(itmp_byte_1d)

   nf_status = nf_INQ_VARID(ncid, 'std_dev_radiance_value_of_valid_pixels', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, rtmp)
   sd = rtmp

   ! qf (DQF, Data Quality Flag)
   ! 0:good, 1:conditionally_usable, 2:out_of_range, 3:no_value
   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   qf(:,:) = imiss
   do j = 1, ny
      do i = 1, nx
         qf(i,j) = itmp_byte_2d(i,j)
      end do
   end do
   deallocate(itmp_byte_2d)

   nf_status = nf_INQ_VARID(ncid, 'planck_fk1', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_fk1)
   nf_status = nf_GET_ATT_REAL(ncid, varid, '_FillValue',  rfill)

   nf_status = nf_INQ_VARID(ncid, 'planck_fk2', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_fk2)

   nf_status = nf_INQ_VARID(ncid, 'planck_bc1', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_bc1)

   nf_status = nf_INQ_VARID(ncid, 'planck_bc2', varid)
   nf_status = nf_GET_VAR_REAL(ncid, varid, planck_bc2)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_short_2d(nx, ny))
   nf_status = nf_INQ_VARID(ncid, 'Rad', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   rad(:,:) = rmiss
   bt(:,:)  = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( itmp_short_2d(i,j) /= ifill ) then
            if ( itmp_short_2d(i,j) .lt. 0_i_short ) STOP 777
            rad(i,j) = offset + itmp_short_2d(i,j) * scalef
            if ( planck_fk1 /= rfill .and. planck_fk2 /= rfill .and. &
                 planck_bc1 /= rfill .and. planck_bc2 /= rfill ) then
              if ( rad(i,j) > 0.0 ) then
               bt(i,j) = (planck_fk2/(log((planck_fk1/rad(i,j))+1.0))-planck_bc1)/planck_bc2
              end if
            end if
         end if
      end do
   end do
   deallocate(itmp_short_2d)

   return
end subroutine read_GRB

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
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
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
   nf_status = nf_INQ_VARID(ncid, 'BCM', varid)
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

subroutine read_L2_TEMP(ncid, nx, ny, ctt, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   real(r_kind),      intent(inout) :: ctt(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_short), allocatable    :: itmp_short_2d(:,:)
   integer(i_kind),  allocatable    :: itmp_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_short)                 :: ifill
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: rmiss = -999.0
   real(r_single)                   :: scalef, offset
   integer(i_kind)                  :: qf(nx,ny)
   character(len=4)                 :: l_unsigned
   integer(i_kind)                  :: xtype
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
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
   allocate(itmp_short_2d(nx, ny))
   nf_status = nf_INQ_VARID(ncid, 'TEMP', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   nf_status = nf_GET_ATT_TEXT(ncid, varid, '_Unsigned', l_unsigned)
   if( nf_status .eq. 0 ) write(*,*) "---- Attribute @_Unsigned = ",l_unsigned
   nf_status = nf_INQ_VARTYPE(ncid, varid, xtype)
   if( nf_status .eq. 0 ) write(*,*) "---- xtype = ", xtype
   if( xtype .eq. nf_ushort ) l_unsigned = "true"  ! nf_ushort is defined in "netcdf.inc"
   ! convert unsigned short to signed short
   if( l_unsigned == "true" ) then
      allocate(itmp_2d(nx,ny))
      itmp_2d = itmp_short_2d
      do j = 1, ny
         do i = 1, nx
            if ( itmp_short_2d(i,j) .lt. 0_i_short ) then
               itmp_2d(i,j) = itmp_2d(i,j) + 65536
            end if
         end do
      end do
   end if
   !write(*,*) "min/max of itmp_short_2d =", minval(itmp_short_2d), maxval(itmp_short_2d)
   !write(*,*) "min/max of itmp_2d =", minval(itmp_2d), maxval(itmp_2d)

   ctt(:,:) = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( itmp_short_2d(i,j) /= ifill ) then
            if (qf(i,j) == 0 ) then ! good quality
               if( l_unsigned == "true" ) then
                  ctt(i,j) = offset + itmp_2d(i,j) * scalef
               else
                  ctt(i,j) = offset + itmp_short_2d(i,j) * scalef
               end if
            end if
         end if
      end do
   end do
   !write(*,*) "min/max of ctt =", minval(ctt), maxval(ctt)
   deallocate(itmp_short_2d)
   if( allocated(itmp_2d) ) deallocate(itmp_2d)

   return
end subroutine read_L2_TEMP

subroutine read_L2_Phase(ncid, nx, ny, ctph, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   real(r_kind),      intent(inout) :: ctph(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: rmiss = -999.0
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
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
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
   nf_status = nf_INQ_VARID(ncid, 'Phase', varid)
   nf_status = nf_GET_VARA_INT1(ncid, varid, istart(1:2), icount(1:2), itmp_byte_2d(:,:))
   ctph(:,:) = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( qf(i,j) == 0 ) then ! good quality
            ctph(i,j) = itmp_byte_2d(i,j)
         end if
      end do
   end do
   deallocate(itmp_byte_2d)

   return
end subroutine read_L2_Phase

subroutine read_L2_HT(ncid, nx, ny, cth, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   real(r_kind),      intent(inout) :: cth(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_short), allocatable    :: itmp_short_2d(:,:)
   integer(i_kind),  allocatable    :: itmp_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_short)                 :: ifill
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: rmiss = -999.0
   real(r_single)                   :: scalef, offset
   integer(i_kind)                  :: qf(nx,ny)
   character(len=4)                 :: l_unsigned
   integer(i_kind)                  :: xtype
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
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
   allocate(itmp_short_2d(nx, ny))
   nf_status = nf_INQ_VARID(ncid, 'HT', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   nf_status = nf_GET_ATT_TEXT(ncid, varid, '_Unsigned', l_unsigned)
   if( nf_status .eq. 0 ) write(*,*) "---- Attribute @_Unsigned = ",l_unsigned 
   nf_status = nf_INQ_VARTYPE(ncid, varid, xtype)
   if( nf_status .eq. 0 ) write(*,*) "---- xtype = ", xtype
   if( xtype .eq. nf_ushort ) l_unsigned = "true"  ! nf_ushort is defined in "netcdf.inc"
   ! convert unsigned short to signed short
   if( l_unsigned == "true" ) then
      allocate(itmp_2d(nx,ny))
      itmp_2d = itmp_short_2d
      do j = 1, ny
         do i = 1, nx
            if ( itmp_short_2d(i,j) .lt. 0_i_short ) then
               itmp_2d(i,j) = itmp_2d(i,j) + 65536
            end if
         end do
      end do
   end if
   !write(*,*) "min/max of itmp_short_2d =", minval(itmp_short_2d), maxval(itmp_short_2d)
   !write(*,*) "min/max of itmp_2d =", minval(itmp_2d), maxval(itmp_2d)

   cth(:,:) = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( itmp_short_2d(i,j) /= ifill ) then
            if (qf(i,j) == 0 ) then ! good quality
               if( l_unsigned == "true" ) then
                  cth(i,j) = offset + itmp_2d(i,j) * scalef
               else
                  cth(i,j) = offset + itmp_short_2d(i,j) * scalef
               end if
            end if
         end if
      end do
   end do
   !write(*,*) "min/max of cth =", minval(cth), maxval(cth)
   deallocate(itmp_short_2d)
   if( allocated(itmp_2d) ) deallocate(itmp_2d)

   return
end subroutine read_L2_HT

subroutine read_L2_PRES(ncid, nx, ny, ctp, time_start)
   implicit none
   integer(i_kind),   intent(in)    :: ncid
   integer(i_kind),   intent(in)    :: nx, ny
   real(r_kind),      intent(inout) :: ctp(nx,ny)
   character(len=22), intent(out)   :: time_start  ! 2017-10-01T18:02:19.6Z
   integer(i_byte),  allocatable    :: itmp_byte_2d(:,:)
   integer(i_short), allocatable    :: itmp_short_2d(:,:)
   integer(i_kind),  allocatable    :: itmp_2d(:,:)
   integer(i_kind)                  :: nf_status
   integer(i_kind)                  :: istart(2), icount(2)
   integer(i_kind)                  :: varid, i, j
   integer(i_short)                 :: ifill
   integer(i_kind)                  :: imiss = -999
   integer(i_kind)                  :: rmiss = -999.0
   real(r_single)                   :: scalef, offset
   integer(i_kind)                  :: qf(nx,ny)
   character(len=4)                 :: l_unsigned
   integer(i_kind)                  :: xtype
   continue

   ! time_start is the same for all bands, but time_end is not
   nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_start', time_start)
   !nf_status = nf_GET_ATT_TEXT(ncid, nf_GLOBAL, 'time_coverage_end',   time_end)

   istart(1) = 1
   icount(1) = nx
   istart(2) = 1
   icount(2) = ny
   allocate(itmp_byte_2d(nx,ny))
   nf_status = nf_INQ_VARID(ncid, 'DQF', varid)
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
   allocate(itmp_short_2d(nx, ny))
   nf_status = nf_INQ_VARID(ncid, 'PRES', varid)
   nf_status = nf_GET_VARA_INT2(ncid, varid, istart(1:2), icount(1:2), itmp_short_2d(:,:))
   nf_status = nf_GET_ATT_INT2(ncid, varid, '_FillValue',  ifill)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'scale_factor', scalef)
   nf_status = nf_GET_ATT_REAL(ncid, varid, 'add_offset', offset)
   nf_status = nf_GET_ATT_TEXT(ncid, varid, '_Unsigned', l_unsigned)
   if( nf_status .eq. 0 ) write(*,*) "---- Attribute @_Unsigned = ",l_unsigned
   nf_status = nf_INQ_VARTYPE(ncid, varid, xtype)
   if( nf_status .eq. 0 ) write(*,*) "---- xtype = ", xtype
   if( xtype .eq. nf_ushort ) l_unsigned = "true"  ! nf_ushort is defined in "netcdf.inc"
   ! convert unsigned short to signed short
   if( l_unsigned == "true" ) then
      allocate(itmp_2d(nx,ny))
      itmp_2d = itmp_short_2d
      do j = 1, ny
         do i = 1, nx
            if ( itmp_short_2d(i,j) .lt. 0_i_short ) then
               itmp_2d(i,j) = itmp_2d(i,j) + 65536
            end if
         end do
      end do
   end if
   !write(*,*) "min/max of itmp_short_2d =", minval(itmp_short_2d), maxval(itmp_short_2d)
   !write(*,*) "min/max of itmp_2d =", minval(itmp_2d), maxval(itmp_2d)

   ctp(:,:) = rmiss
   do j = 1, ny
      do i = 1, nx
         if ( itmp_short_2d(i,j) /= ifill ) then
            if (qf(i,j) == 0 ) then ! good quality
               if( l_unsigned == "true" ) then
                  ctp(i,j) = offset + itmp_2d(i,j) * scalef
               else
                  ctp(i,j) = offset + itmp_short_2d(i,j) * scalef
               end if
            end if
         end if
         if ( itmp_short_2d(i,j) == ifill .and. qf(i,j) == 4 ) then !invalid_due_to_clear_or_probably_clear_sky_qf
            ctp(i,j) = -777.
         end if
      end do
   end do
   !write(*,*) "min/max of ctp =", minval(ctp), maxval(ctp)
   deallocate(itmp_short_2d)
   if( allocated(itmp_2d) ) deallocate(itmp_2d)

   return
end subroutine read_L2_PRES

subroutine decode_nc_fname(fname, finfo, scan_mode, is_BCM, is_TEMP, is_Phase, &
                           is_HT, is_PRES, band_id, sat_id, start_time, jday)
   implicit none
   character(len=*),  intent(in)  :: fname
   character(len=18), intent(out) :: finfo
   character(len=2),  intent(out) :: scan_mode
   logical,           intent(out) :: is_BCM
   logical,           intent(out) :: is_TEMP
   logical,           intent(out) :: is_Phase
   logical,           intent(out) :: is_HT
   logical,           intent(out) :: is_PRES
   integer(i_kind),   intent(out) :: band_id
   character(len=3),  intent(out) :: sat_id
   character(len=22), intent(out) :: start_time
   integer(i_kind),   intent(out) :: jday
   integer(i_kind) :: year, month, day, hour, minute, sec1, sec2

   if ( fname( 1:14) == BCM_id ) then
      is_BCM = .true.
      band_id = -99
   else if ( fname( 1:15) == TEMP_id ) then
      is_TEMP = .true.
      band_id = -99
   else if ( fname( 1:15) == Phase_id ) then
      is_Phase = .true.
      band_id = -99
   else if ( fname( 1:15) == HT_id ) then
      is_HT = .true.
      band_id = -99
   else if ( fname( 1:14) == PRES_id ) then
      is_PRES = .true.
      band_id = -99
   else
      is_BCM = .false.
      is_TEMP = .false.
      is_Phase = .false.
      is_HT = .false.
      is_PRES = .false.
   end if

   !CG_ABI-L2-ACMC-M3_G16_s20180351202275_e20180351205060_c20180351205106.nc
   !OR_ABI-L1b-RadC-M3C16_G16_s20172741802196_e20172741804580_c20172741805015.nc
   !1234567890123456789012345678901234567890123456789012345678901234567890123456
   if ( .not. ( is_BCM .or. is_TEMP .or. is_Phase .or. is_HT .or. is_PRES ) ) then
      read(fname( 1:18), '(a18)') finfo
      read(fname(17:18), '(a2)')  scan_mode
      read(fname(20:21), '(i2)')  band_id
      read(fname(23:25), '(a3)')  sat_id
      read(fname(28:31), '(i4)')  year
      read(fname(32:34), '(i3)')  jday
      read(fname(35:36), '(i2)')  hour
      read(fname(37:38), '(i2)')  minute
      read(fname(39:40), '(i2)')  sec1   ! integer part of second
      read(fname(41:41), '(i1)')  sec2   ! decimal part of second
      ! get month and day from julian day
      call get_date(year, jday, month, day)
      ! 2017-10-01T18:02:19.6Z
      write(start_time,'(i4.4,4(a,i2.2),a,i2.2,a,i1,a)') &
            year, '-', month, '-', day, 'T', hour, ':',  minute, ':', sec1, '.', sec2, 'Z'
   else if ( is_BCM .or. is_PRES ) then
   !OR_ABI-L2-ACMF-M3_G16_s20181050000419_e20181050011186_c20181050011347.nc
   !OR_ABI-L2-CTPF-M3_G16_s20181050000419_e20181050011186_c20181050012223.nc
   !1234567890123456789012345678901234567890123456789012345678901234567890123456
      read(fname( 1:17), '(a17)') finfo
      read(fname(16:17), '(a2)')  scan_mode
      read(fname(19:21), '(a3)')  sat_id
      read(fname(24:27), '(i4)')  year
      read(fname(28:30), '(i3)')  jday
      read(fname(31:32), '(i2)')  hour
      read(fname(33:34), '(i2)')  minute
      read(fname(35:36), '(i2)')  sec1   ! integer part of second
      read(fname(37:37), '(i1)')  sec2   ! decimal part of second
      ! get month and day from julian day
      call get_date(year, jday, month, day)
      ! 2017-10-01T18:02:19.6Z
      write(start_time,'(i4.4,4(a,i2.2),a,i2.2,a,i1,a)') &
            year, '-', month, '-', day, 'T', hour, ':',  minute, ':', sec1, '.', sec2, 'Z'
   else if ( is_TEMP .or. is_Phase .or. is_HT ) then
      !OR_ABI-L2-ACHTF-M3_G16_s20181050000419_e20181050011186_c20181050012223.nc
      !OR_ABI-L2-ACTPF-M3_G16_s20181050000419_e20181050011186_c20181050011460.nc
      !OR_ABI-L2-ACHAF-M3_G16_s20181050000419_e20181050011186_c20181050012223.nc
      !1234567890123456789012345678901234567890123456789012345678901234567890123456
      read(fname( 1:18), '(a18)') finfo
      read(fname(17:18), '(a2)')  scan_mode
      read(fname(20:22), '(a3)')  sat_id
      read(fname(25:28), '(i4)')  year
      read(fname(29:31), '(i3)')  jday
      read(fname(32:33), '(i2)')  hour
      read(fname(34:35), '(i2)')  minute
      read(fname(36:37), '(i2)')  sec1   ! integer part of second
      read(fname(38:38), '(i1)')  sec2   ! decimal part of second
      ! get month and day from julian day
      call get_date(year, jday, month, day)
      ! 2017-10-01T18:02:19.6Z
      write(start_time,'(i4.4,4(a,i2.2),a,i2.2,a,i1,a)') &
            year, '-', month, '-', day, 'T', hour, ':',  minute, ':', sec1, '.', sec2, 'Z'
   else
      write(0,*) 'Error decode_nc_fname'
      stop
   end if
   return
end subroutine decode_nc_fname

 end module mod_goes_abi
