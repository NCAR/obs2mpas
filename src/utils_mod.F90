module utils_mod

   use netcdf_mod, only: open_netcdf_for_write, close_netcdf, &
      def_netcdf_dims, def_netcdf_var, def_netcdf_end, &
      put_netcdf_var

   use control_para !BJJ

   implicit none
   include 'netcdf.inc'

   contains

subroutine interpret_qa_flags(qa_value, clear_flag, cloudy_flag)
  implicit none

  integer, intent(in)  :: qa_value
  integer, intent(out) :: clear_flag
  integer, intent(out) :: cloudy_flag

  integer :: cloud_mask, retrieval_phase

  clear_flag = 0
  cloudy_flag = 0

  ! extract the cloud mask confidence level (bits 4,3) from qa_value, bits 4,3
  cloud_mask = mod(shift_right_logical(qa_value, 3), 4)

  ! extract the cloud retrieval phase flag (bits 6,5) from qa_value, bits 6,5
  retrieval_phase = mod(shift_right_logical(qa_value, 5), 4)

  ! determine if the pixel is confidently clear or confidently cloudy
  ! Clear pixel (Cloud Mask = 00 and Retrieval Phase = 00)
  if (cloud_mask == 0 .and. retrieval_phase == 0) then
     clear_flag = 1  ! Confidently clear

  ! Cloudy pixel (Cloud Mask = 11 and any Retrieval Phase)
  elseif (cloud_mask == 3) then
     cloudy_flag = 1  ! Confidently cloudy
  end if

end subroutine interpret_qa_flags

integer function shift_right_logical(value, n)
  implicit none
  integer, intent(in) :: value, n

  ! perform logical right shift
  shift_right_logical = ishft(value, -n)
end function shift_right_logical


! convert a string to uppercase
function to_upper(s)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: to_upper
    integer :: i

    to_upper = s
    do i = 1, len(s)
        if (iachar(to_upper(i:i)) >= iachar('a') .and. iachar(to_upper(i:i)) <= iachar('z')) then
            to_upper(i:i) = achar(iachar(to_upper(i:i)) - iachar('a') + iachar('A'))
        end if
    end do
end function to_upper

! check if a string is empty (contains only spaces)
logical function is_empty_string(s)
  character(len=*) :: s
  integer :: i

  ! loop through each character and check if any is non-space
  do i = 1, len(s)
    if (s(i:i) /= ' ') then
       is_empty_string = .false.
       return
    end if
  end do

  ! if no non-space characters are found, it is empty
  is_empty_string = .true.
end function is_empty_string

subroutine get_namelist_vars(nfile, fnames, list_files, data_dir, data_id, sat_id, n_subsample, write_iodav1)

   implicit none

   integer(i_kind),                 intent(out) ::  nfile
   character(len=256), allocatable, intent(out) ::  fnames(:)

   integer(i_kind)                 :: nml_unit = 81
   integer(i_kind)                 :: iunit    = 87

   integer(i_kind)                 :: istat, ifile
   logical                         :: isfile

   character(len=256)              :: txtbuf
   character(len=256)              :: list_files  ! the text file that contains a list of netcdf files to process
   character(len=256)              :: data_dir
   character(len=18)               :: data_id
   character(len=3)                :: sat_id
   integer(i_kind)                 :: n_subsample
   logical                         :: write_iodav1

   namelist /data_nml/ list_files, data_dir, data_id, sat_id, n_subsample, write_iodav1

   ! initialize namelist variables from namelist.obs2model
   ! &data_nml
   !    list_files = 'flist.txt'               ! list of file names (no path)
   !    data_dir = '/data/goes/or/himawari',   ! path of the observation files
   !    data_id = '<HS> <OR_ABI-L1b-RadC-M3>', ! prefix of the downloaded files
   !    sat_id = '<H08> <G16> <G18>',          ! satellite ID
   !    n_subsample = 1,                 ! value use for thinning if write_iodav1 = .true.
   !    write_iodav1 = .false.,          ! option to write out an iodav1 file (no superobbing)
   !    /

   ! initialize namelist variable
   ! Note: default values are for ABI observations processing
   list_files        = 'flist.txt'
   data_dir          = '.'
   data_id           = 'OR_ABI-L1b-RadC-M3'
   sat_id            = 'G16'
   n_subsample       = 1
   write_iodav1      = .false.

   ! read namelist
   open(unit=nml_unit, file='namelist.obs2model', status='old', form='formatted')
   read(unit=nml_unit, nml=data_nml, iostat=istat)
   write(0,nml=data_nml)
   if ( istat /= 0 ) then
      write(0,*) 'Error reading namelist data_nml', istat
      stop
   end if

   ! get file names from list_files
   nfile  = 0  ! initialize the number of netcdf files to read
   inquire(file=trim(list_files), exist=isfile)
   if ( .not. isfile ) then
      write(0,*) 'File not found: list_files '//trim(list_files)
      stop 1
   else
      open(unit=iunit, file=trim(list_files), status='old', form='formatted')
      !first find out the number of netcdf files to read
      istat = 0
      do while ( istat == 0 )
         read(unit=iunit, fmt='(a)', iostat=istat) txtbuf
         if ( istat /= 0 ) then
            exit
         else
            nfile = nfile + 1
         end if
      end do
      if ( nfile > 0 ) then
         allocate (fnames(nfile))
         !read the list_files again to get the file names
         rewind(iunit)
         do ifile = 1, nfile
            read(unit=iunit, fmt='(a)', iostat=istat) fnames(ifile)
         end do
      else
         write(0,*) 'File not found from list_files '//trim(list_files)
         stop
      end if
      close(iunit)
   end if !list_files

end subroutine

subroutine get_date(ccyy, jday, month, day)
   implicit none
   integer(i_kind), intent(in)  :: ccyy, jday
   integer(i_kind), intent(out) :: month, day
   integer(i_kind) :: mmday(12) = (/31,28,31,30,31,30,31,31,30,31,30,31/)
   integer(i_kind) :: i, jdtmp
   continue

   if ( MOD(ccyy,4) == 0 ) then
      mmday(2) = 29
      if ( MOD(ccyy,100) == 0 ) then
         mmday(2) = 28
      end if
      if ( MOD(ccyy,400) == 0 ) then
         mmday(2) = 29
      end if
   end if

   jdtmp = 0
   do i = 1, 12
      jdtmp = jdtmp + mmday(i)
      if ( jday <= jdtmp ) then
         month = i
         day = jday - ( jdtmp - mmday(i) )
         exit
      end if
   end do

   return
end subroutine get_date

subroutine read_GRB_dims(ncid, xname, yname, nx, ny)
   implicit none
   integer(i_kind),   intent(in) :: ncid
   character(len=*),  intent(in) :: xname
   character(len=*),  intent(in) :: yname
   integer(i_kind),   intent(out):: nx, ny
   integer(i_kind)               :: dimid
   integer(i_kind)               :: nf_status(4)
   continue
   nf_status(1) = nf_INQ_DIMID(ncid, xname, dimid)
   nf_status(2) = nf_INQ_DIMLEN(ncid, dimid, nx)
   nf_status(3) = nf_INQ_DIMID(ncid, yname, dimid)
   nf_status(4) = nf_INQ_DIMLEN(ncid, dimid, ny)
   if ( any(nf_status /= 0) ) then
      write(0,*) 'Error reading dimensions'
      stop
   end if
   return
end subroutine read_GRB_dims

subroutine output_iodav1( fname, time_start, nx, ny, nband, n_subsample, got_latlon, lat, lon, sat_zen, sun_zen, bt, qf, sdtb, cloudmask )

   implicit none

   character(len=*),   intent(in) :: fname
   character(len=22),  intent(in) :: time_start
   integer(i_kind),    intent(in) :: nx, ny, nband
   integer(i_kind),    intent(in) :: n_subsample
   logical,            intent(in) :: got_latlon(nx,ny)
   real(r_kind),       intent(in) :: lat(nx,ny)
   real(r_kind),       intent(in) :: lon(nx,ny)
   real(r_kind),       intent(in) :: sat_zen(nx,ny)
   real(r_kind),       intent(in) :: sun_zen(nx,ny)
   real(r_kind),       intent(in) :: bt(nband,nx,ny)
   integer(i_kind),    intent(in) :: qf(nband,nx,ny)
   real(r_kind),       intent(in) :: sdtb(nband)
   integer(i_kind),    intent(in), optional :: cloudmask(nx,ny)

   integer(i_kind), parameter :: nstring = 50
   integer(i_kind), parameter :: ndatetime = 20
   integer(i_kind) :: nvars
   integer(i_kind) :: nlocs

   character(len=ndatetime), allocatable  :: datetime(:)   ! ccyy-mm-ddThh:mm:ssZ
   real(r_kind), allocatable :: lat_out(:)
   real(r_kind), allocatable :: lon_out(:)
   real(r_kind), allocatable :: scan_pos_out(:)
   real(r_kind), allocatable :: sat_zen_out(:)
   real(r_kind), allocatable :: sun_zen_out(:)
   real(r_kind), allocatable :: sat_azi_out(:)
   real(r_kind), allocatable :: sun_azi_out(:)
   real(r_kind), allocatable :: bt_out(:,:)
   real(r_kind), allocatable :: err_out(:,:)
   real(r_kind), allocatable :: qf_out(:,:)

   integer(i_kind) :: ncid_nlocs
   integer(i_kind) :: ncid_nvars
   integer(i_kind) :: ncid_nstring
   integer(i_kind) :: ncid_ndatetime
   integer(i_kind) :: ncfileid
   character(len=nstring) :: ncname

   character(len=nstring), allocatable :: name_var_tb(:)
   character(len=4) :: c4

   integer(i_kind) :: iline, isample, iband
   integer(i_kind) :: i, iloc
   integer(i_kind) :: iyear, imonth, iday, ihour, imin, isec

   character(len=60), parameter :: var_tb = "brightness_temperature"

   nvars = nband

   nlocs = 0
   do iline = 1, ny, n_subsample
      do isample = 1, nx, n_subsample
         if ( .not. got_latlon(isample,iline) ) cycle
         if ( sat_zen(isample,iline) > 80.0 ) cycle
         ! qf (DQF, Data Quality Flag)
         ! 0:good, 1:conditionally_usable, 2:out_of_range, 3:no_value
         ! keep only qf=0,1 pixels
         if ( all(qf(:,isample,iline) > 1) ) cycle
         if ( all(bt(:,isample,iline)<0.0) ) cycle
         nlocs = nlocs + 1
      end do
   end do

   write(0,*) 'nlocs = ', nlocs
   if ( nlocs <= 0 ) then
      return
   end if

   allocate (name_var_tb(1:nband))
   allocate (datetime(nlocs))
   allocate (lat_out(nlocs))
   allocate (lon_out(nlocs))
   allocate (scan_pos_out(nlocs))
   allocate (sat_zen_out(nlocs))
   allocate (sat_azi_out(nlocs))
   allocate (sun_zen_out(nlocs))
   allocate (sun_azi_out(nlocs))
   allocate (bt_out(nband,nlocs))
   allocate (err_out(nband,nlocs))
   allocate (qf_out(nband,nlocs))

   read(time_start( 1: 4), '(i4)') iyear
   read(time_start( 6: 7), '(i2)') imonth
   read(time_start( 9:10), '(i2)') iday
   read(time_start(12:13), '(i2)') ihour
   read(time_start(15:16), '(i2)') imin
   read(time_start(18:19), '(i2)') isec

   iloc = 0
   do iline = 1, ny, n_subsample
      do isample = 1, nx, n_subsample
         if ( .not. got_latlon(isample,iline) ) cycle
         if ( sat_zen(isample,iline) > 80.0 ) cycle
         if ( all(qf(:,isample,iline) > 1) ) cycle
         if ( all(bt(:,isample,iline)<0.0) ) cycle
         iloc = iloc + 1
         write(unit=datetime(iloc), fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
               iyear, '-', imonth, '-', iday, 'T', ihour, ':', imin, ':', isec, 'Z'
         lat_out(iloc) = lat(isample,iline)
         lon_out(iloc) = lon(isample,iline)
         sat_zen_out(iloc) = sat_zen(isample,iline)
         sun_zen_out(iloc) = sun_zen(isample,iline)
         bt_out(1:nband,iloc) = bt(1:nband,isample,iline)
         qf_out(1:nband,iloc) = qf(1:nband,isample,iline)
         scan_pos_out(iloc) = isample
         sat_azi_out(iloc) = missing_r
         sun_azi_out(iloc) = missing_r
         err_out(1:nband,iloc) = 1.0 !missing_r
      end do
   end do

   call open_netcdf_for_write(trim(fname),ncfileid)
   call def_netcdf_dims(ncfileid,'nvars',nvars,ncid_nvars)
   call def_netcdf_dims(ncfileid,'nlocs',nlocs,ncid_nlocs)
   call def_netcdf_dims(ncfileid,'nstring',nstring,ncid_nstring)
   call def_netcdf_dims(ncfileid,'ndatetime',ndatetime,ncid_ndatetime)
   do i = 1, nvars
      write(unit=c4, fmt='(i4)') i+6
      name_var_tb(i) = trim(var_tb)//'_'//trim(adjustl(c4))
      ncname = trim(name_var_tb(i))//'@ObsValue'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT,'units','K')
      ncname = trim(name_var_tb(i))//'@ObsError'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
      ncname = trim(name_var_tb(i))//'@PreQC'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_INT)
   end do
   ncname = 'latitude@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'longitude@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'solar_azimuth_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'scan_position@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_azimuth_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'solar_zenith_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_zenith_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_view_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_channel@VarMetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nvars/),NF_INT)
   ncname = 'datetime@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_ndatetime,ncid_nlocs/),NF_CHAR)
   ncname = 'variable_names@VarMetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nstring,ncid_nvars/),NF_CHAR)
   call def_netcdf_end(ncfileid)

   do i = 1, nvars
      ncname = trim(name_var_tb(i))//'@ObsValue'
      call put_netcdf_var(ncfileid,ncname,bt_out(i,:))
      ncname = trim(name_var_tb(i))//'@ObsError'
      call put_netcdf_var(ncfileid,ncname,err_out(i,:))
      ncname = trim(name_var_tb(i))//'@PreQC'
      call put_netcdf_var(ncfileid,ncname,qf_out(i,:))
   end do

   ncname = 'latitude@MetaData'
   call put_netcdf_var(ncfileid,ncname,lat_out)
   ncname = 'longitude@MetaData'
   call put_netcdf_var(ncfileid,ncname,lon_out)
   ncname = 'solar_azimuth_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sun_azi_out)
   ncname = 'scan_position@MetaData'
   call put_netcdf_var(ncfileid,ncname,scan_pos_out)
   ncname = 'sensor_azimuth_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_azi_out)
   ncname = 'solar_zenith_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sun_zen_out)
   ncname = 'sensor_zenith_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_zen_out)
   ncname = 'sensor_view_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_zen_out)
   ncname = 'sensor_channel@VarMetaData'
   call put_netcdf_var(ncfileid,ncname,(/7,8,9,10,11,12,13,14,15,16/))
   ncname = 'datetime@MetaData'
   call put_netcdf_var(ncfileid,ncname,datetime)
   ncname = 'variable_names@VarMetaData'
   call put_netcdf_var(ncfileid,ncname,name_var_tb(1:nband))
   call close_netcdf(trim(fname),ncfileid)

   deallocate (name_var_tb)
   deallocate (datetime)
   deallocate (lat_out)
   deallocate (lon_out)
   deallocate (scan_pos_out)
   deallocate (sat_zen_out)
   deallocate (sat_azi_out)
   deallocate (sun_zen_out)
   deallocate (sun_azi_out)
   deallocate (bt_out)
   deallocate (err_out)
   deallocate (qf_out)

end subroutine output_iodav1

subroutine output_iodav1_o2m(fname, time_start, nC, nband, got_latlon, lat, lon, sat_zen, sun_zen, bt, bt_std)

   implicit none

   character(len=*),   intent(in) :: fname
   character(len=22),  intent(in) :: time_start
   integer(i_kind),    intent(in) :: nC, nband
   logical,            intent(in) :: got_latlon(nC)
   real(r_kind),       intent(in) :: lat(nC)
   real(r_kind),       intent(in) :: lon(nC)
   real(r_kind),       intent(in) :: sat_zen(nC)
   real(r_kind),       intent(in) :: sun_zen(nC)
   real(r_kind),       intent(in) :: bt(nband+3,nC) !BJJ 1:nband for bt, nband+1 for 2d cloud fraction
   real(r_kind),       intent(in) :: bt_std(nband+3,nC) !BJJ 1:nband for bt, nband+1 for # of obs for SO  
! bt from nearest neighbor; bt_std from l_superob
   integer(i_kind), parameter :: nstring = 50
   integer(i_kind), parameter :: ndatetime = 20
   integer(i_kind) :: nvars
   integer(i_kind) :: nlocs

   character(len=ndatetime), allocatable  :: datetime(:)   ! ccyy-mm-ddThh:mm:ssZ
   real(r_kind), allocatable :: lat_out(:)
   real(r_kind), allocatable :: lon_out(:)
   real(r_kind), allocatable :: scan_pos_out(:)
   real(r_kind), allocatable :: sat_zen_out(:)
   real(r_kind), allocatable :: sun_zen_out(:)
   real(r_kind), allocatable :: sat_azi_out(:)
   real(r_kind), allocatable :: sun_azi_out(:)
   real(r_kind), allocatable :: bt_out(:,:)
   real(r_kind), allocatable :: bt_std_out(:,:)   
   real(r_kind), allocatable :: err_out(:,:)
   real(r_kind), allocatable :: qf_out(:,:)
   integer(i_kind), allocatable :: iC_out(:)  !BJJ for cellIndex@MetaData

   integer(i_kind) :: ncid_nlocs
   integer(i_kind) :: ncid_nvars
   integer(i_kind) :: ncid_nstring
   integer(i_kind) :: ncid_ndatetime
   integer(i_kind) :: ncfileid
   character(len=nstring) :: ncname

   character(len=nstring), allocatable :: name_var_tb(:)
   character(len=4) :: c4

   integer(i_kind) :: i, iC, iband
   integer(i_kind) :: iloc
   integer(i_kind) :: iyear, imonth, iday, ihour, imin, isec

   character(len=60), parameter :: var_tb = "brightness_temperature"

   nvars = nband
   nlocs = 0
   do iC = 1, nC
      if ( .not. got_latlon(iC) ) cycle
      if ( sat_zen(iC) >= 70.0 ) cycle !BJJ 80 -> 70 for consistency btw tb & cm
      if ( all(bt(:,iC)<0.0) ) cycle
      nlocs = nlocs + 1
   end do

   write(0,*) 'nlocs = ', nlocs
   if ( nlocs <= 0 ) then
      return
   end if

   allocate (name_var_tb(1:nband))
   allocate (datetime(nlocs))
   allocate (lat_out(nlocs))
   allocate (lon_out(nlocs))
   allocate (scan_pos_out(nlocs))
   allocate (sat_zen_out(nlocs))
   allocate (sat_azi_out(nlocs))
   allocate (sun_zen_out(nlocs))
   allocate (sun_azi_out(nlocs))
   allocate (bt_out(nband+3,nlocs)) 
   allocate (bt_std_out(nband+3,nlocs)) 
   allocate (err_out(nband,nlocs))
   allocate (qf_out(nband,nlocs))
   allocate (iC_out(nlocs))  !BJJ for cellIndex@MetaData

   read(time_start( 1: 4), '(i4)') iyear
   read(time_start( 6: 7), '(i2)') imonth
   read(time_start( 9:10), '(i2)') iday
   read(time_start(12:13), '(i2)') ihour
   read(time_start(15:16), '(i2)') imin
   read(time_start(18:19), '(i2)') isec

   iloc = 0
   do iC = 1, nC
      if ( .not. got_latlon(iC) ) cycle
      if ( sat_zen(iC) >= 70.0 ) cycle !BJJ 80 -> 70 for consistency btw tb & cm
      if ( all(bt(:,iC)<0.0) ) cycle
      iloc = iloc + 1
      write(unit=datetime(iloc), fmt='(i4,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a,i2.2,a)')  &
            iyear, '-', imonth, '-', iday, 'T', ihour, ':', imin, ':', isec, 'Z'
      lat_out(iloc) = lat(iC)
      lon_out(iloc) = lon(iC)
      sat_zen_out(iloc) = sat_zen(iC)
      sun_zen_out(iloc) = sun_zen(iC)
      bt_out(1:nband+1,iloc) = bt(1:nband+1,iC) !nband+1 for 2d clm
      bt_std_out(1:nband+1,iloc) = bt_std(1:nband+1,iC) !nband+1 for 2d clm
      bt_out(1:nband+2,iloc) = bt(1:nband+2,iC) !nband+2 for 2d CloudTopPres
      bt_std_out(1:nband+2,iloc) = bt_std(1:nband+2,iC) !nband+2 for 2d CloudTopPres
      bt_out(1:nband+3,iloc) = bt(1:nband+3,iC) !nband+3 for 2d CloudType
      bt_std_out(1:nband+3,iloc) = bt_std(1:nband+3,iC) !nband+3 for 2d CloudType	  
      qf_out(1:nband,iloc) = 0.0 ! BJJ what this can be for superob/nearest obs ?
      scan_pos_out(iloc) = 0.0   ! BJJ what this can be ?
      sat_azi_out(iloc) = missing_r
      sun_azi_out(iloc) = missing_r
      err_out(1:nband,iloc) = 1.0 !missing_r
      iC_out(iloc) = iC !BJJ
   end do

   call open_netcdf_for_write(trim(fname),ncfileid)
   call def_netcdf_dims(ncfileid,'nvars',nvars,ncid_nvars)
   call def_netcdf_dims(ncfileid,'nlocs',nlocs,ncid_nlocs)
   call def_netcdf_dims(ncfileid,'nstring',nstring,ncid_nstring)
   call def_netcdf_dims(ncfileid,'ndatetime',ndatetime,ncid_ndatetime)
   do i = 1, nvars
      write(unit=c4, fmt='(i4)') i+6
      name_var_tb(i) = trim(var_tb)//'_'//trim(adjustl(c4))
      ncname = trim(name_var_tb(i))//'@ObsValue'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT,'units','K')
      ncname = trim(name_var_tb(i))//'@ObsError'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
      ncname = trim(name_var_tb(i))//'@PreQC'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_INT)
      ncname = trim(name_var_tb(i))//'_so_std@MetaData'
      call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT,'units','K')
   end do
   ncname = 'latitude@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'longitude@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'solar_azimuth_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'scan_position@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_azimuth_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'solar_zenith_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_zenith_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_view_angle@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'sensor_channel@VarMetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nvars/),NF_INT)
   ncname = 'datetime@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_ndatetime,ncid_nlocs/),NF_CHAR)
   ncname = 'variable_names@VarMetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nstring,ncid_nvars/),NF_CHAR)
   ncname = 'cloudAmount@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'cldTopPres@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'cloudType@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)   
   ncname = 'obsNumerForSO@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_FLOAT)
   ncname = 'cellIndex@MetaData'
   call def_netcdf_var(ncfileid,ncname,(/ncid_nlocs/),NF_INT)
   call def_netcdf_end(ncfileid)

   do i = 1, nvars
      ncname = trim(name_var_tb(i))//'@ObsValue'
      call put_netcdf_var(ncfileid,ncname,bt_out(i,:))
      ncname = trim(name_var_tb(i))//'@ObsError'
      call put_netcdf_var(ncfileid,ncname,err_out(i,:))
      ncname = trim(name_var_tb(i))//'@PreQC'
      call put_netcdf_var(ncfileid,ncname,qf_out(i,:))
      ncname = trim(name_var_tb(i))//'_so_std@MetaData'
      call put_netcdf_var(ncfileid,ncname,bt_std_out(i,:))
   end do

   ncname = 'latitude@MetaData'
   call put_netcdf_var(ncfileid,ncname,lat_out)
   ncname = 'longitude@MetaData'
   call put_netcdf_var(ncfileid,ncname,lon_out)
   ncname = 'solar_azimuth_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sun_azi_out)
   ncname = 'scan_position@MetaData'
   call put_netcdf_var(ncfileid,ncname,scan_pos_out)
   ncname = 'sensor_azimuth_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_azi_out)
   ncname = 'solar_zenith_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sun_zen_out)
   ncname = 'sensor_zenith_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_zen_out)
   ncname = 'sensor_view_angle@MetaData'
   call put_netcdf_var(ncfileid,ncname,sat_zen_out)
   ncname = 'sensor_channel@VarMetaData'
   call put_netcdf_var(ncfileid,ncname,(/7,8,9,10,11,12,13,14,15,16/))
   ncname = 'datetime@MetaData'
   call put_netcdf_var(ncfileid,ncname,datetime)
   ncname = 'variable_names@VarMetaData'
   call put_netcdf_var(ncfileid,ncname,name_var_tb(1:nband))
   ncname = 'cloudAmount@MetaData'
   call put_netcdf_var(ncfileid,ncname,bt_out(nband+1,:))
   write(*,*) "check min/max of clm =", minval(bt_out(nband+1,:)), maxval(bt_out(nband+1,:))
   ncname = 'cldTopPres@MetaData'
   call put_netcdf_var(ncfileid,ncname,bt_out(nband+2,:))
   write(*,*) "check min/max of cloud top press =", minval(bt_out(nband+2,:)), maxval(bt_out(nband+2,:))   
   ncname = 'cloudType@MetaData'
   call put_netcdf_var(ncfileid,ncname,bt_out(nband+3,:))   
   write(*,*) "check min/max of cloud type =", minval(bt_out(nband+3,:)), maxval(bt_out(nband+3,:))   
   ncname = 'obsNumerForSO@MetaData'
   call put_netcdf_var(ncfileid,ncname,bt_std_out(nband+1,:))
   ncname = 'cellIndex@MetaData'
   call put_netcdf_var(ncfileid,ncname,iC_out)
   call close_netcdf(trim(fname),ncfileid)

   deallocate (name_var_tb)
   deallocate (datetime)
   deallocate (lat_out)
   deallocate (lon_out)
   deallocate (scan_pos_out)
   deallocate (sat_zen_out)
   deallocate (sat_azi_out)
   deallocate (sun_zen_out)
   deallocate (sun_azi_out)
   deallocate (bt_out)
   deallocate (bt_std_out)
   deallocate (err_out)
   deallocate (qf_out)
   deallocate (iC_out)

end subroutine output_iodav1_o2m

subroutine calc_solar_zenith_angle(nx, ny, xlat, xlon, xtime, julian, solzen, got_latlon)

! the calulcation is adapted from subroutines radconst and calc_coszen in
! WRF phys/module_radiation_driver.F

   implicit none

   integer(i_kind),   intent(in)    :: nx, ny, julian
   real(r_kind),      intent(in)    :: xlat(nx,ny), xlon(nx,ny)
   character(len=22), intent(in)    :: xtime
   real(r_kind),      intent(inout) :: solzen(nx,ny)
   logical,           intent(in)    :: got_latlon(nx,ny)

   real(r_kind) :: obliq = 23.5
   real(r_kind) :: deg_per_day = 360.0/365.0
   real(r_kind) :: slon   ! longitude of the sun
   real(r_kind) :: declin ! declination of the sun
   real(r_kind) :: hrang, da, eot, xt, tloctm, rlat
   integer(i_kind) :: gmt, minute, i, j

   ! calculate longitude of the sun from vernal equinox
   if ( julian >= 80 ) slon = (julian - 80 ) * deg_per_day
   if ( julian <  80 ) slon = (julian + 285) * deg_per_day

   declin = asin(sin(obliq*deg2rad)*sin(slon*deg2rad)) ! in radian

   read(xtime(12:13), '(i2)') gmt
   read(xtime(15:16), '(i2)') minute

   da = 6.2831853071795862*(julian-1)/365.
   eot = (0.000075+0.001868*cos(da)-0.032077*sin(da) &
          -0.014615*cos(2.0*da)-0.04089*sin(2.0*da))*(229.18)
   xt = gmt + (minute + eot)/60.0

   do j = 1, ny
      do i = 1, nx
         if ( .not. got_latlon(i,j) ) cycle
         tloctm = xt + xlon(i,j)/15.0
         hrang = 15.0*(tloctm-12.0) * deg2rad
         rlat = xlat(i,j) * deg2rad
         solzen(i,j) = acos( sin(rlat)*sin(declin) + &
                             cos(rlat)*cos(declin)*cos(hrang) )
         solzen(i,j) = solzen(i,j) * rad2deg
      end do
   end do

   return
end subroutine calc_solar_zenith_angle

subroutine calc_geostationary_satellite_zenith_angle( rlat, rlon, lon_sat, r_eq, h_sat, rzen )
   implicit none
   real(r_kind),   intent(in)  :: rlat    ! in [radian]
   real(r_kind),   intent(in)  :: rlon    ! in [radian]
   real(r_double), intent(in)  :: lon_sat ! satellite longitude, longitude_of_projection_origin
   real(r_double), intent(in)  :: r_eq    ! GRS80 semi-major axis of earth
   real(r_double), intent(in)  :: h_sat   ! satellite height
   real(r_kind),   intent(out) :: rzen    ! in [radian]
   real(r_kind) :: lon_diff, tmp1, theta1, theta2

   lon_diff = abs(rlon-lon_sat)
!   tmp1 = sqrt((2.0*r_eq*sin(lon_diff/2.)-r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2 &
!     +(2.0*r_eq*sin(rlat/2.))**2-(r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2)
   tmp1 = (2.0*r_eq*sin(lon_diff/2.)-r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2 &
     +(2.0*r_eq*sin(rlat/2.))**2-(r_eq*(1.0-cos(rlat))*sin(lon_diff/2.))**2
   if ( tmp1 < 0.0 ) return
   tmp1 = sqrt(tmp1)
   theta1 = 2.0*asin(tmp1/r_eq/2.)
   theta2 = atan(r_eq*sin(theta1)/((h_sat-r_eq)+r_eq*(1.0-sin(theta1))))
   rzen = theta1+theta2
   !gzen(i,j) = 90.0 - atan((cos(lon_diff)*cos(rlat)-0.1512)/(sqrt(1.0-cos(lon_diff)*cos(lon_diff)*cos(rlat)*cos(rlat)))) * rad2deg

   return
end subroutine calc_geostationary_satellite_zenith_angle

!     This subroutine handles errors by printing an error message and
!     exiting with a non-zero status.
subroutine check(errcode)
    use netcdf
    implicit none
    integer, intent(in) :: errcode

    if(errcode /= nf90_noerr) then
       print *, 'Error: ', trim(nf90_strerror(errcode))
       stop 2
    end if
end subroutine check

end module utils_mod
