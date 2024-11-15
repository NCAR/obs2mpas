module  mod_read_write_mpas

   use netcdf
   use control_para !BJJ

   implicit none

   contains

   subroutine read_mpas_latlon (fname, nC, lon, lat, bdymask)

   implicit none
   character(len=256),                  intent( in) :: fname
   integer(i_kind),                     intent(out) :: nC
   real(r_kind), allocatable,           intent(out) :: lon(:), lat(:) 
   real(r_kind), allocatable, optional, intent(out) :: bdymask(:)
   ! loc
   integer(i_kind) :: ncid, nf_status, dimid, varid
   logical :: isfile
     
   inquire(file=trim(fname), exist=isfile)
   if ( .not. isfile ) then
      write(0,*) 'File not found: '//trim(fname)
      stop
   end if

   nf_status = nf90_OPEN(trim(fname), nf90_NOWRITE, ncid)
   if ( nf_status == 0 ) then
      write(0,*) 'Reading '//trim(fname)
   end if

   nf_status = nf90_inq_dimid(ncid, "nCells", dimid)
   nf_status = nf90_inquire_dimension(ncid, dimid, len=nC)
   if ( nf_status /= 0 ) then
      write(0,*) 'Error reading dimensions'
      stop
   end if

   allocate( lon(nC), lat(nC) )
   nf_status = nf90_INQ_VARID(ncid, 'lonCell', varid)
   nf_status = nf90_GET_VAR(ncid, varid, lon(:))
   nf_status = nf90_INQ_VARID(ncid, 'latCell', varid)
   nf_status = nf90_GET_VAR(ncid, varid, lat(:))
   !optionally read 'bdyMaskCell' for regional model domain
   nf_status = nf90_INQ_VARID(ncid, 'bdyMaskCell', varid)
   if ( nf_status == 0 ) then
      allocate( bdymask(nC) )
      nf_status = nf90_GET_VAR(ncid, varid, bdymask(:))
   end if

   nf_status = nf90_CLOSE(ncid)

   end subroutine read_mpas_latlon


   subroutine write_to_mpas (fname, nC, nfield, var, varname, sensor)

   implicit none
   character(len=256), intent(in) :: fname
   integer(i_kind),    intent(in) :: nC
   integer(i_kind),    intent(in) :: nfield
   real(r_kind),       intent(in) :: var(nC,nfield)
   character(len=64),  intent(in) :: varname(nfield)
   character(len=3),   intent(in) :: sensor
   ! loc
   integer(i_kind) :: ncid, nf_status, dimid(2), varid
   integer(i_kind) :: i
   logical :: isfile, l_newfile
     
   l_newfile=.true. !BJJ Set as default for now.

   if ( l_newfile ) then
   !----Create New File----------
      !BJJ https://docs.unidata.ucar.edu/netcdf-fortran/current/f90_datasets.html#f90-nf90_create
      !nf_status = nf90_CREATE(trim(fname), cmode=NF90_NETCDF4, ncid=ncid)
      nf_status = nf90_CREATE(trim(fname), cmode=NF90_64BIT_OFFSET, ncid=ncid) !to make compatible with MPAS-Model NetCDF files
      if ( nf_status == 0 ) then
         write(0,*) 'Writing '//trim(fname)
      else
         write(0,*) "nf90_create:",nf_status
         stop
      end if
      nf_status = nf90_def_dim(ncid,"nCells",nC,dimid(1))
      if ( nf_status /= 0 ) then; write(0,*) "nf90_def_dim:nCells:",nf_status; stop; end if
      nf_status = nf90_def_dim(ncid,"Time",NF90_UNLIMITED,dimid(2))
      if ( nf_status /= 0 ) then; write(0,*) "nf90_def_dim:Time:",nf_status; stop; end if
   else
   !----- add new variable to the existing file
      inquire(file=trim(fname), exist=isfile)
      if ( .not. isfile ) then
         write(0,*) 'File not found: '//trim(fname)
         stop
      end if
      nf_status = nf90_OPEN(trim(fname), mode=nf90_WRITE, ncid=ncid)
      if ( nf_status /= 0 ) then; write(0,*) "nf90_open:",nf_status; stop; end if
      nf_status = nf90_redef(ncid)

      if ( nf_status /= 0 ) then; write(0,*) "nf90_redef:",nf_status; stop; end if
      nf_status = nf90_inq_dimid(ncid, "nCells", dimid(1))
      if ( nf_status /= 0 ) then; write(0,*) "nf90_inq_dimid1:",nf_status;  stop; end if
      nf_status = nf90_inq_dimid(ncid, "Time", dimid(2))
      if ( nf_status /= 0 ) then; write(0,*) "nf90_inq_dimid2:",nf_status; stop; end if
   end if

   !----- applicable to both new & existing file
   do i=1,nfield
      nf_status = nf90_def_var(ncid,trim(varname(i)),NF90_FLOAT,dimid,varid)
      if ( nf_status /= 0 ) then; write(0,*) "nf90_def_var:",nf_status; stop; end if
      select case (varname(i)(1:3))
         case ('BT_')
            !float BT_G16C13(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L1b Brightness Temperature of '//varname(i)(4:9))
            nf_status = nf90_put_att(ncid,varid,'units','K')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('Rad')
            !float Rad_G16C13(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L1b Radiances of '//varname(i)(5:10))
            nf_status = nf90_put_att(ncid,varid,'units','mW m-2 sr-1 (cm-1)-1')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('BCM')
            !float BCM_G16(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L2+ Clear Sky Mask of '//varname(i)(5:7))
            nf_status = nf90_put_att(ncid,varid,'description','0=clear, 1=cloudy')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('TEM')
            !float TEMP_G16(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L2+ Cloud Top Temperature of '//varname(i)(6:8))
            nf_status = nf90_put_att(ncid,varid,'units','K')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('Pha')
            !float Phase_G16(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L2+ Cloud Top Phase of '//varname(i)(7:9))
            nf_status = nf90_put_att(ncid,varid,'description','0=clear_sky, 1=liquid_water, &
                        2=super_cooled_liquid_water, 3=mixed_phase, 4=ice, 5=unknown')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('HT_')
            !float HT_G16(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L2+ Cloud Top Height of '//varname(i)(4:6))
            nf_status = nf90_put_att(ncid,varid,'units','m')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case ('PRE')
            !float PRES_G16(Time, nCells) ;
            nf_status = nf90_put_att(ncid,varid,'long_name',sensor//' L2+ Cloud Top Pressure of '//varname(i)(6:8))
            nf_status = nf90_put_att(ncid,varid,'units','hPa')
            nf_status = nf90_put_att(ncid,varid,'_FillValue',-999.0)
         case default
      end select
   end do
   nf_status = nf90_enddef(ncid)
   if ( nf_status /= 0 ) then; write(0,*) "nf90_enddef:",nf_status; stop; end if

   do i=1,nfield
     nf_status = nf90_inq_varid(ncid,trim(varname(i)),varid)
     if ( nf_status /= 0 ) then; write(0,*) "nf90_inq_varid:",nf_status; stop; end if
     nf_status = nf90_put_var(ncid,varid,reshape(var(:,i), (/nC,1/)) )
     if ( nf_status /= 0 ) then; write(0,*) "nf90_put_var:",nf_status; stop; end if
   end do

   nf_status = nf90_close(ncid)
   if ( nf_status /= 0 ) then; write(0,*) "nf90_close:",nf_status; stop; end if

   end subroutine write_to_mpas

end module mod_read_write_mpas
