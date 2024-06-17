module  mod_read_write_indx

   use netcdf
   use control_para !BJJ

   implicit none

   contains

   subroutine write_indx (fname, nn, nS_valid, nC, interp_indx, cnt_match)

   implicit none
   character(len=256), intent(in) :: fname
   integer(i_kind),    intent(in) :: nn
   integer(i_kind),    intent(in) :: nS_valid
   integer(i_kind),    intent(in) :: nC
   integer(i_kind),    intent(in) :: interp_indx(nn,nS_valid)
   integer(i_kind),    intent(in) :: cnt_match(nC)
   ! loc
   integer(i_kind) :: ncid, nf_status, dimid(2), varid
   logical :: isfile

   inquire(file=trim(fname), exist=isfile)
   if ( isfile ) then
      write(0,*) 'File already exist: '//trim(fname)
      stop
   end if

   nf_status = nf90_CREATE(trim(fname), cmode=NF90_NETCDF4, ncid=ncid)
   if ( nf_status == 0 ) then
      write(0,*) 'Creating '//trim(fname)
   end if

   nf_status = nf90_def_dim(ncid,         "nn",       nn, dimid(1))
   nf_status = nf90_def_dim(ncid,   "nS_valid", nS_valid, dimid(2))
   nf_status = nf90_def_var(ncid,"interp_indx", NF90_INT, dimid(1:2), varid)

   nf_status = nf90_def_dim(ncid,     "nCells",       nC, dimid(1))
   nf_status = nf90_def_var(ncid,  "cnt_match", NF90_INT, dimid(1), varid)

   nf_status = nf90_enddef(ncid)

   nf_status = nf90_inq_varid(ncid,"interp_indx",varid)
   nf_status = nf90_put_var(ncid,varid,interp_indx)

   nf_status = nf90_inq_varid(ncid,  "cnt_match",varid) 
   nf_status = nf90_put_var(ncid,varid,cnt_match)
 
   nf_status = nf90_close(ncid)

   end subroutine write_indx

   subroutine read_indx (fname, nn, nS_valid, nC, interp_indx, cnt_match)

   implicit none
   character(len=256), intent( in) :: fname
   integer(i_kind),    intent( in) :: nn
   integer(i_kind),    intent( in) :: nS_valid
   integer(i_kind),    intent( in) :: nC
   integer(i_kind),    intent(out) :: interp_indx(nn,nS_valid)
   integer(i_kind),    intent(out) :: cnt_match(nC)
   ! loc
   integer(i_kind) :: ncid, nf_status, dimid(2), varid
   logical :: isfile

   inquire(file=trim(fname), exist=isfile)
   if ( .not. isfile ) then
      write(0,*) 'File not Found: '//trim(fname)
      stop
   end if

   nf_status = nf90_OPEN(trim(fname), nf90_NOWRITE, ncid)
   if ( nf_status == 0 ) then
      write(0,*) 'Reading '//trim(fname)
   end if

   nf_status = nf90_INQ_VARID(ncid, "interp_indx", varid)
   nf_status = nf90_GET_VAR(ncid, varid, interp_indx(:,:))
   nf_status = nf90_INQ_VARID(ncid, "cnt_match", varid)
   nf_status = nf90_GET_VAR(ncid, varid, cnt_match(:))

   nf_status = nf90_CLOSE(ncid)

   end subroutine read_indx

end module mod_read_write_indx
