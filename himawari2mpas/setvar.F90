!-----------------------------------------------------
! initially designed and written
! by Yonggang G. Yu (yuxxx135@umn.edu) 31-Oct-2022
!-----------------------------------------------------

module kinds
   implicit none
   save

   integer, parameter, public :: sp = selected_real_kind(6)    ! REAL*4
   integer, parameter, public :: dp = selected_real_kind(15)   ! REAL*8
#ifdef DOUBLEPRECISION
   integer, parameter, public :: rt = dp
#else
   integer, parameter, public :: rt = dp    ! default
#endif
end module kinds


module control_para
   implicit none
   save

   real*8, parameter :: pi      = dacos(-1.d0)
   real*8, parameter :: deg2rad = pi/180.0d0
   real*8, parameter :: rad2deg = 180.0d0/pi
   integer           :: ierr
   integer           :: tval(8)      ! date/time

   !BJJ from goes_abi_mod.F90
   integer, parameter  :: r_single = selected_real_kind(6)  ! single precision
   integer, parameter  :: r_double = selected_real_kind(15) ! double precision
   integer, parameter  :: i_byte   = selected_int_kind(1)   ! byte integer
   integer, parameter  :: i_short  = selected_int_kind(4)   ! short integer
   integer, parameter  :: i_long   = selected_int_kind(8)   ! long integer
   integer, parameter  :: i_kind   = i_long                 ! default integer
   integer, parameter  :: r_kind   = r_single               ! default real
!BJJ   integer, parameter  :: r_kind   = r_double               ! default real

end module control_para
