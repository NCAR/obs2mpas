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
   integer, parameter  :: r_single = selected_real_kind(6)           ! single precision
   integer, parameter  :: r_double = selected_real_kind(15)          ! double precision
   integer, parameter  :: i_byte   = selected_int_kind(1)            ! byte integer
   integer, parameter  :: i_short  = selected_int_kind(4)            ! short integer
   integer, parameter  :: i_long   = selected_int_kind(8)            ! long integer
   integer, parameter  :: i_kind   = i_long                          ! default integer
   integer, parameter  :: r_kind   = r_single                        ! default real
   integer, parameter, private :: llong_t = selected_int_kind(16)    ! MRI -- copied from kinds.f90 [obs2ioda-v2]
   integer, parameter, public  :: i_llong = max( llong_t, i_long )   ! MRI -- copied from kinds.f90 [obs2ioda-v2]
   real(r_kind),    parameter :: missing_r         = -999.0          ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: missing_i         = -999            ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: nstring           = 50              ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: ndatetime         = 20              ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: ninst             = 17              ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: nvar_info         = 9               ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: nsen_info         = 7               ! MRI -- copied from define_mod.f90 [obs2ioda-v2]
   integer(i_kind), parameter :: StrLen            = 512             ! MRI -- copied from define_mod.f90 [obs2ioda-v2]


!BJJ   integer, parameter  :: r_kind   = r_double               ! default real

end module control_para
