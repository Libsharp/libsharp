module sharp
  use iso_c_binding
  implicit none
  ! alm_info flags
  integer, parameter :: SHARP_PACKED = 1

  ! sharp job types
  enum, bind(c)
      enumerator :: SHARP_YtW = 0
      enumerator :: SHARP_Y = 1
      enumerator :: SHARP_Yt = 2
      enumerator :: SHARP_WY = 3
      enumerator :: SHARP_ALM2MAP_DERIV1 = 4
   end enum

  ! sharp job flags
  integer, parameter :: SHARP_DP             = ISHFT(1, 4)
  integer, parameter :: SHARP_ADD            = ISHFT(1, 5)
  integer, parameter :: SHARP_REAL_HARMONICS = ISHFT(1, 6)
  integer, parameter :: SHARP_NO_FFT         = ISHFT(1, 7)

  type sharp_geom_info
     type(c_ptr) :: handle
     integer(c_ptrdiff_t) :: n_local
  end type sharp_geom_info

  type sharp_alm_info
     type(c_ptr) :: handle
     integer(c_ptrdiff_t) :: n_local
  end type sharp_alm_info

  interface

     ! alm_info
     subroutine sharp_make_general_alm_info( &
         lmax, nm, stride, mval, mvstart, flags, alm_info) bind(c)
       use iso_c_binding
       integer(c_int), value, intent(in)    :: lmax, nm, stride, flags
       integer(c_int), intent(in)           :: mval(nm)
       integer(c_ptrdiff_t), intent(in)     :: mvstart(nm)
       type(c_ptr), intent(out)             :: alm_info
     end subroutine sharp_make_general_alm_info

     subroutine c_sharp_make_mmajor_real_packed_alm_info( &
         lmax, stride, nm, ms, alm_info) bind(c, name='sharp_make_mmajor_real_packed_alm_info')
       use iso_c_binding
       integer(c_int), value, intent(in)    :: lmax, nm, stride
       integer(c_int), intent(in), optional :: ms(nm)
       type(c_ptr), intent(out)             :: alm_info
     end subroutine c_sharp_make_mmajor_real_packed_alm_info

     function c_sharp_alm_count(alm_info) bind(c, name='sharp_alm_count')
       use iso_c_binding
       integer(c_ptrdiff_t)           :: c_sharp_alm_count
       type(c_ptr), value, intent(in) :: alm_info
     end function c_sharp_alm_count

     subroutine c_sharp_destroy_alm_info(alm_info) bind(c, name='sharp_destroy_alm_info')
       use iso_c_binding
       type(c_ptr), value                   :: alm_info
     end subroutine c_sharp_destroy_alm_info

     ! geom_info
     subroutine sharp_make_subset_healpix_geom_info ( &
          nside, stride, nrings, rings, weight, geom_info) bind(c)
       use iso_c_binding
       integer(c_int), value, intent(in)    :: nside, stride, nrings
       integer(c_int), intent(in), optional :: rings(nrings)
       real(c_double), intent(in), optional :: weight(2 * nside)
       type(c_ptr), intent(out)             :: geom_info
     end subroutine sharp_make_subset_healpix_geom_info

     subroutine c_sharp_destroy_geom_info(geom_info) bind(c, name='sharp_destroy_geom_info')
       use iso_c_binding
       type(c_ptr), value                   :: geom_info
     end subroutine c_sharp_destroy_geom_info

     function c_sharp_map_size(info) bind(c, name='sharp_map_size')
       use iso_c_binding
       integer(c_ptrdiff_t) :: c_sharp_map_size
       type(c_ptr), value   :: info
     end function c_sharp_map_size


     ! execute
     subroutine c_sharp_execute(type, spin, alm, map, geom_info, alm_info, ntrans, &
                                flags, time, opcnt) bind(c, name='sharp_execute')
       use iso_c_binding
       integer(c_int), value                        :: type, spin, ntrans, flags
       type(c_ptr), value                           :: alm_info, geom_info
       real(c_double), intent(out), optional        :: time
       integer(c_long_long), intent(out), optional  :: opcnt
       type(c_ptr), intent(in)                      :: alm(*), map(*)
     end subroutine c_sharp_execute

     subroutine c_sharp_execute_mpi(comm, type, spin, alm, map, geom_info, alm_info, ntrans, &
                                    flags, time, opcnt) bind(c, name='sharp_execute_mpi_fortran')
       use iso_c_binding
       integer(c_int), value                        :: comm, type, spin, ntrans, flags
       type(c_ptr), value                           :: alm_info, geom_info
       real(c_double), intent(out), optional        :: time
       integer(c_long_long), intent(out), optional  :: opcnt
       type(c_ptr), intent(in)                      :: alm(*), map(*)
     end subroutine c_sharp_execute_mpi

  end interface

  interface sharp_execute
     module procedure sharp_execute_d
  end interface

contains
  ! alm info

  ! if ms is not passed, we default to using m=0..lmax.
  subroutine sharp_make_mmajor_real_packed_alm_info(lmax, ms, alm_info)
    use iso_c_binding
    integer(c_int), value, intent(in)    :: lmax
    integer(c_int), intent(in), optional :: ms(:)
    type(sharp_alm_info), intent(out)    :: alm_info
    !--
    integer(c_int), allocatable          :: ms_copy(:)
    integer(c_int)                       :: nm

    if (present(ms)) then
       nm = size(ms)
       allocate(ms_copy(nm))
       ms_copy = ms
       call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, nm, ms_copy, alm_info=alm_info%handle)
       deallocate(ms_copy)
    else
       call c_sharp_make_mmajor_real_packed_alm_info(lmax, 1, lmax + 1, alm_info=alm_info%handle)
    end if
    alm_info%n_local = c_sharp_alm_count(alm_info%handle)
  end subroutine sharp_make_mmajor_real_packed_alm_info

  subroutine sharp_destroy_alm_info(alm_info)
    use iso_c_binding
    type(sharp_alm_info), intent(inout) :: alm_info
    call c_sharp_destroy_alm_info(alm_info%handle)
    alm_info%handle = c_null_ptr
  end subroutine sharp_destroy_alm_info


  ! geom info
  subroutine sharp_make_healpix_geom_info(nside, rings, weight, geom_info)
    integer(c_int), value                :: nside
    integer(c_int), optional             :: rings(:)
    real(c_double), intent(in), optional :: weight(2 * nside)
    type(sharp_geom_info), intent(out)   :: geom_info
    !--
    integer(c_int) :: nrings
    integer(c_int), allocatable :: rings_copy(:)

    if (present(rings)) then
       nrings = size(rings)
       allocate(rings_copy(nrings))
       rings_copy = rings
       call sharp_make_subset_healpix_geom_info(nside, 1, nrings, rings_copy, &
                                                weight, geom_info%handle)
       deallocate(rings_copy)
    else
       call sharp_make_subset_healpix_geom_info(nside, 1, nrings=4 * nside - 1, &
                                                weight=weight, geom_info=geom_info%handle)
    end if
    geom_info%n_local = c_sharp_map_size(geom_info%handle)
  end subroutine sharp_make_healpix_geom_info

  subroutine sharp_destroy_geom_info(geom_info)
    use iso_c_binding
    type(sharp_geom_info), intent(inout) :: geom_info
    call c_sharp_destroy_geom_info(geom_info%handle)
    geom_info%handle = c_null_ptr
  end subroutine sharp_destroy_geom_info


  ! Currently the only mode supported is stacked (not interleaved) maps.
  !
  ! Note that passing the exact dimension of alm/map is necesarry, it
  ! prevents the caller from doing too crazy slicing prior to pass array
  ! in...
  !
  ! Usage:
  !
  ! The alm array must have shape exactly alm(alm_info%n_local, nmaps)
  ! The maps array must have shape exactly map(map_info%n_local, nmaps).
  subroutine sharp_execute_d(type, spin, nmaps, alm, alm_info, map, geom_info, &
                             add, time, opcnt, comm)
    use iso_c_binding
    use mpi
    implicit none
    integer(c_int), value                        :: type, spin, nmaps
    integer(c_int), optional                     :: comm
    logical, value, optional                     :: add  ! should add instead of replace out

    type(sharp_alm_info)                         :: alm_info
    type(sharp_geom_info)                        :: geom_info
    real(c_double), intent(out), optional        :: time
    integer(c_long_long), intent(out), optional  :: opcnt
    real(c_double), target, intent(inout)        :: alm(0:alm_info%n_local - 1, 1:nmaps)
    real(c_double), target, intent(inout)        :: map(0:geom_info%n_local - 1, 1:nmaps)
    !--
    integer(c_int)         :: mod_flags, ntrans, k
    type(c_ptr), target    :: alm_ptr(nmaps)
    type(c_ptr), target    :: map_ptr(nmaps)

    mod_flags = SHARP_DP
    if (present(add) .and. add) then
       mod_flags = or(mod_flags, SHARP_ADD)
    end if

    if (spin == 0) then
       ntrans = nmaps
    else
       ntrans = nmaps / 2
    end if

    ! Set up pointer table to access maps
    do k = 1, nmaps
       alm_ptr(k) = c_loc(alm(0, k))
       map_ptr(k) = c_loc(map(0, k))
    end do

    if (present(comm)) then
      call c_sharp_execute_mpi(comm, type, spin, alm_ptr, map_ptr, &
          geom_info=geom_info%handle, &
          alm_info=alm_info%handle, &
          ntrans=ntrans, &
          flags=mod_flags, &
          time=time, &
          opcnt=opcnt)
    else
      call c_sharp_execute(type, spin, alm_ptr, map_ptr, &
          geom_info=geom_info%handle, &
          alm_info=alm_info%handle, &
          ntrans=ntrans, &
          flags=mod_flags, &
          time=time, &
          opcnt=opcnt)
   end if
  end subroutine sharp_execute_d



end module
