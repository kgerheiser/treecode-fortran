module fconfig
  use iso_fortran_env, only: real64, real32, int64, int32, iostat_end
  implicit none

  integer, parameter :: MAX_STR_LEN = 256
  private :: value_from_key_str, value_from_key_r4, value_from_key_r8, &
       value_from_key_i4, value_from_key_i8, value_from_key_logical

  type :: string
     character(:), allocatable :: str
  end type string

  type :: config
     character :: separator = ":"
     integer :: num_entries
     type(string), allocatable :: keys(:), values(:)
   contains
     procedure :: find_str_value_with_key
     procedure, private :: find_key_index
     procedure :: read_file
     procedure, private :: value_from_key_str
     procedure, private :: value_from_key_r4
     procedure, private :: value_from_key_r8
     procedure, private :: value_from_key_i4
     procedure, private :: value_from_key_i8
     procedure, private :: value_from_key_logical
     generic :: value_from_key => value_from_key_str, value_from_key_r4, &
          value_from_key_r8, value_from_key_i4, value_from_key_i8, value_from_key_logical
  end type config

contains

  subroutine value_from_key_str(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    character(*), optional, intent(in) :: default_value
    character(:), allocatable, intent(out) :: val
    character(:), allocatable :: str_val

    str_val = conf%find_str_value_with_key(key)

    if (.not. str_val == "") then
       val = str_val    
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_str

  subroutine value_from_key_r4(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    real(real32), intent(out) :: val
    character(:), allocatable :: str_val
    real(real32), optional, intent(in) :: default_value
    integer :: iostat

    str_val = conf%find_str_value_with_key(key)
    
    if (.not. str_val == "") then
       read(str_val, *, iostat=iostat) val
       if (iostat /= 0) call string_conversion_error(str_val, "r4", key)      
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_r4

  subroutine value_from_key_r8(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    real(real64), intent(out) :: val
    real(real64), optional, intent(in) :: default_value
    character(:), allocatable :: str_val
    integer :: iostat

    str_val = conf%find_str_value_with_key(key)

    if (.not. str_val == "") then
       read(str_val, *, iostat=iostat) val
       if (iostat /= 0) call string_conversion_error(str_val, "r8", key)      
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_r8

  subroutine value_from_key_i4(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    integer(int32), intent(out) :: val
    integer(int32), optional, intent(in) :: default_value
    character(:), allocatable :: str_val
    integer :: iostat

    str_val = conf%find_str_value_with_key(key)

    if (.not. str_val == "") then
       read(str_val, *, iostat=iostat) val
       if (iostat /= 0) call string_conversion_error(str_val, "i4", key)      
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_i4

  subroutine value_from_key_i8(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    integer(int64), intent(out) :: val
    integer(int64), optional, intent(in) :: default_value
    character(:), allocatable :: str_val
    integer :: iostat

    str_val = conf%find_str_value_with_key(key)

    if (.not. str_val == "") then
       read(str_val, *, iostat=iostat) val
       if (iostat /= 0) call string_conversion_error(str_val, "i8", key)      
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_i8

  subroutine value_from_key_logical(conf, key, val, default_value)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    logical, intent(out) :: val
    logical, optional, intent(in) :: default_value
    character(:), allocatable :: str_val
    integer :: iostat

    str_val = conf%find_str_value_with_key(key)

    if (.not. str_val == "") then
       read(str_val, *, iostat=iostat) val
       if (iostat /= 0) call string_conversion_error(str_val, "logical", key)      
    else
       if (present(default_value)) then
          val = default_value
       else
          call value_not_found_error(key)
       end if
    end if

  end subroutine value_from_key_logical

  subroutine value_not_found_error(key)
    character(*), intent(in) :: key
    print *, "Value not defined and no default value present for key ", '"', key, '"'
    stop
  end subroutine value_not_found_error

  subroutine string_conversion_error(str_val, type, key)
    character(*), intent(in) :: str_val, type, key
    print *, "Error converting string ", '"', str_val, '" to ', type, " from key ", '"', key, '"'
    stop
  end subroutine string_conversion_error

  integer function find_key_index(conf, key) result(index)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key
    integer :: i, n

    index = 0
    n = conf%num_entries

    do i = 1, n
       if (key == conf%keys(i)%str) index = i
    end do

  end function find_key_index

  function find_str_value_with_key(conf, key) result(str_val)
    class(config), intent(in) :: conf
    character(*), intent(in) :: key

    character(:), allocatable :: str_val
    integer n, key_index

    n = conf%num_entries

    key_index = conf%find_key_index(adjustl(trim(key)))
    str_val = ""

    if (key_index /= 0) then
       str_val = conf%values(key_index)%str
    end if

  end function find_str_value_with_key


  subroutine read_file(conf, file) 
    class(config), intent(inout) :: conf
    character(*), intent(in) :: file

    character(len=MAX_STR_LEN) :: iomsg, line
    integer :: file_unit, iostat, i, num_entries, sub_index

    open(file = file, newunit = file_unit, form = 'formatted', status = 'old', action = 'read',&
         iostat = iostat, iomsg = iomsg)

    if (iostat /= 0) then
       print *, trim(iomsg)
       stop 
    end if

    num_entries = 0

    do 
       read(file_unit, '(a)', iostat = iostat) line
       if (iostat == iostat_end) exit
       if (.not. accept_line(line)) cycle
       num_entries = num_entries + 1
    end do

    conf%num_entries = num_entries
    allocate(conf%keys(num_entries), conf%values(num_entries))
    rewind(file_unit)
    i = 0

    do
       read(file_unit, '(a)', iostat = iostat) line
       if (iostat == iostat_end) exit
       if (.not. accept_line(line)) cycle

       i = i + 1
       sub_index = index(line, ":")

       if (sub_index == 0) then
          print *, "Line missing ':', ", adjustl(trim(line))
          print *, "canceling parse..."
          stop
       else
          conf%keys(i)%str = line(1:sub_index-1)
          conf%values(i)%str = line(sub_index+1:)
       end if
    end do

    do i = 1, size(conf%keys)
       conf%keys(i)%str = adjustl(conf%keys(i)%str)
       conf%keys(i)%str = trim(conf%keys(i)%str)
       
       conf%values(i)%str = adjustl(conf%values(i)%str)
       conf%values(i)%str = trim(conf%values(i)%str)
    end do

  contains

    logical function accept_line(line) result(accept)
      character(*), intent(in) :: line
      character(:), allocatable :: trimmed_line

      trimmed_line = adjustl(trim(line))

      if(len(trimmed_line) == 0) then
         accept = .false. 
      else if(trimmed_line(1:1) == "!" .or. trimmed_line(1:1) == "#") then
         accept = .false.
      else 
         accept = .true.
      end if

    end function accept_line

  end subroutine read_file

end module fconfig

