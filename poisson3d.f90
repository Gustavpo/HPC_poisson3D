program poisson3d

  ! Import different methods
  use precision
  use poisson_methods_m
  use init 

  ! Define options for user
  ! Please note that while we assign on definition, save is implicit.
  ! But please write save to show intent.
  integer(lp), save :: N = 1000
  real(dp), save :: T0 = 0._dp
  integer(lp), save :: itermax = 1000
  real(dp), save :: tolerance = 1.e-4_dp
  integer(lp), save :: output = 4 ! 3 for binary, 4 for VTK
  integer(lp), save :: algorithm = 1 ! 1 for Jacobi, 2 for Gauss-Seidel
  ! don't ask why, this is for historical reasons and to be compatible with the C-version
  ! of the assignment
  real(dp) :: delta
  real(dp) :: count_rate, wtime
  integer(lp) :: count1, count2
  !character(len=25) :: string

  real(dp), allocatable, dimension(:,:,:) :: u, u_old, f
  integer :: iostat

  ! Read in namelist.
  ! Call program namelist
  ! Create a file called: input.dat
  ! The content may be:
  ! &INPUT
  !   N = 100
  !   itermax = 3000
  !   T0 = 10.
  !   tolerance = 1e-5
  !   output = 0
  !   algorithm = 1
  ! /
  ! If one of the keywords are not in the
  ! input.dat namelist, then the default value
  ! will be used.
  call read_namelist()
  delta = 2.0/(N-1)

  allocate(u(N,N,N), stat=iostat)
  call check_iostat(iostat, "Could not allocate 'u' matrix!")
  allocate(f(N,N,N), stat=iostat)
  call check_iostat(iostat, "Could not allocate 'f' matrix!")
  

  call system_clock(count1,count_rate)  
  
  call initialize(u,f,N,T0,delta)
  allocate(u_old(N,N,N), stat=iostat)
  call check_iostat(iostat, "Could not allocate 'u_old' matrix!")
  
  if (algorithm.eq.1) call jacobi(u,u_old,f,N,delta,itermax,tolerance)
  if (algorithm.eq.2) call gauss_seidel(u,u_old,f,N,delta,itermax,tolerance)
  call system_clock(count2,count_rate)
  wtime = (count2-count1)*1.0/count_rate
  
  print*,'wall time: ',wtime, ' seconds'
 

  !write(string,'(A,I3,A)') 'output/timelist',N,'_2.dat'
  !open(12,file=string)
  !write(12,'(A)') 'Time    N'
  !write(12,'(E12.4,I3.2)') wtime, N
  !close(12)
  
  deallocate(f)
  if ( allocated(u_old) ) deallocate(u_old)

  ! Keep u until we have written in out
  select case ( output )
  case ( 0 ) ! pass, valid but not used value
    ! do nothing
  case ( 3 ) ! write to binary file
    call write_binary(u)
  case ( 4 ) ! write to binary file
    call write_vtk(u)
  case default

    write(*,'(a,i0)') 'Unknown output type requested: ', output
    stop

  end select
    
  deallocate(u)

contains

  subroutine read_namelist()
    integer :: unit, io_err
    
    namelist /INPUT/ N, itermax, T0, tolerance, output, algorithm
    
    ! open and read file
    unit = 128273598
    open(unit, FILE="input.dat", action='read', iostat=io_err)
    call check_iostat(io_err, &
        "Could not open file 'input.dat', perhaps it does not exist?")
    read(unit, nml=INPUT, iostat=io_err)
    call check_iostat(io_err, &
        "Error on reading name-list, please correct file content")
    close(unit)

  end subroutine read_namelist

  subroutine check_iostat(iostat, msg)
    integer, intent(in) :: iostat
    character(len=*), intent(in) :: msg

    if ( iostat == 0 ) return

    write(*,'(a,i0,/,tr3,a)') 'ERROR = ', iostat, msg
    stop

  end subroutine check_iostat

  subroutine write_binary(u)
    ! Array to write-out
    real(dp), intent(in) :: u(:,:,:)

    ! Local variables
    character(len=*), parameter :: filename = 'poisson_u.bin'
    integer :: N
    integer :: iostat

    integer, parameter :: unit = 11249

    ! Get size of the array
    N = size(u, 1)

    ! replace == overwrite
    open(unit,file=filename,form='unformatted',access='stream',status='replace', iostat=iostat)
    call check_iostat(iostat, "Could not open file '"//filename//"'!")
    write(unit, iostat=iostat) N
    call check_iostat(iostat, "Could not write variable N")
    write(unit, iostat=iostat) u
    call check_iostat(iostat, "Could not write variable u")
    close(unit)

  end subroutine write_binary

  subroutine write_vtk(u)

    ! This is C-interoperability using bindings
    use, intrinsic :: iso_c_binding

    interface
      subroutine write_vtk_c(n, u) BIND(C, name="write_vtk")
        import c_ptr, c_int
        implicit none
        integer(c_int), value :: n
        type(c_ptr), value :: u
      end subroutine
    end interface

    ! Array to write-out
    real(dp), intent(in), target :: u(:,:,:)

    ! Local variables
    integer(c_int) :: N

    ! Get size of the array
    N = size(u, 1)

    call write_vtk_c(n, c_loc(u(1,1,1)))

  end subroutine write_vtk

end program poisson3d
