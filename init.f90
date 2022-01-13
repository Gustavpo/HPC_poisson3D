module init
use precision
implicit none

contains

subroutine initialize(x,y,N,T0,delta)
  integer(lp) :: i,j,k
  integer(lp), intent(in) :: N
  real(dp), intent(in) :: delta,T0
  real(dp), dimension(N,N,N), intent(inout) :: x,y

  do i = 1,N
    do j = 1,N
      do k = 1,N
        x(i,j,k) = T0
        y(i,j,k) = 0.0
        x(i,1,k) = 0.0
        x(i,N,k) = 20.0
        x(1,j,k) = 20.0
        x(N,j,k) = 20.0 
      enddo
      x(i,j,1) = 20.0
      x(i,j,N) = 20.0
    enddo
  enddo
  
  do i = 1,int(floor(5.0/8.0/delta))
    do j = 1,int(floor(1.0/2.0/delta))
      do k = int(ceiling(1.0/3.0/delta)),int(floor(1.0/delta))
        y(i,j,k) = 200.0
      enddo
    enddo
  enddo
  
end subroutine initialize

end module init
