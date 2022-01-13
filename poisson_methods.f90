module poisson_methods_m

  ! Please do not put any computational variable in this module.
  ! Pass all variables as arguments to the subroutine.
  use precision, only: dp, lp ! use dp as the data-type (and for defining constants!)

  implicit none

  private
  public :: jacobi
  public :: gauss_seidel

contains

  subroutine jacobi(x,x_old,y,N,delta,itermax,tolerance)
    integer(lp), intent(in) :: N, itermax
    real(dp), intent(in) :: delta, tolerance
    integer(lp) :: nits = 0
    integer(lp) :: i,j,k
    real(dp) :: p_error = 1e10
    real(dp), dimension(N,N,N), intent(inout) :: x,x_old
    real(dp), dimension(N,N,N), intent(in) :: y
    x_old = x
    
    do while (p_error>tolerance .AND. nits < itermax)  
        do k = 2,N-1
            do j = 2,N-1
                do i = 2,N-1
                    x(i,j,k) = 1.0/6.0*(x_old(i-1,j,k) + x_old(i+1,j,k) + x_old(i,j-1,k) + &
                    x_old(i,j+1,k) + x_old(i,j,k-1) + x_old(i,j,k+1) + delta**2.0*y(i,j,k))
                enddo
            enddo
        enddo
        p_error = sum(abs(x - x_old))/(N**3.0)
        x_old = x
        nits = nits + 1      
    enddo
    print*,nits    
  end subroutine jacobi

  subroutine gauss_seidel(x,x_old,y,N,delta,itermax,tolerance)
    integer(lp), intent(in) :: N, itermax
    real(dp), intent(in) :: delta, tolerance
    integer(lp) :: nits = 0
    integer(lp) :: i,j,k
    real(dp) :: p_error = 1e10
    real(dp), dimension(N,N,N), intent(inout) :: x,x_old
    real(dp), dimension(N,N,N), intent(in) :: y
    x_old = x
    
    do while (p_error>tolerance .AND. nits < itermax)  
        do k = 2,N-1
            do j = 2,N-1
                do i = 2,N-1
                    x(i,j,k) = 1.0/6.0*(x(i-1,j,k) + x_old(i+1,j,k) + x(i,j-1,k) + &
                    x_old(i,j+1,k) + x(i,j,k-1) + x_old(i,j,k+1) + delta**2.0*y(i,j,k))
                enddo
            enddo
        enddo
        p_error = sum(abs(x - x_old))/(N**3.0)
        x_old = x
        nits = nits + 1      
    enddo
    print*,nits    
    
  end subroutine gauss_seidel

end module poisson_methods_m
