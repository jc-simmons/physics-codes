module funcmod
    
    public 
    real*8, parameter :: pi = 3.1415926535898d0
   


contains


 !factorial function
 real(8) function fact(n)
    implicit none
    !integer*8 :: ans
    real(8) :: ans
    integer :: i, n
    ans = 1d0
    do i=1,n
     ans = ans*i
    end do
    fact = ans
    return

end function fact



! Recursive Hermite polynomial function. Inefficient if many evaluation are required
recursive real(8) function hn(n,x) result(herm)
    implicit none
    integer :: n
    real(8) :: x,b
    b=dble(1)


    if (n==0) then
        herm = dble(1.0)
        return
    else if (n==1) then
        herm = dble(2)*x
        return
    else
        herm = dble(2)*x*hn(n-1,x)-dble(2)*dble(n-1)*hn(n-2,x)
        return
    end if


end function hn



! SOP potential function
real(8) function pot(fexp,fcoef,npot,d,x) result(loc)
    implicit none
    integer :: d, npot,i,fexp(npot,d),j
    real(8)::fcoef(npot),x(d),prod
    loc=dble(0)

    
    do i=1,npot
    prod=1d0
    do j=1,d
    prod=prod*(x(j)**fexp(i,j))
    end do
        loc = loc + prod*fcoef(i)
    end do

    return

    end function pot



SUBROUTINE HERM (N, X, PL)

    IMPLICIT REAL*8 (A-H, O-Z)
    DIMENSION PL(0 : N)

    Y0 = 1.0D0
    Y1 = 2.0D0 * X
    PL(0) = 1.0D0
    PL(1) = 2.0D0 * X
    DO K = 2, N
    C = 2.0D0 * (K - 1.0D0)
    YN = 2.d0 * X * Y1 - C * Y0
    PL(K) = YN
    Y0 = Y1
    Y1 = YN
    ENDDO

    RETURN
END SUBROUTINE HERM




! From http://fortranwiki.org/fortran/show/Matrix+inversion
! Returns the inverse of a matrix calculated by finding the LU
! decomposition.  Depends on LAPACK.
function inv(A) result(Ainv)
    real*8, dimension(:,:), intent(in) :: A
    real*8, dimension(size(A,1),size(A,2)) :: Ainv

    real*8, dimension(size(A,1)) :: work  ! work array for LAPACK
    integer, dimension(size(A,1)) :: ipiv   ! pivot indices
    integer :: n, info

    ! External procedures defined in LAPACK
    external DGETRF
    external DGETRI

    ! Store A in Ainv to prevent it from being overwritten by LAPACK
    Ainv = A
    n = size(A,1)

    ! DGETRF computes an LU factorization of a general M-by-N matrix A
    ! using partial pivoting with row interchanges.
    call DGETRF(n, n, Ainv, n, ipiv, info)

    if (info /= 0) then
    stop 'Matrix is numerically singular!'
    end if

    ! DGETRI computes the inverse of a matrix using the LU factorization
    ! computed by DGETRF.
    call DGETRI(n, Ainv, n, ipiv, work, n, info)

    if (info /= 0) then
    stop 'Matrix inversion failed!'
    end if
end function inv



subroutine lu_nopivot(n,A,L,U)

    implicit none
    integer :: i,j,k,n
    real*8 :: A(n,n),U(n,n),L(n,n)
    
    L=0d0
    
    
    
    do i=1,n
        L(i,i)=1d0
    end do
    
    do k=1,n
        
        if (abs(A(k,k))<0.00000001d0) then
            print*, 'cannot perform LU factorization'
            print*, k , A(k,k)
            stop
        end if
        L(k+1:n,k)=A(k+1:n,k)/A(k,k)
        
        do j=k+1,n
            A(j,:)=A(j,:)-L(j,k)*A(k,:)
        end do
        
    end do
    
    U=A;    
    
    end subroutine lu_nopivot  


end module funcmod