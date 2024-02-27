module buildmatmod

    implicit none
    public 



contains



! Define 1D postion B and KEO B'' matrices
subroutine buildBmats(d, n, pmax, Bmats, kterms, grid)
    use funcmod
    integer, intent(in) :: d, n(d), pmax
    integer, intent(out) :: kterms
    integer :: i,j,k,l
    real*8 :: beta, grid(pmax, d)
    real*8, allocatable, intent(inout) :: Bmats(:,:,:,:,:)

    kterms=d  ! how many terms depends on choice of KEO. Here simple normal coordinate example
    beta = 1d0  ! harmonic oscillator parameter, here dimensionless coordinates
    allocate(Bmats(2,kterms,pmax,pmax+2,d))
    Bmats=0d0



    ! B matrices for nth dimension are stored in the Bmats(1,1,:,:,n) position
    do l=1,d 
    do j=1,pmax
    do k=1,pmax+2
        Bmats(1,1,j,k,l)= 1d0/sqrt(2d0**(k-1)*dble(fact(k-1)))*(1d0/pi)**0.25d0*exp(-(1d0*grid(j,l))**2/2d0)*hn(k-1,1d0*grid(j,l))
    end do
    end do
    end do

    ! print*, Bmats(1,1,:,:,1)
    !print*, ' '


    ! For B'', if now there are also 'm' terms in the keo, these are stored in the Bmats(2,m,:,:,n)
    !
    ! using simple normal coordinate KEO, the ith term acts only on the ith dimension,
    ! for B'', there are d x d one dimensional matrices, d of which are affected by the keo terms:


    ! make all d x d matrices
    do i=1,kterms ! here kterms = d
    do j=1,d

    Bmats(2,i,:,1:n(j),j)=Bmats(1,1,:,1:n(j),j)

    end do
    end do


    ! replace entries corresponding to derivative terms
    do l=1,d
    do j=1,pmax
    do k=1,pmax
    if (k>2) then     
        Bmats(2,l,j,k,l)=-beta**2/2d0*(-(2d0*(k-1)+1d0)*Bmats(1,1,j,k,l)+sqrt(1d0*(k-1)+1d0)*sqrt(1d0*(k-1)+2d0)*Bmats(1,1,j,k+2,l)+&
                sqrt(1d0*(k-1))*sqrt(1d0*(k-1)-1d0)*Bmats(1,1,j,k-2,l))     
    else

        Bmats(2,l,j,k,l)=-beta**2/2d0*(-(2d0*(k-1)+1d0)*Bmats(1,1,j,k,l)+sqrt(1d0*(k-1)+1d0)*sqrt(1d0*(k-1)+2d0)*Bmats(1,1,j,k+2,l))       
    end if
    end do
    end do
    end do


end subroutine buildBmats



subroutine buildLUmats(d, n, pmax, kterms, Bz, Bmats)
    use funcmod
    integer, intent(in) :: d, n(d), pmax, kterms
    real*8, intent(in) :: Bmats(2,kterms,pmax,pmax+2,d)
    integer :: i,j 
    real*8, allocatable, intent(inout) :: Bz(:,:,:,:,:,:)
    real*8, allocatable :: BL(:,:), BU(:,:), A(:,:)

    allocate(Bz(2,2,d,kterms,pmax,pmax))

    Bz = 0d0

      ! LU factorization for B

    do i=1,d

    
        allocate(BL(n(i),n(i)),BU(n(i),n(i)),A(n(i),n(i)))
        BL=0d0
        BU=0d0
        
        A=Bmats(1,1,:,1:n(i),i)

        call lu_nopivot(n(i),A,BL(:,:),BU(:,:))
        
        Bz(1,1,i,1,1:n(i),1:n(i))=BL(:,:)
        Bz(1,2,i,1,1:n(i),1:n(i))=BU(:,:)
    
        deallocate(A,BL,BU)
    
    end do 

    
    ! LU factorization for B''

    do j=1,kterms
    
    do i=1,d
    
    
        allocate(BL(n(i),n(i)),BU(n(i),n(i)),A(n(i),n(i)))
        BL=0d0
        BU=0d0
        
        A=Bmats(2,j,:,1:n(i),i)


        call lu_nopivot(n(i),A,BL(:,:),BU(:,:))
        
        Bz(2,1,i,j,1:n(i),1:n(i))=BL(:,:)
        Bz(2,2,i,j,1:n(i),1:n(i))=BU(:,:)
    
        deallocate(A,BL,BU)
    
    end do 

    end do


end subroutine buildLUmats


subroutine buildINVmats(d, n, pmax, kterms, Bz, Binv)

    integer, intent(in) :: d, n(d), pmax, kterms
    real*8, intent(in) :: Bz(2,2,d,kterms,pmax,pmax)
    integer :: i, info, lworkl
    integer, allocatable :: ipiv(:), work(:)
    real*8, allocatable ::  A(:,:)
    real*8, allocatable, intent(inout) :: Binv(:,:,:,:)
    
    allocate(Binv(d,2,pmax,pmax))

    Binv=0d0
 
    do i=1,d
    
        allocate(A(n(i),n(i)),ipiv(n(i)))
        lworkl=n(i)*2
        allocate(work(lworkl))
        
        A=0d0
        ipiv=0
        info=0
        
        ipiv=0
        info=0  
    
        A=Bz(1,1,i,1,1:n(i),1:n(i))
        
        call dgetrf(n(i),n(i),A,n(i),ipiv,info)
        call dgetri(n(i),A,n(i),ipiv,work,lworkl,info)
        
        Binv(i,1,1:n(i),1:n(i))=A
        
        
        
        ipiv=0
        info=0
        A=0d0
        
        
        A=Bz(1,2,i,1,1:n(i),1:n(i))
        
        call dgetrf(n(i),n(i),A,n(i),ipiv,info)
        call dgetri(n(i),A,n(i),ipiv,work,lworkl,info)
        
        Binv(i,2,1:n(i),1:n(i))=A
        
    
        
        deallocate(A,ipiv,work)
    
    
    end do


end subroutine buildINVmats

end module buildmatmod
