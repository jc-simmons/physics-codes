module indmod
    implicit none
    public indgen
    private

contains

    subroutine indgen(d,b,m,ind,imod,aflag)

        integer :: d,b,m,imod(d),aflag
        integer, allocatable, intent(inout) :: ind(:,:)

        if (d==3) then
            call indgen3D(d,b,m,ind,imod,aflag)
        else if (d==9) then
            call indgen9D(d,b,m,ind,imod)
        end if

    end subroutine indgen



    subroutine indgen3D(d,b,m,ind,imod,aflag)

        implicit none
        
        integer :: d,b,m,aflag
        integer :: i1,i2,i3,imod(d)
        integer, allocatable, intent(inout) :: ind(:,:)
    
    
        m=0
        
        do i1=0,b+imod(1)
        do i2=0,b-i1+imod(2)
        do i3=0,b-i1-i2+imod(3)
        
        
        m=m+1
        
        end do
        end do
        end do
    
        if (aflag ==1) then
            allocate(ind(d,m))
        end if

        
        ind=0
        m=0
    
        do i1=0,b+imod(1)
        do i2=0,b-i1+imod(2)
        do i3=0,b-i1-i2+imod(3)
    
        
        
        m=m+1
        
        ind(1,m)=i1+1
        ind(2,m)=i2+1
        ind(3,m)=i3+1
    
        
        end do
        end do
        end do
    

        
    end subroutine indgen3D


    subroutine indgen9D(d,b,m,ind,imod)

        implicit none
        
        integer :: d,b,m
        integer :: i1,i2,i3,i4,i5,i6,i7,i8,i9,imod(d)
        integer, allocatable, intent(inout) :: ind(:,:)


        m=0
        
        do i1=0,b+imod(1)
        do i2=0,b-i1+imod(2)
        do i3=0,b-i1-i2+imod(3)
        do i4=0,b-i1-i2-i3+imod(4)
        do i5=0,b-i1-i2-i3-i4+imod(5)
        do i6=0,b-i1-i2-i3-i4-i5+imod(6)
        do i7=0,b-i1-i2-i3-i4-i5-i6+imod(7)
        do i8=0,b-i1-i2-i3-i4-i5-i6-i7+imod(8)
        do i9=0,b-i1-i2-i3-i4-i5-i6-i7-i8+imod(9)
        
        
        m=m+1
        
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        
        allocate(ind(d,m))
        
        m=0
        ind=0
        
        do i1=0,b+imod(1)
        do i2=0,b-i1+imod(2)
        do i3=0,b-i1-i2+imod(3)
        do i4=0,b-i1-i2-i3+imod(4)
        do i5=0,b-i1-i2-i3-i4+imod(5)
        do i6=0,b-i1-i2-i3-i4-i5+imod(6)
        do i7=0,b-i1-i2-i3-i4-i5-i6+imod(7)
        do i8=0,b-i1-i2-i3-i4-i5-i6-i7+imod(8)
        do i9=0,b-i1-i2-i3-i4-i5-i6-i7-i8+imod(9)
        
        
        m=m+1
        
        ind(1,m)=i1+1
        ind(2,m)=i2+1
        ind(3,m)=i3+1
        ind(4,m)=i4+1
        ind(5,m)=i5+1
        ind(6,m)=i6+1
        ind(7,m)=i7+1
        ind(8,m)=i8+1
        ind(9,m)=i9+1
        
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        end do
        

        
    end subroutine indgen9D


end module indmod