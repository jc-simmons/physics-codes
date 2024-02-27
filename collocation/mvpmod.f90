module mvpmod

    implicit none
    public

contains


subroutine mvp(u,init,np,nb,nx,d,terms,Bz,Binv,potp,sortlens,sortlims,usort,dsort,ssort,fcoef)
 
    integer :: i,j,d,k,l,tind(2,2),np,nb,terms,nx,p
    integer :: deriv(d),idim(d,d)
    integer :: sortlens(4,2,2,d),sortlims(4,2,d,np),usort(2,d,np,nx),dsort(2,d,np,nx),ssort(2,d,np,nx)
    real*8 :: u(nb),v(np),Bz(2,2,d,terms,nx,nx),Binv(d,2,nx,nx)
    real*8 :: res(np+1),uu(np)
    real*8 :: init(nb),fcoef(d)
    real*8 :: potp(np)
    !integer :: time1,time2,rate
    !call system_clock(count_rate=rate)
    
    
 
    idim=0   
    do i=1,d
     idim(1,i)=d+1-i
    end do
  
    do i=2,d
      idim(i,1)=idim(1,i)
      idim(i,2:i)=idim(1,1:i-1)
       
 
      if (i<d) then
 
        idim(i,i+1:d)=idim(1,i+1:d)
      
      end if
 
    end do
    
  
    idim=1
    do i=1,d
     idim(i,i)=2
    end do
 
    uu=0d0
    v=0d0
    u=0d0
    deriv=1
    deriv(1)=2
 
    tind=1
    tind(1,1)=2
    tind(2,2)=2    
 
    res=0d0
 !  indb=0
    res(2:nb+1)=init  
 
    
 
 
 !----multiply by BU-----
    do i=1,d
     
         v=0d0
 
 
     !$omp parallel do
     do j=1,sortlens(4,1,2,i)
       do k=1,sortlims(4,2,i,j)
         do l=1,sortlims(4,1,i,j)
           v(ssort(1,i,j,k))=v(ssort(1,i,j,k))+Bz(1,2,d-i+1,1,k,l)*res(ssort(1,i,j,l)+1)
         end do
       end do
     end do
 
     
 
     res(2:)=v
 
 
 
    end do
 
 
 !----multiply by BL-----
 
    do i=1,d
     
         v=0d0
 
     !$omp parallel do
     do j=1,sortlens(1,1,2,i)
       do k=1,sortlims(1,2,i,j)
         do l=1,sortlims(1,1,i,j)
           v(usort(2,i,j,k))=v(usort(2,i,j,k))+Bz(1,1,d-i+1,1,k,l)*res(usort(1,i,j,l)+1)
         end do
       end do
     end do
 
     
 
     res(2:)=v
 
 
    end do
    
    
   
 
 !-----------multiply by V-------------
 
    !$omp parallel do
    do i=1,np
       uu(i)=uu(i)+res(i+1)*potp(i)
    end do    
 
 !-----------multiply by B''-------------
   
 
  
   
    
   do l=1,terms
   
 
     res=0d0
     res(2:nb+1)=init
     
     
    do i=1,d
   
    v=0d0  
  
    !$omp parallel do
    do j=1,sortlens(4,1,2,i)
      do k=1,sortlims(4,2,i,j)
        do p=1,sortlims(4,1,i,j)
          v(ssort(1,i,j,k))=v(ssort(1,i,j,k))+Bz(2,2,d-i+1,terms-l+1,k,p)*res(ssort(1,i,j,p)+1)
        end do
      end do
    end do
 
    
    res(2:)=v
 
 
    end do
      
      
      do i=1,d
       
       v=0d0  
    
    !$omp parallel do
    do j=1,sortlens(1,1,2,i)
      do k=1,sortlims(1,2,i,j)
        do p=1,sortlims(1,1,i,j)
          v(usort(2,i,j,k))=v(usort(2,i,j,k))+Bz(2,1,d-i+1,terms-l+1,k,p)*res(usort(1,i,j,p)+1)
        end do
      end do
    end do
     
   res(2:)=v
 
   end do
      
   uu=uu+res(2:)*fcoef(d-l+1)
 
   end do  
     
   
   
  
 !----------- Multiply by inv(BL) ---------  
     res=0d0
     res(2:)=uu 
     
   do i=1,d
     
    v=0d0  
    
    !$omp parallel do
    do j=1,sortlens(3,1,2,i)
      do k=1,sortlims(3,2,i,j)
        do p=1,sortlims(3,1,i,j)
          v(ssort(2,i,j,k))=v(ssort(2,i,j,k))+ Binv(d-i+1,1,k,p)*res(ssort(2,i,j,p)+1)
        end do
      end do
    end do
     
   res(2:)=v
    
 
   end do
   
   
    
    uu=res(2:)
    res=0d0
    res(2:)=uu 
    
    !----------- Multiply by inv(BU) ---------  
     
   do i=1,d
     
    v=0d0  
    !$omp parallel do
    do j=1,sortlens(2,1,2,i)
      do k=1,sortlims(2,2,i,j)
        do p=1,sortlims(2,1,i,j)
          v(dsort(2,i,j,k))=v(dsort(2,i,j,k))+ Binv(d-i+1,2,k,p)*res(dsort(1,i,j,p)+1)
        end do
      end do
    end do
     
   res(2:)=v
    
 
   end do
    
    
 
    u=v
 
end subroutine mvp

end module mvpmod