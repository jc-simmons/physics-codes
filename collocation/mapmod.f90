module mapmod

    implicit none
    public

contains

subroutine indmap(d, nb, np, bmax, pmax, binds, pinds, usort, dsort, ssort, sortlens, sortlims)
    use sortmod
    use indmod

    integer, intent(in) :: d, np, nb, bmax, pmax, binds(d,nb), pinds(d,np)
    integer :: i, l, imod(d), perml(np), indb(d,np)
    integer, allocatable, intent(inout) :: usort(:,:,:,:), dsort(:,:,:,:), ssort(:,:,:,:)
    integer, allocatable, intent(inout) :: sortlens(:,:,:,:), sortlims(:,:,:,:)
    integer, allocatable :: indp(:,:)


    allocate(usort(2,d,np,pmax))
    allocate(ssort(2,d,np,pmax))
    allocate(dsort(2,d,np,pmax))
    allocate(sortlens(4,2,2,d))
    allocate(sortlims(4,2,d,np))    
    allocate(indp(d,np))

    sortlims=0
    usort=0
    dsort=0
    ssort=0
    sortlens=0
    sortlens(1,1,:,1)=nb
    perml=0
    indp=0
    indp(:,1:nb)=binds



     ! basis -> points 
    do l=1,d
    
        indb=indp
     
        indb(d,:)=indp(d-l+1,:)
        indb(d-l+1,:)=indp(d,:)
        call hpsort2(indb,d,sortlens(1,1,1,l),np,perml)
        call sortinds(indb,usort(1,l,:,:),sortlims(1,1,l,:),d,np,pmax,sortlens(1,1,:,l),perml)
       
      
        !This generates the intermediate indices after each 1D sum
        ! sortlens(1,2,1,l) is the number of index combinations in this list
        imod(d-l+1)=pmax-bmax
        !deallocate(indp)
        call indgen(d,bmax-1,sortlens(1,2,1,l),indp,imod,0)
        sortlens(1,1,1,2:d)=sortlens(1,2,1,1:d-1)
        
        
        
        indb(1:d,1:sortlens(1,2,1,l))=indp
        indb(d,1:sortlens(1,2,1,l))=indp(d-l+1,:)
        indb(d-l+1,1:sortlens(1,2,1,l))=indp(d,:)
        call hpsort2(indb,d,sortlens(1,2,1,l),np,perml)
        call sortinds(indb,usort(2,l,:,:),sortlims(1,2,l,:),d,np,pmax,sortlens(1,2,:,l),perml)
  
  
      end do
    
    !deallocate(indp)
    !allocate(indp(d,np))
     
    if (sortlens(1,2,1,d).ne.np) then
      print*, '1sort failed'
      stop
    end if
     
    indp=pinds
    sortlens(2,1,1,1)=np
    
    imod=0
    
    
    ! points -> basis
    do l=1,d
    
    
        indb=indp
        
        indb(d,:)=indp(d-l+1,:)
        indb(d-l+1,:)=indp(d,:)
        call hpsort2(indb,d,sortlens(2,1,1,l),np,perml)
        call sortinds(indb,dsort(1,l,:,:),sortlims(2,1,l,:),d,np,pmax,sortlens(2,1,:,l),perml)
      
      
        !This generates the intermediate indices after each 1D sum
        ! sortlens(2,2,1,l) is the number of index combinations in this list
        imod(d-l+1)=bmax-pmax  
        !deallocate(indp)
        call indgen(d,pmax-1,sortlens(2,2,1,l),indp,imod,0)
        sortlens(2,1,1,2:d)=sortlens(2,2,1,1:d-1)
        
        do i=1,d/2
          
        indb(i,:)=indp(d-i+1,:)
        indb(d-i+1,:)=indp(i,:)
        
        end do
        
        
        indb(1:d,1:sortlens(2,2,1,l))=indp
        indb(d,:)=indp(d-l+1,:)
        indb(d-l+1,:)=indp(d,:)
        call hpsort2(indb,d,sortlens(2,2,1,l),np,perml)
        call sortinds(indb,dsort(2,l,:,:),sortlims(2,2,l,:),d,np,pmax,sortlens(2,2,:,l),perml)
    
    
    end do 


     
    if (sortlens(2,2,1,d).ne.nb) then
      print*, '2sort failed'
      stop
    end if 


    !deallocate(indp)
    !allocate(indp(d,np))
  
  
    sortlens(3,1,1,1)=np
    indp=pinds
    sortlens(3,1,1,:)=np
    sortlens(3,2,1,:)=np
    
    imod=0
    
    
    
    ! points -> points 
    do l=1,d
    
       
     indb=indp
     indb(d,:)=indp(d-l+1,:)
     indb(d-l+1,:)=indp(d,:)
     call hpsort2(indb,d,sortlens(3,1,1,l),np,perml)
     
     
     call sortinds(indb,ssort(1,l,:,:),sortlims(3,1,l,:),d,np,pmax,sortlens(3,1,:,l),perml)
          
     call sortinds(indb,ssort(2,l,:,:),sortlims(3,2,l,:),d,np,pmax,sortlens(3,2,:,l),perml)
    
    
    end do 


    if (sortlens(3,2,1,d).ne.np) then
        print*, '3sort failed'
        stop
      end if 
    
    
  
    sortlens(4,1,1,1)=nb
    indp(:,1:nb)=binds
    sortlens(4,1,1,:)=nb
    sortlens(4,2,1,:)=nb
    
    imod=0
     
     
    ! basis -> basis  
    do l=1,d
     
     
     indb=indp
     indb(d,:)=indp(d-l+1,:)
     indb(d-l+1,:)=indp(d,:)
     call hpsort2(indb,d,sortlens(4,1,1,l),np,perml)
     
     
     call sortinds(indb,ssort(1,l,:,:),sortlims(4,1,l,:),d,np,pmax,sortlens(4,1,:,l),perml)
     sortlims(4,2,l,:)=sortlims(4,1,l,:)
   
          
    end do 
  
  
  
    if (sortlens(4,2,1,d).ne.nb) then
     print*, '4rectangular failed'
     stop
    end if 


end subroutine indmap

end module mapmod