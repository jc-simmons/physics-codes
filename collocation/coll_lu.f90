program main
  use sortmod
  use indmod
  use funcmod
  use buildmatmod
  use mapmod
  use mvpmod
  implicit none 
  !$use omp_lib
  integer :: i,j,d,info,nb,np,i2,kterms,pmax,bmax
  integer :: nx,nev,ncv,iparam(11),ido,ipntr(14),maxiter,lworkl,npot
  integer, allocatable :: iselect(:),fexp(:,:),pinds(:,:),binds(:,:),imod(:),n(:)
  integer, allocatable :: usort(:,:,:,:),dsort(:,:,:,:),sortlens(:,:,:,:),sortlims(:,:,:,:),ssort(:,:,:,:)
  integer :: time1,time2,rate,nthreads
  real*8, allocatable :: Bmats(:,:,:,:,:),workd(:),workl(:),resid(:)
  real*8, allocatable :: Bz(:,:,:,:,:,:),dr(:),di(:),zeig(:,:),workev(:)
  real*8 , allocatable :: nmvp(:),v(:,:),dr2(:),potp(:),point(:),fcoef(:),Binv(:,:,:,:),grid(:,:)
  real*8 :: ran,tol
  real*8 :: sigmar,sigmai
  logical :: rvec, add
  character(len=1) :: bmat, howmny
  character(len=2) :: which
  character(len=64):: out_file, pot_file,pstring,bstring
  
  
  ! # processes. Only parallel section is MVP evaluation
  nthreads=1
  CALL OMP_SET_NUM_THREADS(nthreads)
  
  call system_clock(time1)
  call system_clock(count_rate=rate)
  
  call random_seed()
  call random_number(ran)



  open(unit=12,file='params.dat')
    read(12,*) pot_file   ! name of potential file
    read(12,*) d          ! dimension
    read(12,*) bmax       ! basis pruning parameter
    read(12,*) pmax       ! point pruning parameter
    read(12,*) npot       ! # terms in potential
    read(12,*) out_file   ! name of output file = outfile_(bmax+1)_(pmax+1).dat
  close(unit=12)
  

  
  !------------- basis and point pruning indices ------------------
  
  
  
  allocate(imod(d),n(d))
  !generate basis/grid indices 
  imod=0

  bmax=bmax+1
  pmax=pmax+1
  nb=0
  np=0

  ! get basis function indices
  call indgen(d,bmax-1,nb,binds,imod,1)
  ! get grid indices
  call indgen(d,pmax-1,np,pinds,imod,1)


  
  !set max indices for all dimension = pmax, can vary for general case based on pruning scheme
  n=pmax
  
  print*, 'basis/grid size: ' , nb,np
  print*, bmax,pmax
  
  
  !-------------------- point selection and potential energy surface --------------------------
  

  allocate(fexp(npot,d),fcoef(npot))

  ! read potential file (here SOP)
  
  open(unit=12,file=pot_file)
  do i=1,npot
      read(12,*) fexp(i,:), fcoef(i)
  end do
  close(unit=12)
    

  ! point selection 
  allocate(grid(pmax,d))

  grid=0d0

  !choice of 1D ordered points
  open(unit=12,file='hobleja.dat')
    read(12,*) grid(:,1)
  close(unit=12)    
  
  ! set same point selection for all dimensions
  do i=2,d
   grid(:,i)=grid(:,1)
  end do



  ! generate potential at np points stored in potp vector
  allocate(potp(np),point(d))
  point=0d0
  potp=0d0

  do i=1,np
  
    do j=1,d
   
      ! point corresponding to the indices in position i in pinds array
      point(j)=grid(pinds(j,i),j) 
  
    end do 

    !call potential routine here, just a simple SOP
    potp(i)=pot(fexp,fcoef,npot,d,point)
    
  end do  


  !------------------- Matrices and Mapping --------------------------------------



  ! generate B and B'' matrices
  call buildBmats(d,n, pmax, Bmats, kterms, grid )

  ! generate LU decompositions of B and B'' matrices
  call buildLUmats(d,n,pmax,kterms,Bz,Bmats)
  
  ! generate inverses of the LU decompositions
  call buildINVmats(d,n,pmax,kterms,Bz,Binv)

  ! generates mapping arrays for the different index transformations in the MVPs
  call indmap(d,nb,np,bmax,pmax,binds,pinds,usort,dsort,ssort,sortlens,sortlims)


   
  !---------------------------- ARPACK section ------------------------------------
 
 
  print*, 'starting ARPACK..'
  nx=nb
  allocate(resid(nb))
  allocate(nmvp(nb))
  resid=1d0
  nmvp=0d0
  resid(1)=1d0


  do i=1,nx
    call random_number(ran)
    resid(i)=ran
  end do


  maxiter=150
  resid=0d0  
  iparam=0
  iparam(1)=1
  iparam(3)=2*maxiter
  iparam(4)=1
  iparam(7)=1
  nev=min(nx-2,50)
  ncv=min(nx,120)
  tol=0d0
  ido=0 
  which='SR'
  bmat='I'
  howmny='A'
  
  ipntr(1)=1
  info=0
  allocate(workd(3*nx))
  lworkl= 3*ncv**2+8*ncv
  allocate(workl(lworkl))
  workd(1:nx)=resid
  workd(nx+1:2*nx)=resid
  workd(2*nx+1:3*nx)=resid  
  
  allocate(v(nx,ncv))
  v=0d0
  add=.false.
 
  i2=0

  print*, 'starting dnaupd...'  
  
  call system_clock(time1)
  do while(ido < 99)

    call dnaupd(ido,bmat,nx,which,nev,tol,resid,ncv,v,nx,iparam,ipntr,workd,&
         workl,lworkl,info) 
         

  
    if (info <0) then
      print*, 'Error of type info = ', info
    end if

    !print*, info
   
    if (ido.eq.-1.or.ido.eq.1) then
      nmvp=0d0

      call mvp(nmvp,workd(ipntr(1):ipntr(1)+nx-1),np,nb,pmax,d,kterms,Bz,Binv,potp,sortlens,sortlims,usort,dsort,ssort,fcoef)

      workd(ipntr(2):ipntr(2)+nx-1) = nmvp
      i2=i2+1
    end if
  end do
  
  if (info==0) then
  
    print*, 'normal exit'

  end if

  print*, 'iterations done:',iparam(3)
  print*, '# converged: ', iparam(5)
  print*, 'ido ,' ,ido
  print*, 'info, ' ,info


  
  rvec=.True.
  howmny='A'
  allocate(dr(nev+1),dr2(nev+1),di(nev+1),zeig(nx,nev+1),workev(3*ncv))
  allocate(iselect(ncv))
  iselect=1 !.True.
  zeig=v(:,1:nev+1)
  dr=0d0
  di=0d0
 

  print*, 'starting dneupd...' 

  call dneupd2(0,howmny,iselect,dr,di,zeig,nx,sigmar,sigmai,workev,bmat,nx,which, &
       nev,tol,resid,ncv,v,nx,iparam,ipntr,workd,workl,lworkl,info)


!---------------------------- sorting and output ------------------------------------

  dr2=dr
  call sort(dr2(1:nev),nev)

  call system_clock(time2)

  write(bstring,*) bmax
  write(pstring,*) pmax

 
  
  !display file name, basis prune, point prune
  out_file=trim(out_file)//'_'//trim(adjustl(bstring))//'_'//trim(adjustl(pstring))//'.dat'
  

      
  open(unit=12,file=out_file)
  do i=1,nev
   write(12,*) dr2(i)
  end do
  write(12,*) ' '
  write(12,*)  1d0*(time2-time1)/(rate*1d0)
  close(unit=12)





 end program main  