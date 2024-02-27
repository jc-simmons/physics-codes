module sortmod

    public sort, sortinds
    private

contains

   !  used in the sort() subroutine 
   integer function FindMinimum(x, Start, sto)
      IMPLICIT  NONE
      INTEGER, INTENT(IN)                :: Start, sto
      real*8 :: x(start-sto+1)
      real*8                             :: Minimum
      INTEGER                            :: Location
      INTEGER                            :: i

      Minimum  = x(Start)       ! assume the first is the min
      Location = Start          ! record its position
      DO i = Start+1, sto       ! start with next elements
         IF (x(i) < Minimum) THEN   !   if x(i) less than the min?
            Minimum  = x(i)     !      Yes, a new minimum found
            Location = i                !      record its position
         END IF
      END DO
      FindMinimum = Location            ! return the position
   END FUNCTION  FindMinimum



   ! --------------------------------------------------------------------
   ! SUBROUTINE  Swap():
   !    This subroutine swaps the values of its two formal arguments.
   ! --------------------------------------------------------------------

   SUBROUTINE  Swap(a, b)
      IMPLICIT  NONE
      real*8, INTENT(INOUT) :: a, b
      real*8                :: Temp

      Temp = a
      a    = b
      b    = Temp
   END SUBROUTINE  Swap

   ! --------------------------------------------------------------------
   ! SUBROUTINE  Sort():
   !    This subroutine receives an array x() and sorts it into ascending
   ! order.
   ! --------------------------------------------------------------------


   SUBROUTINE  Sort(x, len)
      IMPLICIT  NONE
      real*8 :: x(len)
      INTEGER, INTENT(IN)                   :: len
      INTEGER                               :: i
      INTEGER                               :: Location
      !INTEGER :: FindMinimum


      DO i = 1, len-1           ! except for the last
         Location = FindMinimum(x, i, len)  ! find min from this to last
         CALL  Swap(x(i), x(Location))  ! swap this and the minimum
      END DO
   END SUBROUTINE  Sort




   subroutine sortinds(presort,postsort,sortlims,d,np,n,k,perm)
   
      integer :: n,d,m,np
      integer :: i,j,k(2)
      integer :: presort(d,np), postsort(np,n),inter(d-1),sortlims(np),perm(np)
      
      
   
   
      postsort=0
      
   
   
      inter=0
   
      k(2)=m
   
   
      sortlims=0
      postsort=0
      
      inter(:)=presort(1:d-1,1)
      
      j=1
      
      
      do i=1,k(1)
      
      if (all(inter(:)==presort(1:d-1,i))) then
      
      
            sortlims(j)=sortlims(j)+1
            postsort(j,sortlims(j))=perm(i)
   
         
      else
      
         j=j+1
         
         sortlims(j)=sortlims(j)+1
         postsort(j,sortlims(j))=perm(i)
         inter(:)=presort(1:d-1,i)
      
   
      end if
                              
      end do
      
      k(2)=j
      
   end subroutine sortinds


end module sortmod