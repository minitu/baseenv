!
! This program explores array relaxation as a function of size and number of
! cores/processes running independent sweeps.
 program blksweep
   use mpi
   implicit none
      double precision mysecond, t, tmin, tprev, rate, tminb, tmin2d, tminlx, &
           tminlx2, rateb, rate2d, ratelx, ratelx2
      character*32 argv
      double precision, allocatable :: a(:,:,:), b(:,:,:)
      double precision :: a1, a2, a3, b1, b2, b3
      integer maxn, i, j, k, ktest, n, bsize
      integer ii, jj, kk
      integer provided, wsize, wrank, color, ierr, nc, comm
      logical debug
      interface
         subroutine initarrays(a,b,n)
           double precision, allocatable :: a(:,:,:), b(:,:,:)
           integer n
         end subroutine initarrays
         subroutine sdummy(n,b,a)
           double precision, allocatable :: a(:,:,:), b(:,:,:)
           integer n
         end subroutine sdummy
      end interface

      debug = .false.
      maxn = 128
      call get_command_argument(1,argv)
      if (len_trim(argv) > 0) then
         read(argv,*) maxn
      endif

      call MPI_Init_thread(MPI_THREAD_SINGLE, provided, ierr)
      call MPI_Comm_size(MPI_COMM_WORLD, wsize, ierr)
      call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierr)

      ! For each number of processes (cores), starting at 1, to wsize
      do nc=1,wsize
         if (debug .and. wrank == 0) print *, "nc = ", nc
         call MPI_Barrier(MPI_COMM_WORLD, ierr)
         if (wrank < nc) then
            color = 1
         else
            color = 0
         endif
         call MPI_Comm_split(MPI_COMM_WORLD, color, wrank, comm, ierr)
         if (debug .and. wrank == 0) print *, "split complete"
         if (wrank >= nc) then
            call MPI_Comm_free(comm, ierr)
            cycle
         endif

         n = 8
         bsize = 32
         tprev = 0.0
         do while (n <= maxn)
            allocate(a(0:n+1,0:n+1,0:n+1))
            allocate(b(0:n+1,0:n+1,0:n+1))
!
            call initarrays(a,b,n)
            tmin = 1.0e10
            do ktest=1, 10
               call MPI_Barrier(comm, ierr)
               t = mysecond()
               do k=1,n
                  do j=1,n
                     do i=1,n
                        b(i,j,k) = a(i+1,j,k)+a(i-1,j,k)+a(i,j+1,k)+          &
                             a(i,j-1,k)+a(i,j,k+1)+a(i,j,k-1) -6.0*a(i,j,k)
                     enddo
                  enddo
               enddo
               t = mysecond() - t
               if (t < tmin) tmin = t
               call MPI_Allreduce(MPI_IN_PLACE, tmin, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, comm, ierr)
               call sdummy(n,b,a)
            enddo
!
            call initarrays(a,b,n)
            tminb = 1.0e10
            bsize = n
            if (n > 32) bsize = n/2
            do ktest=1, 10
               call MPI_Barrier(comm, ierr)
               t = mysecond()
               do kk=1,n,bsize
                  do jj=1,n,bsize
                     do ii=1,n,bsize
                        do k=kk,min(kk+bsize-1,n)
                           do j=jj,min(jj+bsize-1,n)
                              do i=ii,min(ii+bsize-1,n)
                                 b(i,j,k) = a(i+1,j,k)+a(i-1,j,k)+a(i,j+1,k)+  &
                                      a(i,j-1,k)+a(i,j,k+1)+a(i,j,k-1) -       &
                                      6.0*a(i,j,k)
                              enddo
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               t = mysecond() - t
               if (t < tminb) tminb = t
               call MPI_Allreduce(MPI_IN_PLACE, tminb, 1, MPI_DOUBLE_PRECISION,&
                    MPI_MIN, comm, ierr)
               call sdummy(n,b,a)
            enddo
!
            call initarrays(a,b,n)
            tmin2d = 1.0e10
            bsize  = n
            if (n > 32) bsize = n/2
            do ktest=1, 10
               call MPI_Barrier(comm, ierr)
               t = mysecond()
               do kk=1,n,bsize
                  do jj=1,n,bsize
                     do k=kk,min(kk+bsize-1,n)
                        do j=jj,min(jj+bsize-1,n)
                           do i=1,n
                              b(i,j,k) = a(i+1,j,k)+a(i-1,j,k)+a(i,j+1,k)+  &
                                   a(i,j-1,k)+a(i,j,k+1)+a(i,j,k-1) -       &
                                   6.0*a(i,j,k)
                           enddo
                        enddo
                     enddo
                  enddo
               enddo
               t = mysecond() - t
               if (t < tmin2d) tmin2d = t
               call MPI_Allreduce(MPI_IN_PLACE, tmin2d, 1,                  &
                    MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
               call sdummy(n,b,a)
            enddo
!
            call initarrays(a,b,n)
            tminlx = 1.0e10
            do ktest=1, 10
               call MPI_Barrier(comm, ierr)
               t = mysecond()
               do k=1,n
                  do j=1,n
                     a1 =a(0,j,k)
                     a2 =a(1,j,k)
                     do i=1,n
                        a3 = a(i+1,j,k)
                        b(i,j,k) = a3+a1+a(i,j+1,k)+        &
                             a(i,j-1,k)+a(i,j,k+1)+a(i,j,k-1) -6.0*a2
                        a1=a2
                        a2=a3
                     enddo
                  enddo
               enddo
               t = mysecond() - t
               if (t < tminlx) tminlx = t
               call MPI_Allreduce(MPI_IN_PLACE, tminlx, 1,                  &
                    MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
               call sdummy(n,b,a)
            enddo
!
            call initarrays(a,b,n)
            tminlx2 = 1.0e10
            do ktest=1, 10
               call MPI_Barrier(comm, ierr)
               t = mysecond()
               do k=1,n
                  do j=1,n,2
                     a1 =a(0,j,k)
                     a2 =a(1,j,k)
                     b1 =a(0,j+1,k)
                     b2 =a(1,j+1,k)
                     do i=1,n
                        a3 = a(i+1,j,k)
                        b3 = a(i+1,j+1,k)
                        b(i,j,k) = a3+a1+b2+                       &
                             a(i,j-1,k)+a(i,j,k+1)+a(i,j,k-1) -6.0*a2
                        b(i,j+1,k)= b3+b1+a(i,j+2,k)+                      &
                             a2+a(i,j+1,k+1)+a(i,j,k-1)-6.0*b2
                        a1=a2
                        a2=a3
                        b1=b2
                        b2=b3
                     enddo
                  enddo
               enddo
               t = mysecond() - t
               if (t < tminlx2) tminlx2 = t
               call MPI_Allreduce(MPI_IN_PLACE, tminlx2, 1,                  &
                    MPI_DOUBLE_PRECISION, MPI_MIN, comm, ierr)
               call sdummy(n,b,a)
            enddo
!
            deallocate(a)
            deallocate(b)
            ! Rate in MF
            if (wrank == 0) then
               rate    = 7e-6*(n/tmin)*(1.0*n)*(1.0*n)
               rateb   = 7e-6*(n/tminb)*(1.0*n)*(1.0*n)
               rate2d  = 7e-6*(n/tmin2d)*(1.0*n)*(1.0*n)
               ratelx  = 7e-6*(n/tminlx)*(1.0*n)*(1.0*n)
               ratelx2 = 7e-6*(n/tminlx2)*(1.0*n)*(1.0*n)
               if (tprev > 0) then
                  print *, nc, n, tmin, rate, tminb, rateb, tmin2d, rate2d, &
                       tminlx, ratelx, tminlx2, ratelx2, tmin/tprev
               else
                  print *, nc, n, tmin, rate, tminb, rateb, tmin2d, rate2d, &
                       tminlx, ratelx, tminlx2, ratelx2
               endif
               tprev = tmin
            endif
            n = n * 2
         enddo
         call MPI_Comm_free(comm, ierr)
      enddo
 
 end program blksweep

 subroutine initarrays(a, b, n)
   integer n
   double precision, allocatable :: a(:,:,:), b(:,:,:)
   integer i, j, k
   do k=0,n+1
      do j=0,n+1
         do i=0,n+1
            a(i,j,k) = i + j + k
            b(i,j,k) = 0
         enddo
      enddo
   enddo
 end subroutine initarrays
