!
! Example of sparse matrix-vector multiply, using CSR (compressed sparse
! row format).  Because this code is in C, the index values in the ia and ja
! arrays are zero origin rather than one origin as they would be in Fortran.
!
! While this code does time the sparse matrix-vector product, as noted below,
! these timings are too simple to be used for serious analysis.
!
      program main
      integer max_row, max_nnz
      parameter (max_row=10000000)
      parameter (max_nnz=50000000)
      integer ia(max_row), ja(max_nnz)
      double precision a(max_nnz)
      double precision x(max_row), y(max_row)
      integer row, i, j, k, idx, n
      double precision mysecond, t, tmin, sum
      character*32 argv

      n = 10
      call get_command_argument(1,argv)
      if (len_trim(argv) > 0) then
         read(argv,*) n
      endif

! Warning: if n > sqrt(2^31), you may get integer overflow

      nrows  = n * n
      nnzMax = nrows * 5
      if (nrows > max_row .or. nnzMax > max_nnz) then
         print *, 'n = ', n, ' is too large'
         stop
      endif

! Create a pentadiagonal matrix, representing very roughly a finite
! difference approximation to the Laplacian on a square n x n mesh
      row = 1
      nnz = 0
      do i=1,n
         do j=1,n
            ia(row) = nnz
            if (i>1) then
               ja(nnz) = row - n
               a(nnz) = -1.0
               nnz = nnz + 1
            endif
            if (j>1) then
               ja(nnz) = row - 1
               a(nnz) = -1.0
               nnz = nnz + 1
            endif
            ja(nnz) = row
            a(nnz) = 4.0
            nnz = nnz + 1
            if (j<n) then
               ja(nnz) = row + 1
               a(nnz) = -1.0
               nnz = nnz + 1
            endif
            if (i<n) then
               ja(nnz) = row + n
               a(nnz) = -1.0
               nnz = nnz + 1
            endif
            row = row + 1
         enddo
      enddo
      ia(row) = nnz

! Create the source (x) vector
      do i=1, nrows
         x(i) = 1.0d0
      enddo

! Perform a matrix-vector multiply: y = A*x
! Warning: To use this for timing, you need to (a) handle cold start
! (b) perform enough tests to make timing quantum relatively small
      tmin = 1.0e12
      do k=1, 10
         t = mysecond()
         do row=1, nrows
            sum = 0.0
            do idx=ia(row), ia(row+1)-1
               sum = sum + a(idx) * x(ja(idx))
            enddo
            y(row) = sum
         enddo
         t = mysecond() - t
         if (t < tmin) then
            tmin = t
         endif
! Compute with the result to keep the compiler for marking the
! code as dead
         do row=1, nrows
            if (y(row) .lt. 0.0) then
               print *, 'y(',row,') fails consistency test'
            endif
         enddo
      enddo

      print *, 'Time for Sparse Ax, nrows=', nrows, ', nnz=', nnz,        &
     &     ', T = ', tmin

      end
