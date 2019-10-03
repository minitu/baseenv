! This is a test of the nodecart shim file.
	program main
	use mpi
	integer wrank, wsize, dims(2), periods(2), i, ierror
	integer west, east, north, south, ndim, crank, coords(2)
	integer createcalls, subcalls
        ! Type(MPI_Comm) :: cartcomm
	integer cartcomm

	call mpi_init(ierror)


	call MPI_Comm_rank(MPI_COMM_WORLD, wrank, ierror)
	call MPI_Comm_size(MPI_COMM_WORLD, wsize, ierror)

        ndim = 2
	do i=1,ndim
	   dims(i)    = 0
	   periods(i) = 0
	enddo

	call MPI_Dims_create(wsize, ndim, dims, ierror)
	call MPI_Cart_create(MPI_COMM_WORLD, ndim, dims, periods, 1, &
                             cartcomm, ierror)
        call MPI_Cart_shift(cartcomm, 0, 1, west, east, ierror)
	call MPI_Cart_shift(cartcomm, 1, 1, north, south, ierror)
	call MPI_Comm_rank(cartcomm, crank, ierror)
	call MPI_Cart_coords(cartcomm, crank, ndim, coords, ierror)

        ! This is a special routine provided by the shim library to make it
        ! easier to confirm that the shim library has been correctly applied
        ! and that the routine interception has occured
	call MPIX_CartShim_info(createcalls, subcalls)
	if (wrank .eq. 0) then
	   if (createcalls .ne. 1 .and. subcalls .ne. 0) then
	        print *, "Test FAILED: createcalls = ", createcalls, &
                   "(should be 1) and subcalls = ", subcalls, &
                   "(should be 0)"
	    else
	        print *, "Test PASSED"
            endif
	endif

	call MPI_Comm_free(cartcomm, ierror)
	call MPI_Finalize(ierror)
end
