
!! Compile with (gnu fortran): mpif90 -cpp -O3 chunk_cubep3m_v1.5.f90 -o chunk -DGFORTRAN
!! Compile with (intel fortran): mpif90 -fpp -O3 chunk_cubep3m_v1.5.f90 -o chunk -DBINARY



program chunk_cubep3m 

  implicit none
  include 'mpif.h'

  integer(4) :: PID_FLAG
  real :: z

  integer(4) :: chunk_x
  integer(4) :: chunk_y
  integer(4) :: chunk_z

  real(4) :: buffer_size

  ! Cubep3m file variables

  character(400) :: file, ifile, redshift, particle_data_path, chunk_output_path
  integer(4) :: np_local, nts, cur_checkpoint, cur_projection, cur_halofind
  real(4) :: a, t, tau, dt_f_acc, dt_pp_acc, dt_c_acc, mass_p
  real(4), allocatable, dimension(:,:) :: xv
  integer(8), allocatable, dimension(:) :: PID  
  integer(4), parameter :: blocksize = 32*1024*1024/24 !32*1024*1024/24 ! reading blocksize in bytes default is 32Mb
  integer(4) :: read_offset, num_reads
  real(4) :: ncc, box
  integer(4) :: nodes, nc, nodes_dim
  real(4) cubep3m_x_min, cubep3m_x_max, cubep3m_y_min, cubep3m_y_max, cubep3m_z_min, cubep3m_z_max
  integer(4), allocatable :: node_coords(:,:)


  ! MPI

  integer nodes_returned, ierr, rank

  ! chunk coordinate bounds

  real(4),  allocatable :: x_min(:), x_max(:)
  real(4), allocatable :: y_min(:), y_max(:)
  real(4), allocatable :: z_min(:), z_max(:)
  real(4) x_chunk_length, y_chunk_length, z_chunk_length

  integer(4), allocatable :: chunk_coords(:,:)
  integer(4) :: no_of_chunks

  integer(4), allocatable :: chunk_flag(:,:,:) ! array of flags (1 or 0) if chunk file is required
  integer(8), allocatable :: chunk_particles(:,:,:,:)
  integer(8), allocatable :: chunk_buffer_particles(:,:,:,:)
  integer(8), allocatable :: chunk_tot_particles(:,:,:,:)
  integer(4) :: no_chunk_files_req, chunk_file_int
  integer, allocatable :: file_number(:)
  integer, allocatable :: chunk_flag_global(:,:,:,:), chunk_flag_local(:,:,:,:)

  ! buffer coodinate bounds

  real(4), allocatable :: buf_x_min(:), buf_x_max(:)
  real(4), allocatable :: buf_y_min(:), buf_y_max(:)
  real(4), allocatable :: buf_z_min(:), buf_z_max(:)
  real(4) x_chunk_buf_length, y_chunk_buf_length, z_chunk_buf_length

  ! Particle numbers for checking

  integer(8) np_total, np_check, np_local_check, np_local_buffer, np_local_total
  integer(8), allocatable :: tot_np_in_chunk(:)
  integer(8), allocatable :: buf_np_in_chunk(:)  
  integer(8), allocatable :: check_np_in_chunk(:)  
  integer(8), allocatable :: tot_np_for_each_chunk_by_rank(:)
  integer(8), allocatable :: buf_np_for_each_chunk_by_rank(:)
  integer(8), allocatable :: check_np_for_each_chunk_by_rank(:)

  real(4), allocatable :: max_x(:), max_y(:), max_z(:), min_x(:), min_y(:), min_z(:)
  real(4), allocatable :: max_x_global(:), max_y_global(:), max_z_global(:)
  real(4), allocatable :: min_x_global(:), min_y_global(:), min_z_global(:)

  ! utility variables

  integer(4) :: i, j, k, n, p, dummy
  real(4) :: dummy_real

  !------------------------------------------------------------------------------------------!

  call mpi_initialise

  if(rank .eq. 0) then

    call begin_print_screen

  endif

  if(rank .eq. 0) call read_parameters

  call broadcast_parameters

  call variable_initialise


  call mpi_barrier(mpi_comm_world,ierr)

  call node_coords_initialise
  call initialise_chunks
  call chunk_file_calc ! calculate which chunk files this node overlaps with
  call chunk_file_calc_buf ! calculate which chunk file buffer zones this node overlaps with
  call open_chunk_files

  call mpi_barrier(mpi_comm_world,ierr)

  call open_cubep3m_file
  call read_cubep3m_header

  if(PID_FLAG .eq. 1) call open_cubep3m_PID_file
  if(PID_FLAG .eq. 1) call read_cubep3m_PID_header

  call mpi_barrier(mpi_comm_world,ierr)

  call total_particle_calculation    
  call write_chunk_file_header_draft

  ! Now let's process the particles

  read_offset = 0 ! variable for handling large files in smaller sections given by blocksize

  num_reads = np_local/blocksize+1

  do n = 1, num_reads

    if(rank .eq. 0) print*, 'Reading ',n,'of ',num_reads,' on task 0'

    if(np_local - read_offset .gt. blocksize) then ! if we still have more than the blocksize data remaining

      allocate(xv(6,blocksize))
      if(PID_FLAG .eq. 1) allocate(PID(blocksize))
      read(10) xv
      if(PID_FLAG .eq. 1) read(11) PID           
      read_offset = read_offset + blocksize

      call process_particles(blocksize)

      deallocate(xv)
      if(PID_FLAG .eq. 1) deallocate(PID)

    else ! if we are writing the last block (smaller than blocksize)

      allocate(xv(6,np_local - read_offset))
      if(PID_FLAG .eq. 1) allocate(PID(np_local - read_offset))
      read(10) xv
      if(PID_FLAG .eq. 1) read(11) PID

      call process_particles(np_local - read_offset)

      deallocate(xv)
      if(PID_FLAG .eq. 1) deallocate(PID)        

    endif    

  enddo

  call close_cubep3m_file

  call particle_number_check
  call write_chunk_file_header_final

  call close_chunk_files

  call particle_bound_check



  call finalise_chunker



contains

  !------------------------------------------------------------------------------------------!

  subroutine mpi_initialise

      call mpi_init(ierr)
      if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

      call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
      if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)


      call mpi_comm_rank(mpi_comm_world,rank,ierr)
      if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

  end subroutine mpi_initialise

  !------------------------------------------------------------------------------------------!


  subroutine variable_initialise

      nodes = nodes_dim*nodes_dim*nodes_dim
      no_of_chunks = chunk_x * chunk_y * chunk_z

      ncc = real(nc)/real(nodes_dim)

      buffer_size = buffer_size * ncc / box

      if (nodes_returned /= nodes ) then
	write(*,*) 'chunker compiled for a different number of nodes'
	write(*,*) 'mpirun nodes=',nodes_returned,'chunker nodes=',nodes 
	call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

      allocate(node_coords(0:nodes-1,3))
      allocate(x_min(chunk_x))
      allocate(x_max(chunk_x))
      allocate(y_min(chunk_y))
      allocate(y_max(chunk_y))
      allocate(z_min(chunk_z))
      allocate(z_max(chunk_z))
      allocate(chunk_coords(0:chunk_x * chunk_y * chunk_z,3))
      allocate(chunk_flag(chunk_x, chunk_y, chunk_z))
      allocate(chunk_particles(0:no_of_chunks-1,chunk_x, chunk_y, chunk_z))
      allocate(chunk_buffer_particles(0:no_of_chunks-1,chunk_x, chunk_y, chunk_z))
      allocate(chunk_tot_particles(0:no_of_chunks-1,chunk_x, chunk_y, chunk_z))
      allocate(file_number(0:no_of_chunks-1))
      allocate(chunk_flag_global(0:nodes-1, chunk_x, chunk_y, chunk_z))
      allocate(chunk_flag_local(0:nodes-1, chunk_x, chunk_y, chunk_z))
      allocate(buf_x_min(chunk_x))
      allocate(buf_x_max(chunk_x))
      allocate(buf_y_min(chunk_y))
      allocate(buf_y_max(chunk_y))
      allocate(buf_z_min(chunk_z))
      allocate(buf_z_max(chunk_z))
      allocate(tot_np_in_chunk(0:no_of_chunks-1))
      allocate(buf_np_in_chunk(0:no_of_chunks-1))
      allocate(check_np_in_chunk(0:no_of_chunks-1))
      allocate(tot_np_for_each_chunk_by_rank(0:no_of_chunks-1))
      allocate(buf_np_for_each_chunk_by_rank(0:no_of_chunks-1))
      allocate(check_np_for_each_chunk_by_rank(0:no_of_chunks-1))
      allocate(max_x(0:no_of_chunks-1))
      allocate(min_x(0:no_of_chunks-1))
      allocate(max_y(0:no_of_chunks-1))
      allocate(min_y(0:no_of_chunks-1))
      allocate(max_z(0:no_of_chunks-1))
      allocate(min_z(0:no_of_chunks-1))
      allocate(max_x_global(0:no_of_chunks-1))
      allocate(min_x_global(0:no_of_chunks-1))
      allocate(max_y_global(0:no_of_chunks-1))
      allocate(min_y_global(0:no_of_chunks-1))
      allocate(max_z_global(0:no_of_chunks-1))
      allocate(min_z_global(0:no_of_chunks-1))




      np_total = 0
      chunk_particles = 0 
      chunk_buffer_particles = 0
      np_check = 0
      np_local_check = 0
      np_local_buffer = 0
      np_local_total = 0
      chunk_flag = 0
      tot_np_for_each_chunk_by_rank = 0
      buf_np_for_each_chunk_by_rank = 0  
      check_np_for_each_chunk_by_rank = 0  
      tot_np_in_chunk = 0
      buf_np_in_chunk = 0
      check_np_in_chunk = 0

      write(redshift,'(f8.3)') z
      redshift=adjustl(redshift)	

      max_x = 0
      max_y = 0
      max_z = 0
      min_x = 1000
      min_y = 1000
      min_z = 1000

      max_x_global = 0
      max_y_global = 0
      max_z_global = 0
      min_x_global = 1000
      min_y_global = 1000
      min_z_global = 1000







  end subroutine variable_initialise


  !------------------------------------------------------------------------------------------!


  subroutine begin_print_screen


      print*
      print*
      print*, '       |--------------------------------|'
      print*, '       |                                |'
      print*, '       |         CUBEP3M CHUNKER        |'
      print*, '       |                                |'
      print*, '       |                                |'
      print*, '       |              by                |'   
      print*, '       |                                |'
      print*, '       |             The                |'
      print*, '       |                                |'
      print*, '       |            Sirius              |'
      print*, '       |                                |'                      
      print*, '       |          Cybernetics           |'
      print*, '       |                                |'
      print*, '       |        Corporation (TM)        |'
      print*, '       |                                |'                
      print*, '       |--------------------------------|'
      print*
      print*
      print*
  end subroutine begin_print_screen

  !------------------------------------------------------------------------------------------!

  subroutine node_coords_initialise

      integer i, j, k, n

      !! Calculate node_coords
      do n=0,nodes_dim**3
	do k=1,nodes_dim
	  do j=1,nodes_dim
	    do i=1,nodes_dim
	      if (n == (i-1)+(j-1)*nodes_dim+(k-1)*nodes_dim**2)  then
		node_coords(n,:)=(/(i-1),(j-1),(k-1)/)
	      endif 
	    enddo
	  enddo
	enddo

      enddo

  end subroutine node_coords_initialise

  !------------------------------------------------------------------------------------------!

  subroutine initialise_chunks

      integer i, j, k, n

      ! Set up coordinates of chunks

      x_chunk_length = real(nc) / real(chunk_x)
      y_chunk_length = real(nc) / real(chunk_y)
      z_chunk_length = real(nc) / real(chunk_z)
      x_chunk_buf_length = real(nc)*(1+2*buffer_size) / real(chunk_x)
      y_chunk_buf_length = real(nc)*(1+2*buffer_size) / real(chunk_y)
      z_chunk_buf_length = real(nc)*(1+2*buffer_size) / real(chunk_z)

      do i = 1 , chunk_x
	x_min(i) = (i-1) * real(nc) / real(chunk_x)
	x_max(i) = i * real(nc) / real(chunk_x)
      enddo
      do i = 1 , chunk_y
	y_min(i) = (i-1) * real(nc) / real(chunk_y)
	y_max(i) = i * real(nc) / real(chunk_y)
      enddo
      do i = 1 , chunk_z
	z_min(i) = (i-1) * real(nc) / real(chunk_z)
	z_max(i) = i * real(nc) / real(chunk_z)
      enddo

      ! label chunk files (x then y then z):

      do n = 0, no_of_chunks
	do k=1,chunk_z
	  do j=1,chunk_y
	    do i=1,chunk_x
	      if (n == (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y)  then
		chunk_coords(n,:)=(/(i-1),(j-1),(k-1)/)
	      endif
	    enddo
	  enddo
	enddo		
      enddo


#ifdef DEBUG

      if(rank .eq. 0) then
	print*, 'chunk coords are:'

	do k = 1, chunk_z
	  do j = 1, chunk_y
	    do i = 1, chunk_x

	      print*,'x', x_min(i), x_max(i), 'y', y_min(j), y_max(j),&
		  & 'z', z_min(k), z_max(k)

	    enddo
	  enddo
	enddo
	print*

      endif
#endif


      ! set up buffers for chunks

      do i = 1, chunk_x
	buf_x_min(i) = x_min(i) - buffer_size 
	buf_x_max(i) = x_max(i) + buffer_size 
      enddo
      do i = 1, chunk_y
	buf_y_min(i) = y_min(i) - buffer_size 
	buf_y_max(i) = y_max(i) + buffer_size
      enddo
      do i = 1, chunk_z
	buf_z_min(i) = z_min(i) - buffer_size
	buf_z_max(i) = z_max(i) + buffer_size
      enddo

      ! set the new bounds to the buffer bounds (NB at this point the max and min
      ! buffer bounds are larger than the box size)



#ifdef DEBUG
      if(rank .eq. 0) then
	print*, 'chunk buffer coords are:'

	do k = 1, chunk_z
	  do j = 1, chunk_y
	    do i = 1, chunk_x

	      print*,'x', buf_x_min(i), buf_x_max(i), 'y', buf_y_min(j), buf_y_max(j),&
		  & 'z', buf_z_min(k), buf_z_max(k)

	    enddo
	  enddo
	enddo
      endif
#endif

  end subroutine initialise_chunks

  !------------------------------------------------------------------------------------------!

  subroutine chunk_file_calc

      integer i, j, k

      ! How many chunk files do we need to open?

      do i = 1, chunk_x 	
	if(x_chunk_length .ge. ncc) then

	  if((node_coords(rank,1)*ncc .ge. buf_x_min(i) .and. &			!if the cubep3m and chunk overlap in x
	      & node_coords(rank,1)*ncc .lt. buf_x_max(i)) .or. &
	      & (node_coords(rank,1)*ncc+ncc .gt. buf_x_min(i) .and. &
	      & node_coords(rank,1)*ncc+ncc .lt. buf_x_max(i))) then

	  do j = 1, chunk_y		
	    if(y_chunk_length .ge. ncc) then
	      if((node_coords(rank,2)*ncc .ge. buf_y_min(j) .and. &			!if the cubep3m and chunk overlap in y
		  & node_coords(rank,2)*ncc .lt. buf_y_max(j)) .or. &
		  & (node_coords(rank,2)*ncc+ncc .gt. buf_y_min(j) .and. &
		  & node_coords(rank,2)*ncc+ncc .lt. buf_y_max(j))) then

	      do k = 1, chunk_z		
		if(z_chunk_length .ge. ncc) then
		  if((node_coords(rank,3)*ncc .ge. buf_z_min(k) .and. &			!if the cubep3m and chunk overlap in z
		      & node_coords(rank,3)*ncc .lt. buf_z_max(k)) .or. &
		      & (node_coords(rank,3)*ncc+ncc .gt. buf_z_min(k) .and. &
		      & node_coords(rank,3)*ncc+ncc .lt. buf_z_max(k))) then

		  chunk_flag(i,j,k) = 1
		endif
	      else
		if((buf_z_min(k) .ge. node_coords(rank,3)*ncc .and. &			!if the cubep3m and chunk overlap in z
		    & buf_z_min(k) .lt. node_coords(rank,3)*ncc+ncc) .or. &
		    & (buf_z_max(k) .gt. node_coords(rank,3)*ncc .and. &
		    &  buf_z_max(k) .lt. node_coords(rank,3)*ncc+ncc)) then

		chunk_flag(i,j,k) = 1

	      endif
	    endif
	  enddo	

	endif
      else
	if((buf_y_min(j) .ge. node_coords(rank,2)*ncc .and. &			!if the cubep3m and chunk overlap in y
	    & buf_y_min(j) .lt. node_coords(rank,2)*ncc+ncc) .or. &
	    & (buf_y_max(j) .gt. node_coords(rank,2)*ncc .and. &
	    &  buf_y_max(j) .lt. node_coords(rank,2)*ncc+ncc)) then

	do k = 1, chunk_z		
	  if(z_chunk_length .ge. ncc) then
	    if((node_coords(rank,3)*ncc .ge. buf_z_min(k) .and. &			!if the cubep3m and chunk overlap in z
		& node_coords(rank,3)*ncc .lt. buf_z_max(k)) .or. &
		& (node_coords(rank,3)*ncc+ncc .gt. buf_z_min(k) .and. &
		& node_coords(rank,3)*ncc+ncc .lt. buf_z_max(k))) then

	    chunk_flag(i,j,k) = 1

	  endif
	else
	  if((buf_z_min(k) .ge. node_coords(rank,3)*ncc .and. &			!if the cubep3m and chunk overlap in z
	      & buf_z_min(k) .lt. node_coords(rank,3)*ncc+ncc) .or. &
	      & (buf_z_max(k) .gt. node_coords(rank,3)*ncc .and. &
	      &  buf_z_max(k) .lt. node_coords(rank,3)*ncc+ncc)) then

	  chunk_flag(i,j,k) = 1

	endif
      endif
    enddo

  endif
endif
				enddo
			      endif
			    else 
			      if((buf_x_min(i) .ge. node_coords(rank,1)*ncc .and. &			!if the cubep3m and chunk overlap in x
				  & buf_x_min(i) .lt. node_coords(rank,1)*ncc+ncc) .or. &
				  & (buf_x_max(i) .gt. node_coords(rank,1)*ncc .and. &
				  &  buf_x_max(i) .lt. node_coords(rank,1)*ncc+ncc)) then

			      do j = 1, chunk_y	
				if(y_chunk_length .ge. ncc) then
				  if((node_coords(rank,2)*ncc .ge. buf_y_min(j) .and. &			!if the cubep3m and chunk overlap in y
				      & node_coords(rank,2)*ncc .lt. buf_y_max(j)) .or. &
				      & (node_coords(rank,2)*ncc+ncc .gt. buf_y_min(j) .and. &
				      & node_coords(rank,2)*ncc+ncc .lt. buf_y_max(j))) then

				  do k = 1, chunk_z		
				    if(z_chunk_length .ge. ncc) then
				      if((node_coords(rank,3)*ncc .ge. buf_z_min(k) .and. &			!if the cubep3m and chunk overlap in z
					  & node_coords(rank,3)*ncc .lt. buf_z_max(k)) .or. &
					  & (node_coords(rank,3)*ncc+ncc .gt. buf_z_min(k) .and. &
					  & node_coords(rank,3)*ncc+ncc .lt. buf_z_max(k))) then

				      chunk_flag(i,j,k) = 1

				    endif
				  else
				    if((buf_z_min(k) .ge. node_coords(rank,3)*ncc .and. &			!if the cubep3m and chunk overlap in z
					& buf_z_min(k) .lt. node_coords(rank,3)*ncc+ncc) .or. &
					& (buf_z_max(k) .gt. node_coords(rank,3)*ncc .and. &
					&  buf_z_max(k) .lt. node_coords(rank,3)*ncc+ncc)) then

				    chunk_flag(i,j,k) = 1

				  endif
				endif
			      enddo	

			    endif
			  else
			    if((buf_y_min(j) .ge. node_coords(rank,2)*ncc .and. &			!if the cubep3m and chunk overlap in y
				& buf_y_min(j) .lt. node_coords(rank,2)*ncc+ncc) .or. &
				& (buf_y_max(j) .gt. node_coords(rank,2)*ncc .and. &
				&  buf_y_max(j) .lt. node_coords(rank,2)*ncc+ncc)) then

			    do k = 1, chunk_z		
			      if(z_chunk_length .ge. ncc) then
				if((node_coords(rank,3)*ncc .ge. buf_z_min(k) .and. &			!if the cubep3m and chunk overlap in z
				    & node_coords(rank,3)*ncc .lt. buf_z_max(k)) .or. &
				    & (node_coords(rank,3)*ncc+ncc .gt. buf_z_min(k) .and. &
				    & node_coords(rank,3)*ncc+ncc .lt. buf_z_max(k))) then

				chunk_flag(i,j,k) = 1										
			      endif
			    else
			      if((buf_z_min(k) .ge. node_coords(rank,3)*ncc .and. &			!if the cubep3m and chunk overlap in z
				  & buf_z_min(k) .lt. node_coords(rank,3)*ncc+ncc) .or. &
				  & (buf_z_max(k) .gt. node_coords(rank,3)*ncc .and. &
				  &  buf_z_max(k) .lt. node_coords(rank,3)*ncc+ncc)) then

			      chunk_flag(i,j,k) = 1

			    endif
			  endif
			enddo	

		      endif
		    endif
		  enddo	
		endif
	      endif
	    enddo

	end subroutine chunk_file_calc

	!------------------------------------------------------------------------------------------!



	subroutine chunk_file_calc_buf

	    integer i, j, k

	    do i = 1, chunk_x ! +x direction
	      if(buf_x_min(i) .lt. 0 ) then 

		buf_x_min(i) = buf_x_min(i) + nc
		buf_x_max(i) = buf_x_max(i) + nc


		call chunk_file_calc

		! put the bounds back

		buf_x_min(i) = buf_x_min(i) - nc
		buf_x_max(i) = buf_x_max(i) - nc

	      endif  
	    enddo

	    do i = 1, chunk_x ! -x direction
	      if(buf_x_max(i) .gt. nc  ) then

		buf_x_min(i) = buf_x_min(i) - nc
		buf_x_max(i) = buf_x_max(i) - nc

		call chunk_file_calc

		! put the bounds back

		buf_x_min(i) = buf_x_min(i) + nc
		buf_x_max(i) = buf_x_max(i) + nc

	      endif  
	    enddo

	    do i = 1, chunk_y ! +y direction
	      if(buf_y_min(i) .lt. 0 ) then 

		buf_y_min(i) = buf_y_min(i) + nc
		buf_y_max(i) = buf_y_max(i) + nc

		call chunk_file_calc

		! put the bounds back

		buf_y_min(i) = buf_y_min(i) - nc
		buf_y_max(i) = buf_y_max(i) - nc

	      endif  
	    enddo

	    do i = 1, chunk_y ! -y direction
	      if(buf_y_max(i) .gt. nc  ) then

		buf_y_min(i) = buf_y_min(i) - nc
		buf_y_max(i) = buf_y_max(i) - nc

		call chunk_file_calc

		! put the bounds back

		buf_y_min(i) = buf_y_min(i) + nc
		buf_y_max(i) = buf_y_max(i) + nc

	      endif  
	    enddo

	    do i = 1, chunk_z ! +z direction
	      if(buf_z_min(i) .lt. 0 ) then 

		buf_z_min(i) = buf_z_min(i) + nc
		buf_z_max(i) = buf_z_max(i) + nc

		call chunk_file_calc

		! put the bounds back

		buf_z_min(i) = buf_z_min(i) - nc
		buf_z_max(i) = buf_z_max(i) - nc

	      endif  
	    enddo

	    do i = 1, chunk_z ! -z direction
	      if(buf_z_max(i) .gt. nc  ) then

		buf_z_min(i) = buf_z_min(i) - nc
		buf_z_max(i) = buf_z_max(i) - nc

		call chunk_file_calc

		! put the bounds back

		buf_z_min(i) = buf_z_min(i) + nc
		buf_z_max(i) = buf_z_max(i) + nc

	      endif  
	    enddo

	    do i = 1, chunk_x ! +x+y direction
	      do j = 1, chunk_y
		if(buf_x_min(i) .lt. 0 .and. buf_y_min(j) .lt. 0) then 

		  buf_x_min(i) = buf_x_min(i) + nc
		  buf_x_max(i) = buf_x_max(i) + nc

		  buf_y_min(j) = buf_y_min(j) + nc
		  buf_y_max(j) = buf_y_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_x_min(i) = buf_x_min(i) - nc
		  buf_x_max(i) = buf_x_max(i) - nc
		  buf_y_min(j) = buf_y_min(j) - nc
		  buf_y_max(j) = buf_y_max(j) - nc

		endif  
	      enddo
	    enddo


	    do i = 1, chunk_x ! -x+y direction
	      do j = 1, chunk_y
		if(buf_x_max(i) .gt. nc .and. buf_y_min(j) .lt. 0) then 

		  buf_x_min(i) = buf_x_min(i) - nc
		  buf_x_max(i) = buf_x_max(i) - nc

		  buf_y_min(j) = buf_y_min(j) + nc
		  buf_y_max(j) = buf_y_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_x_min(i) = buf_x_min(i) + nc
		  buf_x_max(i) = buf_x_max(i) + nc
		  buf_y_min(j) = buf_y_min(j) - nc
		  buf_y_max(j) = buf_y_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_x ! +x-y direction
	      do j = 1, chunk_y
		if(buf_y_max(i) .gt. nc .and. buf_x_min(j) .lt. 0) then 

		  buf_y_min(i) = buf_y_min(i) - nc
		  buf_y_max(i) = buf_y_max(i) - nc

		  buf_x_min(j) = buf_x_min(j) + nc
		  buf_x_max(j) = buf_x_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_y_min(i) = buf_y_min(i) + nc
		  buf_y_max(i) = buf_y_max(i) + nc
		  buf_x_min(j) = buf_x_min(j) - nc
		  buf_x_max(j) = buf_x_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_x ! -x-y direction
	      do j = 1, chunk_y
		if(buf_y_max(i) .gt. nc .and. buf_x_max(j) .gt. nc) then 

		  buf_y_min(i) = buf_y_min(i) - nc
		  buf_y_max(i) = buf_y_max(i) - nc

		  buf_x_min(j) = buf_x_min(j) - nc
		  buf_x_max(j) = buf_x_max(j) - nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_y_min(i) = buf_y_min(i) + nc
		  buf_y_max(i) = buf_y_max(i) + nc
		  buf_x_min(j) = buf_x_min(j) + nc
		  buf_x_max(j) = buf_x_max(j) + nc

		endif  
	      enddo
	    enddo  

	    do i = 1, chunk_x ! +x+z direction
	      do j = 1, chunk_z
		if(buf_x_min(i) .lt. 0 .and. buf_z_min(j) .lt. 0) then 

		  buf_x_min(i) = buf_x_min(i) + nc
		  buf_x_max(i) = buf_x_max(i) + nc

		  buf_z_min(j) = buf_z_min(j) + nc
		  buf_z_max(j) = buf_z_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_x_min(i) = buf_x_min(i) - nc
		  buf_x_max(i) = buf_x_max(i) - nc
		  buf_z_min(j) = buf_z_min(j) - nc
		  buf_z_max(j) = buf_z_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_x ! -x+z direction
	      do j = 1, chunk_z
		if(buf_x_max(i) .gt. nc .and. buf_z_min(j) .lt. 0) then 

		  buf_x_min(i) = buf_x_min(i) - nc
		  buf_x_max(i) = buf_x_max(i) - nc

		  buf_z_min(j) = buf_z_min(j) + nc
		  buf_z_max(j) = buf_z_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_x_min(i) = buf_x_min(i) + nc
		  buf_x_max(i) = buf_x_max(i) + nc
		  buf_z_min(j) = buf_z_min(j) - nc
		  buf_z_max(j) = buf_z_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_x ! +x-z direction
	      do j = 1, chunk_z
		if(buf_z_max(i) .gt. nc .and. buf_x_min(j) .lt. 0) then 

		  buf_z_min(i) = buf_z_min(i) - nc
		  buf_z_max(i) = buf_z_max(i) - nc

		  buf_x_min(j) = buf_x_min(j) + nc
		  buf_x_max(j) = buf_x_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_z_min(i) = buf_z_min(i) + nc
		  buf_z_max(i) = buf_z_max(i) + nc
		  buf_x_min(j) = buf_x_min(j) - nc
		  buf_x_max(j) = buf_x_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_x ! -x-z direction
	      do j = 1, chunk_z
		if(buf_z_max(i) .gt. nc .and. buf_x_max(j) .gt. nc) then 

		  buf_z_min(i) = buf_z_min(i) - nc
		  buf_z_max(i) = buf_z_max(i) - nc

		  buf_x_min(j) = buf_x_min(j) - nc
		  buf_x_max(j) = buf_x_max(j) - nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_z_min(i) = buf_z_min(i) + nc
		  buf_z_max(i) = buf_z_max(i) + nc
		  buf_x_min(j) = buf_x_min(j) + nc
		  buf_x_max(j) = buf_x_max(j) + nc

		endif  
	      enddo
	    enddo  

	    do i = 1, chunk_y ! +y+z direction
	      do j = 1, chunk_z
		if(buf_y_min(i) .lt. 0 .and. buf_z_min(j) .lt. 0) then  

		  buf_y_min(i) = buf_y_min(i) + nc
		  buf_y_max(i) = buf_y_max(i) + nc

		  buf_z_min(j) = buf_z_min(j) + nc
		  buf_z_max(j) = buf_z_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_y_min(i) = buf_y_min(i) - nc
		  buf_y_max(i) = buf_y_max(i) - nc
		  buf_z_min(j) = buf_z_min(j) - nc
		  buf_z_max(j) = buf_z_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_y ! -y+z direction
	      do j = 1, chunk_z
		if(buf_y_max(i) .gt. nc .and. buf_z_min(j) .lt. 0) then 

		  buf_y_min(i) = buf_y_min(i) - nc
		  buf_y_max(i) = buf_y_max(i) - nc

		  buf_z_min(j) = buf_z_min(j) + nc
		  buf_z_max(j) = buf_z_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_y_min(i) = buf_y_min(i) + nc
		  buf_y_max(i) = buf_y_max(i) + nc
		  buf_z_min(j) = buf_z_min(j) - nc
		  buf_z_max(j) = buf_z_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_y ! +y-z direction
	      do j = 1, chunk_z
		if(buf_z_max(i) .gt. nc .and. buf_y_min(j) .lt. 0) then 

		  buf_z_min(i) = buf_z_min(i) - nc
		  buf_z_max(i) = buf_z_max(i) - nc

		  buf_y_min(j) = buf_y_min(j) + nc
		  buf_y_max(j) = buf_y_max(j) + nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_z_min(i) = buf_z_min(i) + nc
		  buf_z_max(i) = buf_z_max(i) + nc
		  buf_y_min(j) = buf_y_min(j) - nc
		  buf_y_max(j) = buf_y_max(j) - nc

		endif  
	      enddo
	    enddo

	    do i = 1, chunk_y ! -y-z direction
	      do j = 1, chunk_z
		if(buf_z_max(i) .gt. nc .and. buf_y_max(j) .gt. nc) then 

		  buf_z_min(i) = buf_z_min(i) - nc
		  buf_z_max(i) = buf_z_max(i) - nc

		  buf_y_min(j) = buf_y_min(j) - nc
		  buf_y_max(j) = buf_y_max(j) - nc

		  call chunk_file_calc

		  ! put the bounds back

		  buf_z_min(i) = buf_z_min(i) + nc
		  buf_z_max(i) = buf_z_max(i) + nc
		  buf_y_min(j) = buf_y_min(j) + nc
		  buf_y_max(j) = buf_y_max(j) + nc

		endif  
	      enddo
	    enddo  

	    do i = 1, chunk_x ! +x+y+z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_min(i) .lt. 0 .and. buf_y_min(j) .lt. 0 .and. buf_z_min(k) .lt. 0) then 

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc

		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc

		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc
		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc
		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		  endif  
		enddo
	      enddo
	    enddo

	    do i = 1, chunk_x ! -x+y+z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_max(i) .gt. nc .and. buf_y_min(j) .lt. 0 .and. buf_z_min(k) .lt. 0) then 

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc

		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc

		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc
		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc
		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		  endif  
		enddo
	      enddo
	    enddo  

	    do i = 1, chunk_x ! +x-y+z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_min(i) .lt. 0 .and. buf_y_max(j) .gt. nc .and. buf_z_min(k) .lt. 0) then 

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc

		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc

		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc
		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc
		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		  endif  
		enddo
	      enddo
	    enddo    

	    do i = 1, chunk_x ! +x+y-z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_min(i) .lt. 0 .and. buf_y_min(j) .lt. 0 .and. buf_z_max(k) .gt. nc) then 

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc

		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc

		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc
		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc
		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		  endif  
		enddo
	      enddo
	    enddo  

	    do i = 1, chunk_x ! -x-y+z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_max(i) .gt. nc .and. buf_y_max(j) .gt. nc .and. buf_z_min(k) .lt. 0) then 

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc

		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc

		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc
		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc
		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		  endif  
		enddo
	      enddo
	    enddo  

	    do i = 1, chunk_x ! -x+y-z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_max(i) .gt. nc .and. buf_y_min(j) .lt. 0 .and. buf_z_max(k) .gt. nc) then 

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc

		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc

		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc
		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc
		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		  endif  
		enddo
	      enddo
	    enddo  

	    do i = 1, chunk_x ! +x-y-z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_min(i) .lt. 0 .and. buf_y_max(j) .gt. nc .and. buf_z_max(k) .gt. nc) then 

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc

		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc

		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc
		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc
		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		  endif  
		enddo
	      enddo
	    enddo    

	    do i = 1, chunk_x ! -x-y-z direction
	      do j = 1, chunk_y
		do k = 1, chunk_z
		  if(buf_x_max(i) .gt. nc .and. buf_y_max(j) .gt. nc .and. buf_z_max(k) .gt. nc) then 

		    buf_x_min(i) = buf_x_min(i) - nc
		    buf_x_max(i) = buf_x_max(i) - nc

		    buf_y_min(j) = buf_y_min(j) - nc
		    buf_y_max(j) = buf_y_max(j) - nc

		    buf_z_min(k) = buf_z_min(k) - nc
		    buf_z_max(k) = buf_z_max(k) - nc

		    call chunk_file_calc

		    ! put the bounds back

		    buf_x_min(i) = buf_x_min(i) + nc
		    buf_x_max(i) = buf_x_max(i) + nc
		    buf_y_min(j) = buf_y_min(j) + nc
		    buf_y_max(j) = buf_y_max(j) + nc
		    buf_z_min(k) = buf_z_min(k) + nc
		    buf_z_max(k) = buf_z_max(k) + nc

		  endif  
		enddo
	      enddo
	    enddo

#ifdef DEBUG
	    if(rank .eq. 1) then
	      print*, 'chunk_flag matrix:'
	      do k = 1, chunk_z
		do j = 1, chunk_y
		  print*, chunk_flag(:,j,k)
		enddo
		print*
	      enddo
	    endif
#endif

	    no_chunk_files_req = sum(chunk_flag)

#ifdef DEBUG  
	    print*, 'number of chunk files required for rank',rank,'=', no_chunk_files_req
#endif  

	end subroutine chunk_file_calc_buf

	!------------------------------------------------------------------------------------------!

	subroutine read_cubep3m_header

	    ! read header

	    read(10) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
		cur_projection,cur_halofind,mass_p

	end subroutine read_cubep3m_header

	!------------------------------------------------------------------------------------------!

	subroutine read_cubep3m_PID_header

	    ! read header

	    read(11) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint, &
		cur_projection,cur_halofind,mass_p

	end subroutine read_cubep3m_PID_header

	!------------------------------------------------------------------------------------------!


	subroutine total_particle_calculation

	    real*8 :: dummy_real

	    call mpi_allreduce(real(np_local,kind=8),dummy_real,1,mpi_real8,mpi_sum,mpi_comm_world,ierr)
	    np_total = int(dummy_real,8)

	    if(rank .eq. 0) then

	      print*
	      print*, 'Total particles in cubep3m simulation = ', np_total

	    endif

	    call mpi_barrier(mpi_comm_world,ierr)

	end subroutine total_particle_calculation

	!------------------------------------------------------------------------------------------!

	subroutine write_chunk_file_header_draft

	    integer(4) i, j, k, n, dummy
	    integer(8) dummy_long_int
	    real(8), dimension(3) :: position_offset, position_offset_low, position_offset_high
	    real(8), dimension(3) :: position_offset_low_buf, position_offset_high_buf
	    character(16) :: header_char = 'CHUNKED CUBEPM  '

	    dummy = 1

	    ! write the header information into chunk files

	    do k=1,chunk_z
	      do j=1,chunk_y
		do i=1,chunk_x
		  if(chunk_flag(i,j,k) == 1) then

		    chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

		    ! position offset to translate chunk coords back to global coords (for position of origin)

		    position_offset(1) = chunk_coords(chunk_file_int,1)*x_chunk_length - 2*buffer_size          
		    position_offset(2) = chunk_coords(chunk_file_int,2)*y_chunk_length - 2*buffer_size          
		    position_offset(3) = chunk_coords(chunk_file_int,3)*z_chunk_length - 2*buffer_size  

		    ! position of buffer zone bottom corner in chunk coords

		    position_offset_low_buf(1) =  buffer_size
		    position_offset_low_buf(2) =  buffer_size
		    position_offset_low_buf(3) =  buffer_size                    

		    ! position of buffer zone top corner in chunk coords

		    position_offset_high_buf(1) = x_chunk_length + 4*buffer_size
		    position_offset_high_buf(2) = y_chunk_length + 4*buffer_size
		    position_offset_high_buf(3) = z_chunk_length + 4*buffer_size          

		    ! position of particle zone bottom corner in chunk coords

		    position_offset_low(1) = 2*buffer_size
		    position_offset_low(2) = 2*buffer_size
		    position_offset_low(3) = 2*buffer_size          

		    ! position of particle zone top corner in chunk coords

		    position_offset_high(1) = x_chunk_length + 2*buffer_size
		    position_offset_high(2) = y_chunk_length + 2*buffer_size
		    position_offset_high(3) = z_chunk_length + 2*buffer_size


		    write(100+chunk_file_int) dummy                
		    write(100+chunk_file_int) header_char
		    write(100+chunk_file_int) dummy_long_int ! will eventually be total particle numbers in chunk
		    write(100+chunk_file_int) dummy_long_int ! will eventually be particle number in file
		    write(100+chunk_file_int) dummy ! will eventually be the number of files to read for this chunk
		    write(100+chunk_file_int) position_offset
		    write(100+chunk_file_int) position_offset_low_buf
		    write(100+chunk_file_int) position_offset_high_buf
		    write(100+chunk_file_int) position_offset_low
		    write(100+chunk_file_int) position_offset_high

		    write(100+chunk_file_int) np_local,a,t,tau,nts,dt_f_acc,dt_pp_acc,dt_c_acc,cur_checkpoint,&
			&cur_projection,cur_halofind,mass_p


		  endif


		enddo
	      enddo
	    enddo

	end subroutine write_chunk_file_header_draft

	!------------------------------------------------------------------------------------------!

	subroutine write_chunk_file_header_final

	    integer(4) :: i, j, k, num_files
	    integer(8) :: total_particles
	    integer(8) :: particles_in_this_file

	    ! write the header information into chunk files
	    ! the header is 1 * long integer for particles in file, then 1 * real for the buffer size
	    ! then 4 * integers and 7 * floats = 1*8+1*4+4*4+7*4=56 bytes

	    do k=1,chunk_z
	      do j=1,chunk_y
		do i=1,chunk_x
		  if(chunk_flag(i,j,k) == 1) then

		    chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

		    total_particles = tot_np_in_chunk(chunk_file_int)
		    particles_in_this_file = tot_np_for_each_chunk_by_rank(chunk_file_int)
		    num_files = sum(chunk_flag_global(:,i,j,k))


		    write(100+chunk_file_int, POS=21) total_particles
		    write(100+chunk_file_int) particles_in_this_file
		    write(100+chunk_file_int) num_files

		  endif

		enddo
	      enddo
	    enddo

	end subroutine write_chunk_file_header_final

	!------------------------------------------------------------------------------------------!

	subroutine process_particles(size)

	    integer i, j, k, p, size

	    ! place particle into full box coordinates and write them to their
	    ! respective chunk files

	    do p = 1, size

	      xv(1,p) = xv(1,p) + node_coords(rank,1)*ncc
	      xv(2,p) = xv(2,p) + node_coords(rank,2)*ncc
	      xv(3,p) = xv(3,p) + node_coords(rank,3)*ncc                  

	      if(PID_FLAG .eq. 0) call particle_assignation_no_ids(xv(:,p))
	      if(PID_FLAG .eq. 1) call particle_assignation_ids(xv(:,p),PID(p))

	    enddo

	    ! Now make sure that particles that are on the other side of the box are written into the buffers for the files 
	    ! on the other side of the box

	    do p = 1, size ! +x direction

	      if(xv(1,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p)
		xv(3,p) = xv(3,p)

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            
		! put the particle back

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)
		xv(3,p) = xv(3,p)        


	      endif
	    enddo

	    do p = 1, size ! +y direction

	      if(xv(2,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)         
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

	      endif
	    enddo

	    do p = 1, size ! +z direction

	      if(xv(3,p) .lt. buffer_size) then


		xv(1,p) = xv(1,p)         
		xv(2,p) = xv(2,p)
		xv(3,p) = xv(3,p) + nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)         
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) - nc

	      endif
	    enddo

	    do p = 1, size ! -x direction

	      if(xv(1,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p)
		xv(3,p) = xv(3,p)        


	      endif
	    enddo

	    do p = 1, size ! -y direction

	      if(xv(2,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p)
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)        

	      endif
	    enddo

	    do p = 1, size ! -z direction

	      if(xv(3,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) - nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) + nc

	      endif
	    enddo

	    do p = 1, size ! +x+y direction

	      if(xv(1,p) .lt. buffer_size .and. xv(2,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

	      endif
	    enddo

	    do p = 1, size ! +x+z direction

	      if(xv(1,p) .lt. buffer_size .and. xv(3,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) + nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) - nc

	      endif
	    enddo

	    do p = 1, size ! +y+z direction

	      if(xv(2,p) .lt. buffer_size .and. xv(3,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p)
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p) + nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p) - nc

	      endif
	    enddo

	    do p = 1, size ! -x+y direction

	      if(xv(1,p) + buffer_size .gt. nc .and. xv(2,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

	      endif
	    enddo

	    do p = 1, size ! +x-y direction

	      if(xv(2,p) + buffer_size .gt. nc .and. xv(1,p) .lt. buffer_size) then


		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back


		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)        

	      endif
	    enddo

	    do p = 1, size ! -x-y direction

	      if(xv(1,p) + buffer_size .gt. nc .and. xv(2,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p)        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) + nc        
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p)

	      endif
	    enddo

	    do p = 1, size ! -x+z direction

	      if(xv(1,p) + buffer_size .gt. nc .and. xv(3,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) + nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) - nc

	      endif
	    enddo        


	    do p = 1, size ! +x-z direction

	      if(xv(3,p) + buffer_size .gt. nc .and. xv(1,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p) + nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) - nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)        
		xv(3,p) = xv(3,p) + nc                 

	      endif
	    enddo

	    do p = 1, size ! -x-z direction

	      if(xv(1,p) + buffer_size .gt. nc .and. xv(3,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p) - nc
		xv(2,p) = xv(2,p)         
		xv(3,p) = xv(3,p) - nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p) + nc 
		xv(2,p) = xv(2,p)                
		xv(3,p) = xv(3,p) + nc

	      endif
	    enddo

	    do p = 1, size ! -y+z direction

	      if(xv(2,p) + buffer_size .gt. nc .and. xv(3,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p) + nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)       
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p) - nc

	      endif
	    enddo          

	    do p = 1, size ! +y-z direction

	      if(xv(3,p) + buffer_size .gt. nc .and. xv(2,p) .lt. buffer_size) then

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) + nc
		xv(3,p) = xv(3,p) - nc        

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p) + nc


	      endif
	    enddo


	    do p = 1, size ! -y-z direction

	      if(xv(2,p) + buffer_size .gt. nc .and. xv(3,p) + buffer_size .gt. nc) then

		xv(1,p) = xv(1,p)        
		xv(2,p) = xv(2,p) - nc
		xv(3,p) = xv(3,p) - nc

		if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
		if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

		! put the particle back

		xv(1,p) = xv(1,p)
		xv(2,p) = xv(2,p) + nc        
		xv(3,p) = xv(3,p) + nc

	      endif
	    enddo


	    do p = 1, size ! +x+y+z direction

	      if(xv(1,p) .lt. buffer_size .and. xv(2,p) .lt. buffer_size&
		  & .and. xv(3,p) .lt. buffer_size) then

	      xv(1,p) = xv(1,p) + nc
	      xv(2,p) = xv(2,p) + nc
	      xv(3,p) = xv(3,p) + nc  

	      if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
	      if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

	      ! put the particle back

	      xv(1,p) = xv(1,p) - nc
	      xv(2,p) = xv(2,p) - nc
	      xv(3,p) = xv(3,p) - nc  

	    endif
	  enddo


	  do p = 1, size ! -x+y+z direction

	    if(xv(1,p) + buffer_size .gt. nc .and. xv(2,p) .lt. buffer_size&
		& .and. xv(3,p) .lt. buffer_size) then

	    xv(1,p) = xv(1,p) - nc
	    xv(2,p) = xv(2,p) + nc
	    xv(3,p) = xv(3,p) + nc  

	    if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
	    if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

	    ! put the particle back

	    xv(1,p) = xv(1,p) + nc
	    xv(2,p) = xv(2,p) - nc
	    xv(3,p) = xv(3,p) - nc  

	  endif
	enddo

	do p = 1, size ! +x-y+z direction

	  if(xv(2,p) + buffer_size .gt. nc .and. xv(1,p) .lt. buffer_size&
	      & .and. xv(3,p) .lt. buffer_size) then

	  xv(1,p) = xv(1,p) + nc
	  xv(2,p) = xv(2,p) - nc
	  xv(3,p) = xv(3,p) + nc  

	  if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
	  if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

	  ! put the particle back


	  xv(1,p) = xv(1,p) - nc
	  xv(2,p) = xv(2,p) + nc
	  xv(3,p) = xv(3,p) - nc  

	endif
      enddo

      do p = 1, size ! +x+y-z direction

	if(xv(3,p) + buffer_size .gt. nc .and. xv(2,p) .lt. buffer_size&
	    & .and. xv(1,p) .lt. buffer_size) then

	xv(1,p) = xv(1,p) + nc         
	xv(2,p) = xv(2,p) + nc
	xv(3,p) = xv(3,p) - nc

	if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
	if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

	! put the particle back

	xv(1,p) = xv(1,p) - nc  
	xv(2,p) = xv(2,p) - nc
	xv(3,p) = xv(3,p) + nc

      endif
    enddo

    do p = 1, size ! -x-y+z direction

      if(xv(1,p) + buffer_size .gt. nc .and. xv(2,p)+buffer_size .gt. nc&
	  & .and. xv(3,p) .lt. buffer_size) then

      xv(1,p) = xv(1,p) - nc
      xv(2,p) = xv(2,p) - nc
      xv(3,p) = xv(3,p) + nc  

      if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
      if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

      ! put the particle back

      xv(1,p) = xv(1,p) + nc
      xv(2,p) = xv(2,p) + nc
      xv(3,p) = xv(3,p) - nc  

    endif
  enddo

  do p = 1, size ! -x+y-z direction

    if(xv(1,p) + buffer_size .gt. nc .and. xv(3,p)+buffer_size .gt. nc&
	& .and. xv(2,p) .lt. buffer_size) then

    xv(1,p) = xv(1,p) - nc
    xv(2,p) = xv(2,p) + nc  
    xv(3,p) = xv(3,p) - nc        

    if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
    if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

    ! put the particle back

    xv(1,p) = xv(1,p) + nc
    xv(2,p) = xv(2,p) - nc  
    xv(3,p) = xv(3,p) + nc

  endif
enddo

do p = 1, size ! +x-y-z direction

  if(xv(3,p) + buffer_size .gt. nc .and. xv(2,p)+buffer_size .gt. nc&
      & .and. xv(1,p) .lt. buffer_size) then

  xv(1,p) = xv(1,p) + nc  
  xv(2,p) = xv(2,p) - nc
  xv(3,p) = xv(3,p) - nc

  if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
  if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

  ! put the particle back


  xv(1,p) = xv(1,p) - nc  
  xv(2,p) = xv(2,p) + nc        
  xv(3,p) = xv(3,p) + nc

endif
	enddo

	do p = 1, size ! -x-y-z direction

	  if(xv(1,p) + buffer_size .gt. nc .and. xv(2,p)+buffer_size .gt. nc&
	      & .and. xv(3,p) + buffer_size .gt. nc) then

	  xv(1,p) = xv(1,p) - nc
	  xv(2,p) = xv(2,p) - nc
	  xv(3,p) = xv(3,p) - nc  

	  if(PID_FLAG .eq. 0) call buffer_particle_assignation_no_ids(xv(:,p))
	  if(PID_FLAG .eq. 1) call buffer_particle_assignation_ids(xv(:,p),PID(p))            

	  ! put the particle back

	  xv(1,p) = xv(1,p) + nc
	  xv(2,p) = xv(2,p) + nc
	  xv(3,p) = xv(3,p) + nc  

	endif
      enddo

  end subroutine process_particles

  !------------------------------------------------------------------------------------------!

  subroutine particle_assignation_no_ids(xv)

      integer i, j, k, chunk_file_int
      real(4), dimension(6) :: xv

      do k = 1, chunk_z
	do j = 1, chunk_y
	  do i = 1, chunk_x

	    if(xv(1) .ge. buf_x_min(i) .and. xv(1) .lt. buf_x_max(i) .and. &
		& xv(2) .ge. buf_y_min(j) .and. xv(2) .lt. buf_y_max(j) .and. &
		& xv(3) .ge. buf_z_min(k) .and. xv(3) .lt. buf_z_max(k)) then

	    chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	    if(xv(1) .ge. x_min(i) .and. xv(1) .lt. x_max(i) .and. & ! count the particles into either 
		& xv(2) .ge. y_min(j) .and. xv(2) .lt. y_max(j) .and. &  ! the buffer or internal volume
		& xv(3) .ge. z_min(k) .and. xv(3) .lt. z_max(k)) then

	    chunk_particles(chunk_file_int,i,j,k) = chunk_particles(chunk_file_int,i,j,k) + 1

	  else

	    chunk_buffer_particles(chunk_file_int,i,j,k) = chunk_buffer_particles(chunk_file_int,i,j,k) + 1

	  endif  

	  ! Put the particle into the required coordinates for the chunk file:

	  xv(1) = xv(1) - chunk_coords(chunk_file_int,1)*x_chunk_length + 2*buffer_size
	  xv(2) = xv(2) - chunk_coords(chunk_file_int,2)*y_chunk_length + 2*buffer_size
	  xv(3) = xv(3) - chunk_coords(chunk_file_int,3)*z_chunk_length + 2*buffer_size

	  ! write the particle into its relevant chunk file

	  write(100+chunk_file_int) xv

	  if(xv(1) .gt. max_x(chunk_file_int)) max_x(chunk_file_int) = xv(1)
	  if(xv(2) .gt. max_y(chunk_file_int)) max_y(chunk_file_int) = xv(2)
	  if(xv(3) .gt. max_x(chunk_file_int)) max_z(chunk_file_int) = xv(3)

	  if(xv(1) .lt. min_x(chunk_file_int)) min_x(chunk_file_int) = xv(1)
	  if(xv(2) .lt. min_y(chunk_file_int)) min_y(chunk_file_int) = xv(2)
	  if(xv(3) .lt. min_x(chunk_file_int)) min_z(chunk_file_int) = xv(3)

	  ! Put the particle back into global coordinates

	  xv(1) = xv(1) + chunk_coords(chunk_file_int,1)*x_chunk_length - 2*buffer_size
	  xv(2) = xv(2) + chunk_coords(chunk_file_int,2)*y_chunk_length - 2*buffer_size
	  xv(3) = xv(3) + chunk_coords(chunk_file_int,3)*z_chunk_length - 2*buffer_size

	endif

      enddo
    enddo
  enddo

end subroutine particle_assignation_no_ids

!------------------------------------------------------------------------------------------!

subroutine buffer_particle_assignation_no_ids(xv)

    integer i, j, k, chunk_file_int
    real(4), dimension(6) :: xv

    do k = 1, chunk_z
      do j = 1, chunk_y
	do i = 1, chunk_x

	  if(xv(1) .gt. buf_x_min(i) .and. xv(1) .le. buf_x_max(i) .and. &
	      & xv(2) .gt. buf_y_min(j) .and. xv(2) .le. buf_y_max(j) .and. &
	      & xv(3) .gt. buf_z_min(k) .and. xv(3) .le. buf_z_max(k)) then

	  chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	  ! Put the particle into the required coordinates for the chunk file:

	  xv(1) = xv(1) - chunk_coords(chunk_file_int,1)*x_chunk_length + 2*buffer_size
	  xv(2) = xv(2) - chunk_coords(chunk_file_int,2)*y_chunk_length + 2*buffer_size
	  xv(3) = xv(3) - chunk_coords(chunk_file_int,3)*z_chunk_length + 2*buffer_size

	  ! write the particle into its relevant chunk file

	  write(100+chunk_file_int) xv

	  if(xv(1) .gt. max_x(chunk_file_int)) max_x(chunk_file_int) = xv(1)
	  if(xv(2) .gt. max_y(chunk_file_int)) max_y(chunk_file_int) = xv(2)
	  if(xv(3) .gt. max_x(chunk_file_int)) max_z(chunk_file_int) = xv(3)

	  if(xv(1) .lt. min_x(chunk_file_int)) min_x(chunk_file_int) = xv(1)
	  if(xv(2) .lt. min_y(chunk_file_int)) min_y(chunk_file_int) = xv(2)
	  if(xv(3) .lt. min_x(chunk_file_int)) min_z(chunk_file_int) = xv(3)


	  ! Put the particle back into global coordinates

	  xv(1) = xv(1) + chunk_coords(chunk_file_int,1)*x_chunk_length - 2*buffer_size
	  xv(2) = xv(2) + chunk_coords(chunk_file_int,2)*y_chunk_length - 2*buffer_size
	  xv(3) = xv(3) + chunk_coords(chunk_file_int,3)*z_chunk_length - 2*buffer_size

	  chunk_buffer_particles(chunk_file_int,i,j,k) = chunk_buffer_particles(chunk_file_int,i,j,k) + 1

	endif

      enddo
    enddo
  enddo

end subroutine buffer_particle_assignation_no_ids


!------------------------------------------------------------------------------------------!

subroutine particle_assignation_ids(xv,PID)

    integer i, j, k, chunk_file_int
    real(4), dimension(6) :: xv
    integer(8) :: PID

    do k = 1, chunk_z
      do j = 1, chunk_y
	do i = 1, chunk_x

	  if(xv(1) .ge. buf_x_min(i) .and. xv(1) .lt. buf_x_max(i) .and. &
	      & xv(2) .ge. buf_y_min(j) .and. xv(2) .lt. buf_y_max(j) .and. &
	      & xv(3) .ge. buf_z_min(k) .and. xv(3) .lt. buf_z_max(k)) then

	  chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	  if(xv(1) .ge. x_min(i) .and. xv(1) .lt. x_max(i) .and. & ! count the particles into either 
	      & xv(2) .ge. y_min(j) .and. xv(2) .lt. y_max(j) .and. &  ! the buffer or internal volume
	      & xv(3) .ge. z_min(k) .and. xv(3) .lt. z_max(k)) then

	  chunk_particles(chunk_file_int,i,j,k) = chunk_particles(chunk_file_int,i,j,k) + 1

	else

	  chunk_buffer_particles(chunk_file_int,i,j,k) = chunk_buffer_particles(chunk_file_int,i,j,k) + 1

	endif  

	! Put the particle into the required coordinates for the chunk file:

	xv(1) = xv(1) - chunk_coords(chunk_file_int,1)*x_chunk_length + 2*buffer_size
	xv(2) = xv(2) - chunk_coords(chunk_file_int,2)*y_chunk_length + 2*buffer_size
	xv(3) = xv(3) - chunk_coords(chunk_file_int,3)*z_chunk_length + 2*buffer_size

	! write the particle into its relevant chunk file

	write(100+chunk_file_int) xv,PID

	if(xv(1) .gt. max_x(chunk_file_int)) max_x(chunk_file_int) = xv(1)
	if(xv(2) .gt. max_y(chunk_file_int)) max_y(chunk_file_int) = xv(2)
	if(xv(3) .gt. max_x(chunk_file_int)) max_z(chunk_file_int) = xv(3)

	if(xv(1) .lt. min_x(chunk_file_int)) min_x(chunk_file_int) = xv(1)
	if(xv(2) .lt. min_y(chunk_file_int)) min_y(chunk_file_int) = xv(2)
	if(xv(3) .lt. min_x(chunk_file_int)) min_z(chunk_file_int) = xv(3)

	! Put the particle back into global coordinates

	xv(1) = xv(1) + chunk_coords(chunk_file_int,1)*x_chunk_length - 2*buffer_size
	xv(2) = xv(2) + chunk_coords(chunk_file_int,2)*y_chunk_length - 2*buffer_size
	xv(3) = xv(3) + chunk_coords(chunk_file_int,3)*z_chunk_length - 2*buffer_size

      endif

    enddo
  enddo
enddo

end subroutine particle_assignation_ids

!------------------------------------------------------------------------------------------!

subroutine buffer_particle_assignation_ids(xv,PID)

    integer i, j, k, chunk_file_int
    real(4), dimension(6) :: xv
    integer(8) :: PID

    do k = 1, chunk_z
      do j = 1, chunk_y
	do i = 1, chunk_x

	  if(xv(1) .gt. buf_x_min(i) .and. xv(1) .le. buf_x_max(i) .and. &
	      & xv(2) .gt. buf_y_min(j) .and. xv(2) .le. buf_y_max(j) .and. &
	      & xv(3) .gt. buf_z_min(k) .and. xv(3) .le. buf_z_max(k)) then

	  chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	  ! Put the particle into the required coordinates for the chunk file:

	  xv(1) = xv(1) - chunk_coords(chunk_file_int,1)*x_chunk_length + 2*buffer_size
	  xv(2) = xv(2) - chunk_coords(chunk_file_int,2)*y_chunk_length + 2*buffer_size
	  xv(3) = xv(3) - chunk_coords(chunk_file_int,3)*z_chunk_length + 2*buffer_size

	  ! write the particle into its relevant chunk file

	  write(100+chunk_file_int) xv,PID

	  if(xv(1) .gt. max_x(chunk_file_int)) max_x(chunk_file_int) = xv(1)
	  if(xv(2) .gt. max_y(chunk_file_int)) max_y(chunk_file_int) = xv(2)
	  if(xv(3) .gt. max_x(chunk_file_int)) max_z(chunk_file_int) = xv(3)

	  if(xv(1) .lt. min_x(chunk_file_int)) min_x(chunk_file_int) = xv(1)
	  if(xv(2) .lt. min_y(chunk_file_int)) min_y(chunk_file_int) = xv(2)
	  if(xv(3) .lt. min_x(chunk_file_int)) min_z(chunk_file_int) = xv(3)


	  ! Put the particle back into global coordinates

	  xv(1) = xv(1) + chunk_coords(chunk_file_int,1)*x_chunk_length - 2*buffer_size
	  xv(2) = xv(2) + chunk_coords(chunk_file_int,2)*y_chunk_length - 2*buffer_size
	  xv(3) = xv(3) + chunk_coords(chunk_file_int,3)*z_chunk_length - 2*buffer_size

	  chunk_buffer_particles(chunk_file_int,i,j,k) = chunk_buffer_particles(chunk_file_int,i,j,k) + 1

	endif

      enddo
    enddo
  enddo

end subroutine buffer_particle_assignation_ids

!------------------------------------------------------------------------------------------!

subroutine open_cubep3m_file

    write(file,'(i4)') rank
    file=adjustl(file)

    ifile=particle_data_path(1:len_trim(particle_data_path))//'/'//&
	&redshift(1:len_trim(redshift))//'xv'//file(1:len_trim(file))//'.dat'	

#ifdef GFORTRAN		           
    open (unit=10,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=10,file=ifile,form='binary')
#endif        

    print*, 'successfully opened:',ifile(1:len_trim(ifile)), ' on rank', rank

end subroutine open_cubep3m_file

!------------------------------------------------------------------------------------------!

subroutine close_cubep3m_file

    close(10)

end subroutine close_cubep3m_file 

!------------------------------------------------------------------------------------------!


subroutine open_chunk_files

    integer i,j,k
    character(4) dummy_c

    file_number = 0

    call create_file_numbers

    ! open the chunk files:

    do k=1,chunk_z
      do j=1,chunk_y
	do i=1,chunk_x
	  if(chunk_flag(i,j,k) .eq. 1) then

	    chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	    write(file,'(i4)') file_number(chunk_file_int)
	    file=adjustl(file)

	    write(dummy_c,'(i4)') chunk_file_int
	    dummy_c=adjustl(dummy_c)


	    ifile=chunk_output_path(1:len_trim(chunk_output_path))//'/'//&
		&'z_'//redshift(1:len_trim(redshift))&
		&//'/chunk_'//dummy_c(1:len_trim(dummy_c))//&
		&'/'//redshift(1:len_trim(redshift))//&
		&'xv_chunk_'//dummy_c(1:len_trim(dummy_c))//&
		&'_'//file(1:len_trim(file))//'.dat'


#ifdef GFORTRAN		           
	    open (unit=100+chunk_file_int,file=ifile,access='stream')
#endif
#ifdef BINARY
	    open (unit=100+chunk_file_int,file=ifile,form='binary',access='stream')
#endif              

	  endif
	enddo
      enddo
    enddo		

    call mpi_barrier(mpi_comm_world,ierr)

end subroutine open_chunk_files

!------------------------------------------------------------------------------------------!

subroutine create_file_numbers

    ! ensures that each chunk is produced with a file number that is sequential

    integer i, j, k, n, l

    chunk_flag_local = 0
    chunk_flag_global = 0

    do k=1,chunk_z
      do j=1,chunk_y
	do i=1,chunk_x

	  if(chunk_flag(i,j,k) .eq. 1) then

	    chunk_flag_local(rank,i,j,k) = 1

	  endif

	enddo
      enddo
    enddo

    do n = 0,nodes-1
      do k=1,chunk_z
	do j=1,chunk_y
	  do i=1,chunk_x

	    call mpi_allreduce(chunk_flag_local(n,i,j,k),chunk_flag_global(n,i,j,k),1,mpi_integer,&
		&mpi_max, mpi_comm_world,ierr)

	  enddo
	enddo
      enddo
    enddo

    do k=1,chunk_z
      do j=1,chunk_y
	do i=1,chunk_x
	  if(chunk_flag_global(rank,i,j,k) .eq. 1) then

	    n = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y
	    file_number(n) = file_number(n) + sum(chunk_flag_global(0:rank-1,i,j,k))
	  endif    
	enddo     
      enddo                     
    enddo        

end subroutine create_file_numbers

!------------------------------------------------------------------------------------------!

subroutine close_chunk_files

    integer i, j, k

    do k=1,chunk_z
      do j=1,chunk_y
	do i=1,chunk_x
	  if(chunk_flag(i,j,k) == 1) then

	    chunk_file_int = (i-1)+(j-1)*chunk_x+(k-1)*chunk_x*chunk_y

	    close(100+chunk_file_int)

	  endif

	enddo
      enddo
    enddo

end subroutine close_chunk_files  

!------------------------------------------------------------------------------------------!

subroutine particle_number_check

    integer n, i, j, k
    integer(8) dummy_long_int

    do n = 0, no_of_chunks-1  
      do k = 1, chunk_z
	do j = 1, chunk_y
	  do i = 1, chunk_x

	    chunk_tot_particles(n,i,j,k) = chunk_buffer_particles(n,i,j,k) + chunk_particles(n,i,j,k)

	  enddo        
	enddo  
      enddo 

    enddo

    do n = 0, no_of_chunks-1

      do k = 1, chunk_z
	do j = 1, chunk_y
	  do i = 1, chunk_x

	    np_local_check = chunk_particles(n,i,j,k)+np_local_check

	    np_local_buffer = chunk_buffer_particles(n,i,j,k)+np_local_buffer

	    np_local_total = chunk_tot_particles(n,i,j,k)+np_local_total

	    tot_np_for_each_chunk_by_rank(n) = chunk_tot_particles(n,i,j,k)+tot_np_for_each_chunk_by_rank(n)
	    buf_np_for_each_chunk_by_rank(n) = chunk_buffer_particles(n,i,j,k)+buf_np_for_each_chunk_by_rank(n)
	    check_np_for_each_chunk_by_rank(n) = chunk_particles(n,i,j,k)+check_np_for_each_chunk_by_rank(n)


	  enddo        
	enddo  
      enddo  

    enddo

    do n = 0, no_of_chunks-1

      dummy_long_int = 0

      call mpi_allreduce(tot_np_for_each_chunk_by_rank(n),dummy_long_int,1,mpi_integer8, &
	  & mpi_sum,mpi_comm_world,ierr)

      tot_np_in_chunk(n) = tot_np_in_chunk(n) + dummy_long_int

      dummy_long_int = 0

      call mpi_allreduce(buf_np_for_each_chunk_by_rank(n),dummy_long_int,1,mpi_integer8, &
	  & mpi_sum,mpi_comm_world,ierr)

      buf_np_in_chunk(n) = buf_np_in_chunk(n) + dummy_long_int    

      dummy_long_int = 0

      call mpi_allreduce(check_np_for_each_chunk_by_rank(n),dummy_long_int,1,mpi_integer8, &
	  & mpi_sum,mpi_comm_world,ierr)

      check_np_in_chunk(n) = check_np_in_chunk(n) + dummy_long_int

    enddo

    call mpi_barrier(mpi_comm_world,ierr)

    if(rank .eq. 0) then 

      print*
      print*,'Particles on nodes being written into chunk files:'
      print*,'Rank                      Total Particles     Buffer Particles       &
	  &Non-Buffer Particles'

    endif

    call mpi_barrier(mpi_comm_world,ierr)

    do n = 0, nodes-1

      if(rank .eq. n) print*, rank, np_local_total, np_local_buffer, np_local_check

      call mpi_barrier(mpi_comm_world, ierr)

    enddo



    call mpi_allreduce(np_local_check,np_check,1,mpi_integer8, &
	mpi_sum,mpi_comm_world,ierr)


    if(rank .eq. 0 ) then

      print*
      print*, 'Particle number check:'
      print*, 'np_total = ',np_total
      print*, 'np_check = ',np_check

      if(np_total .eq. np_check) then
	print*
	print*, '     - Passed'    
      else
	print*
	print*, 'WARNING - PARTICLE MISMATCH! This is only a problem if'    
	print*, 'the number of particles is significantly off. Otherwise'
	print*, 'all that has occurred is that some particles that are on' 
	print*, 'the border of a cubep3m boundary have been counted as'
	print*, 'belonging to the buffer rather than the original volume'
      endif

      print*
      print*, 'Particles in chunks:        Total              Buffer             Non-Buffer'

      do n = 0, no_of_chunks-1

	print*, n, '=', tot_np_in_chunk(n), buf_np_in_chunk(n), check_np_in_chunk(n)

      enddo

    endif


end subroutine particle_number_check

!------------------------------------------------------------------------------------------!

subroutine particle_bound_check

    integer :: i

    do i = 0, no_of_chunks-1

      call mpi_allreduce(max_x(i),max_x_global(i),1,mpi_real,&
	  &mpi_max, mpi_comm_world,ierr)
      call mpi_allreduce(max_y(i),max_y_global(i),1,mpi_real,&
	  &mpi_max, mpi_comm_world,ierr)
      call mpi_allreduce(max_z(i),max_z_global(i),1,mpi_real,&
	  &mpi_max, mpi_comm_world,ierr)
      call mpi_allreduce(min_x(i),min_x_global(i),1,mpi_real,&
	  &mpi_min, mpi_comm_world,ierr)
      call mpi_allreduce(min_y(i),min_y_global(i),1,mpi_real,&
	  &mpi_min, mpi_comm_world,ierr)
      call mpi_allreduce(min_z(i),min_z_global(i),1,mpi_real,&
	  &mpi_min, mpi_comm_world,ierr)

    enddo


    if(rank .eq. 0) then

      print*  
      print*, 'particle boundary check:'

      do i = 0, no_of_chunks-1

	print*,'max x, y, z, min x, y, z, on rank',i,' are:',&
	    &max_x_global(i), max_y_global(i), max_z_global(i), &
	    &min_x_global(i), min_y_global(i), min_z_global(i)


      enddo


    endif


end subroutine

!------------------------------------------------------------------------------------------!
subroutine finalise_chunker

    if(rank .eq. 0) then

      print*, 'Chunking complete.' 
      print*
      print*,'"Share and Enjoy" -TM the Sirius Cybernetics Corp'

    endif


    call mpi_finalize(ierr)

end subroutine finalise_chunker

!------------------------------------------------------------------------------------------!


subroutine open_cubep3m_PID_file

    write(file,'(i4)') rank
    file=adjustl(file)

    ifile=particle_data_path(1:len_trim(particle_data_path))//'/'//&
	&redshift(1:len_trim(redshift))//'PID'//file(1:len_trim(file))//'.dat'	

#ifdef GFORTRAN		           
    open (unit=11,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=11,file=ifile,form='binary')
#endif        

    print*, 'successfully opened:',ifile(1:len_trim(ifile)), ' on rank', rank

end subroutine open_cubep3m_PID_file

!------------------------------------------------------------------------------------------!

subroutine close_cubep3m_PID_file

    close(11)

end subroutine close_cubep3m_PID_file 


!------------------------------------------------------------------------------------------------------------------------ 

subroutine read_parameters

    character(300)  :: line1
    integer(4)	  :: Num_args, IO, param_count
    character(400)  :: arg1, parameterfile

    ! input parameter variables:

    real(4)	  :: redshift_param
    character(100)  :: redshift_file_param
    character(300)  :: particle_data_path_param, chunk_output_path_param
    real(4)	  :: boxlength_param, nc_param, buffer_size_param
    integer(4)	  :: nodes_per_dimension_param, pid_flag_param
    integer(4)	  :: chunk_x_param, chunk_y_param, chunk_z_param



    Num_args = IARGC()

    if (Num_args .ne. 1) then
      write(*,'(a)') "COMMAND LINE ERROR! : no parameterfile provided&
      &... include paramterfile name after executable."
      call mpi_abort(mpi_comm_world,ierr,ierr)

    else

      CALL GETARG(1,arg1)
      read(arg1,*) parameterfile
      write(*,*) "Reading parameter file: ",&
	  &parameterfile(1:len_trim(parameterfile))

    endif

    open(unit=82,file=parameterfile(1:len_trim(parameterfile)),form='formatted',status='old')

    param_count = 1

    do

      read(82,'(a)',iostat = IO) line1

      if ( IO < 0 ) exit

      if(line1(1:1) .ne. '#') then

	if(param_count .eq. 1) read(line1,*) redshift_param
	if(param_count .eq. 2) read(line1,'(a100)') redshift_file_param
	if(param_count .eq. 3) read(line1,'(a300)') particle_data_path_param
	if(param_count .eq. 4) read(line1,'(a300)') chunk_output_path_param
	if(param_count .eq. 5) read(line1,*) boxlength_param
	if(param_count .eq. 6) read(line1,*) nodes_per_dimension_param
	if(param_count .eq. 7) read(line1,*) nc_param
	if(param_count .eq. 8) read(line1,*) pid_flag_param
	if(param_count .eq. 9) read(line1,*) buffer_size_param
	if(param_count .eq. 10) read(line1,*) chunk_x_param
	if(param_count .eq. 11) read(line1,*) chunk_y_param
	if(param_count .eq. 12) read(line1,*) chunk_z_param


	param_count = param_count + 1

      endif

    enddo

    if(redshift_param .eq. -1.0) then

      redshift_file_param = adjustl(redshift_file_param)


      print*,"REDSHIFT SET TO -1... Reading redshift from file..."
      print*,"File to read = ",redshift_file_param(1:len_trim(redshift_file_param))

      open(unit=83,file=redshift_file_param,form='formatted',status='old',iostat=IO)

      if ( IO < 0 ) then

	print*,"Error - redshift file specified (redshift in parameters set to -1)"
	print*,"but the file could not be opened! Filename: "
	print*, redshift_file_param(1:len_trim(redshift_file_param))

	call exit(0)

      endif

      read(83,*) redshift_param

      print*,"Redshift read, value = ",redshift_param

      close(83)

    endif  

    print*
    print*, "finished reading ",parameterfile(1:len_trim(parameterfile))
    print*
    print*, "Parameters are as follows:"
    print*
    if(redshift_param .eq. -1.0) then
      print*,"Redshift specified from file: ",redshift_file_param(1:len_trim(redshift_file_param))
    endif
    print*, "Redshift: ",redshift_param
    print*, "Particle data path = ",particle_data_path_param(1:len_trim(particle_data_path_param))
    print*, "Chunk output path = ",chunk_output_path_param(1:len_trim(chunk_output_path_param))
    print*, "Box length (Mpc/h) = ", boxlength_param
    print*, "Cubep3m nodes per dimension = ", nodes_per_dimension_param
    print*, "Cubep3m fine cells per dimension = ", nc_param
    print*, "PID flag (0 for off, 1 for on) = ", pid_flag_param
    print*, "Buffer size (Mpc/h) = ", buffer_size_param
    print*, "Chunks in x dimension = ", chunk_x_param
    print*, "Chunks in y dimension = ", chunk_y_param
    print*, "Chunks in z dimension = ", chunk_z_param

    ! initialise parameters on task 1:

    z = redshift_param
    particle_data_path = adjustl(particle_data_path_param) 
    chunk_output_path = adjustl(chunk_output_path_param)
    box = boxlength_param
    nodes_dim = nodes_per_dimension_param
    nc = nc_param
    PID_FLAG = pid_flag_param
    buffer_size = buffer_size_param
    chunk_x = chunk_x_param
    chunk_y = chunk_y_param
    chunk_z = chunk_z_param


end subroutine
!------------------------------------------------------------------------------------------------------------------------ 

subroutine broadcast_parameters

    call MPI_Bcast(z, 1, MPI_REAL, 0, mpi_comm_world,ierr)
    call MPI_Bcast(nc, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(particle_data_path, 400, MPI_CHARACTER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(chunk_output_path, 400, MPI_CHARACTER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(box, 1, MPI_REAL, 0, mpi_comm_world,ierr)
    call MPI_Bcast(nodes_dim, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(nc, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(PID_FLAG, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(buffer_size, 1, MPI_REAL, 0, mpi_comm_world,ierr)
    call MPI_Bcast(chunk_x, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(chunk_y, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
    call MPI_Bcast(chunk_z, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)

end subroutine


end program chunk_cubep3m
