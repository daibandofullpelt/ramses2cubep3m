
! Converts Gadget ICs to Cubep3m ICs. THIS PROGRAM NEEDS TO BE RUN ON THE SAME NUMBER OF MPI TASKS AS THE CUBEP3M
! SIMULATION. Tasks read the gadget files (the number of gadget files does not need to be the same as the number of 
! Cubep3m nodes) and then distributes the particles to the appropriate Cubep3m task. The cubep3m tasks then write the 
! particles they recieve. The size of the particle batch that is read and sent is set by the 'read_block_size' parameter.
! If this is set low then there will be more reads and sends but the code will have a much lighter memory footprint.
!
! To setup the run edit the Cubep3m file details and Gadget file details sections below to match the simulation requirements
!
! email: ww60@sussex.ac.uk
!
! compile with:
! (willmb): mpif90 -cpp GADGET2CUBEP3M_v0.11.f90 -o G2CPM -DGFORTRAN -fdefault-integer-8
! (superMUC): mpif90 -cpp GADGET2CUBEP3M_v0.11.f90 -o G2CPM -DBINARY


program GADGET2CUBEP3M

  implicit none

  include 'mpif.h'

  logical :: debug = .false.

  integer*4, parameter :: read_block_size = 50000*24    ! tunable parameter -- how large a block size of data
  ! to read on each node. MAKE SURE IT STRIDES AN INTEGER
  ! NUMBER OF PARTICLES (i.e. that read_length is an integer, so
  ! keep the factor of 24 in the sum and it should be fine) 
  ! Currently there are two arrays
  ! that are 6 x read_length x number of nodes. As the nodes increase in size
  ! this will become an issue

  ! At the moment this is also the size of the MPI buffer for 
  ! sending particles so set with care.
  ! 
  ! Note that for reading the Intel compiler may require this 
  ! to be less than 24Mb depending on system architecture

  integer*4, parameter :: read_length = read_block_size / 24      ! i.e. divided by 24 bytes for 3 positions and velocities



  ! Gadget file details:

  character(*), parameter :: gadget_path ='/Users/ww60/Documents/work/code/LgenIC/ICs/'
  character(*), parameter :: gadget_root ='ic'  
  integer*4, parameter :: num_files_gadget = 32
  integer*4, parameter :: gadget_length_flag = 0 ! Flag for the units of the gadget boxlength. 0 = Mpc/h, 1 = kpc/h



  ! Cubep3m file details:

  integer*4, parameter :: nodes_dim = 2 ! number of cubep3m nodes per dimension
  character(*), parameter :: cpm_path ='/Users/ww60/Documents/work/code_output/cubep3m/gadget_IC_converter_test/'


  ! Cubep3m Variables

  integer*4 :: nc
  real*4 :: ncc
  integer*4, parameter :: nodes = (nodes_dim*nodes_dim*nodes_dim)
  integer*4 :: node_coords(0:nodes-1,3)


  ! Gadget Header Block:

  integer*4 :: header_length_ghead      ! header begins with an integer sizeof(header)
  integer*4 :: npart_ghead(6)           ! number of particles of each type in this file
  real*8 :: mass_ghead(6)               ! mass of particles of each type
  real*8 :: time_ghead                  ! expansion factor
  real*8 :: redshift_ghead              ! redshift of snapshot file
  integer*4 :: flag_sfr_ghead           ! flags whether the simulation was including star formation 
  integer*4 :: flag_feedback_ghead      ! flags whether feedback was included (obsolete) 
  integer*4 :: npartTotal_ghead(6)      ! total number of particles of each type in this snapshot. 
  integer*4 :: npartTotal_ghead_8byte(6)
  integer*4 :: flag_cooling_ghead       ! flags whether cooling was included  
  integer*4 :: num_files_ghead          ! number of files in multi-file snapshot 
  real*8 :: BoxSize_ghead               ! box-size of simulation in case periodic boundaries were used. Could be kpc/h or Mpc/h set with gadget_length_flag
  real*8 :: Omega0_ghead                ! matter density in units of critical density 
  real*8 :: OmegaLambda_ghead           ! cosmological constant parameter 
  real*8 :: HubbleParam_ghead           ! Hubble parameter in units of 100 km/sec/Mpc 
  integer*4 :: flag_stellarage_ghead    ! flags whether the file contains formation times of star particles 
  integer*4 :: flag_metals_ghead                ! flags whether the file contains metallicity values for gas and star particles
  integer*4 :: npartTotalHighWord_ghead(6)      ! High word of the total number of particles of each type !!!! CHECK THIS - unsigned!!!!!!
  integer*8 :: npartTotalHighWord_ghead_8byte(6)
  integer*4 :: flag_entropy_instead_u_ghead     ! flags that IC-file contains entropy instead of u 
  character(60) :: fill_ghead           ! fill rest of header to 256 bytes
  integer*4 :: part_data_size

  ! Particle data handling:

  integer*8 :: total_number_of_particles, particles_on_rank
  integer*8, allocatable :: particles_on_rank_and_file(:)
  real*4 :: xv_read(6),xv(6,read_length),xv_to_send(6,read_length,0:nodes-1)
  real*4 :: xv_to_recv(6,read_length+1,0:nodes-1),send_xv_buf(6,read_length),recv_xv_buf(6,read_length)
  integer*8 :: coord_read_position, vel_read_position, particle_count, file_size
  integer*4, allocatable :: particle_counts_for_sending(:,:) ! count of particles to be sent (rank,num_read)
  integer*4, allocatable :: particle_counts_for_receiving(:,:) ! reciprical array for count of particles to be sent (rank,num_read)
  integer*4 :: particle_stack_count_send(0:nodes-1) ! records how many particles a rank has stacked to send to other ranks
  integer*4 :: particle_stack_count_recv(0:nodes-1) ! records how many particles a rank has to prepare to receive from other ranks
  integer*4 :: particle_stack_count_send_max(0:nodes-1) ! records how many particles a rank has stacked to send to other ranks
  real*4 :: min_x_gadget, min_y_gadget, min_z_gadget, max_x_gadget, max_y_gadget, max_z_gadget  
  real*4 :: min_x_cubep3m, min_y_cubep3m, min_z_cubep3m, max_x_cubep3m, max_y_cubep3m, max_z_cubep3m  
  integer*4 :: particle_write_count, particle_read_count
  integer*8 :: particle_read_count_global,particle_write_count_global


  ! file handling

  integer*4 :: read_count, num_reads_local,num_reads_max,num_reads_check
  character(400) :: file, ifile
  integer*4, allocatable :: file_list(:)
  integer*4 :: file_loop_counter

  ! MPI variables

  integer*4 :: nodes_returned, ierr, rank,count
  integer*4 :: tag, rtag, send_rank,recv_rank
  integer*4 :: requests(2*nodes)
  integer*4, dimension(MPI_STATUS_SIZE,2*nodes) :: wait_status      
  integer*4, dimension(MPI_STATUS_SIZE) :: status

  ! Utility Variables       

  integer*4 :: i, j, k, n, p, dummy




  !------------------------------------------------------------------------------------------------------------!      


  call mpi_initialise
  call node_coords_initialise        ! sets up cubep3m node coords  

  call file_distribution_calculation  ! works out file reading mapping based on number of gadget files and 
  ! number of cubep3m nodes (= number of cubep3m files)

  particle_write_count = 0
  particle_read_count = 0
  particle_write_count_global = 0
  particle_read_count_global = 0
  num_reads_check = -1       

  ! first let us inspect the gadget file headers and make sure the simulation parameters make sense:

  do file_loop_counter = 1, num_files_gadget/nodes + 1  ! if we're reading multiple gadget files
    ! per cubep3m node we iterate in this loop
    call open_gadget_file(file_list(file_loop_counter)) ! file list of gadget files to open on each node 
    ! is set in the file_distribution subroutine
    call read_gadget_header                       

    if(file_list(file_loop_counter) .ne. -1) print*, 'Particles in file ',file_list(file_loop_counter),&
      &' = ',npart_ghead(2)

    if(num_reads_check .lt. npart_ghead(2) / read_length + 1) then
      num_reads_check = npart_ghead(2) / read_length + 1
    endif              

  enddo

  if(debug) call print_gadget_header(0) ! call print gadget header on rank 0

  call calculate_total_particles

  call mpi_allreduce(num_reads_check,num_reads_max,1,mpi_integer, &
    &mpi_max,mpi_comm_world,ierr)

  if(rank .eq. 0 .and. debug) print*, 'MAX NUM READS = ',num_reads_max

  allocate(particle_counts_for_sending(0:nodes-1,num_reads_max*(num_files_gadget/nodes + 1)))
  allocate(particle_counts_for_receiving(0:nodes-1,num_reads_max*(num_files_gadget/nodes + 1)))
  particle_counts_for_sending = 0
  particle_counts_for_receiving = 0      



  call check_particle_counts      

  if(rank == 0)  then
    call check_simulation_parameters
  endif

  call mpi_barrier(mpi_comm_world,ierr)

  ! broadcast simulation info

  call MPI_Bcast(nc, 1, MPI_INTEGER, 0, mpi_comm_world,ierr)
  call MPI_Bcast(ncc, 1, MPI_REAL, 0, mpi_comm_world,ierr)

  call mpi_barrier(mpi_comm_world,ierr)

  ! open cubep3m files

  call open_cubep3m_files

  ! write a dummy integer for the cubp3m header -- to be updated with the full particle count as a final step

  call write_dummy_cubep3m_header

  ! now assess particle data to configure memory for MPI sends and receives

  do file_loop_counter = 1, num_files_gadget/nodes + 1  ! if we're reading multiple gadget files
    ! per cubep3m node we iterate in this loop              

    min_x_gadget = 1.0e12
    min_y_gadget = 1.0e12
    min_z_gadget = 1.0e12     
    max_x_gadget = -1.0e12
    max_y_gadget = -1.0e12
    max_z_gadget = -1.0e12    

    call open_gadget_file(file_list(file_loop_counter))
    call read_gadget_header
    call assess_particle_data

    if(debug) call show_min_max_gadget


  enddo ! end of file_loop_counter iteration      



  ! if(debug) call write_particle_map_to_file ! produces files (the number of them
  ! depends on the number of nodes) 
  ! for debugging. The files contain the 
  ! particle send counts arrays, i.e. which 
  ! rank is sending how many particles where on 
  ! each read iteration



  call check_particle_map ! checks number of particles to send is correct      

  call send_particle_map  ! sends the particle_counts_for_sending array info so 
  ! that all tasks can assign particle_counts_for_recieving data



  call check_particle_map_recv ! checks number of particles received by tasks is correct



  ! now read the particle data, distribute it to the correct cubep3m nodes, and write it to disk                  

  do file_loop_counter = 1, num_files_gadget/nodes + 1  ! if we're reading multiple gadget files
    ! per cubep3m node we iterate in this loop                         

    call open_gadget_file(file_list(file_loop_counter))
    call read_gadget_header
    call read_and_send_particle_data
    call mpi_barrier(mpi_comm_world,ierr)      

  enddo ! end of file_loop_counter iteration      

  call check_particle_numbers

  if(debug) call show_min_max_cubep3m  

  call finalise_cubep3m_header

  call mpi_finalize(ierr)

contains

  !------------------------------------------------------------------------------------------!  

  subroutine open_gadget_file(num)

    integer*4 :: num

    if(file_list(file_loop_counter) .ne. -1) then

      write(file,'(i4)') num
      file=adjustl(file)

      ifile=gadget_path//gadget_root//'.'//file(1:len_trim(file))

#ifdef GFORTRAN		           
      open (unit=10,file=ifile,access='stream',status='old')
#endif
#ifdef BINARY
      open (unit=10,file=ifile,access='stream',form='binary',status='old')
#endif        

      if(debug) print*, 'successfully opened: ',ifile(1:len_trim(ifile)), ' on rank', rank

    endif

  end subroutine 

  !------------------------------------------------------------------------------------------!

  subroutine print_gadget_header(num)

    integer*8 :: num

    if(rank .eq. num) then

      print*, 'Gadget Header Size =    ', header_length_ghead
      print*, 'Particles in this file (rank = ',rank,'):'
      print*, 'Gas =                   ', npart_ghead(1)
      print*, 'DM =                    ', npart_ghead(2)  
      print*, 'Disk =                  ', npart_ghead(3)
      print*, 'Bulge =                 ', npart_ghead(4)  
      print*, 'Stars =                 ', npart_ghead(5)
      print*, 'Bndry =                 ', npart_ghead(6)  
      print*
      print*, 'Particle Masses (in Gadget internal units):'
      print*, 'Gas =                   ', mass_ghead(1)
      print*, 'DM =                    ', mass_ghead(2)  
      print*, 'Disk =                  ', mass_ghead(3)
      print*, 'Bulge =                 ', mass_ghead(4)  
      print*, 'Stars =                 ', mass_ghead(5)
      print*, 'Bndry =                 ', mass_ghead(6)  
      print*
      print*, 'Expansion Factor:       ', time_ghead
      print*, 'Redshift:               ', redshift_ghead
      print*, 'Omega M0:               ', Omega0_ghead
      print*, 'Omega Lambda0:          ', OmegaLambda_ghead
      if(gadget_length_flag .eq. 0) print*, 'Boxsize (Mpc/h):        ', BoxSize_ghead
      if(gadget_length_flag .eq. 1) print*, 'Boxsize (kpc/h):        ', BoxSize_ghead
      print*, 'H_0 (km/s):             ', HubbleParam_ghead
      print* 
      print*, 'Number of files for ICs:', num_files_ghead
      print* 
      print*, 'Total Particles in Simulation'
      print*, 'Gas =                   ', npartTotal_ghead_8byte(1)
      print*, 'DM =                    ', npartTotal_ghead_8byte(2)  
      print*, 'Disk =                  ', npartTotal_ghead_8byte(3)
      print*, 'Bulge =                 ', npartTotal_ghead_8byte(4)  
      print*, 'Stars =                 ', npartTotal_ghead_8byte(5)
      print*, 'Bndry =                 ', npartTotal_ghead_8byte(6)  
      print*
      print*, 'Star Formation Flag:    ', flag_sfr_ghead
      print*, 'Feedback Flag:          ', flag_feedback_ghead
      print*, 'Cooling Flag:           ', flag_cooling_ghead
      print*, 'Steller Age Flag:       ', flag_stellarage_ghead
      print*, 'Flag Metals:            ', flag_metals_ghead
      print*, 'Entropy Flag:           ', flag_entropy_instead_u_ghead
      print* 
      print*, 'Total Particles High Word'
      print*, 'Gas =                   ', npartTotalHighWord_ghead_8byte(1)
      print*, 'DM =                    ', npartTotalHighWord_ghead_8byte(2)  
      print*, 'Disk =                  ', npartTotalHighWord_ghead_8byte(3)
      print*, 'Bulge =                 ', npartTotalHighWord_ghead_8byte(4)  
      print*, 'Stars =                 ', npartTotalHighWord_ghead_8byte(5)
      print*, 'Bndry =                 ', npartTotalHighWord_ghead_8byte(6)  
      print*

    endif

  end subroutine print_gadget_header

  !------------------------------------------------------------------------------------------!

  subroutine read_gadget_header

    integer*4 :: i, int_4byte
    integer*8 :: int_8byte

    if(file_list(file_loop_counter) .ne. -1) then

      read(10) header_length_ghead,npart_ghead,mass_ghead,time_ghead,&
        &redshift_ghead,flag_sfr_ghead,flag_feedback_ghead,npartTotal_ghead,&
        &flag_cooling_ghead,num_files_ghead,BoxSize_ghead,Omega0_ghead,OmegaLambda_ghead,&
        &HubbleParam_ghead,flag_stellarage_ghead,flag_metals_ghead,&
        &npartTotalHighWord_ghead,flag_entropy_instead_u_ghead, fill_ghead


      ! convert unsigned integers into 8 byte fortran integers:

      do i = 1, 6

        int_4byte = npartTotal_ghead(i)
        int_8byte = int(int_4byte,8)

        if(int_8byte .lt. 0) int_8byte = 4294967296 - abs(int_8byte)

        npartTotal_ghead_8byte(i) = int_8byte

      enddo

      do i = 1, 6

        int_4byte = npartTotalHighWord_ghead(i)

        int_8byte = int(int_4byte,8)

        if(int_8byte .lt. 0) int_8byte = 4294967296 - abs(int_8byte)

        npartTotalHighWord_ghead_8byte(i) = int_8byte

      enddo


      ! check gadget file count agrees:


      if(rank .eq. 0) then

        if(num_files_ghead .ne. num_files_gadget) then

          print*, 'GADGET FILE COUNT INCORRECT!!!!'
          print*, 'Gadget header file count = ',num_files_ghead
          print*, 'Input file count = ', num_files_gadget
          call mpi_abort(mpi_comm_world,ierr,ierr)

        endif

      endif

      particles_on_rank_and_file(file_loop_counter-1) = int(npart_ghead(2),8)

    endif

  end subroutine


  !------------------------------------------------------------------------------------------!


  subroutine assess_particle_data


    integer*4 target_rank

    ! this subroutine reads all the particle data in order to work out which node each particle
    ! needs to be sent to. This allows arrays of particle data to be allocated on each node to send
    ! and receive the appropriate number of particles.

    if(file_list(file_loop_counter) .ne. -1) then ! if this is a reading node

      !      if(debug .and. file_loop_counter .eq. 1 .and. rank .eq. 0)&
      !      & print*, 'Read length = ',read_length

      xv = 0

      particle_count = 0


      num_reads_local = npart_ghead(2) / read_length + 1


      if(debug .and. file_loop_counter .eq. 1) print*, &
        &'Number of reads = ', num_reads_local,' on rank ',rank

      read(10) dummy,part_data_size
      if(debug .and. file_loop_counter .eq. 1) print*, &
        &'Total bytes of position data in this file: ',part_data_size,' on rank ',rank

      coord_read_position = 4+256+4+4+1    ! i.e. header length (4) plus 
      ! header (256) plus 
      ! header length (4) plus part_data_size (4)
      ! plus 1 to reach the start of the coord data
      do read_count = 1, num_reads_local
        if(read_count .lt. num_reads_local) then
          do i = 1, read_length
            read(10,POS=coord_read_position) xv_read(1:3)

            ! convert to cubep3m coordinates:

            xv_read(1:3) = real((real(xv_read(1:3),8)*real(nc,8)/real(BoxSize_ghead,8)),4)

            ! calculate where this particle needs to be in the cubic decomposition

	if(int(xv_read(1)/ncc).lt.nodes_dim .and.&
	    & int(xv_read(2)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(3)/ncc).lt.nodes_dim) then

	target_rank = int(xv_read(1)/ncc)+nodes_dim*int(xv_read(2)/ncc)+(nodes_dim**2.0)*int(xv_read(3)/ncc)

      elseif(int(xv_read(1)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(2)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(3)/ncc).lt.nodes_dim) then

	target_rank = (nodes_dim-1)+nodes_dim*int(xv_read(2)/ncc)+(nodes_dim**2.0)*int(xv_read(3)/ncc)

      elseif(int(xv_read(1)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(2)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(3)/ncc).lt.nodes_dim) then

	target_rank = int(xv_read(1)/ncc)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*int(xv_read(3)/ncc)

      elseif(int(xv_read(1)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(2)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(3)/ncc).ge.nodes_dim) then

	target_rank = int(xv_read(1)/ncc)+nodes_dim*int(xv_read(2)/ncc)+(nodes_dim**2.0)*(nodes_dim-1)


      elseif(int(xv_read(1)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(2)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(3)/ncc).lt.nodes_dim) then

	target_rank = (nodes_dim-1)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*int(xv_read(3)/ncc)

      elseif(int(xv_read(1)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(2)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(3)/ncc).ge.nodes_dim) then

	target_rank = (nodes_dim-1)+nodes_dim*int(xv_read(2)/ncc)+(nodes_dim**2.0)*(nodes_dim-1)

      elseif(int(xv_read(1)/ncc).lt.nodes_dim.and.&
	    & int(xv_read(2)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(3)/ncc).ge.nodes_dim) then

	target_rank = int(xv_read(1)/ncc)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)

      elseif(int(xv_read(1)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(2)/ncc).ge.nodes_dim.and.&
	    & int(xv_read(3)/ncc).ge.nodes_dim) then

	target_rank = (nodes_dim-1)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)


      endif



            if(target_rank .gt. nodes-1 .or. target_rank .lt. 0) then

              print*,"ERROR - TARGET RANK NOT ALLOWED FOR PARTICLE ",xv_read(1:3),&
                      &" ON RANK", rank, "INTENDED TARGET = ", target_rank

                      call mpi_abort(mpi_comm_world,ierr,ierr)

            endif

            particle_counts_for_sending(target_rank,read_count*file_loop_counter) = &
              &particle_counts_for_sending(target_rank,read_count*file_loop_counter) + 1

            if(max_x_gadget .lt. xv_read(1)) max_x_gadget = xv_read(1)
            if(min_x_gadget .gt. xv_read(1)) min_x_gadget = xv_read(1)
            if(max_y_gadget .lt. xv_read(2)) max_y_gadget = xv_read(2)
            if(min_y_gadget .gt. xv_read(2)) min_y_gadget = xv_read(2)
            if(max_z_gadget .lt. xv_read(3)) max_z_gadget = xv_read(3)
            if(min_z_gadget .gt. xv_read(3)) min_z_gadget = xv_read(3)

            particle_count = particle_count + 1        
            coord_read_position = coord_read_position + 12 ! i.e. plus 3 * 4 bytes of pos data

          enddo
        else ! read_count .ge. num_reads_local i.e this is the last read iteration
          do i = 1, npart_ghead(2) - (read_count-1)*read_length

            read(10,POS=coord_read_position) xv_read(1:3)

            ! convert to cubep3m coordinates:
            xv_read(1:3) = real((real(xv_read(1:3),8)*real(nc,8)/real(BoxSize_ghead,8)/1000.0),4)

            ! calculate where this particle needs to be in the cubic decomposition
            target_rank = int(xv_read(1)/(ncc+0.001))+nodes_dim*int(xv_read(2)/(ncc+0.001))+&
              &(nodes_dim**2.0)*int(xv_read(3)/(ncc+0.001))

            if(target_rank .gt. nodes-1 .or. target_rank .lt. 0) then

              print*,"ERROR - TARGET RANK NOT ALLOWED FOR PARTICLE ",xv_read(1:3),&
                      &" ON RANK", rank, "INTENDED TARGET = ", target_rank

                    call mpi_abort(mpi_comm_world,ierr,ierr)
            endif


            particle_counts_for_sending(target_rank,read_count*file_loop_counter) = &
              &particle_counts_for_sending(target_rank,read_count*file_loop_counter) + 1

            if(max_x_gadget .lt. xv_read(1)) max_x_gadget = xv_read(1)
            if(min_x_gadget .gt. xv_read(1)) min_x_gadget = xv_read(1)
            if(max_y_gadget .lt. xv_read(2)) max_y_gadget = xv_read(2)
            if(min_y_gadget .gt. xv_read(2)) min_y_gadget = xv_read(2)
            if(max_z_gadget .lt. xv_read(3)) max_z_gadget = xv_read(3)
            if(min_z_gadget .gt. xv_read(3)) min_z_gadget = xv_read(3)

            particle_count = particle_count + 1        
            coord_read_position = coord_read_position + 12 ! i.e. plus 3 * 4 bytes of pos data

          enddo
        endif
      enddo

      !      if(debug) print*, 'Particles read check sum (on assessment) = ', &
      !                              &particle_count-npart_ghead(2),' on rank ',rank

      if(particle_count-npart_ghead(2) .ne. 0) then
        print*, 'ERROR -- PARTICLE COUNT (on read) INCORRECT!!!'
        call mpi_abort(mpi_comm_world,ierr,ierr)
      endif

    endif    ! if this is a reading node


  end subroutine assess_particle_data

  !------------------------------------------------------------------------------------------!

  subroutine read_and_send_particle_data

    ! this subroutine reads all the particle data (again) and sends it in batches to the correct desintation
    ! nodes. The batch size is controlled by the 'read_block_size' parameter


    integer*4 :: target_rank, i_max

    ! first clear any left over send/receive info

    xv_to_send = 0
    particle_stack_count_send_max = 0
    particle_stack_count_send = 0
    xv = 0      
    particle_count = 0
    particle_counts_for_sending = 0
    i_max = 0

    if(file_list(file_loop_counter) .ne. -1) then  ! if this is a reading node

      !        if(debug) print*, 'Read length = ',read_length,' on rank ',rank

      num_reads_local = npart_ghead(2) / read_length + 1        

      read(10) dummy,part_data_size

      coord_read_position = 4+256+4+4+1    ! i.e. header length (4) plus 
      ! header (256) plus 
      ! header length (4) plus part_data_size (4)
      ! plus 1 to reach the start of the coord data


      vel_read_position = 4+256+4+4+part_data_size+4+4+1
      ! i.e. header length (4) plus 
      ! header (256) plus 
      ! header length (4) plus part_data_size (4)
      ! plus part data (part_data_size) +
      ! part(vel)_data_size (4) 
      ! plus 1 to reach the start of the vel data

      file_size = vel_read_position - 1 + part_data_size + 4

      read(10,POS=vel_read_position-4) part_data_size

    endif  ! if this is a reading node


    do read_count = 1, num_reads_max

      particle_stack_count_recv = 0
      particle_stack_count_send = 0


      if(read_count .le. num_reads_local .and. file_list(file_loop_counter) .ne. -1) then

        particle_stack_count_send = 0  

        if(read_count .lt. num_reads_local) i_max = read_length
        if(read_count .eq. num_reads_local) i_max = npart_ghead(2) - (read_count-1)*read_length

        do i = 1, i_max
          read(10,POS=coord_read_position) xv_read(1:3)

          particle_read_count = particle_read_count + 1

          coord_read_position = coord_read_position + 12 ! i.e. plus 3 * 4 bytes of pos data                                

          ! convert to cubep3m coordinates:
          xv_read(1:3) = real((real(xv_read(1:3),8)*real(nc,8)/real(BoxSize_ghead,8)),4)
          particle_count = particle_count + 1

            target_rank = int(xv_read(1)/(ncc+0.001))+nodes_dim*int(xv_read(2)/(ncc+0.001))+&
              &(nodes_dim**2.0)*int(xv_read(3)/(ncc+0.001))

          particle_stack_count_send(target_rank) = particle_stack_count_send(target_rank) + 1

          ! add to particle data array:                
          xv(1:3,i) = xv_read(1:3)
        enddo

        do i = 1, i_max
          read(10,POS=vel_read_position) xv_read(4:6)                      
          vel_read_position = vel_read_position + 12 ! i.e. plus 3 * 4 bytes of pos data

          ! convert to cubep3m velocities:

          xv_read(4:6) = xv_read(4:6)*sqrt(1.0/(1.0+redshift_ghead)) ! converts from gadget to km/s


          if(gadget_length_flag .eq. 0) then

            xv_read(4:6) = xv_read(4:6)/(BoxSize_ghead*150.0*Omega0_ghead**0.5/real(nc,8))/(1.0+redshift_ghead)

          elseif(gadget_length_flag .eq. 1) then

            xv_read(4:6) = xv_read(4:6)/((BoxSize_ghead/1000.0)*150.0*Omega0_ghead**0.5/real(nc,8))/(1.0+redshift_ghead)

          endif

          ! add to particle data array:                
          xv(4:6,i) = xv_read(4:6)
        enddo

        ! Now prepare array for sending particles:

        particle_stack_count_send_max = particle_stack_count_send
        particle_stack_count_send = 0 

        do i = 1, i_max

          ! calculate where this particle needs to be in the cubic decomposition
            target_rank = int(xv(1,i)/(ncc+0.001))+nodes_dim*int(xv(2,i)/(ncc+0.001))+&
              &(nodes_dim**2.0)*int(xv(3,i)/(ncc+0.001))

          particle_stack_count_send(target_rank) = particle_stack_count_send(target_rank) + 1                
          xv_to_send(1:6,particle_stack_count_send(target_rank),target_rank) = xv(1:6,i)

        enddo  

        ! check numbers of particles match:

        do i = 0, nodes_dim-1
          if(particle_stack_count_send(i) .ne. particle_stack_count_send_max(i)) then
            print*, 'ABORT - PARTICLE STACK COUNT MIS-MATCH IN SEND ARRAY on rank ',rank
            print*, 'expected particle counts:    actual particle counts:'
            do k = 0, nodes-1
              print*, particle_stack_count_send_max(k), particle_stack_count_send(k)
            enddo
            print*
            call mpi_abort(mpi_comm_world,ierr,ierr)            
          endif
        enddo

      endif ! if this is a reading node


      ! now send and receive particles. Each rank knows how many particles is it expecting
      ! from the entries in the particle_counts_for_receiving array

      call send_particle_stack_counts
      call send_receive_particle_data
      call write_cubep3m_data          
      call mpi_barrier(mpi_comm_world,ierr)          

    enddo


    !     if(debug) print*, 'Particles read check sum = (on read+send)',&   
    !                            &particle_count-npart_ghead(2),' on rank ',rank

    if(file_list(file_loop_counter) .ne. -1 .and. particle_count-npart_ghead(2) .ne. 0) then
      print*, 'ERROR -- PARTICLE COUNT (on read) INCORRECT!!!'
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif



  end subroutine read_and_send_particle_data

  !------------------------------------------------------------------------------------------!

  subroutine calculate_total_particles

    total_number_of_particles = npartTotal_ghead_8byte(2) + &
      & npartTotalHighWord_ghead_8byte(2)*4294967296

  end subroutine calculate_total_particles

  !------------------------------------------------------------------------------------------!

  subroutine check_particle_counts

    integer*8 npart_total_check, i

    particles_on_rank = sum(particles_on_rank_and_file)

    call mpi_allreduce(particles_on_rank,npart_total_check,1,mpi_integer8, &
      mpi_sum,mpi_comm_world,ierr)

    do i = 1, nodes

      if(rank .eq. i .and. npart_total_check .ne. total_number_of_particles) then

        print*, 'TOTAL PARTICLES INCORRECT ON RANK', rank
        print*, 'Total (assumed) = ', total_number_of_particles
        print*, 'Total (MPI_allreduce) = ', npart_total_check
        call mpi_abort(mpi_comm_world,ierr,ierr)

      endif

      call mpi_barrier(mpi_comm_world,ierr)


    enddo

  end subroutine check_particle_counts

  !------------------------------------------------------------------------------------------!

  subroutine check_simulation_parameters

    real*16 :: x,y,z

    x = 1.0
    y = 3.0
    z = 2.0

    ! check total particles is a cube of an integer

    if(abs((real(total_number_of_particles,16)**(x/y))**y - &
      &real(total_number_of_particles,16)) .lt. 0.001) then

      print*    
      print*, 'TOTAL PARTICLES = ',total_number_of_particles,&
        ' = ',nint((real(total_number_of_particles,16)**(x/y))),'^3'
      print*

    else

      print*    
      print*, 'TOTAL PARTICLES NOT EQUAL TO A CUBIC NUMBER!!!'
      print*, '           Particles:    Cube root:  '
      print*, total_number_of_particles,'   ',&
        &(real(total_number_of_particles,8)**(x/y))
      print*
      call mpi_abort(mpi_comm_world,ierr,ierr)

    endif


    ! check nodes_dim divides into npart^(1/3)

    if(real(total_number_of_particles,16)**(x/y) /&
      & real(nodes_dim,16) - real(nint(real(total_number_of_particles,16)**(x/y) / &
      & real(nodes_dim,16)),16) .lt.0.0001) then

      print*    
      print*, 'Particles per dimension = ',nint((real(total_number_of_particles,16)**(x/y)))
      print*, 'Nodes per dimension = ', nodes_dim
      print*, 'Particles per node = ', nint((real(total_number_of_particles,16)**(x/y)))/nodes_dim
      print*

    else  

      print*    
      print*, 'ERROR: PARTICLE COUNT DOES NOT MATCH NODE COUNT!!!'
      print*, 'Particles per dimension = ',nint((real(total_number_of_particles,16)**(x/y)))
      print*, 'Nodes per dimension = ', nodes_dim
      print*, 'Particles per node = ', (real(total_number_of_particles,16)**(x/y))/real(nodes_dim,16)
      print*

      call mpi_abort(mpi_comm_world,ierr,ierr)

    endif

    ! check nodes_dim divides into fine_cells^(1/3)

    if((z*real(total_number_of_particles,16))**(x/y) /&
      & real(nodes_dim,16) - real(nint(z*(real(total_number_of_particles,16))**(x/y) / &
      & real(nodes_dim,16)),16) .lt.0.0001) then

      print*    
      print*, 'Fine cells per dimension = ',nint(z*(real(total_number_of_particles,16))**(x/y))
      print*, 'Nodes per dimension = ', nodes_dim
      print*, 'Fine cells per node = ', nint(z*(real(total_number_of_particles,16))**(x/y))/nodes_dim
      print*

      ncc = nint(z*(real(total_number_of_particles,16))**(x/y))/nodes_dim
      nc = ncc*nodes_dim

    else

      print*    
      print*, 'ERROR: FINE CELL COUNT DOES NOT MATCH NODE COUNT!!!'
      print*, 'Fine cells per dimension = ',nint(z*(real(total_number_of_particles,16))**(x/y))
      print*, 'Nodes per dimension = ', nodes_dim
      print*, 'Fine cells per node = ', nint(z*(real(total_number_of_particles,16))**(x/y))/nodes_dim
      print*

      call mpi_abort(mpi_comm_world,ierr,ierr)

    endif

    ! Print simulation attributes

    print*
    print*        
    print*, 'SIMULATION ATTRIBUTES:'
    print*
    print*, 'NP = ',nint((real(total_number_of_particles,16)**(x/y))),'^3'
    print*, 'NODES = ', nodes_dim
    print*
    print*, 'INITIAL REDSHIFT = ', redshift_ghead
    print*
    if(gadget_length_flag .eq. 0) print*, 'BOXSIZE (Mpc/h) = ', BoxSize_ghead
    if(gadget_length_flag .eq. 1) print*, 'BOXSIZE (kpc/h) = ', BoxSize_ghead
    print*
    print*, 'NB, check your choice of units for the box size!'
    print*    
    print*, 'COSMOLOGY:'
    print*
    print*, 'OMEGA_M0 = ', Omega0_ghead
    print*, 'OMEGA_LAMBDA0 = ', OmegaLambda_ghead
    print*, 'H0 (km/s) = ', HubbleParam_ghead
    print*
    print*        


  end subroutine check_simulation_parameters

  !------------------------------------------------------------------------------------------!

  subroutine mpi_initialise

    call mpi_init(ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    call mpi_comm_size(mpi_comm_world,nodes_returned,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    if (nodes_returned /= nodes ) then
      write(*,*) 'IC converter compiled for a different number of nodes'
      write(*,*) 'mpirun nodes=',nodes_returned,'converter nodes=',nodes 
      call mpi_abort(mpi_comm_world,ierr,ierr)
    endif

    call mpi_comm_rank(mpi_comm_world,rank,ierr)
    if (ierr /= mpi_success) call mpi_abort(mpi_comm_world,ierr,ierr)

    requests = 0

  end subroutine mpi_initialise

  !------------------------------------------------------------------------------------------!

  subroutine file_distribution_calculation

    integer*4 :: i,j,k

    ! number of gadget files = num_files_gadget
    ! number of cubep3m files = nodes

    if(num_files_gadget .eq. nodes) then ! this is the ideal scenario

      allocate(file_list(1)) ! one-to-one correspondence between files
      allocate(particles_on_rank_and_file(1))

      do i = 0, nodes-1

        if(rank .eq. i .and. i .lt. num_files_gadget) file_list(1) = i

      enddo


    elseif(num_files_gadget .lt. nodes) then ! we have some redundant nodes for reading

      allocate(file_list(1))
      allocate(particles_on_rank_and_file(1))

      do i = 0, nodes-1

        if(rank .eq. i .and. i .lt. num_files_gadget) file_list(1) = i
        if(rank .eq. i .and. i .ge. num_files_gadget) file_list(1) = -1 ! redundant nodes

      enddo


    elseif(num_files_gadget .gt. nodes) then ! we have to read multiple gadget files on some/all nodes

      k = 0

      allocate(file_list(num_files_gadget/nodes + 1))
      allocate(particles_on_rank_and_file(0:num_files_gadget/nodes))        


      do j = 1, num_files_gadget/nodes+1

        do i = 0, nodes-1

          if(rank .eq. i) file_list(j) = k ! files to be read on this node

          if(k .ge. num_files_gadget-1 .and. (rank+nodes*(j-1)) .ge. num_files_gadget) file_list(j) = -1


          k = k + 1                                            

        enddo

      enddo



    endif



      do i = 0, nodes

        if(rank .eq. i) print*, 'Filelist on rank ',rank,':',file_list
        call mpi_barrier(mpi_comm_world,ierr)
      enddo



    particles_on_rank_and_file = 0

  end subroutine file_distribution_calculation

  !------------------------------------------------------------------------------------------!

  subroutine node_coords_initialise

    integer i, j, k, n

    !! Calculate node_coords
    do n=0,nodes_dim**3-1
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

  subroutine show_min_max_gadget

    if(rank.eq.0) print*,'           MIN MAX PARTICLE POSITIONS ON EACH RANK (GADGET):'
    if(rank.eq.0) print*,'           (RANK, min_x,y,z, max_x,y,z)'      

    do i = 0, nodes-1
      if(rank .eq. i) print*, rank, min_x_gadget, min_y_gadget,min_z_gadget,max_x_gadget,max_y_gadget,max_z_gadget            
      call mpi_barrier(mpi_comm_world,ierr)
    enddo

  end subroutine show_min_max_gadget

  !------------------------------------------------------------------------------------------!

  subroutine write_particle_map_to_file

    integer*4 :: j
    character(20) :: format_string, num_ints

    write(file,'(i4)') rank
    file=adjustl(file)
    write(num_ints,'(i4)') size(particle_counts_for_sending(1,:))+1
    num_ints=adjustl(num_ints)



    format_string=num_ints(1:len_trim(num_ints))//'i10'


    ifile='./particle_send_map_'//file(1:len_trim(file))//'.dat'

    open (unit=99,file=ifile,form='formatted')



    do j = 0, nodes-1
      write(99,'('//format_string(1:len_trim(format_string))//')') j,particle_counts_for_sending(j,:)
    enddo


    close(99)


  end subroutine write_particle_map_to_file

  !------------------------------------------------------------------------------------------!

  subroutine check_particle_map

    integer*8 :: total_particle_check_map, local_particle_check_map

    local_particle_check_map = int(sum(particle_counts_for_sending),8)

    call mpi_allreduce(local_particle_check_map,total_particle_check_map,1,mpi_integer8, &
      mpi_sum,mpi_comm_world,ierr)

    if(rank .eq. 0 .and.total_particle_check_map.eq.total_number_of_particles) then

      if(debug) print*, "TOTAL PARTICLE MAP SUM CHECK = OK"

    else
      if(rank .eq. 0) then
        print*, "TOTAL PARTICLE MAP SUM CHECK = INCORRECT"
        print*, 'Total calculated from maps = ',total_particle_check_map
        print*, 'Total particles = ',total_number_of_particles
        call  mpi_abort(mpi_comm_world,ierr,ierr)
      endif
    endif

  end subroutine check_particle_map

  !------------------------------------------------------------------------------------------!

  subroutine send_particle_map

    integer*4 :: i, j, count, send,recv
    integer*4, allocatable :: send_data(:),recv_data(:)
    integer*4 :: array_row_count

    array_row_count = num_reads_max*(num_files_gadget/nodes+1)

    allocate(send_data(array_row_count))
    allocate(recv_data(array_row_count))

    send_data = 0
    recv_data = 0

    particle_counts_for_receiving(rank,:) = particle_counts_for_sending(rank,:)

    do i = 0, nodes-1 


      rtag = modulo(rank-i,nodes)**2
      tag = rank**2       

      send_rank = modulo(i+rank,nodes) ! rank this rank is sending to
      recv_rank = modulo(rank-i,nodes) ! rank this rank is receiving from

      send_data(:) = particle_counts_for_sending(send_rank,:)

      if(MOD(rank,2) .ne. 0) then ! odd ranks

        call mpi_isend(send_data(1),array_row_count,mpi_integer,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

        call mpi_irecv(recv_data(1),array_row_count,mpi_integer,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

      endif

      if(MOD(rank,2) .eq. 0) then ! even ranks

        call mpi_irecv(recv_data(1),array_row_count,mpi_integer,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

        call mpi_isend(send_data(1),array_row_count,mpi_integer,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

      endif


      call mpi_waitall(2*nodes,requests, wait_status, ierr)

      particle_counts_for_receiving(recv_rank,:) = recv_data(:)

    enddo


    deallocate(send_data)
    deallocate(recv_data)



  end subroutine send_particle_map

  !------------------------------------------------------------------------------------------!

  subroutine check_particle_map_recv

    integer*8 :: total_particle_check_map, local_particle_check_map

    local_particle_check_map = int(sum(particle_counts_for_receiving),8)

    call mpi_allreduce(local_particle_check_map,total_particle_check_map,1,mpi_integer8, &
      mpi_sum,mpi_comm_world,ierr)

    if(rank .eq. 0 .and.total_particle_check_map.eq.total_number_of_particles) then

      if(debug) print*, "TOTAL PARTICLE MAP SUM CHECK (RECEIVED) = OK"

    else
      if(rank .eq. 0) then
        print*, "TOTAL PARTICLE MAP SUM CHECK (RECV) = INCORRECT"
        print*, 'Total calculated from maps = ',total_particle_check_map
        print*, 'Total particles = ',total_number_of_particles
        call  mpi_abort(mpi_comm_world,ierr,ierr)
      endif
    endif

  end subroutine check_particle_map_recv

  !------------------------------------------------------------------------------------------!

  subroutine send_particle_stack_counts

    integer*4 :: i, j, count, send,recv
    integer*4, allocatable :: send_data(:),recv_data(:)
    integer*4 :: array_row_count

    array_row_count = 1 ! we are just sending one number (the relevent particle stack count) to each node from each node

    allocate(send_data(array_row_count))
    allocate(recv_data(array_row_count))

    send_data = 0
    recv_data = 0
    particle_stack_count_recv = 0

    particle_stack_count_recv(rank) = particle_stack_count_send(rank) ! although we are sending zero particles from 
    ! a rank to itself we still need the memory to 
    ! be free to accommodate the particle data

    do i = 1, nodes-1 

      rtag = modulo(rank-i,nodes)**2
      tag = rank**2       

      send_rank = modulo(i+rank,nodes) ! rank this rank is sending to
      recv_rank = modulo(rank-i,nodes) ! rank this rank is receiving from

      send_data(1) = particle_stack_count_send(send_rank)

      if(MOD(rank,2) .ne. 0) then ! odd ranks

        call mpi_isend(send_data(1),array_row_count,mpi_integer,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

        call mpi_irecv(recv_data(1),array_row_count,mpi_integer,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

      endif

      if(MOD(rank,2) .eq. 0) then ! even ranks

        call mpi_irecv(recv_data(1),array_row_count,mpi_integer,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

        call mpi_isend(send_data(1),array_row_count,mpi_integer,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

      endif

      call mpi_waitall(2*nodes,requests, wait_status, ierr)

      particle_stack_count_recv(recv_rank) = recv_data(1)

    enddo


    deallocate(send_data)
    deallocate(recv_data)

  end subroutine send_particle_stack_counts

  !------------------------------------------------------------------------------------------!

  subroutine send_receive_particle_data

    integer*4 :: i, j, count, send,recv
    real*4, allocatable :: send_data(:),recv_data(:)
    integer*4 :: array_row_count_send,array_row_count_recv
    integer*4 :: dummy = -1.0 ! I've added this into the send/recv data to handle
    ! the case where there is zero data to be sent (it was 
    ! creating issues with arrays of zero size)


    xv_to_recv(:,1:read_length,rank) = xv_to_send(:,:,rank) ! no need to send data from a rank to itself

    do i = 1, nodes-1 

      rtag = modulo(rank-i,nodes)**2
      tag = rank**2       

      send_rank = modulo(i+rank,nodes) ! rank this rank is sending to
      recv_rank = modulo(rank-i,nodes) ! rank this rank is receiving from

      array_row_count_send = 6*particle_stack_count_send(send_rank)+1 ! the extra +1 is for the dummy
      array_row_count_recv = 6*particle_stack_count_recv(recv_rank)+1 ! the extra +1 is for the dummy

      allocate(send_data(array_row_count_send))
      allocate(recv_data(array_row_count_recv))

      call mpi_barrier(mpi_comm_world,ierr)


      send_data = 0
      recv_data = 0

      ! put data in send buffer if there is data to send

      if(array_row_count_send .ne. 1) then

        send_data(1) = dummy ! send a -1.0 in the first block as the dummy

        do j = 1,particle_stack_count_send(send_rank)
          do k = 1, 6
            send_data(6*(j-1)+k+1) = xv_to_send(k,j,send_rank)
          enddo
        enddo
      else

        send_data(1) = dummy ! send a -1.0 in the first block as the dummy

      endif

      if(MOD(rank,2) .ne. 0) then ! odd ranks

        call mpi_isend(send_data(1),array_row_count_send,mpi_real,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

        call mpi_irecv(recv_data(1),array_row_count_recv,mpi_real,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

      endif

      if(MOD(rank,2) .eq. 0) then ! even ranks

        call mpi_irecv(recv_data(1),array_row_count_recv,mpi_real,&
          &recv_rank,rtag,mpi_comm_world,requests(rank+1+nodes),ierr)  

        call mpi_isend(send_data(1),array_row_count_send,mpi_real,&
          &send_rank,tag,mpi_comm_world,requests(rank+1),ierr)

      endif

      call mpi_waitall(2*nodes,requests, wait_status, ierr)

      ! unpack recv_data if there is data to receive

      if(array_row_count_recv .ne. 1) then
        ! check the dummy value is correct
        if(recv_data(1) .ne. dummy) then
          print*,'SEND ERROR - DUMMY ARGUMENT NOT RECEIVED on rank ',rank
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif

        do j = 1,particle_stack_count_recv(recv_rank)
          do k = 1, 6
            xv_to_recv(k,j,recv_rank) = recv_data(6*(j-1)+k+1)
          enddo
        enddo
      else
        ! check the dummy value is correct
        if(recv_data(1) .ne. dummy) then
          print*,'SEND ERROR - DUMMY ARGUMENT NOT RECEIVED on rank ',rank
          call mpi_abort(mpi_comm_world,ierr,ierr)
        endif
      endif

      deallocate(send_data)
      deallocate(recv_data)

    enddo

    call mpi_barrier(mpi_comm_world,ierr)     

    ! we now have the data to write into the cubep3m files


  end subroutine send_receive_particle_data

  !------------------------------------------------------------------------------------------!

  subroutine write_cubep3m_data

    integer :: j,k,rank_number
    real*4 :: xv(6)



    do rank_number = 0, nodes-1
      do j = 1,particle_stack_count_recv(rank_number)

        !       write(89) xv_to_recv(:,j,rank_number)

        write(89) xv_to_recv(1,j,rank_number)-node_coords(rank,1)*ncc, &
          &xv_to_recv(2,j,rank_number)-node_coords(rank,2)*ncc, &
          &xv_to_recv(3,j,rank_number)-node_coords(rank,3)*ncc, &
          &xv_to_recv(4:6,j,rank_number)

        particle_write_count = particle_write_count + 1
      enddo
    enddo



  end subroutine write_cubep3m_data

  !------------------------------------------------------------------------------------------!

  subroutine open_cubep3m_files

    write(file,'(i4)') rank
    file=adjustl(file)

    ifile=cpm_path//'xv'//file(1:len_trim(file))//'.ic'

#ifdef GFORTRAN		           
    open (unit=89,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=89,file=ifile,form='binary',access='stream')
#endif        

    print*, 'successfully opened:',ifile(1:len_trim(ifile)), ' on rank', rank

  end subroutine open_cubep3m_files

  !------------------------------------------------------------------------------------------!

  subroutine write_dummy_cubep3m_header

    integer(4) :: dummy_particle_count = 0

    write(89) dummy_particle_count

  end subroutine write_dummy_cubep3m_header

  !------------------------------------------------------------------------------------------!

  subroutine check_particle_numbers

    integer :: rank_count

    if(debug) then

      if(rank .eq. 0) then
        print*
        print*,'Particle Read Counts on all ranks:'
        print*
      endif


      do rank_count = 0, nodes-1

        if(rank .eq. rank_count) print*,'Rank = ',rank,' Particles read = ',particle_read_count

        call mpi_barrier(mpi_comm_world,ierr)    

      enddo

      if(rank .eq. 0) then
        print*
        print*,'Particle Write Counts on all ranks:'
        print*
      endif



      do rank_count = 0, nodes-1

        if(rank .eq. rank_count) print*,'Rank = ',rank,' Particles written = ',particle_write_count

        call mpi_barrier(mpi_comm_world,ierr)    

      enddo

    endif

    call mpi_allreduce(int(particle_read_count,8),particle_read_count_global,1,mpi_integer8, &
      mpi_sum,mpi_comm_world,ierr)

    call mpi_allreduce(int(particle_write_count,8),particle_write_count_global,1,mpi_integer8, &
      mpi_sum,mpi_comm_world,ierr)


    if(rank .eq. 0) then
      if(particle_read_count_global .ne. particle_write_count_global) then

        print*,'MAJOR ERROR - GLOBAL READ COUNT NOT EQUAL TO GLOBAL WRITE COUNT!'
        print*,'PARTICLES READ =',particle_read_count_global
        print*,'PARTICLES WRITTEN =',particle_write_count_global
        call mpi_abort(mpi_comm_world,ierr)

      else

        print*,'PARTICLE READ/WRITE CHECK OK:'
        print*,'PARTICLES READ =',particle_read_count_global
        print*,'PARTICLES WRITTEN =',particle_write_count_global

      endif
    endif

  end subroutine check_particle_numbers

  !------------------------------------------------------------------------------------------!

  subroutine finalise_cubep3m_header



    write(89,pos=1) particle_write_count

    close(89)

  end subroutine

  !------------------------------------------------------------------------------------------!

  subroutine show_min_max_cubep3m

    if(rank.eq.0) print*,'           MIN MAX PARTICLE POSITIONS ON EACH RANK (CUBEP3M):'
    if(rank.eq.0) print*,'           (RANK, min_x,y,z, max_x,y,z)'  



    !  do i = 0, nodes-1
    !      if(rank .eq. i) print*, rank, min_x_cubep3m, min_y_cubep3m,&
    !                        &min_z_cubep3m,max_x_cubep3m,max_y_cubep3m,max_z_cubep3m            
    !  call mpi_barrier(mpi_comm_world,ierr)
    !  enddo

  end subroutine show_min_max_cubep3m

  !------------------------------------------------------------------------------------------!


end program

