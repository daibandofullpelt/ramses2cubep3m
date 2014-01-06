
! Takes output from AHF in chunks and creates
! a full catalogue of halos

! intel fortran: ifort -fpp AHF_chunk_catalogue_maker_v2.0.f90 -o AHF_chunk_cat -DBINARY
! gnu fortran: gfortran -cpp AHF_chunk_catalogue_maker_v2.0.f90 -o AHF_chunk_cat -DGFORTRAN

! options: -DPROFILES         processes profile files too
!          -DAHFBINARY        handles AHF binary format
!          -DPIDS             process PID files too
!          -DCUBEPM_DOMAINS   output halos in cubep3m subdomains
!          -DOUTPUT_BINARY    outputs binary data  
!          -DWITH_BUFFER      writes halos even if they are in the buffer zone of the chunk.

program AHF_chunk_catalogue_maker

  implicit none


  real(4) :: z, box
  integer(4) :: nc, nodes_dim
  integer :: number_of_chunks 
  integer :: num_AHF_files 
  integer(4) :: chunk_min
  integer(4) :: chunk_max
  character(400) :: AHF_data_path, output_path

  ! add any random offset in here (shouldn't need this):

  real, parameter :: x_offset = 0.0
  real, parameter :: y_offset = 0.0
  real, parameter :: z_offset = 0.0

  character(10) :: redshift

  integer(8) :: total_particles_in_chunk, total_halos
  integer(8) :: particles_in_chunk_file
  integer(4) :: number_chunk_files

  integer(4) :: halo_write_flag

  real :: nc_chunk
  real(8) :: box_in_kpc

  real(4) :: buffer_size
  real(8), dimension(3) :: pos_lowest_value, pos_highest_value

  real(4) :: ncc 
  real(4) :: ncc_kpc 

  integer i,j,k,chunk_number,AHF_file, number_of_halos_written

  ! AHF halo variables

  integer(8) :: ID, HostHalo, ID_count
  integer(4) :: NumSubStruct,npart,nbins
  real(4) :: pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir
  real(4) :: R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE
  real(4) :: Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx, M_vir
  real(4) :: Ecy,Ecz,ovdens,fMhires,Ekin,Epot,SurfP,Phi0,cNFW

  ! AHF profile variables

  real(4) :: radius,m_in_r,ovdense,dens,vcirc,sigv,PL_x,P_Ly,P_Lz,PEax,PEay,PEaz
  real(4) :: b_axis,PEbx,PEby,PEbz,c_axis,PEcx,PEcy,PEcz,P_Ekin,P_Epot,Vesc
  integer(4) :: npart_prof, read_check,nbins_prof
  real(4) :: mass_last

  ! AHF PID variables

  integer(8) :: npart_pid, pid, nhalo_pid
  real(4) :: ptype_pid

  ! cubep3m domain variables

  integer :: nodes
  real(4) :: cubep3m_x_min, cubep3m_x_max, cubep3m_y_min, cubep3m_y_max, cubep3m_z_min, cubep3m_z_max
  integer(4), allocatable :: node_coords(:,:)
  integer(4) :: cubep3m_target_file, cubep3m_target_file_previous

  ! chunk details

  real(4) x_chunk_length, y_chunk_length, z_chunk_length
  real(8), dimension(3) :: chunk_offset
  integer(4), allocatable :: chunk_coords(:,:)

  !------------------------------------------------------------------------------------------------------------------------  

  call read_parameters

  ncc = nc/nodes_dim
  ncc_kpc = box*1000.0/nodes_dim
  nodes = nodes_dim*nodes_dim*nodes_dim

  allocate(node_coords(0:nodes-1,3))

  call initialise_chunks

  write(redshift,'(f8.3)') z
  redshift=adjustl(redshift)

  ! open first chunk file to find out how many chunks there are:

  total_halos = 0
  chunk_number = 0

  nc_chunk = nc / number_of_chunks**(1./3.)
  print*, 'nc_chunk = ',nc_chunk
  print*, 'nc = ', nc
  print*, 'boxsize (Mpc) = ', box

  box_in_kpc = box * 1000.0
  print*, 'boxsize (kpc) = ', box_in_kpc
  print*

  pos_lowest_value = 1000000
  pos_highest_value = 0

  print*, "Buffer_size (Mpc/h) = ", buffer_size
  print*, "Giving a buffer length in code units of: ",buffer_size * real(ncc) / box

  buffer_size=2.0*buffer_size*real(ncc)/box  ! note that in this code the buffer is 2 times the buffer in the chunking code

  print*
  print*, 'bounds are:'
  print*, buffer_size, ' code units - ', (nc_chunk + buffer_size),'code units'
  print*, buffer_size*box_in_kpc/float(nc), ' kpc/h - ', (nc_chunk + buffer_size)*box_in_kpc/float(nc),'kpc/h'
  print*, "Note that the bounds should be twice the buffer!"

  ID_count = 0

#ifdef CUBEPM_DOMAINS
  call node_coords_initialise
#endif

  do chunk_number = chunk_min, chunk_max

    chunk_offset(1) = chunk_coords(chunk_number,1)*x_chunk_length - buffer_size          
    chunk_offset(2) = chunk_coords(chunk_number,2)*y_chunk_length - buffer_size          
    chunk_offset(3) = chunk_coords(chunk_number,3)*z_chunk_length - buffer_size  


    number_of_halos_written = 0

    print*, 'offset:', chunk_offset  

#ifndef CUBEPM_DOMAINS
#ifdef WITH_BUFFER
    call open_chunk_catalogue_file_w_buf
#else
    call open_chunk_catalogue_file
#endif

#ifdef PROFILES
    call open_chunk_profiles_file
#endif

#ifdef PIDS
    call open_chunk_pids_file
#endif

#endif ! CUBEPM_DOMAINS


#ifndef OUTPUT_BINARY    
    call write_AHF_halos_header
#endif

    do AHF_file = 0, num_AHF_files - 1

      call open_AHF_file(AHF_file)

#ifdef PROFILES
      call open_AHF_profiles_file(AHF_file)
#endif

#ifdef PIDS
      call open_AHF_pids_file(AHF_file)
#endif

      call process_AHF_halo_data

      call close_AHF_file

#ifdef PROFILES
      call close_AHF_profiles_file
#endif    

#ifdef PIDS
      call close_AHF_pids_file
#endif    


    enddo ! AHF FILE DO LOOP

#ifndef CUBEPM_DOMAINS

#ifdef WITH_BUFFER
    call close_chunk_catalogue_file_w_buf
#else
    call close_chunk_catalogue_file
#endif


#ifdef PROFILES
    call close_chunk_profiles_file
#endif    

#ifdef PIDS
    call close_chunk_pids_file
#endif    

#endif !CUBEPM_DOMAINS

    print*, 'Chunk catalogue written for chunk', chunk_number, 'halos written = ', number_of_halos_written&
	&,'ID count = ',ID_count

    total_halos = total_halos + number_of_halos_written


  enddo ! CHUNK FILE DO LOOP

  print*, "Total halos = ", total_halos
  print*, "Min (x,y,z) = ",pos_lowest_value*float(nc)/box_in_kpc, ' (code units)'
  print*, "Max (x,y,z) = ",pos_highest_value*float(nc)/box_in_kpc, ' (code units)'  
  print*, "Min (x,y,z) = ",pos_lowest_value, ' kpc/h'
  print*, "Max (x,y,z) = ",pos_highest_value, ' kpc/h'

  !  enddo ! LOOP OVER ALL CHECKPOINTS


contains


  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine close_chunk_catalogue_file

      close(20)

  end subroutine close_chunk_catalogue_file

  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine open_chunk_catalogue_file_w_buf

      character(10) :: chunk_string
      character(400) :: ifile

      write(chunk_string,'(i4)') chunk_number
      chunk_string=adjustl(chunk_string)

#ifdef OUTPUT_BINARY 

      ifile=output_path(1:len_trim(output_path))//'/'//&
	  &redshift(1:len_trim(redshift))//'_AHF_halos_chunk_w_buf_'//&
	  &chunk_string(1:len_trim(chunk_string))//'.dat_bin'

#ifdef GFORTRAN		           
      open (unit=20,file=ifile,access='stream')
#endif
#ifdef BINARY
      open (unit=20,file=ifile,form='binary')
#endif                 

#else

      ifile=output_path(1:len_trim(output_path))//'/'//&
	  &redshift(1:len_trim(redshift))//'_AHF_halos_chunk_w_buf_'//&
	  &chunk_string(1:len_trim(chunk_string))//'.dat'

      open (unit=20,file=ifile,form='formatted')

#endif


  end subroutine open_chunk_catalogue_file_w_buf



  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine close_chunk_catalogue_file_w_buf

      close(20)

  end subroutine close_chunk_catalogue_file_w_buf

  !------------------------------------------------------------------------------------------------------------------------ 


  subroutine open_AHF_file(file_number)

      integer file_number, i
      character(10) :: file_string, chunk_string
      character(400) :: ifile

      write(file_string,'(i4)') file_number

      do i=1,4
	if (file_string(i:i) .eq. ' ') then 
	  file_string(i:i)='0'
	end if
      enddo

      write(chunk_string,'(i4)') chunk_number
      chunk_string=adjustl(chunk_string)

#ifdef AHFBINARY
      ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	  &'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	  &//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	  &'.z'//redshift(1:len_trim(redshift))//'.AHF_halos_bin' 
#else
      ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	  &'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	  &//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	  &'.z'//redshift(1:len_trim(redshift))//'.AHF_halos'
#endif


#ifdef AHFBINARY
#ifdef GFORTRAN		           
      open (unit=30,file=ifile,access='stream')
#endif
#ifdef BINARY
      open (unit=30,file=ifile,form='binary')
#endif                 

      print*,'Opened AHF binary file: ',ifile(1:len_trim(ifile))

#else
      open (unit=30,file=ifile,form='formatted')
#endif

  end subroutine open_AHF_file

  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine close_AHF_file

      close(30)

  end subroutine close_AHF_file

  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine read_AHF_header

      integer*4:: number_of_columns ! technically this is an unsigned int but I don't image the column number will
      ! ever be higher than 2147483647

      integer*8 :: number_of_haloes ! technically this is an unsigned int but I don't image the number of haloes in
      ! one AHF output file will ever be over 9223372036854775807

      integer*4 :: one
      character(382) :: line


#ifdef AHFBINARY
      read(30) one
      read(30) number_of_haloes
      read(30) number_of_columns



      print*, 'Number of haloes and columns (HALOES) =', number_of_haloes, number_of_columns

#else
      read(30,'(382a)') line ! read first data line
#endif    

  end subroutine read_AHF_header

  !------------------------------------------------------------------------------------------------------------------------

  subroutine read_AHF_profiles_file_header

      integer*4:: number_of_columns ! technically this is an unsigned int but I don't image the column number will
      ! ever be higher than 2147483647

      integer*8 :: number_of_haloes ! technically this is an unsigned int but I don't image the number of haloes in
      ! one AHF output file will ever be over 9223372036854775807

      integer*4 :: one
      character(342) :: line


#ifdef AHFBINARY
      read(52) one
      read(52) number_of_haloes
      read(52) number_of_columns

      print*, 'Number of haloes and columns (PROFILES) =', number_of_haloes, number_of_columns

#else
      read(52,'(342a)') line ! read first data line
#endif    


  end subroutine read_AHF_profiles_file_header

  !------------------------------------------------------------------------------------------------------------------------ 

  subroutine read_AHF_pids_file_header

      integer*4:: number_of_columns ! technically this is an unsigned int but I don't image the column number will
      ! ever be higher than 2147483647

      integer*8 :: number_of_haloes ! technically this is an unsigned int but I don't image the number of haloes in
      ! one AHF output file will ever be over 9223372036854775807

      integer*4 :: one



#ifdef AHFBINARY
      read(62) one
      read(62) number_of_haloes
      read(62) number_of_columns

      print*, 'Number of haloes and columns (PIDS) =', number_of_haloes, number_of_columns

#else
      read(62,*) number_of_haloes ! read first data line
#endif    

      nhalo_pid = number_of_haloes

  end subroutine read_AHF_pids_file_header

  !------------------------------------------------------------------------------------------------------------------------   

  subroutine process_AHF_halo_data

      integer i,j,k
      logical lopen

      100 format(3i26,es20.5,i16,3f20.8,5f14.2,3f18.6,f14.2,f18.6,f14.2,2f18.6,f16.4,f18.6,&
	  &f16.4,11f18.6,f14.2,i16,f18.6,4es20.5,f18.6)

      101 format(f20.4,i16,es20.5,3f16.2,f20.6,f16.2,3es20.3,11f20.6,2es20.3)

#ifdef AHFBINARY
      call read_AHF_header
#ifdef PROFILES  
      call read_AHF_profiles_file_header
#endif
#ifdef PIDS
      call read_AHF_pids_file_header
#endif
#else
      if(AHF_file .eq. 0) call read_AHF_header
#ifdef PROFILES 
      if(AHF_file .eq. 0) call read_AHF_profiles_file_header
#endif
#ifdef PIDS 
      call read_AHF_pids_file_header
#endif
#endif

      k = 0
      read_check = 1
      mass_last = -1
      cubep3m_target_file_previous = -1

      do ! This is unbounded as it should exit when the file has been read	

	halo_write_flag = 0

	! read halo data

#ifdef AHFBINARY  
	read(30,end=123) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	    &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	    &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW 
#else
	read(30,*,end=123) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	    &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	    &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW 
#endif


	if(pos_x .ge. buffer_size*box_in_kpc/float(nc) .and. pos_x .lt. (nc_chunk + buffer_size)*box_in_kpc/float(nc) .and. &
	    & pos_y .ge. buffer_size*box_in_kpc/float(nc) .and. pos_y .lt. (nc_chunk + buffer_size)*box_in_kpc/float(nc) .and. &
	    & pos_z .ge. buffer_size*box_in_kpc/float(nc) .and. pos_z .lt. (nc_chunk + buffer_size)*box_in_kpc/float(nc)) then

	halo_write_flag = 1 ! this halo is not in the buffer and will be written

	! put halo positions back into global coordinates:

	pos_x = pos_x + (chunk_offset(1)-x_offset)*box_in_kpc/float(nc)
	pos_y = pos_y + (chunk_offset(2)-y_offset)*box_in_kpc/float(nc)
	pos_z = pos_z + (chunk_offset(3)-z_offset)*box_in_kpc/float(nc)      


	! wrap the halos to make sure they are in the box boundarys

	if(pos_x .lt. 0) pos_x = pos_x + box_in_kpc
	if(pos_y .lt. 0) pos_y = pos_y + box_in_kpc
	if(pos_z .lt. 0) pos_z = pos_z + box_in_kpc

	if(pos_x .gt. box_in_kpc) pos_x = pos_x - box_in_kpc
	if(pos_y .gt. box_in_kpc) pos_y = pos_y - box_in_kpc 
	if(pos_z .gt. box_in_kpc) pos_z = pos_z - box_in_kpc

#ifdef CUBEPM_DOMAINS
	! find target file for halo:

	if(int(pos_x/ncc_kpc).lt.nodes_dim .and.&
	    & int(pos_y/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).lt.nodes_dim) then

	cubep3m_target_file = int(pos_x/ncc_kpc)+nodes_dim*int(pos_y/ncc_kpc)+(nodes_dim**2.0)*int(pos_z/ncc_kpc)

      elseif(int(pos_x/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).lt.nodes_dim) then

	cubep3m_target_file = (nodes_dim-1)+nodes_dim*int(pos_y/ncc_kpc)+(nodes_dim**2.0)*int(pos_z/ncc_kpc)

      elseif(int(pos_x/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).lt.nodes_dim) then

	cubep3m_target_file = int(pos_x/ncc_kpc)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*int(pos_z/ncc_kpc)

      elseif(int(pos_x/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).ge.nodes_dim) then

	cubep3m_target_file = int(pos_x/ncc_kpc)+nodes_dim*int(pos_y/ncc_kpc)+(nodes_dim**2.0)*(nodes_dim-1)


      elseif(int(pos_x/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).lt.nodes_dim) then

	cubep3m_target_file = (nodes_dim-1)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*int(pos_z/ncc_kpc)

      elseif(int(pos_x/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).ge.nodes_dim) then

	cubep3m_target_file = (nodes_dim-1)+nodes_dim*int(pos_y/ncc_kpc)+(nodes_dim**2.0)*(nodes_dim-1)

      elseif(int(pos_x/ncc_kpc).lt.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).ge.nodes_dim) then

	cubep3m_target_file = int(pos_x/ncc_kpc)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)

      elseif(int(pos_x/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_y/ncc_kpc).ge.nodes_dim.and.&
	    & int(pos_z/ncc_kpc).ge.nodes_dim) then

	cubep3m_target_file = (nodes_dim-1)+nodes_dim*(nodes_dim-1)+(nodes_dim**2.0)*(nodes_dim-1)


      endif

#endif


    endif


#ifdef PROFILES
    ! read the profile data:


#ifdef CUBEPM_DOMAINS ! open the relevent cubepm_domain profiles file for this halo if it is not already open...

    if(halo_write_flag .eq. 1) then
      inquire(unit=72,opened=lopen)    

      if(cubep3m_target_file_previous .ne. cubep3m_target_file) then

	if(lopen .eqv. .true.) call close_cubepm_domain_profiles_file

	call open_cubepm_domain_profiles_file

      else

	if(lopen .eqv. .false.) call open_cubepm_domain_profiles_file

      endif

    endif  

#endif  

#ifdef AHFBINARY

    read(52,end=125) nbins_prof 

#ifdef OUTPUT_BINARY              
    if(halo_write_flag .eq. 1) then

#ifdef CUBEPM_DOMAINS
      write(72) nbins_prof
#else
      write(53) nbins_prof
#endif
    endif
#endif

    do i = 1, nbins_prof

      read(52,end=125) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	  &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot

      if(halo_write_flag .eq. 1) then

#ifdef OUTPUT_BINARY
#ifdef CUBEPM_DOMAINS
	write(72) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#else          
	write(53) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#endif
#else
#ifdef CUBEPM_DOMAINS
	write(72,101) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#else
	write(53,101) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#endif
#endif
      endif

    enddo



#else ! AHFBINARY ELSE

    do

      if(read_check .eq. 1) then

	read(52,*,end=125) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot

      endif

      if(npart_prof .lt. mass_last) goto 124

      ! if the halo is to be included write the profiles to the profiles file:

      if(halo_write_flag .eq. 1) then

#ifdef OUTPUT_BINARY
#ifdef CUBEPM_DOMAINS
	write(72) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#else
	write(53) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#endif
#else
#ifdef CUBEPM_DOMAINS
	write(72,101) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#else
	write(53,101) radius,npart_prof,m_in_r,ovdense,dens,vcirc,Vesc,sigv,PL_x,P_Ly,P_Lz,b_axis,c_axis,PEax,PEay,PEaz,&
	    &PEbx,PEby,PEbz,PEcx,PEcy,PEcz,P_Ekin,P_Epot
#endif
#endif


      endif

      if(mass_last .eq. npart_prof) read_check = 1

      mass_last=npart_prof

    enddo

#endif !AHFBINARY ELSE    

#endif !AHF PROFILES


    124 read_check = 0 

    125 mass_last=npart_prof


#ifdef PIDS


#ifdef CUBEPM_DOMAINS ! open the relevent cubepm_domain profiles file for this halo if it is not already open...

    if(halo_write_flag .eq. 1) then
      inquire(unit=73,opened=lopen)    

      if(cubep3m_target_file_previous .ne. cubep3m_target_file) then

	if(lopen .eqv. .true.) call close_cubepm_domain_pids_file

	call open_cubepm_domain_pids_file

      else

	if(lopen .eqv. .false.) call open_cubepm_domain_pids_file

      endif

    endif  

#endif  



#ifdef AHFBINARY
    read(62,end=126) npart_pid
#else
    read(62,*,end=126) npart_pid
#endif


    ! if the halo is to be included write the pid to the pids file:

    if(halo_write_flag .eq. 1) then

#ifdef OUTPUT_BINARY
#ifdef CUBEPM_DOMAINS
      write(73) npart_pid
#else
      write(63) npart_pid
#endif
#else
#ifdef CUBEPM_DOMAINS
      write(73,*) npart_pid
#else
      write(63,*) npart_pid
#endif
#endif

    endif


    do i = 1, npart_pid

#ifdef AHFBINARY
      read(62,end=126) pid, ptype_pid
#else
      read(62,*,end=126) pid, ptype_pid
#endif

      ! if the halo is to be included write the pid to the pids file:

      if(halo_write_flag .eq. 1) then

#ifdef OUTPUT_BINARY
#ifdef CUBEPM_DOMAINS
	write(73) pid, ptype_pid
#else
	write(63) pid, ptype_pid
#endif
#else
#ifdef CUBEPM_DOMAINS
	write(73,*) pid, ptype_pid
#else
	write(63,*) pid, ptype_pid
#endif
#endif
      endif

    enddo

#endif !PIDS

    126 if(1 .gt. 2) print*,'Blah' ! just a dummy label for flow control (couldn't think how else to do it)

    ! is the halo in the buffer? if so exclude it...
    ! there is the buffer size around 
    ! the halos that are not in the buffer


#ifndef WITH_BUFFER
    if(halo_write_flag .eq. 1) then
#endif      


      ! Alter the IDs to account for multple files:

      ID = ID + ID_count
      if(HostHalo .ne. 0) HostHalo = HostHalo + ID_count

      ! write the halo into the catalogue

#ifdef CUBEPM_DOMAINS

      inquire(unit=71,opened=lopen)    

      if(cubep3m_target_file_previous .ne. cubep3m_target_file) then

	if(lopen .eqv. .true.) call close_cubepm_domain_halo_file

	call open_cubepm_domain_halo_file

      else

	if(lopen .eqv. .false.) call open_cubepm_domain_halo_file

      endif

#ifdef OUTPUT_BINARY
      write(71) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	  &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	  &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW         
#else      
      write(71,100) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	  &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	  &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW         
#endif        

#else ! CUBEPM_DOMAINS ELSE



#ifdef OUTPUT_BINARY
      write(20) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	  &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	  &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW         
#else      
      write(20,100) ID,HostHalo,NumSubStruct,M_vir,npart,pos_x,pos_y,pos_z,vel_x,vel_y,vel_z,R_vir,  &  
	  &  R_max,r2,mbp_offset,com_offset,V_max,v_esc,sig_V,lambda,lambdaE,Lx,Ly,Lz,b,c,Eax,Eay,Eaz,Ebx,Eby,Ebz,Ecx,& 
	  &  Ecy,Ecz,ovdens,nbins,fMhires,Ekin,Epot,SurfP,Phi0,cNFW         
#endif        

#endif

      number_of_halos_written = number_of_halos_written + 1

      if(pos_x .lt. pos_lowest_value(1)) pos_lowest_value(1) = pos_x
      if(pos_y .lt. pos_lowest_value(2)) pos_lowest_value(2) = pos_y
      if(pos_z .lt. pos_lowest_value(3)) pos_lowest_value(3) = pos_z
      if(pos_x .gt. pos_highest_value(1)) pos_highest_value(1) = pos_x
      if(pos_y .gt. pos_highest_value(2)) pos_highest_value(2) = pos_y
      if(pos_z .gt. pos_highest_value(3)) pos_highest_value(3) = pos_z

#ifndef WITH_BUFFER    
    endif
#endif

    k = k + 1

    cubep3m_target_file_previous = cubep3m_target_file

  enddo

  123 close(30)

  ID_count = number_of_halos_written + ID_count


end subroutine process_AHF_halo_data

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_AHF_profiles_file(file_number)

    integer file_number, i
    character(10) :: file_string, chunk_string
    character(400) :: ifile

    write(file_string,'(i4)') file_number

    do i=1,4
      if (file_string(i:i) .eq. ' ') then 
	file_string(i:i)='0'
      end if
    enddo

    write(chunk_string,'(i4)') chunk_number
    chunk_string=adjustl(chunk_string)

#ifdef AHFBINARY

    ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	&'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	&//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	&'.z'//redshift(1:len_trim(redshift))//'.AHF_profiles_bin'

#ifdef GFORTRAN		           
    open (unit=52,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=52,file=ifile,form='binary')
#endif                 

    print*,'Opened AHF binary file: ',ifile(1:len_trim(ifile))

#else

    ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	&'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	&//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	&'.z'//redshift(1:len_trim(redshift))//'.AHF_profiles'


    open (unit=52,file=ifile,form='formatted')

#endif


end subroutine open_AHF_profiles_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_AHF_pids_file(file_number)

    integer file_number, i
    character(10) :: file_string, chunk_string
    character(400) :: ifile

    write(file_string,'(i4)') file_number

    do i=1,4
      if (file_string(i:i) .eq. ' ') then 
	file_string(i:i)='0'
      end if
    enddo

    write(chunk_string,'(i4)') chunk_number
    chunk_string=adjustl(chunk_string)

#ifdef AHFBINARY

    ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	&'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	&//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	&'.z'//redshift(1:len_trim(redshift))//'.AHF_particles_bin'

#ifdef GFORTRAN		           
    open (unit=62,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=62,file=ifile,form='binary')
#endif                 

    print*,'Opened AHF binary file: ',ifile(1:len_trim(ifile))

#else

    ifile=AHF_data_path(1:len_trim(AHF_data_path))//'/'//'z_'//redshift(1:len_trim(redshift))//&
	&'/chunk_'//chunk_string(1:len_trim(chunk_string))//'/'&
	&//redshift(1:len_trim(redshift))//'xv.'//file_string(1:len_trim(file_string))//&
	&'.z'//redshift(1:len_trim(redshift))//'.AHF_particles'

    open (unit=62,file=ifile,form='formatted')

#endif


end subroutine open_AHF_pids_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_chunk_catalogue_file

    character(10) :: chunk_string
    character(400) :: ifile

    write(chunk_string,'(i4)') chunk_number
    chunk_string=adjustl(chunk_string)



#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//&
	&'_AHF_halos_chunk_'//chunk_string(1:len_trim(chunk_string))//'_halos.dat_bin'


#ifdef GFORTRAN		           
    open (unit=20,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=20,file=ifile,form='binary')
#endif                 


#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//&
	&'_AHF_halos_chunk_'//chunk_string(1:len_trim(chunk_string))//'_halos.dat'

    open (unit=20,file=ifile,form='formatted')

#endif

end subroutine open_chunk_catalogue_file 


!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_chunk_profiles_file

    character(10) :: chunk_string
    character(400) :: ifile

    write(chunk_string,'(i4)') chunk_number
    chunk_string=adjustl(chunk_string)


#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
	&//chunk_string(1:len_trim(chunk_string))//'_profiles.dat_bin'


#ifdef GFORTRAN		           
    open (unit=53,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=53,file=ifile,form='binary')
#endif                 

#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
	&//chunk_string(1:len_trim(chunk_string))//'_profiles.dat'

    open (unit=53,file=ifile,form='formatted')

#endif

end subroutine open_chunk_profiles_file 

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_chunk_pids_file

    character(10) :: chunk_string
    character(400) :: ifile

    write(chunk_string,'(i4)') chunk_number
    chunk_string=adjustl(chunk_string)


#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
	&//chunk_string(1:len_trim(chunk_string))//'_PIDs.dat_bin'


#ifdef GFORTRAN		           
    open (unit=63,file=ifile,access='stream')
#endif
#ifdef BINARY
    open (unit=63,file=ifile,form='binary')
#endif                 

#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_chunk_'&
	&//chunk_string(1:len_trim(chunk_string))//'_PIDs.dat'

    open (unit=63,file=ifile,form='formatted')

#endif

end subroutine open_chunk_pids_file 


!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_AHF_profiles_file

    close(52)

end subroutine close_AHF_profiles_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_chunk_profiles_file

    close(53)

end subroutine close_chunk_profiles_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_AHF_pids_file

    close(62)

end subroutine close_AHF_pids_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_chunk_pids_file

    close(63)

end subroutine close_chunk_pids_file

!------------------------------------------------------------------------------------------------------------------------ 


subroutine write_AHF_halos_header

    write(20,'(804a)') '#ID(1)       & 
	&  hostHalo(2)     &
	&  numSubStruct(3) & 
	&  Mvir(4)         & 
	&  npart(5)        & 
	&  Xc(6)           & 
	&  Yc(7)           & 
	&  Zc(8)           & 
	&  VXc(9)          & 
	&  VYc(10)          & 
	&  VZc(11)          & 
	&  Rvir(12)         & 
	&  Rmax(13)         & 
	&  r2(14)           & 
	&  mbp_offset(15)   & 
	&  com_offset(16)   & 
	&  Vmax(17)         & 
	&  v_esc(18)        & 
	&  sigV(19)         & 
	&  lambda(20)       & 
	&  lambdaE(21)      & 
	&  Lx(22)           & 
	&  Ly(23)           & 
	&  Lz(24)           & 
	&  b(25)            & 
	&  c(26)            & 
	&  Eax(27)          & 
	&  Eay(28)          & 
	&  Eaz(29)          & 
	&  Ebx(30)          & 
	&  Eby(31)          & 
	&  Ebz(32)          & 
	&  Ecx(33)          & 
	&  Ecy(34)          & 
	&  Ecz(35)          & 
	&  ovdens(36)       & 
	&  nbins(37)        & 
	&  fMhires(38)      & 
	&  Ekin(39)         & 
	&  Epot(40)         & 
	&  SurfP(41)        & 
	&  Phi0(42)         & 
	&  cNFW(43)'



end subroutine write_AHF_halos_header

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

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_cubepm_domain_halo_file

    character(10) :: file_string
    character(400) :: ifile

    write(file_string,'(i4)') cubep3m_target_file
    file_string=adjustl(file_string)

#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_halos.dat_bin'

#ifdef GFORTRAN		           
    open (unit=71,file=ifile,access='stream',POSITION = 'APPEND')
#endif
#ifdef BINARY
    open (unit=71,file=ifile,form='binary',POSITION = 'APPEND')
#endif                 

#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_halos.dat'

    open (unit=71,file=ifile,form='formatted',POSITION = 'APPEND')

#endif


end subroutine open_cubepm_domain_halo_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_cubepm_domain_halo_file


    close(71)


end subroutine close_cubepm_domain_halo_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_cubepm_domain_profiles_file

    character(10) :: file_string
    character(400) :: ifile

    write(file_string,'(i4)') cubep3m_target_file
    file_string=adjustl(file_string)

#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_profiles.dat_bin'

#ifdef GFORTRAN		           
    open (unit=72,file=ifile,access='stream',POSITION = 'APPEND')
#endif
#ifdef BINARY
    open (unit=72,file=ifile,form='binary',POSITION = 'APPEND')
#endif                 

#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_profiles.dat'

    open (unit=72,file=ifile,form='formatted',POSITION = 'APPEND')

#endif


end subroutine open_cubepm_domain_profiles_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_cubepm_domain_profiles_file


    close(72)


end subroutine close_cubepm_domain_profiles_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine open_cubepm_domain_pids_file

    character(10) :: file_string
    character(400) :: ifile

    write(file_string,'(i4)') cubep3m_target_file
    file_string=adjustl(file_string)

#ifdef OUTPUT_BINARY          

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_pids.dat_bin'

#ifdef GFORTRAN		           
    open (unit=73,file=ifile,access='stream',POSITION = 'APPEND')
#endif
#ifdef BINARY
    open (unit=73,file=ifile,form='binary',POSITION = 'APPEND')
#endif                 

#else	           

    ifile=output_path(1:len_trim(output_path))//'/'//redshift(1:len_trim(redshift))//'_AHF_halos_cubepm_domain_'&
	&//file_string(1:len_trim(file_string))//'_pids.dat'

    open (unit=73,file=ifile,form='formatted',POSITION = 'APPEND')

#endif


end subroutine open_cubepm_domain_pids_file

!------------------------------------------------------------------------------------------------------------------------ 

subroutine close_cubepm_domain_pids_file


    close(73)


end subroutine close_cubepm_domain_pids_file


!------------------------------------------------------------------------------------------------------------------------ 


subroutine initialise_chunks

    integer i, j, k, n

    allocate(chunk_coords(0:number_of_chunks+1,3))

    ! Set up coordinates of chunks

    x_chunk_length = real(nc) / real(number_of_chunks)**(1./3.)
    y_chunk_length = real(nc) / real(number_of_chunks)**(1./3.)
    z_chunk_length = real(nc) / real(number_of_chunks)**(1./3.)

    do n = 0,number_of_chunks
      do k=1,int((number_of_chunks+0.1)**(1./3.))
	do j=1,int((number_of_chunks+0.1)**(1./3.))
	  do i=1,int((number_of_chunks+0.1)**(1./3.))
	    if (n == (i-1)+(j-1)*int((number_of_chunks+0.1)**(1./3.))+&
		&(k-1)*int((number_of_chunks+0.1)**(1./3.))**2)  then
	    chunk_coords(n,:)=(/(i-1),(j-1),(k-1)/)
	  endif
	enddo
      enddo
    enddo		
  enddo


end subroutine initialise_chunks

!------------------------------------------------------------------------------------------------------------------------ 

subroutine read_parameters

    character(300)  :: line1
    integer(4)	  :: Num_args, IO, param_count
    character(400)  :: arg1, parameterfile

    ! input parameter variables:

    real(4)	  :: redshift_param
    character(100)  :: redshift_file_param
    character(300)  :: AHF_data_path_param, output_path_param
    real(4)	  :: boxlength_param, nc_param, buffer_size_param
    integer(4)	  :: nodes_per_dimension_param, pid_flag_param
    integer(4)	  :: ahf_tasks_param, number_of_chunks_param
    integer(4) :: chunk_start_param, chunk_end_param



    Num_args = IARGC()

    if (Num_args .ne. 1) then
      write(*,'(a)') "COMMAND LINE ERROR! : no parameterfile provided&
      &... include paramterfile name after executable."
      call exit(0)

    else

      CALL GETARG(1,arg1)
      read(arg1,*) parameterfile
      write(*,*) "Reading parameter file: ",&
	  &parameterfile(1:len_trim(parameterfile))

    endif

    open(unit=82,file=parameterfile,form='formatted')

    param_count = 1

    do

      read(82,'(a)',iostat = IO) line1

      if ( IO < 0 ) exit

      if(line1(1:1) .ne. '#') then

	if(param_count .eq. 1) read(line1,*) chunk_start_param
	if(param_count .eq. 2) read(line1,*) chunk_end_param	
	if(param_count .eq. 3) read(line1,*) redshift_param
	if(param_count .eq. 4) read(line1,'(a100)') redshift_file_param
	if(param_count .eq. 5) read(line1,'(a300)') AHF_data_path_param
	if(param_count .eq. 6) read(line1,'(a300)') output_path_param
	if(param_count .eq. 7) read(line1,*) boxlength_param
	if(param_count .eq. 8) read(line1,*) nodes_per_dimension_param
	if(param_count .eq. 9) read(line1,*) nc_param
	if(param_count .eq. 10) read(line1,*) buffer_size_param
	if(param_count .eq. 11) read(line1,*) number_of_chunks_param
	if(param_count .eq. 12) read(line1,*) ahf_tasks_param



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
    print*, "Catalogue code will run on chunks ",chunk_start_param, " to ",chunk_end_param
    print*
    if(redshift_param .eq. -1.0) then
      print*,"Redshift specified from file: ",redshift_file_param(1:len_trim(redshift_file_param))
    endif
    print*, "Redshift: ",redshift_param
    print*, "Particle data path = ",AHF_data_path_param(1:len_trim(AHF_data_path_param))
    print*, "Chunk output path = ",output_path_param(1:len_trim(output_path_param))
    print*, "Box length (Mpc/h) = ", boxlength_param
    print*, "Cubep3m nodes per dimension = ", nodes_per_dimension_param
    print*, "Cubep3m fine cells per dimension = ", nc_param
    print*, "Buffer size (Mpc/h) = ", buffer_size_param
    print*, "Number of chunks in total = ", number_of_chunks_param
    print*, "Number of AHF files per chunk = ",ahf_tasks_param

    ! initialise parameters on task 1:
    chunk_min = chunk_start_param
    chunk_max = chunk_end_param
    z = redshift_param
    AHF_data_path = adjustl(AHF_data_path_param) 
    output_path = adjustl(output_path_param)
    box = boxlength_param
    nodes_dim = nodes_per_dimension_param
    nc = nc_param
    buffer_size = buffer_size_param
    number_of_chunks = number_of_chunks_param
    num_AHF_files = ahf_tasks_param

end subroutine


!------------------------------------------------------------------------------------------------------------------------ 

end program	
