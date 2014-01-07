!Converts RAMSES formatted data to Cubep3m format.
!First converts to Gadget format (following ramses2gadget.f90 (c) 2009 by Timur Doumler) then to Cubep3m (following GADGET2CUBEP3M_v0.11.f90 email: ww60@sussex.ac.uk)
!See respective codes for documentation

!v0.01; Will convert all ramses to gadget, then all gadget to cubep3m (will be inefficient but easy to debug)
!Furute versions: perform ramses -> gadget -> cubep3m file by file

!! Compile with TODO
!! Must use MPI for whole program as needed for GADGET to Cubep3m

!! TODO Test that outputs are the same as vanilla ramses2gadget.f90

!! Dir number hard-coded for now

program ramses2cube3m

	!use mpi

	implicit none

	include 'mpif.h'

	logical :: debug = .true.

		!! Start memory allocation for ramses -> gadget
	! counters
	integer :: i, j, ilevel, idim, idm, istar, ifile, icpu, ivar, ind,col
	integer :: nfiles, ncpus, ndim, twotondim, firstfile, lastfile

	! flags
	logical :: gas, stars, interpolate, refined(8), allrefined

	! cosmological and physical parameters
	real(kind=8)            :: omega_m, omega_l, H0, h, a, gamma, boxsize, dbldummy

	! unit conversion factors and parameters needed for that
	real(kind=8)            :: ramses_length_unit, ramses_dens_unit, ramses_time_unit
	real(kind=8)            :: gadget_length_unit, gadget_mass_unit, gadget_velocity_unit, &
		& gadget_etherm_unit
	real(kind=8)            :: xfactor, vfactor, mfactor, ufactor
	real(kind=8), parameter :: kpc_in_cm = 3.08d21 ! the factor is from ramses: dont use more digits!
	real(kind=8), parameter :: Msun_in_cm = 1.9891d33

	! variables needed for file management
	character(len=128) :: dir_name, part_filename, amr_filename, hydro_filename, info_filename
	character(len=128) :: output_filename
	character(len=5)   :: dir_number,suffix
	character(len=5)   :: output_file_suffix
	integer :: part_file, amr_file, hydro_file, info_file, output_file, scratch_file !i/o unit numbers

	! arrays that will hold the particle data (change integer kind to 8 for more than 2^31 particles!)
	integer (kind=4) :: npart, ngaspart, nramsespart, ndmpart, nstarpart, nsink, intdummy
	integer (kind=4) :: ngaspart_total, ndmpart_total, nstarpart_total
	real(kind=8), dimension(:,:), allocatable :: gaspart_pos, ramsespart_pos, dmpart_pos, starpart_pos
	real(kind=8), dimension(:,:), allocatable :: gaspart_vel, ramsespart_vel, dmpart_vel, starpart_vel
	real(kind=8), dimension(:),   allocatable :: gaspart_m,   ramsespart_m,   dmpart_m,   starpart_m
	integer(kind=4),dimension(:), allocatable :: gaspart_id,  ramsespart_id,  dmpart_id,  starpart_id
	real(kind=8), dimension(:),   allocatable :: gaspart_u
	real(kind=8), dimension(:),   allocatable :: ramsespart_age

	! temporary variables for the particle data
	real(kind=8) :: x, y, z, vx, vy, vz, m, u, m_group
	integer(kind=4) :: id

	! variables needed to read AMR data  
	integer           :: nboundary,nx,ny,nz,nlevelmax,ngrida,nvarh,ix,iy,iz
	real(kind=8)      :: dx
	character(len=80) :: ordering
	integer,      dimension(:,:),   allocatable :: son,ngridfile,ngridlevel,ngridbound
	real(kind=8), dimension(1:8,1:3)            :: xc
	real(kind=8), dimension(1:3)                :: xbound
	real(kind=8), dimension(:,:),   allocatable :: xg
	real(kind=8), dimension(:,:,:), allocatable :: var  

	! little helpers for the gadget header
	integer(kind=4),parameter :: gadget_fillheader(15) = 0
	integer(kind=4)           :: flagsfr, ngaspart_dummy, ndmpart_dummy, nstarpart_dummy

	! mpi variables
	integer         :: mpi_ierr
	integer(kind=4) :: mpi_reduce_buffer, files_per_cpu

	! Initialise MPI
	call mpi_initialise

	! Run ramses2gadget
	call ramses2gadget

	contains

		subroutine ramses2gadget

			! ramses2gadget version 1.0 !
			! ========================= !
			! (c) 2009 by Timur Doumler !

			!! end memory allocation for ramses -> gadget

			if (debug) write(*,*) 'Entering ramses2gadget subroutine'

			! Display hello message
			if (icpu == 0) call hello

			! Read arguments and set directory name and flags for gas/interpolation mode
			call read_args(icpu, dir_name, gas, interpolate)

			! Display status message
			if (icpu == 0) then
				if (.not.gas) write(*,*) 'Running in dark matter mode.'
				if (gas.and.(.not.interpolate)) write (*,*) 'Running with hydro in one-particle-per-cell mode.'
				if (gas.and.interpolate) write (*,*) 'Running with hydro in interpolation mode.'
			endif

			! Open the info_....txt file in order to read cosmological parameters

			! Construct filename, verify existence and assign I/O unit number
			! dir_number=dir_name( index(dir_name,'output_')+7 : index(dir_name,'output_')+13 ) 
			dir_number='00001'
			info_filename   = trim(dir_name) // '/info_'   // trim(dir_number) // '.txt'
			call inquire_file (info_filename)
			info_file   = 988

			! in case we have gas, we need to open some more files and gather some more parameters...
			if (gas) then 
				!Open cpu1 file to verify existence and assign I/O unit number
				amr_filename    = trim(dir_name) // '/amr_'    // trim(dir_number) // '.out00001'
				hydro_filename  = trim(dir_name) // '/hydro_'  // trim(dir_number) // '.out00001'
				call inquire_file (amr_filename)
				call inquire_file (hydro_filename)  
				amr_file    = 989
				hydro_file  = 990  
			end if

			!Open and read info....txt file
			open(unit=info_file, file=info_filename, form='formatted', status='old', action='read', position='rewind')
			read(info_file,'(13X,I11)')    nfiles
			read(info_file,'(13X,I11)')    ndim
			read(info_file,*)              ! skip levelmin
			read(info_file,*)              ! skip levelmax
			read(info_file,*)              ! skip ngridmax
			read(info_file,*)              ! skip nstep_coarse
			read(info_file,*)              ! skip empty line
			read(info_file,'(13X,E23.15)') boxsize ! called 'boxlen' in the file
			read(info_file,*)              ! skip time
			read(info_file,'(13X,E23.15)') a ! called 'aexp' in the file
			read(info_file,'(13X,E23.15)') H0
			read(info_file,'(13X,E23.15)') omega_m
			read(info_file,'(13X,E23.15)') omega_l
			read(info_file,*)              ! skip omega_k
			read(info_file,*)              ! skip omega_b
			read(info_file,'(13X,E23.15)') ramses_length_unit
			read(info_file,'(13X,E23.15)') ramses_dens_unit
			read(info_file,'(13X,E23.15)') ramses_time_unit
			read(info_file,*)              ! skip empty line
			read(info_file,'(14X,A80)')    ordering
			close(info_file)

			if (gas) then
				! open first amr file (.out00001) and read some required quantities from the header
				open(unit=amr_file, file=amr_filename,status='old', form='unformatted', action='read', position='rewind')
				read(amr_file)  ! skip nfiles
				read(amr_file)  ! skip ndim
				read(amr_file)  nx, ny, nz
				read(amr_file)  nlevelmax
				read(amr_file)  ! skip ngridmax
				read(amr_file)  nboundary
				close(amr_file)

				! open first hydro file (.out00001) and read gamma, this is needed for unit conversion
				open(unit=hydro_file, file=hydro_filename, status='old', form='unformatted', action='read', position='rewind')
				read(hydro_file) ! skip nfiles
				read(hydro_file) ! skip nvar
				read(hydro_file) ! skip ndim
				read(hydro_file) ! skip nlevelmax
				read(hydro_file) ! skip nboundary
				read(hydro_file) gamma
				close(hydro_file)
			end if

			! display status message
			if (icpu==0) write(*,'(A,I4,A,I4,A)') ' Working with', ncpus, ' cpu(s) on', nfiles, ' fileset(s).'

		end subroutine ramses2gadget

end program ramses2cube3m

!! Start ramses2gadget subroutines

!------------------------------------------------------------------------------------------------------------!

subroutine mpi_initialise

	!! ramses2gadget version for the time being

	!call MPI_INIT(mpi_ierr)
	!call MPI_COMM_RANK(MPI_COMM_WORLD, icpu, mpi_ierr)
	!call MPI_COMM_SIZE(MPI_COMM_WORLD, ncpus, mpi_ierr)

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

end subroutine mpi_initialise

!------------------------------------------------------------------------------------------------------------!

subroutine hello

	write(*,*) 'ramses2gadget version 1.0 -- (c) 2009 by Timur Doumler'
	write(*,*) '======================================================'
	write(*,*) 'Compiled with MPI.'

end subroutine hello

!------------------------------------------------------------------------------------------------------------!

subroutine read_args(icpu, dir_name, gas, interpolate)

	! This subroutine handles program arguments for ramses2gadget and sets mode flags

	implicit none

	integer,           intent(in)  :: icpu
	integer                        :: n, iargc
	character(len=2)               :: arg1
	character(len=128)             :: arg2
	character(len=128),intent(out) :: dir_name
	logical,           intent(out) :: gas, interpolate

	! retrieve arguments    
	n = iargc()
	if (n /= 2) then
		call wrongcall(icpu)
	else
		! first argument tells us in which mode to run the code
		call getarg(1,arg1)
		if (arg1 == '-g') then
			gas = .true.
			interpolate = .false.
		else if (arg1 == '-i') then
			gas = .true.
			interpolate = .true.
		else if (arg1 == '-d') then
			gas = .false.
			interpolate = .false.
		else
			call wrongcall
		end if

		! second argument tells us the input directory path
		call getarg(2,arg2)
		dir_name = trim(arg2)  
	end if

	return  

end subroutine read_args

!------------------------------------------------------------------------------------------------------------!

subroutine inquire_file(filename)

	! this checks if a file is present

	implicit none

	character(LEN=128), intent(in)::filename
	logical :: ok

	inquire(file=filename, exist=ok) 
	if ( .not. ok ) then
		 write (*,*) 'Error in input files: ' // trim(filename)// ' not found.'
		 call terminate(6)
	endif

	return
end subroutine inquire_file

!------------------------------------------------------------------------------------------------------------!

!! End ramses2gadget subroutines

