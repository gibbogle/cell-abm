module cells_mod
use global
use boundary
use chemokine
use cellstate
use winsock  
use csim_abm

IMPLICIT NONE

contains 

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine CellML_load_file(file_array, buflen, nvar, ncomp) BIND(C,NAME='CellML_load_file')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_load_file
character(c_char),dimension(*) :: file_array
integer(c_int),value :: buflen
integer(c_int) :: nvar, ncomp
integer :: k

cellmlfile = ' '
do k = 1,buflen
	cellmlfile(k:k) = file_array(k)
enddo
call loadCellML
nvar = nvariables
ncomp = ncomponents
end subroutine

subroutine loadCellML
character(64) :: ID, comp_name
logical :: hit
integer :: nlen, k, i, clen

if (allocated(IDlist)) deallocate(IDlist)
if (allocated(csim_state)) deallocate(csim_state)
if (allocated(csim_savestate)) deallocate(csim_savestate)
if (allocated(component)) deallocate(component)
!cellmlfile = "file:///E:/testing/CSim/VPH-MIP_Case_Study/experiments/periodic-stimulus.xml"C
call setupCellml(cellmlfile)
call getNvariables(nvariables)
write(nflog,*) 'nvariables: ',nvariables
allocate(IDlist(0:nvariables-1))
allocate(csim_state(0:nvariables-1))
allocate(csim_savestate(0:nvariables-1))
allocate(component(0:nvariables-1))

ncomponents = 0
do k = 0,nvariables-1
	call getID(k,ID,clen)
!	IDlist(k)%name = ID(1:clen)
	i = index(ID,'.')
	comp_name = ID(1:i-1)
	IDlist(k)%comp_name = comp_name
	IDlist(k)%var_name = ID(i+1:clen)
	if (ncomponents == 0) then
		component(ncomponents) = comp_name
		IDlist(k)%comp_index = ncomponents
		ncomponents = ncomponents + 1
	else
		hit = .false.
		do i = 0,ncomponents-1
			if (comp_name == component(i)) then
				hit = .true.
				IDlist(k)%comp_index = i
				exit
			endif
		enddo
		if (.not.hit) then
			component(ncomponents) = comp_name
			IDlist(k)%comp_index = ncomponents
			ncomponents = ncomponents + 1
		endif
	endif		
enddo
write(nflog,*) 'Components: '
write(nflog,'(a)') component(1:ncomponents)
write(nflog,*) 'Variables: '
do k = 0,nvariables-1
	write(nflog,'(i4,2x,2a30,i3)') k,IDlist(k)%comp_name(1:30),IDlist(k)%var_name(1:30),IDlist(k)%comp_index
enddo

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine setCellStateVectors
integer :: k
do k = 1,nlist
	call getState(cell_list(k)%cellml_state)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine freeCellStateVectors
integer :: k
do k = 1,nlist
	deallocate(cell_list(k)%cellml_state)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! This subroutine is called to initialize a simulation run.
! ncpu = the number of processors to use 
! infile = file with the input data
! outfile = file to hold the output
!-----------------------------------------------------------------------------------------
subroutine Setup(ncpu,infile,outfile,ok)
integer :: ncpu
character*(*) :: infile, outfile
logical :: ok
character*(64) :: msg
integer :: error

ok = .true.
initialized = .false.
par_zig_init = .false.

inputfile = infile
outputfile = outfile
call logger("ReadCellParams")
call ReadCellParams(ok)
if (.not.ok) return
call logger("did ReadCellParams")

if (ncpu == 0) then
	ncpu = ncpu_input
endif
Mnodes = ncpu
write(logmsg,*) 'ncpu: ',ncpu 
call logger(logmsg)

#if defined(OPENMP) || defined(_OPENMP)
    call logger("OPENMP defined")
    call omp_initialisation(ok)
    if (.not.ok) return
#else
    call logger("OPENMP NOT defined")
    if (Mnodes > 1) then
        write(logmsg,'(a)') 'No OpenMP, using one thread only'
        call logger(logmsg)
        Mnodes = 1
    endif
#endif

call ArrayInitialisation(ok)
if (.not.ok) return
call logger('did ArrayInitialisation')

call SetupChemo

if (.not.use_TCP) then
	call loadCellML
endif

call PlaceCells(ok)
call SetRadius(Nsites)
write(logmsg,*) 'did PlaceCells: Ncells: ',Ncells,Radius
call logger(logmsg)
if (.not.ok) return

!call test_get_path

call CreateBdryList

call SetupODEdiff
call InitConcs
call SetupMedium
call AdjustMM

!call setCellStateVectors

t_simulation = 0
istep = 0
write(logmsg,'(a,i6)') 'Startup procedures have been executed: initial T cell count: ',Ncells0
call logger(logmsg)

end subroutine

!----------------------------------------------------------------------------------------- 
! Initialise medium concentrations, etc.
!-----------------------------------------------------------------------------------------
subroutine SetupMedium
integer :: ichemo
real(REAL_KIND) :: V, V0, R1

call SetRadius(Nsites)
R1 = Radius*DELTA_X			! cm
V0 = medium_volume0			! cm3
V = V0 - (4./3.)*PI*R1**3	! cm3
!write(*,'(a,3f10.3)') 'SetupMedium: ',R1,V0,V
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_Cext = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_Cbnd = chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_M = V*chemo(ichemo)%bdry_conc
	chemo(ichemo)%medium_U = 0
!	write(*,'(i4,2e12.4)') ichemo,chemo(ichemo)%medium_Cext,chemo(ichemo)%medium_M
enddo
end subroutine

!----------------------------------------------------------------------------------------- 
!-----------------------------------------------------------------------------------------
subroutine omp_initialisation(ok)
logical :: ok
integer :: npr, nth

ok = .true.
!if (Mnodes == 1) return
#if defined(OPENMP) || defined(_OPENMP)
write(logmsg,'(a,i2)') 'Requested Mnodes: ',Mnodes
call logger(logmsg)
npr = omp_get_num_procs()
write(logmsg,'(a,i2)') 'Machine processors: ',npr
call logger(logmsg)

nth = omp_get_max_threads()
write(logmsg,'(a,i2)') 'Max threads available: ',nth
call logger(logmsg)
if (nth < Mnodes) then
    Mnodes = nth
    write(logmsg,'(a,i2)') 'Setting Mnodes = max thread count: ',nth
	call logger(logmsg)
endif

call omp_set_num_threads(Mnodes)
!$omp parallel
nth = omp_get_num_threads()
write(logmsg,*) 'Threads, max: ',nth,omp_get_max_threads()
call logger(logmsg)
!$omp end parallel
#endif

call logger('did omp_initialisation')
!call test_omp1

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp
integer, parameter :: n = 10
integer :: i

integer :: sum1, sum2
integer, allocatable :: y1(:)
integer :: y2(n)

allocate(y1(n))
y1 = 1
y2 = 1

sum1 = 0
sum2 = 0
!$omp parallel do
do i = 1,n
	sum1 = sum1 + y1(i)
	sum2 = sum2 + y2(i)
enddo
!$omp end parallel do
write(*,*) 'sum1: ',sum1
write(*,*) 'sum2: ',sum2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_omp1
integer :: n = 1000000
!integer :: ncpu = 2
integer :: i, k, ith

real(REAL_KIND) :: x, y, z

real(REAL_KIND) :: sum1, sum2
real(REAL_KIND), allocatable :: y1(:), ysum(:)
!integer :: y2(100)

allocate(omp_x(n))
allocate(omp_y(n))
allocate(omp_z(n))
omp_x = 1.23456
omp_y = 9.87654

!call omp_set_num_threads(ncpu)

allocate(y1(n))
allocate(ysum(0:n-1))
y1 = 1
!y2 = 1
ysum = 0

sum1 = 0
sum2 = 0
!$omp parallel do private(x, y, z, k)
do i = 1,n
	z = 0
	omp_x(i) = omp_x(i)/i
	omp_y(i) = omp_y(i)/i + 1.0
	if (i < 4) then
		x = omp_x(i)**i
	elseif (i < 8) then
		x = 1.0/omp_x(i)**i
	else 
		x = 1./(1 + omp_x(i))
	endif
	omp_x(i) = x
	y = sqrt(i*omp_y(i))
	y = omp_y(i)
	do k = 1,1000
		y = 0.5*y
		call omp_sub(i,x,y)
		z = z + 1/omp_z(i)
	enddo
!	sum1 = sum1 + z
	ith = omp_get_thread_num()
	ysum(ith) = ysum(ith) + z
enddo
!$omp end parallel do
!write(*,*) 'sum1: ',sum1
write(*,*) 'ysum: ',sum(ysum)
stop
end subroutine

subroutine omp_sub(i,x,y)
integer :: i
real(REAL_KIND) :: x, y
omp_z(i) = x/y + omp_x(i)/omp_y(i)
omp_y(i) = y
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine ArrayInitialisation(ok)
logical :: ok
integer :: x,y,z,k, ichemo
integer :: MAXX, z1, z2, nc0, inflow
integer :: cog_size
real(REAL_KIND) :: d, rr(3)

ok = .false.
call RngInitialisation
call logger("did RngInitialisation")

! These are deallocated here instead of in subroutine wrapup so that when a simulation run ends 
! it will still be possible to view the cell distributions and chemokine concentration fields.
if (allocated(occupancy)) deallocate(occupancy)
if (allocated(cell_list)) deallocate(cell_list)
if (allocated(allstate)) deallocate(allstate)
if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
call logger('did deallocation')

!nsteps_per_min = 1.0/DELTA_T
NY = NX
NZ = NX
ngaps = 0
max_ngaps = NY*NZ
nlist = 0

allocate(zoffset(0:2*Mnodes))
allocate(zdomain(NZ))
x0 = (NX + 1.0)/2.        ! global value
y0 = (NY + 1.0)/2.
z0 = (NZ + 1.0)/2.
Centre = (/x0,y0,z0/)   ! now, actually the global centre (units = grids)
call SetRadius(initial_count)
write(logmsg,*) 'Initial radius, nc0, max_nlist: ',Radius, initial_count, max_nlist
call logger(logmsg)

allocate(cell_list(max_nlist))
allocate(occupancy(NX,NY,NZ))
!allocate(gaplist(max_ngaps))

call make_jumpvec

lastNcells = 0
nadd_sites = 0

ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RngInitialisation
integer, allocatable :: zig_seed(:)
integer :: i
integer :: npar, grainsize = 32

npar = Mnodes
write(logmsg,*) 'npar = ',npar,seed
call logger(logmsg)
allocate(zig_seed(0:npar-1))
do i = 0,npar-1
    zig_seed(i) = seed(1)*seed(2)*(i+1)
enddo
call par_zigset(npar,zig_seed,grainsize)
par_zig_init = .true.
end subroutine

!----------------------------------------------------------------------------------------1123
!----------------------------------------------------------------------------------------
subroutine ReadCellParams(ok)
logical :: ok
integer :: i, imetab, itestcase, Nmm3, ichemo, iuse_extra, iuse_relax, iuse_par_relax
integer :: iuse_oxygen, iuse_glucose, iuse_tracer, iV_depend, iV_random
integer :: ictype, idisplay
real(REAL_KIND) :: days, bdry_conc, percent
real(REAL_KIND) :: sigma, DXmm
character*(128) :: testfile

ok = .true.
chemo(:)%used = .false.

open(nfcell,file=inputfile,status='old')

read(nfcell,*) NX							! size of grid
read(nfcell,*) initial_count				! initial number of tumour cells
read(nfcell,*) divide_time_median
read(nfcell,*) divide_time_shape
read(nfcell,*) iV_depend
read(nfcell,*) iV_random
read(nfcell,*) days							! number of days to simulate
read(nfcell,*) DELTA_T						! time step size (sec)
read(nfcell,*) NT_CONC						! number of subdivisions of DELTA_T for diffusion computation
read(nfcell,*) Nmm3							! number of cells/mm^3
read(nfcell,*) fluid_fraction				! fraction of the (non-necrotic) tumour that is fluid
read(nfcell,*) medium_volume0				! initial total volume (medium + spheroid) (cm^3)
read(nfcell,*) d_layer						! thickness of the unstirred layer around the spheroid (cm)
read(nfcell,*) Vdivide0						! nominal cell volume multiple for division
read(nfcell,*) dVdivide						! variation about nominal divide volume
read(nfcell,*) MM_THRESHOLD					! O2 concentration threshold Michaelis-Menten "soft-landing" (uM)
read(nfcell,*) ANOXIA_THRESHOLD			    ! O2 threshold for anoxia (uM)
read(nfcell,*) anoxia_tag_hours				! hypoxic time leading to tagging to die by anoxia (h)
read(nfcell,*) anoxia_death_hours			! time after tagging to death by anoxia (h)
read(nfcell,*) itestcase                    ! test case to simulate
read(nfcell,*) seed(1)						! seed vector(1) for the RNGs
read(nfcell,*) seed(2)						! seed vector(2) for the RNGs
read(nfcell,*) ncpu_input					! for GUI just a placeholder for ncpu, used only when execute parameter ncpu = 0
read(nfcell,*) Ncelltypes					! maximum number of cell types in the spheroid
do ictype = 1,Ncelltypes
	read(nfcell,*) percent
	celltype_fraction(ictype) = percent/100
!	read(nfcell,*) idisplay
!	celltype_display(ictype) = (idisplay == 1)
enddo
read(nfcell,*) NT_GUI_OUT					! interval between GUI outputs (timesteps)
read(nfcell,*) show_progeny                 ! if != 0, the number of the cell to show descendents of
read(nfcell,*) iuse_oxygen		! chemo(OXYGEN)%used
read(nfcell,*) chemo(OXYGEN)%diff_coef
read(nfcell,*) chemo(OXYGEN)%medium_diff_coef
read(nfcell,*) chemo(OXYGEN)%membrane_diff
read(nfcell,*) chemo(OXYGEN)%bdry_conc
read(nfcell,*) chemo(OXYGEN)%max_cell_rate
read(nfcell,*) chemo(OXYGEN)%MM_C0
read(nfcell,*) chemo(OXYGEN)%Hill_N
read(nfcell,*) iuse_glucose		!chemo(GLUCOSE)%used
read(nfcell,*) chemo(GLUCOSE)%diff_coef
read(nfcell,*) chemo(GLUCOSE)%medium_diff_coef
read(nfcell,*) chemo(GLUCOSE)%membrane_diff
read(nfcell,*) chemo(GLUCOSE)%bdry_conc
read(nfcell,*) chemo(GLUCOSE)%max_cell_rate
read(nfcell,*) chemo(GLUCOSE)%MM_C0
read(nfcell,*) chemo(GLUCOSE)%Hill_N
read(nfcell,*) iuse_tracer		!chemo(TRACER)%used
read(nfcell,*) chemo(TRACER)%diff_coef
read(nfcell,*) chemo(TRACER)%medium_diff_coef
read(nfcell,*) chemo(TRACER)%membrane_diff
read(nfcell,*) chemo(TRACER)%bdry_conc
read(nfcell,*) chemo(TRACER)%max_cell_rate
read(nfcell,*) chemo(TRACER)%MM_C0
read(nfcell,*) chemo(TRACER)%Hill_N

read(nfcell,*) O2cutoff(1)
read(nfcell,*) O2cutoff(2)
read(nfcell,*) O2cutoff(3)
read(nfcell,*) growthcutoff(1)
read(nfcell,*) growthcutoff(2)
read(nfcell,*) growthcutoff(3)
read(nfcell,*) spcrad_value
read(nfcell,*) iuse_extra
read(nfcell,*) iuse_relax
read(nfcell,*) iuse_par_relax
read(nfcell,*) signal_decay_coef
read(nfcell,'(a)') cellmlfile
close(nfcell)

! Put string-terminating NULL at first " "
cellmlfile = adjustl(cellmlfile)
i = index(cellmlfile,' ')
cellmlfile(i:i) = char(0)

if (chemo(OXYGEN)%Hill_N /= 1 .and. chemo(OXYGEN)%Hill_N /= 2) then
	call logger('Error: OXYGEN_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
if (chemo(GLUCOSE)%Hill_N /= 1 .and. chemo(GLUCOSE)%Hill_N /= 2) then
	call logger('Error: GLUCOSE_HILL_N must be 1 or 2')
	ok = .false.
	return
endif
MM_THRESHOLD = MM_THRESHOLD/1000					! uM -> mM
ANOXIA_THRESHOLD = ANOXIA_THRESHOLD/1000			! uM -> mM
O2cutoff = O2cutoff/1000							! uM -> mM
relax = (iuse_relax == 1)
use_parallel = (iuse_par_relax == 1)
chemo(OXYGEN)%used = (iuse_oxygen == 1)
chemo(GLUCOSE)%used = (iuse_glucose == 1)
chemo(TRACER)%used = (iuse_tracer == 1)
chemo(OXYGEN)%MM_C0 = chemo(OXYGEN)%MM_C0/1000		! uM -> mM
chemo(GLUCOSE)%MM_C0 = chemo(GLUCOSE)%MM_C0/1000	! uM -> mM
divide_dist%class = LOGNORMAL_DIST
divide_time_median = 60*60*divide_time_median		! hours -> seconds
sigma = log(divide_time_shape)
!divide_dist%p1 = log(divide_time_mean/exp(sigma*sigma/2))	
divide_dist%p1 = log(divide_time_median)	
divide_dist%p2 = sigma
divide_time_mean = exp(divide_dist%p1 + 0.5*divide_dist%p2**2)	! mean = median.exp(sigma^2/2)
write(logmsg,'(a,2e12.4)') 'shape, sigma: ',divide_time_shape,sigma
call logger(logmsg)
write(logmsg,'(a,2e12.4)') 'Median, mean divide time: ',divide_time_median/3600,divide_time_mean/3600
call logger(logmsg)
use_V_dependence = (iV_depend == 1)
randomise_initial_volume = (iV_random == 1)
use_extracellular_O2 = (iuse_extra == 1)
t_anoxic_limit = 60*60*anoxia_tag_hours				! hours -> seconds
anoxia_death_delay = 60*60*anoxia_death_hours		! hours -> seconds
DXmm = 1.0/(Nmm3**(1./3))
DELTA_X = DXmm/10									! mm -> cm
Vsite_cm3 = DELTA_X*DELTA_X*DELTA_X					! total site volume (cm^3)
Vextra_cm3 = fluid_fraction*Vsite_cm3				! extracellular volume in a site (cm^3)
cell_radius = (3*(1-fluid_fraction)*Vsite_cm3/(4*PI))**(1./3.)
! In a well-oxygenated tumour the average cell fractional volume is intermediate between vdivide0/2 and vdivide0.
! We assume that 0.75*vdivide0*Vcell_cm3 = (1 - fluid_fraction)*Vsite_cm3
Vcell_cm3 = (1 - fluid_fraction)*Vsite_cm3/(0.75*vdivide0)	! nominal cell volume (cm^3) (corresponds to %volume = 1)
															! fractional volume: cell%volume = (cell volume in cm^3)/Vcell_cm3, 
															! actual volume (cm^3) = Vcell_cm3*%cellvolume

write(logmsg,'(a,3e12.4)') 'DELTA_X, cell_radius: ',DELTA_X,cell_radius
call logger(logmsg)
write(logmsg,'(a,4e12.4)') 'Volumes: site, extra, cell (average, base): ',Vsite_cm3, Vextra_cm3, Vsite_cm3-Vextra_cm3, Vcell_cm3
call logger(logmsg)

! Setup test_case
test_case = .false.
if (itestcase /= 0) then
    test_case(itestcase) = .true.
endif

if (mod(NX,2) /= 0) NX = NX+1					! ensure that NX is even
open(nfout,file=outputfile,status='replace')
write(logmsg,*) 'Opened nfout: ',outputfile
call logger(logmsg)

Nsteps = days*24*60*60/DELTA_T		! DELTA_T in seconds
write(logmsg,'(a,2i6,f6.0)') 'nsteps, NT_CONC, DELTA_T: ',nsteps,NT_CONC,DELTA_T
call logger(logmsg)

ok = .true.

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine PlaceCells(ok)
logical :: ok
integer :: x, y, z, k, site(3), ichemo, kpar=0
real(REAL_KIND) :: r2lim,r2,rad(3)
real(REAL_KIND) :: R, tpast, tdiv, c_rate, r_mean
real(REAL_KIND),allocatable :: state(:)

allocate(state(0:nvariables-1))
call getState(state)
occupancy(:,:,:)%indx(1) = 0
occupancy(:,:,:)%indx(2) = 0
Radius = max(Radius,1.)
r2lim = Radius*Radius
lastID = 0
k = 0
do x = 1,NX
	do y = 1,NY
		do z = 1,NZ
			rad = (/x-x0,y-y0,z-z0/)
			r2 = dot_product(rad,rad)
			if (r2 < r2lim .and. k < initial_count) then
				k = k+1
				lastID = lastID + 1
				site = (/x,y,z/)
				cell_list(k)%ID = lastID
				cell_list(k)%celltype = random_choice(celltype_fraction,Ncelltypes,kpar)
				cell_list(k)%site = site
				cell_list(k)%state = 1
				cell_list(k)%exists = .true.
				cell_list(k)%active = .true.
!				do
!					R = par_uni(kpar)
!					tpast = -R*divide_time_median
!					tdiv = DivideTime()
!					if (tdiv + tpast > 0) exit
!				enddo
!				cell_list(k)%divide_volume = Vdivide0
				R = par_uni(kpar)
				cell_list(k)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
				R = par_uni(kpar)
				if (randomise_initial_volume) then
					cell_list(k)%volume = cell_list(k)%divide_volume*0.5*(1 + R)
				else
					cell_list(k)%volume = 1.0
				endif
				call getGrowthRateParameters(c_rate,r_mean)
				cell_list(k)%c_rate = c_rate
				cell_list(k)%r_mean = r_mean
				cell_list(k)%t_divide_last = 0		! not used
				cell_list(k)%t_hypoxic = 0
				cell_list(k)%conc = 0
				cell_list(k)%conc(OXYGEN) = chemo(OXYGEN)%bdry_conc
				cell_list(k)%conc(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
				cell_list(k)%conc(TRACER) = chemo(TRACER)%bdry_conc
				cell_list(k)%M = 0
				allocate(cell_list(k)%cellml_state(0:nvariables-1))
				cell_list(k)%cellml_state = state
				! Now randomise the growth stage
				R = par_uni(kpar)
				cell_list(k)%cellml_state(1) = 1.0 + (2.2 - 1.0)*R	! growth
				cell_list(k)%cellml_state(2) = 1.0 + (2.0 - 1.0)*R	! size
!				cell_list(k)%cellml_state(1) = 1.5 ! growth
!				cell_list(k)%cellml_state(2) = 1.5 ! size
				occupancy(x,y,z)%indx(1) = k
				if (x == NX/2 .and. y == NY/2 .and. z == NZ/2) then
					idbug = k
					write(nfout,*) 'Mid-blob cell: idbug: ',idbug, x,y,z
				endif
			else
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		enddo
	enddo
enddo
do ichemo = 1,MAX_CHEMO
    occupancy(:,:,:)%C(ichemo) = 0
enddo
occupancy(:,:,:)%C(OXYGEN) = chemo(OXYGEN)%bdry_conc
occupancy(:,:,:)%C(GLUCOSE) = chemo(GLUCOSE)%bdry_conc
occupancy(:,:,:)%C(TRACER) = chemo(TRACER)%bdry_conc
nlist = k
Nsites = k
Ncells = k
Ncells0 = Ncells
Nreuse = 0	
deallocate(state)
ok = .true.
write(logmsg,*) 'idbug: ',idbug
call logger(logmsg)
end subroutine


!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine get_dimensions(NX_dim, NY_dim, NZ_dim, nsteps_dim, deltat, maxchemo, cused, dfraction) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_dimensions
use, intrinsic :: iso_c_binding
integer(c_int) :: NX_dim,NY_dim,NZ_dim,nsteps_dim, maxchemo
real(c_double) :: deltat, dfraction
logical(c_bool) :: cused(*)
integer :: ichemo

NX_dim = NX
NY_dim = NY
NZ_dim = NZ
nsteps_dim = nsteps
deltat = DELTA_T
maxchemo = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
	cused(ichemo) = chemo(ichemo)%used
enddo
dfraction = 2*cell_radius/DELTA_X
end subroutine

!----------------------------------------------------------------------------------
!----------------------------------------------------------------------------------
subroutine InitConcs
integer :: nextra, ic, ichemo, kcell, site(3), ntvars

write(logmsg,*) 'InitConcs: ',nchemo
call logger(logmsg)
allocate(allstate(MAX_VARS,MAX_CHEMO))
allocate(work_rkc(8+5*MAX_VARS))

do kcell = 1,Ncells
    site = cell_list(kcell)%site
    do ic = 1,nchemo
	    ichemo = chemomap(ic)
        occupancy(site(1),site(2),site(3))%C(ic) = chemo(ichemo)%bdry_conc
    enddo
enddo
ntvars = ODEdiff%nextra + ODEdiff%nintra
do ic = 1,nchemo
	ichemo = chemomap(ic)
	allstate(1:ntvars,ichemo) = chemo(ichemo)%bdry_conc
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Advance simulation through one big time step (DELTA_T)
! The concentration fields are first solved through NT_CONC subdivisions of DELTA_T,
! then the cell states are updated, which includes cell death and division.
! On either death or division cell positions are adjusted, and site concentrations
! (if necessary), and the ODE solver mappings in ODEdiff.
!-----------------------------------------------------------------------------------------
subroutine simulate_step(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: simulate_step  
use, intrinsic :: iso_c_binding
integer(c_int) :: res
integer :: kcell, site(3), hour, nthour, kpar=0
real(REAL_KIND) :: r(3), rmax, tstart, dt
!integer, parameter :: NT_CONC = 6
integer :: nchemo, i
logical :: ok
!real(REAL_KIND) :: interval_dt, cellml_dt

if (Ncells == 0) then
    res = 2
	write(logmsg,'(a,i2)') 'Nells = 0: ',res
	call logger(logmsg)
    return
endif
nthour = 3600/DELTA_T
dt = DELTA_T/NT_CONC
if (mod(istep,nthour) == 0) then
	write(logmsg,*) 'istep, hour: ',istep,istep/nthour,nlist,ncells,nsites-ncells
	call logger(logmsg)
	if (bdry_changed) then
		call UpdateBdrylist
	endif
!	write(logmsg,'(a,2e12.3,i6)') 'Oxygen U, Cbnd, Nbnd: ', chemo(OXYGEN)%medium_U,chemo(OXYGEN)%medium_Cbnd, Nbnd
!	call logger(logmsg)
endif
istep = istep + 1
t_simulation = (istep-1)*DELTA_T	! seconds

!cellml_dt = 0.05
!interval_dt = DELTA_T/3600.	! Cell_cycle uses time unit = hour
!do kcell = 1,nlist
!	write(*,*) 'kcell: ',kcell
!	if (cell_list(kcell)%state == DEAD) cycle
!	call setState(cell_list(kcell)%cellml_state)
!	call resetIntegrator
!	call multiStep(0.0, interval_dt, cellml_dt);
!	call getState(cell_list(kcell)%cellml_state)
!enddo

ok = .true.
call GrowCells(DELTA_T,ok)
if (.not.ok) then
	res = 3
	write(logmsg,'(a,i2)') 'Error: GrowCells: ',res
	call logger(logmsg)
	return
endif
	
	! Update Cbnd using current M, R1 and previous U, Cext
	call UpdateCbnd

	call SetupODEdiff
	call SiteCellToState
!	do it_solve = 1,NT_CONC
!		tstart = (it_solve-1)*dt
!		t_simulation = (istep-1)*DELTA_T + tstart
!		call Solver(it_solve,tstart,dt,Ncells)
!	enddo
!	if (idbug /= 0) then
!		i = cell_list(idbug)%iv
!		write(nfout,'(i6,2f10.4)') istep,allstate(i-1,OXYGEN),allstate(i,OXYGEN)
!	endif
	call StateToSiteCell
!	res = 0
!
!	! Compute U and update M, Cext, Cbnd
!	call UpdateMedium(DELTA_T)

if (mod(istep,60) == -1) then
	rmax = 0
	do kcell = 1,nlist
		r = cell_list(kcell)%site - Centre
		rmax = max(rmax,norm(r))
	enddo
	hour = istep/60
	write(logmsg,'(3i6,2f6.1)') istep, hour, Ncells, Radius, rmax
	call logger(logmsg)
	call CheckBdryList
!	call ShowConcs
	call check_bdry
endif
end subroutine

!--------------------------------------------------------------------------------
! TC = tumour cell
!--------------------------------------------------------------------------------
subroutine get_scene(nTC_list,TC_list) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_scene
use, intrinsic :: iso_c_binding
integer(c_int) :: nTC_list, TC_list(*)
integer :: k, kc, kcell, site(3), j, jb, colour, nlive
integer :: col(3)
integer :: x, y, z
integer :: itcstate, ctype, stage, region, highlight
integer :: last_id1, last_id2
logical :: ok, highlighting
integer, parameter :: axis_centre = -2	! identifies the spheroid centre
integer, parameter :: axis_end    = -3	! identifies the spheroid extent in 5 directions
integer, parameter :: axis_bottom = -4	! identifies the spheroid extent in the -Y direction, i.e. bottom surface
integer, parameter :: ninfo = 6			! the size of the info package for a cell (number of integers)
integer, parameter :: nax = 6			! number of points used to delineate the spheroid
real(REAL_KIND) :: cell_g, cell_size, cell_r
real(REAL_KIND), parameter :: growth_max = 2.2
real(REAL_KIND), parameter :: cell_size_max = 2.0

highlighting = (show_progeny /= 0)
nTC_list = 0

k = 0
last_id1 = 0
if (1 == 0) then
	! Need some markers to delineate the follicle extent.  These nax "cells" are used to convey (the follicle centre
	! and) the approximate ellipsoidal blob limits in the 3 axis directions.
	do k = 1,nax
		select case (k)
	!	case (1)
	!		x = Centre(1) + 0.5
	!		y = Centre(2) + 0.5
	!		z = Centre(3) + 0.5
	!		site = (/x, y, z/)
	!		itcstate = axis_centre
		case (1)
	!		x = Centre(1) - Radius%x - 2
			x = blobrange(1,1) - 1
			y = Centre(2) + 0.5
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (2)
	!		x = Centre(1) + Radius%x + 2
			x = blobrange(1,2) + 1
			y = Centre(2) + 0.5
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (3)
			x = Centre(1) + 0.5
	!		y = Centre(2) - Radius%y - 2
			y = blobrange(2,1) - 1
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_bottom
		case (4)
			x = Centre(1) + 0.5
	!		y = Centre(2) + Radius%y + 2
			y = blobrange(2,2) + 1
			z = Centre(3) + 0.5
			site = (/x, y, z/)
			itcstate = axis_end
		case (5)
			x = Centre(1) + 0.5
			y = Centre(2) + 0.5
	!		z = Centre(3) - Radius%z - 2
			z = blobrange(3,1) - 1
			site = (/x, y, z/)
			itcstate = axis_end
		case (6)
			x = Centre(1) + 0.5
			y = Centre(2) + 0.5
	!		z = Centre(3) + Radius%z + 2
			z = blobrange(3,2) + 1
			site = (/x, y, z/)
			itcstate = axis_end
		end select

		j = ninfo*(k-1)
		TC_list(j+1) = k-1
		TC_list(j+2:j+4) = site
		TC_list(j+5) = itcstate
		TC_list(j+6) = 1
		last_id1 = k-1
	enddo
	k = last_id1 + 1
endif

! Cells
nlive = 0
do kcell = 1,nlist
!	if (idbug /= 0 .and. cell_list(kcell)%ID /= idbug) cycle
	if (cell_list(kcell)%exists) then
		k = k+1
		j = ninfo*(k-1)
		site = cell_list(kcell)%site
!	    if (highlighting .and. cell_list(kcell)%ID == show_progeny) then
!	        highlight = 1
!	    else
!	        highlight = 0
!	    endif
!	    if (use_celltype_colour) then
!			colour = cell_list(kcell)%celltype
!		else
!			call CellColour(kcell,highlight,col)
!			colour = rgb(col)
!		endif
		TC_list(j+1) = kcell + last_id1
		TC_list(j+2:j+4) = site
		cell_g = cell_list(kcell)%cellml_state(1)
		cell_size = cell_list(kcell)%cellml_state(2)
		cell_r = (cell_size/cell_size_max)**(1./3.)
		TC_list(j+5) = cell_g*100			! 100*growth_stage
		TC_list(j+6) = min(1.0,cell_r)*100	! 100*fractional radius
!		TC_list(j+5) = colour
!		TC_list(j+6) = highlight
		last_id2 = kcell + last_id1
		nlive = nlive + 1
	endif
enddo
!nTC_list = last_id2
nTC_list = k
!write(logmsg,*) 'nlive: ',nlive
!call logger(logmsg)
end subroutine

!-----------------------------------------------------------------------------------------
! Rendered cell colour may depend on stage, state, receptor expression level.
! col(:) = (r,g,b)
!-----------------------------------------------------------------------------------------
subroutine CellColour(kcell,highlight,col)
integer :: kcell, highlight, col(3)
integer :: stage, status
integer, parameter :: WHITE(3) = (/255,255,255/)
integer, parameter :: RED(3) = (/255,0,0/)
integer, parameter :: GREEN(3) = (/0,255,0/)
integer, parameter :: BLUE(3) = (/0,0,255/)
integer, parameter :: DEEPRED(3) = (/200,0,0/)
integer, parameter :: DEEPBLUE(3) = (/30,20,255/)
integer, parameter :: DEEPGREEN(3) = (/0,150,0/)
integer, parameter :: LIGHTRED(3) = (/255,70,90/)
integer, parameter :: LIGHTBLUE(3) = (/0,200,255/)
integer, parameter :: LIGHTGREEN(3) = (/50,255,150/)
integer, parameter :: DEEPORANGE(3) = (/240,70,0/)
integer, parameter :: LIGHTORANGE(3) = (/255,130,0/)
integer, parameter :: YELLOW(3) = (/255,255,0/)
integer, parameter :: DEEPPURPLE(3) = (/180,180,30/)
integer, parameter :: LIGHTPURPLE(3) = (/230,230,100/)
integer, parameter :: DEEPBROWN(3) = (/130,70,0/)
integer, parameter :: LIGHTBROWN(3) = (/200,100,0/)
integer, parameter :: GRAY(3) = (/128,128,128/)

integer, parameter :: Qt_white = 3
integer, parameter :: Qt_black = 2
integer, parameter :: Qt_red = 7
integer, parameter :: Qt_darkRed = 13
integer, parameter :: Qt_green = 8
integer, parameter :: Qt_darkGreen = 14
integer, parameter :: Qt_blue = 9
integer, parameter :: Qt_darkBlue = 15
integer, parameter :: Qt_cyan = 10
integer, parameter :: Qt_darkCyan = 16
integer, parameter :: Qt_magenta = 11
integer, parameter :: Qt_darkMagenta = 17
integer, parameter :: Qt_yellow = 12
integer, parameter :: Qt_darkYellow = 18
integer, parameter :: Qt_gray = 5
integer, parameter :: Qt_darkGray = 4
integer, parameter :: Qt_lightGray = 6

if (highlight == 0) then
    col = LIGHTORANGE
else
    col = LIGHTRED
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Pack the colours (r,g,b) into an integer.
!-----------------------------------------------------------------------------------------
integer function rgb(col)
integer :: col(3)

rgb = ishft(col(1),16) + ishft(col(2),8) + col(3)
end function

!-----------------------------------------------------------------------------------------
! Live cells = Ncells
! Total deaths = Ndead
! Current tagged = Ntodie - Ntagdead 
!-----------------------------------------------------------------------------------------
subroutine get_summary(summaryData, i_growth_cutoff) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_summary
use, intrinsic :: iso_c_binding
integer(c_int) :: summaryData(*), i_growth_cutoff
integer :: Ndead, Ntagged, Ntodie, Ntagdead, diam_um, vol_mm3_100000
integer :: ngrowth(3), growth_percent_10, necrotic_percent_10,Ntagged_anoxia
real(REAL_KIND) :: vol_cm3, vol_mm3, hour

hour = istep*DELTA_T/3600.
vol_cm3 = Vsite_cm3*Nsites			! total volume in cm^3
vol_mm3 = vol_cm3*1000				! volume in mm^3
vol_mm3_100000 = vol_mm3*100000		! 100000 * volume in mm^3
diam_um = 2*DELTA_X*Radius*10000
Ntagged_anoxia = Nanoxia_tag - Nanoxia_dead				! number currently tagged by anoxia
call getGrowthCount(ngrowth)
growth_percent_10 = (1000*ngrowth(i_growth_cutoff))/Ncells
necrotic_percent_10 = (1000*(Nsites-Ncells))/Nsites
summaryData(1:8) = (/ istep, Ncells, Nanoxia_dead, Ntagged_anoxia, &
	diam_um, vol_mm3_100000, growth_percent_10, necrotic_percent_10 /)
write(nfres,'(i8,f8.2,f8.4,8i7,7f7.3)') istep, hour, vol_mm3, diam_um, Ncells, &
	ngrowth(:)/real(Ncells), (Nsites-Ncells)/real(Nsites)
end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getHypoxicCount(nhypoxic)
integer :: nhypoxic(3)
integer :: kcell, i

nhypoxic = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%conc(OXYGEN) < O2cutoff(i)) nhypoxic(i) = nhypoxic(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Need to compare growth rate with a fraction of average growth rate
!--------------------------------------------------------------------------------
subroutine getGrowthCount(ngrowth)
integer :: ngrowth(3)
integer :: kcell, i
real(REAL_KIND) :: r_mean

r_mean = Vdivide0/(2*divide_time_mean)
ngrowth = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	do i = 1,3
		if (cell_list(kcell)%dVdt < growthcutoff(i)*r_mean) ngrowth(i) = ngrowth(i) + 1
	enddo
enddo
end subroutine

!--------------------------------------------------------------------------------
! Note: axis = 0,1,2
!--------------------------------------------------------------------------------
subroutine get_fieldinfo(nxx, axis, fraction, ns, nc, cused, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fieldinfo
use, intrinsic :: iso_c_binding
integer(c_int) :: nxx, axis, ns, nc, cused(*), res
real(c_double) :: fraction
integer rng(3,2), ichemo, kcell, x, y, z

!call logger('get_fieldinfo')
res = 0
nxx = NX
nc = MAX_CHEMO
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cused(ichemo) = 1
    else
        cused(ichemo) = 0
    endif
enddo
cused(MAX_CHEMO+1) = 1		! Growth rate
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
!write(logmsg,*) 'Centre, Radius, axis, fraction: ',Centre, Radius, axis, fraction
!call logger(logmsg)
!write(logmsg,*) 'rng: ',rng
!call logger(logmsg)
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell == OUTSIDE_TAG) cycle
            ns = ns + 1
        enddo
    enddo
enddo
!write(logmsg,*) 'get_fieldinfo: ns: ',ns
!call logger(logmsg)
end subroutine

!--------------------------------------------------------------------------------
! Need to transmit medium concentration data.  This could be a separate subroutine.
!--------------------------------------------------------------------------------
subroutine get_fielddata(axis, fraction, nfdata, nc, fdata, res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_fielddata
use, intrinsic :: iso_c_binding
real(c_double) :: fraction
integer(c_int) :: axis, nc, nfdata, res
type(FIELD_DATA) :: fdata(*)
integer rng(3,2), kcell, x, y, z, i, ns
real(REAL_KIND) :: growthrate

!write(logmsg,*) 'get_fielddata: nfdata, nc: ',nfdata, nc, MAX_CHEMO
!call logger(logmsg)
res = 0
if (nc > MAX_CHEMO) then
	write(logmsg,*) 'Error: get_fielddata: dimension of conc(MAX_CHEMO) not big enough!'
	call logger(logmsg)
	res = 1
	return
endif
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
rng(axis,:) = Centre(axis) + fraction*Radius
ns = 0
do z = rng(3,1),rng(3,2)
    do y = rng(2,1),rng(2,2)
        do x = rng(1,1),rng(1,2)
            kcell = occupancy(x,y,z)%indx(1)
            if (kcell == OUTSIDE_TAG) cycle
            ns = ns + 1
	        i = ODEdiff%ivar(x,y,z)
            fdata(ns)%site = (/x,y,z/)
            fdata(ns)%state = 1
            if (kcell > 0) then
                fdata(ns)%volume = cell_list(kcell)%volume
                call getGrowthrate(kcell,growthrate)
            else
                fdata(ns)%volume = 0
                growthrate = 0
            endif
            fdata(ns)%dVdt = growthrate
            fdata(ns)%conc(1:MAX_CHEMO) = allstate(i,1:MAX_CHEMO)	! cell_list(kcell)%conc(1:MAX_CHEMO)
        enddo
    enddo
enddo
if (ns /= nfdata) then
    write(logmsg,*) 'Error: inconsistent nsites: ',nfdata, ns
    call logger(logmsg)
    res = 2
    return
endif

end subroutine

!--------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
subroutine getGrowthrate(kcell,growthrate)
integer :: kcell
real(REAL_KIND) :: growthrate

growthrate = cell_list(kcell)%dVdt
end subroutine

!--------------------------------------------------------------------------------
! Returns all the extracellular concentrations along a line through the blob centre.
! Together with the growth rate dVdt
!--------------------------------------------------------------------------------
subroutine get_concdata(ns, dx, ex_conc) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_concdata
use, intrinsic :: iso_c_binding
integer(c_int) :: ns
real(c_double) :: dx, ex_conc(*)
real(REAL_KIND) :: cbnd, cmin = 1.0e-6
integer rng(3,2), i, k, ichemo, kcell, x, y, z

!call logger('get_concdata')
dx = DELTA_X
rng(:,1) = Centre(:) - (Radius + 2)
rng(:,2) = Centre(:) + (Radius + 2)
!rng(axis,:) = Centre(axis) + fraction*Radius
y = Centre(2) + 0.5
z = Centre(3) + 0.5
ns = 1
do x = rng(1,1),rng(1,2)
    kcell = occupancy(x,y,z)%indx(1)
    if (kcell == OUTSIDE_TAG) cycle
    ns = ns + 1
    do ichemo = 1,MAX_CHEMO+1
        i = ODEdiff%ivar(x,y,z)
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
        if (ichemo <= MAX_CHEMO) then
			if (chemo(ichemo)%used) then
				if (i > 0) then
					ex_conc(k) = allstate(i,ichemo)
				else
					ex_conc(k) = 0
					call logger('get_concdata: ODEdiff%ivar = 0')
				endif
			else
				ex_conc(k) = 0
			endif
        elseif (ichemo == MAX_CHEMO+1) then	! growth rate
			if (kcell > 0) then
				ex_conc(k) = cell_list(kcell)%dVdt
			else
				ex_conc(k) = 0
			endif
		endif
    enddo
enddo
! Add concentrations at the two boundaries 
! At ns=1, at at ns=ns+1
ns = ns+1
do ichemo = 1,MAX_CHEMO
    if (chemo(ichemo)%used) then
        cbnd = BdryConc(ichemo,t_simulation)
        k = ichemo
        ex_conc(k) = cbnd
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
        ex_conc(k) = cbnd
    else
        k = ichemo
        ex_conc(k) = 0
        k = (ns-1)*(MAX_CHEMO+1) + ichemo
        ex_conc(k) = 0
    endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of cell volume.
! nv is passed from the GUI
! Min divide volume = Vdivide0 - dVdivide
! therefore Vmin = (Vdivide0 - dVdivide)/2
! Max divide volume = Vmax = Vdivide0 + dVdivide
! dv = (Vmax - Vmin)/nv
! v0 = Vmin + dv/2
!--------------------------------------------------------------------------------
subroutine get_volprob(nv, v0, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_volprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax

!call logger('get_volprob')
Vmin = (Vdivide0 - dVdivide)/2
Vmax = Vdivide0 + dVdivide
dv = (Vmax - Vmin)/nv
v0 = Vmin + dv/2
prob(1:nv) = 0
n = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%volume
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!--------------------------------------------------------------------------------
! Get component name and number of variables
!--------------------------------------------------------------------------------
subroutine CellML_get_component_info(icomp,component_name, nvars) BIND(C,name='CellML_get_component_info')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_component_info
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp
integer(c_int) :: nvars
character(c_char),dimension(1) :: component_name
integer :: nlen, i

nlen = index(component(icomp),' ')
do i = 1,nlen-1
	component_name(i) = component(icomp)(i:i)
enddo
component_name(nlen) = char(0)
nvars = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		nvars = nvars + 1
	endif
enddo
write(nflog,*) icomp,nvars,component(icomp)(1:nlen-1)
end subroutine

!--------------------------------------------------------------------------------
! Get variable name
!--------------------------------------------------------------------------------
subroutine CellML_get_variable_name(icomp,ivar,variable_name) BIND(C,name='CellML_get_variable_name')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_variable_name
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
character(c_char),dimension(1) :: variable_name
character(64) :: varname
integer :: nv, i, k, nlen

nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			nlen = index(IDlist(i)%var_name,' ')
			do k = 1,nlen-1
				variable_name(k) = IDlist(i)%var_name(k:k)
			enddo
			variable_name(nlen) = char(0)
			write(nflog,*) icomp,ivar,IDlist(i)%var_name(1:nlen-1)
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Get variable value
!--------------------------------------------------------------------------------
subroutine CellML_get_variable_value(icomp,ivar,variable_value) BIND(C,name='CellML_get_variable_value')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_get_variable_value
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
real(c_double) :: variable_value
integer :: nv, i

write(nflog,*) 'CellML_get_variable_value: ',icomp,ivar
call getState(csim_state)
!write(nflog,'(10f7.3)') csim_state
nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			variable_value = csim_state(i)
			write(nflog,*) icomp,ivar,variable_value
!			if (icomp == 2 .and. ivar == 5) stop
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Set variable value
!--------------------------------------------------------------------------------
subroutine CellML_set_variable_value(icomp,ivar,variable_value) BIND(C,name='CellML_set_variable_value')
!DEC$ ATTRIBUTES DLLEXPORT :: CellML_set_variable_value
use, intrinsic :: iso_c_binding
integer(c_int),value :: icomp, ivar
real(c_double),value :: variable_value
integer :: nv, i

nv = 0
do i = 0,nvariables-1
	if (IDlist(i)%comp_index == icomp) then
		if (ivar == nv) then
			call setStateValue(i,variable_value)
			write(nflog,*) icomp,ivar,variable_value
			return
		endif
		nv = nv + 1
	endif
enddo
end subroutine

!--------------------------------------------------------------------------------
! Returns the distribution of intracellular O2 level
! nv is passed from the GUI
!--------------------------------------------------------------------------------
subroutine get_oxyprob(nv, dv, prob) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: get_oxyprob
use, intrinsic :: iso_c_binding
integer(c_int) :: nv
real(c_double) :: v0, dv, prob(*)
integer :: n, kcell, k
real(REAL_KIND) :: v, Vmin, Vmax, O2max

!call logger('get_oxyprob')
Vmin = 0
Vmax = chemo(OXYGEN)%bdry_conc
dv = (Vmax - Vmin)/nv
prob(1:nv) = 0
n = 0
O2max = 0
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	v = cell_list(kcell)%conc(OXYGEN)
	O2max = max(O2max,v)
	k = (v - Vmin)/dv + 1
	k = min(k,nv)
	k = max(k,1)
	prob(k) = prob(k) + 1
	n = n+1
enddo
prob(1:nv) = prob(1:nv)/n
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Execute(ncpu,infile_array,inbuflen,outfile_array,outbuflen) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: execute
use, intrinsic :: iso_c_binding
character(c_char) :: infile_array(128), outfile_array(128)
integer(c_int) :: ncpu, inbuflen, outbuflen
character*(128) :: infile, outfile
logical :: ok, success
integer :: i, res

infile = ''
do i = 1,inbuflen
	infile(i:i) = infile_array(i)
enddo
outfile = ''
do i = 1,outbuflen
	outfile(i:i) = outfile_array(i)
enddo

open(nflog,file='spheroid.log',status='replace')
open(nfres,file='spheroid_ts.out',status='replace')
write(nfres,'(a)') 'istep hour vol_mm3 diam_um Ncells &
 f_growth_1 f_growth_2 f_growth_3'

#ifdef GFORTRAN
    write(logmsg,'(a)') 'Built with GFORTRAN'
	call logger(logmsg)
#endif

logmsg = 'OS??'
#ifdef LINUX
    write(logmsg,'(a)') 'OS is Linux'
#endif
#ifdef OSX
    write(logmsg,'(a)') 'OS is OS-X'
#endif
#ifdef _WIN32
    write(logmsg,'(a)') 'OS is Windows'
#endif
#ifdef WINDOWS
    write(logmsg,'(a)') 'OS is Windows'
#endif
call logger(logmsg)

!#ifdef OPENMP
#if defined(OPENMP) || defined(_OPENMP)
    write(logmsg,'(a)') 'Executing with OpenMP'
	call logger(logmsg)
#endif

write(logmsg,*) 'inputfile:  ', infile
call logger(logmsg)
write(logmsg,*) 'outputfile: ', outfile 
call logger(logmsg)
if (use_TCP) then
	call connecter(ok)
	if (.not.ok) then
		call logger('Failed to make TCP connections')
		return
	endif
endif

istep = 0
NX=100
NY=100
NZ=100
DELTA_T = 600
nsteps = 100
res=0

call Setup(ncpu,infile,outfile,ok)
if (ok) then
!	clear_to_send = .true.
!	simulation_start = .true.
else
	call logger('=== Setup failed ===')
endif
if (ok) then
	res = 0
else
	res = 1
endif
execute_t1 = wtime()

end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine DisableTCP
!DEC$ ATTRIBUTES DLLEXPORT :: disableTCP
!DEC$ ATTRIBUTES STDCALL, REFERENCE, MIXED_STR_LEN_ARG, ALIAS:"DISABLETCP" :: disableTCP

use_TCP = .false.   ! because this is called from spheroid_main()	
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connection(awp,port,error)
TYPE(winsockport) :: awp
integer :: port, error
integer :: address = 0
!!!character*(64) :: ip_address = "127.0.0.1"C      ! need a portable way to make a null-terminated C string
character*(64) :: host_name = "localhost"

if (.not.winsock_init(1)) then
    call logger("winsock_init failed")
    stop
endif

awp%handle = 0
awp%host_name = host_name
awp%ip_port = port
awp%protocol = IPPROTO_TCP
call Set_Winsock_Port (awp,error)

if (.not.awp%is_open) then
    write(nflog,*) 'Error: connection: awp not open: ',port
else
    write(nflog,*) 'connection: awp open: ',port, error
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Connecter(ok)
logical :: ok
integer :: error

! Main connection
ok = .true.
error = 0
call Connection(awp_0,TCP_PORT_0,error)
if (awp_0%handle < 0 .or. error /= 0) then
    write(logmsg,'(a)') 'TCP connection to TCP_PORT_0 failed'
    call logger(logmsg)
    ok = .false.
    return
endif
if (.not.awp_0%is_open) then
	write(logmsg,'(a)') 'No connection to TCP_PORT_0'
    call logger(logmsg)
    ok = .false.
    return
endif
write(logmsg,'(a)') 'Connected to TCP_PORT_0  '
call logger(logmsg)

if (use_CPORT1) then
	call connection(awp_1,TCP_PORT_1,error)
	if (awp_1%handle < 0 .or. error /= 0) then
		write(logmsg,'(a)') 'TCP connection to TCP_PORT_1 failed'
		call logger(logmsg)
		ok = .false.
		return
	endif
	if (.not.awp_1%is_open) then
		write(logmsg,'(a)') 'No connection to TCP_PORT_1'
		call logger(logmsg)
		ok = .false.
		return
	endif
	write(logmsg,'(a)') 'Connected to TCP_PORT_1  '
	call logger(logmsg)
endif
! Allow time for completion of the connection
call sleeper(2)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine terminate_run(res) BIND(C)
!DEC$ ATTRIBUTES DLLEXPORT :: terminate_run 
use, intrinsic :: iso_c_binding
integer(c_int) :: res
character*(8), parameter :: quit = '__EXIT__'
integer :: error, i

call logger('terminate_run')
call Wrapup

if (res == 0) then
	call logger(' Execution successful!')
elseif (res == -1) then
	call logger(' Execution stopped')
else
	call logger('  === Execution failed ===')
endif
write(logmsg,'(a,f10.2)') 'Execution time (min): ',(wtime() - execute_t1)/60
call logger(logmsg)

close(nflog)

if (use_TCP) then
	if (stopped) then
	    call winsock_close(awp_0)
	    if (use_CPORT1) call winsock_close(awp_1)
	else
	    call winsock_send(awp_0,quit,8,error)
	    call winsock_close(awp_0)
		if (use_CPORT1) then
			call winsock_send(awp_1,quit,8,error)
			call winsock_close(awp_1)
		endif
	endif
endif
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine Wrapup
integer :: ierr, ichemo
logical :: isopen

call logger('doing wrapup ...')
ierr = 0
if (allocated(zoffset)) deallocate(zoffset)
if (allocated(zdomain)) deallocate(zdomain)
!if (allocated(occupancy)) deallocate(occupancy)
!if (allocated(cell_list)) deallocate(cell_list)
!if (allocated(allstate)) deallocate(allstate)
if (allocated(allstatep)) deallocate(allstatep)
if (allocated(work_rkc)) deallocate(work_rkc)
do ichemo = 1,MAX_CHEMO
	if (allocated(chemo(ichemo)%coef)) deallocate(chemo(ichemo)%coef)
	if (allocated(chemo(ichemo)%conc)) deallocate(chemo(ichemo)%conc)
	if (allocated(chemo(ichemo)%grad)) deallocate(chemo(ichemo)%grad)
enddo
!if (allocated(ODEdiff%ivar)) deallocate(ODEdiff%ivar)
if (allocated(ODEdiff%varsite)) deallocate(ODEdiff%varsite)
if (allocated(ODEdiff%icoef)) deallocate(ODEdiff%icoef)
call freeCellStateVectors
call logger('deallocated all arrays')

! Close all open files
inquire(unit=nfout,OPENED=isopen)
if (isopen) then
	close(nfout)
	call logger('closed nfout')
endif
inquire(nfres,OPENED=isopen)
if (isopen) close(nfres)
call logger('closed files')

if (par_zig_init) then
	call par_zigfree
endif
call logger('freed par_zig')
end subroutine

end module
