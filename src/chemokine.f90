! Chemokine data

module chemokine

use global

implicit none

!type receptor_type
!	character(8) :: name
!	logical :: used
!	integer :: chemokine
!	integer :: sign
!	real(REAL_KIND) :: strength
!	real(REAL_KIND) :: level(5)
!	real(REAL_KIND) :: saturation_threshold
!	real(REAL_KIND) :: refractory_time
!end type
!type(receptor_type), target :: receptor(MAX_RECEPTOR)

type chemokine_type
	character(16) :: name
	logical :: used
	logical :: use_secretion
	real(REAL_KIND) :: bdry_rate
	real(REAL_KIND) :: bdry_conc
	real(REAL_KIND) :: diff_coef
	real(REAL_KIND) :: membrane_diff
	logical :: decay
	real(REAL_KIND) :: halflife
	real(REAL_KIND) :: decay_rate
!	real(REAL_KIND) :: bdry_halflife
!	real(REAL_KIND) :: bdry_decay_rate
	real(REAL_KIND) :: max_cell_rate	! Vmax
	real(REAL_KIND) :: MM_C0			! Km
	real(REAL_KIND) :: Hill_N
	real(REAL_KIND) :: medium_diff_coef	! diffusion coefficient in the medium
	real(REAL_KIND) :: medium_dlayer	! unstirred layer thickness
	real(REAL_KIND) :: medium_M			! mass of constituent
	real(REAL_KIND) :: medium_U			! total blob uptake rate
	real(REAL_KIND) :: medium_Cext		! far-field concentration
	real(REAL_KIND) :: medium_Cbnd		! boundary concentration
	real(REAL_KIND), allocatable :: coef(:,:)
	real(REAL_KIND), allocatable :: conc(:,:,:)
	real(REAL_KIND), allocatable :: grad(:,:,:,:)
end type
type(chemokine_type), target :: chemo(MAX_CHEMO)

type ODEdiff_type
	integer :: ichemo
	integer :: nextra
	integer :: nintra
	integer :: nvars
	integer, allocatable :: ivar(:,:,:)
	integer, allocatable :: varsite(:,:)
	integer, allocatable :: icoef(:,:)
	integer, allocatable :: iexcoef(:,:)
	integer, allocatable :: vartype(:)
	integer, allocatable :: cell_index(:)
	integer, allocatable :: isite_extra(:)
	integer, allocatable :: isite_intra(:)
	integer, allocatable :: extra_isite(:)
!	integer, allocatable :: ncoef(:)
	real(REAL_KIND) :: deltaC_soft
	real(REAL_KIND) :: k_soft
	real(REAL_KIND) :: C1_soft
end type
type(ODEdiff_type) :: ODEdiff

contains

!----------------------------------------------------------------------------------------
! Convert halflife in hours to a decay rate /sec
!----------------------------------------------------------------------------------------
real(REAL_KIND) function DecayRate(halflife)
real(REAL_KIND) :: halflife

if (halflife == 0) then		! No decay
	DecayRate = 0
else
	DecayRate = log(2.0)/(halflife*60*60)
endif
end function

!----------------------------------------------------------------------------------------
! Note units:
! distance		cm
! volume		cm^3
! time			s
! diff coeff	cm^2.s^-1
! mass			mol
! concentration	mM where 1 mM = 1.0e-6 mol.cm^-3
! consumption	mol.cell^-1.s^-1
! production	mol.cell^-1.s^-1
!----------------------------------------------------------------------------------------
subroutine SetupChemo
integer :: ichemo

chemo(OXYGEN)%name = 'Oxygen'
chemo(GLUCOSE)%name = 'Glucose'
chemo(TRACER)%name = 'Tracer'
chemo(OXYGEN)%decay_rate = 0
chemo(GLUCOSE)%decay_rate = 0
chemo(TRACER)%decay_rate = 0

do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
	chemo(ichemo)%medium_dlayer = d_layer	! for now, using the same unstirred layer thickness for all constituents
enddo
call AllocateConcArrays
!call SetMMParameters
end subroutine

!----------------------------------------------------------------------------------------
! Consumption rate = Rmax*C/(MM_C0 + C)  (Michaelis-Menten)
! where Km = MM_C0 >= Rmax*T*10^6/Vextra_cm3 is the O2 concentration at which cell uptake is halved,
! and T is the maximum expected time step
! This is to ensure that all the O2 in the site is not depleted in a time step, making
! O2 conc go negative.  For now it seems reasonable to assume that glucose uptake
! varies in proportion to O2 uptake, i.e. we only need Oxygen M-M
!----------------------------------------------------------------------------------------
subroutine SetMMParameters

!chemo(OXYGEN)%MM_C0 = 0.00133		! 1 mmHg = 1.33 uM (Kevin suggests 1 - 2 uM)
!chemo(GLUCOSE)%MM_C0 = chemo(OXYGEN)%MM_C0*chemo(GLUCOSE)%max_cell_rate/chemo(OXYGEN)%max_cell_rate
write(logmsg,'(a,e12.4)') 'Oxygen MM_C0: ',chemo(OXYGEN)%MM_C0
call logger(logmsg)
write(logmsg,'(a,e12.4)') 'Glucose MM_C0: ',chemo(GLUCOSE)%MM_C0
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine AllocateConcArrays
integer :: ic

write(logmsg,*) 'AllocateConcArrays'
call logger(logmsg)
do ic = 1,MAX_CHEMO
	if (chemo(ic)%used) then
		if (allocated(chemo(ic)%conc)) then
			call logger("chemo(ic)%conc already allocated")
			deallocate(chemo(ic)%conc)
!			stop
		endif
		allocate(chemo(ic)%conc(NX,NY,NZ))
		if (allocated(chemo(ic)%grad)) then
			call logger("chemo(ic)%grad already allocated")
			deallocate(chemo(ic)%grad)
!			stop
		endif
		allocate(chemo(ic)%grad(3,NX,NY,NZ))
		chemo(ic)%conc = chemo(ic)%bdry_conc	
	endif
enddo
write(logmsg,*) 'did AllocateConcArrays'
call logger(logmsg)
end subroutine

!----------------------------------------------------------------------------------------
! Set the concentrations at a site when it is on the boundary
!----------------------------------------------------------------------------------------
subroutine SetBdryConcs(site)
integer :: site(3)
integer :: ichemo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		chemo(ichemo)%conc(site(1),site(2),site(3)) = chemo(ichemo)%bdry_conc
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
! Reset the concentrations at a site that is no longer on the boundary
!----------------------------------------------------------------------------------------
subroutine ResetConcs(site0)
integer :: site0(3)
integer :: site(3), ichemo, k, n
real(REAL_KIND) :: csum(MAX_CHEMO)

csum = 0
n = 0
do k = 1,27
	if (k == 14) cycle
	site = site0 + jumpvec(:,k)
	if (occupancy(site(1),site(2),site(3))%indx(1) == OUTSIDE_TAG) cycle
!	if (associated(occupancy(site(1),site(2),site(3))%bdry)) cycle
	n = n+1
	do ichemo = 1,MAX_CHEMO
		if (chemo(ichemo)%used) then
			csum(ichemo) = csum(ichemo) + chemo(ichemo)%conc(site(1),site(2),site(3))
		endif
	enddo
enddo

do ichemo = 1,MAX_CHEMO
	if (chemo(ichemo)%used) then
		chemo(ichemo)%conc(site0(1),site0(2),site0(3)) = csum(ichemo)/n
	endif
enddo
end subroutine

!----------------------------------------------------------------------------------------
!----------------------------------------------------------------------------------------
subroutine ShowConcs
integer :: x, y, z, k, x1, site(3), ichemo
real(REAL_KIND) :: C(MAX_CHEMO,100)
logical :: start

y = Centre(2) 
z = Centre(3)
start = .false.
k = 0
do x = blobrange(1,1),blobrange(1,2)
	if (occupancy(x,y,z)%indx(1) == OUTSIDE_TAG) then
		if (start) then
			exit
		else
			cycle
		endif
	endif
	if (.not.start) x1 = x
	start = .true.
	k = k+1
	do ichemo = 1,MAX_CHEMO
		if (.not.chemo(ichemo)%used) cycle
	    C(ichemo,k) = chemo(ichemo)%conc(x,y,z)
	enddo
enddo
do ichemo = 1,MAX_CHEMO
	if (.not.chemo(ichemo)%used) cycle
    write(*,'(i4,2x,10f7.2)') istep/60,C(ichemo,1:10)
enddo
!if (C(1) < 100) then
!	write(logmsg,*) 'First conc < 100: istep: ',istep, C(1),occupancy(x1,y,z)%indx(1)
!	call logger(logmsg)
!	site = (/x1,y,z/)
!	write(logmsg,'(a,L)') 'isbdry: ',isbdry(site)
!	call logger(logmsg)
!	stop
!endif
end subroutine

end module