! Cancer cell state development 

module cellstate
use global
use boundary
use chemokine
use ode_diffuse
use csim_abm
implicit none

!real(REAL_KIND), parameter :: Vdivide0 = 1.6
!real(REAL_KIND), parameter :: dVdivide = 0.05
integer :: kcell_dividing = 0

contains

!-----------------------------------------------------------------------------------------
! Need to initialize site and cell concentrations when a cell divides and when there is
! cell death.
!-----------------------------------------------------------------------------------------
subroutine GrowCells(dt,ok)
real(REAL_KIND) :: dt
logical :: ok

ok = .true.
call CellGrowth(dt,ok)
!if (.not.ok) return
!if (use_death) then
!	call CellDeath(dt,ok)
!	if (.not.ok) return
!endif
!if (use_migration) then
!	call CellMigration(ok)
!	if (.not.ok) return
!endif
end subroutine

!-----------------------------------------------------------------------------------------
! Cells move to preferable nearby sites.
! For now this is turned off - need to formulate a sensible interpretation of "preferable"
!-----------------------------------------------------------------------------------------
subroutine CellMigration(ok)
logical :: ok
integer :: kcell, j, indx, site0(3), site(3), jmax
real(REAL_KIND) :: C0(MAX_CHEMO), C(MAX_CHEMO), v0, v, vmax, d0, d

call logger('CellMigration is not yet implemented')
ok = .false.
return

ok = .true.
do kcell = 1,nlist
	if (cell_list(kcell)%state == DEAD) cycle
	site0 = cell_list(kcell)%site
	C0 = occupancy(site0(1),site0(2),site0(3))%C(:)
	v0 = SiteValue(C0)
	d0 = cdistance(site0)
	jmax = 0
	vmax = -1.0e10
	do j = 1,27
		if (j == 14) cycle
		site = site0 + jumpvec(:,j)
		indx = occupancy(site(1),site(2),site(3))%indx(1)
!		if (indx < -100) then	! necrotic site
		if (indx == 0) then	!	vacant site
			C = occupancy(site(1),site(2),site(3))%C(:)
			v = SiteValue(C)
			d = cdistance(site)
			if (d > d0 .and. v > v0) then
				vmax = v
				jmax = j
			endif
		endif
	enddo
	if (jmax > 0) then
		site = site0 + jumpvec(:,jmax)
		indx = occupancy(site(1),site(2),site(3))%indx(1)
		cell_list(kcell)%site = site
		occupancy(site(1),site(2),site(3))%indx(1) = kcell
		occupancy(site0(1),site0(2),site0(3))%indx(1) = indx
!		write(*,'(i2,2f8.4,3i4,f6.1,4x,3i4,f6.1)') jmax,v0,vmax,site0,d0,site,d
	endif
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A measure of the attractiveness of a site with concentrations C(:)
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function SiteValue(C)
real(REAL_KIND) :: C(:)

SiteValue = C(OXYGEN)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDies
integer :: kcell, i, kpar=0

do i = 1,10
	kcell = random_int(1,nlist,kpar)
	call CellDies(kcell)
enddo
stop
end subroutine

!-----------------------------------------------------------------------------------------
! If the dying cell site is less than a specified fraction f_migrate of the blob radius,
! the site migrates towards the blob centre.
! %indx -> 0
! If the site is on the boundary, it is removed from the boundary list, and %indx -> OUTSIDE_TAG
! The cell contents should be released into the site.
!-----------------------------------------------------------------------------------------
subroutine CellDies(kcell)
integer :: kcell
integer :: site(3)
real(REAL_KIND) :: V

cell_list(kcell)%state = DEAD
cell_list(kcell)%exists = .false.
Ncells = Ncells - 1
site = cell_list(kcell)%site
occupancy(site(1),site(2),site(3))%indx(1) = 0
if (associated(occupancy(site(1),site(2),site(3))%bdry)) then
	call bdrylist_delete(site,bdrylist)
    nullify(occupancy(site(1),site(2),site(3))%bdry)
	occupancy(site(1),site(2),site(3))%indx(1) = OUTSIDE_TAG
	Nsites = Nsites - 1
	bdry_changed = .true.
	call OutsideNeighbours(site)
	call AddToMedium(kcell,site)
elseif (vary_conc) then
	V = cell_list(kcell)%volume*Vcell_cm3
	occupancy(site(1),site(2),site(3))%C = ((Vsite_cm3 - V)*occupancy(site(1),site(2),site(3))%C + V*cell_list(kcell)%conc)/Vsite_cm3
endif
!call NecroticMigration(site)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AddToMedium(kcell,site)
integer :: kcell, site(3)
integer :: ic
real(REAL_KIND) :: V, Cex(MAX_CHEMO), Cin(MAX_CHEMO)

if (.not.vary_conc) return
Cex = occupancy(site(1),site(2),site(3))%C
Cin = cell_list(kcell)%conc
V = cell_list(kcell)%volume*Vcell_cm3
do ic = 1,MAX_CHEMO
	if (.not.chemo(ic)%used) cycle
	chemo(ic)%medium_M = chemo(ic)%medium_M + V*Cin(ic) + (Vsite_cm3 - V)*Cex(ic)
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine RemoveFromMedium
integer :: ic

do ic = 1,MAX_CHEMO
	if (.not.chemo(ic)%used) cycle
	chemo(ic)%medium_M = chemo(ic)%medium_M - Vsite_cm3*chemo(ic)%medium_Cbnd
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! When a bdry cell dies, need to check the neighbours to see if a site is now made outside.
! To be made outside a site must have %indx(1) = 0 and at least one Neumann neighbour outside.
!-----------------------------------------------------------------------------------------
subroutine OutsideNeighbours(site0)
integer :: site0(3)
integer :: k, j, x, y, z, xx, yy, zz
logical :: isout, done

done = .false.
do while(.not.done)
	done = .true.
	do k = 1,27
		if (k == 14) cycle
		x = site0(1) + jumpvec(1,k)
		y = site0(2) + jumpvec(2,k)
		z = site0(3) + jumpvec(3,k)
		if (outside_xyz(x,y,z)) cycle
		if (occupancy(x,y,z)%indx(1) == 0) then
			isout = .false.
			do j = 1,6
				xx = x + neumann(1,j)
				yy = y + neumann(2,j)
				zz = z + neumann(3,j)
				if (outside_xyz(xx,yy,zz)) cycle
				if (occupancy(xx,yy,zz)%indx(1) == OUTSIDE_TAG) then
					isout = .true.
					exit
				endif
			enddo
			if (isout) then
				done = .false.
				occupancy(x,y,z)%indx(1) = OUTSIDE_TAG
			endif
		endif
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! A necrotic site migrates towards the blob centre, stopping when another necrotic 
! site is reached
!-----------------------------------------------------------------------------------------
subroutine NecroticMigration(site0)
integer :: site0(3)
integer :: site1(3), site2(3), site(3), j, jmin, kcell, tmp_indx
real(REAL_KIND) :: d1, d2, dmin

!write(logmsg,*) 'NecroticMigration: site0: ',site0
!call logger(logmsg)
site1 = site0
do
	d1 = cdistance(site1)
	dmin = 1.0e10
	jmin = 0
	do j = 1,27
		if (j == 14) cycle
		site = site1 + jumpvec(:,j)
!		if (occupancy(site(1),site(2),site(3))%indx(1) < 0) cycle	! do not swap with another necrotic site
		if (occupancy(site(1),site(2),site(3))%indx(1) == 0) cycle	! do not swap with another necrotic site
		d2 = cdistance(site)
		if (d2 < dmin) then
			dmin = d2
			jmin = j
		endif
	enddo
!	write(*,*) 'd1, dmin, jmin: ',d1,dmin, jmin
	if (dmin >= d1) then
		site2 = site1
		exit
	endif
	if (jmin <= 0) then
		write(*,*) 'Error: jmin: ',jmin
		stop
	endif
!	write(*,*) site1,jmin,jumpvec(:,jmin)
	site2 = site1 + jumpvec(:,jmin)
!	write(*,*) site2
	! Now swap site1 and site2
	kcell = occupancy(site2(1),site2(2),site2(3))%indx(1)
	tmp_indx = occupancy(site1(1),site1(2),site1(3))%indx(1)
	occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell
	occupancy(site2(1),site2(2),site2(3))%indx(1) = tmp_indx
	cell_list(kcell)%site = site1
	site1 = site2
enddo
!write(*,*) 'NecroticMigration: site2: ',site2
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_CellDivision(ok)
logical :: ok
integer :: kcell, kpar=0
kcell = random_int(1,nlist,kpar)
call CellDivider(kcell,ok)
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine UpdateCellCycle(kcell,tcycle)
integer :: kcell
real(REAL_KIND) :: tcycle
integer :: i

if (cell_list(kcell)%cell_cycle_phase == CELL_CYCLE_G0) then
	if (cell_list(kcell)%cellml_state(58) >= cyclin_threshold) then
		cell_list(kcell)%cell_cycle_phase = CELL_CYCLE_G1
	endif
else
	do i = 1,4
		if (tcycle <= cell_cycle_endtime(i)) then
			cell_list(kcell)%cell_cycle_phase = i
			if (i == 1) then
				cell_list(kcell)%volume = 1 + cell_cycle_rate*tcycle
			elseif (i == 2) then
				cell_list(kcell)%volume = 1 + cell_cycle_rate*cell_cycle_endtime(1)
			elseif (i == 3) then
				cell_list(kcell)%volume = 1 + cell_cycle_rate*(cell_cycle_endtime(1) + tcycle)
			elseif (i == 4) then
				cell_list(kcell)%volume = 2
			endif
			exit
		endif
		cell_list(kcell)%cell_cycle_phase = CELL_CYCLE_D
	enddo
endif
end subroutine

!-----------------------------------------------------------------------------------------
! Cell growth, death and division are handled here.  Division occurs when cell volume 
! exceeds the divide volume. 
! As the cell grows we need to adjust both Cin and Cex to maintain mass conservation.
!
! Nazanin's CELLML_CELL_CYCLE model:
! 0 = time
! 1 = growth
! 2 = size
! 3 = colour
! 4 = signal
! 5 = R = growth rate
! 6 = R_max
! 7 = K_m = Hill K
! 8 = n = Hill n
!-----------------------------------------------------------------------------------------
subroutine CellGrowth(dt,ok)
real(REAL_KIND) :: dt
logical :: ok
integer :: kcell, nlist0, site(3)
integer :: divide_list(10000), ndivide, i, xmax, xmin, dx
real(REAL_KIND) :: interval_dt, cellml_dt
real(REAL_KIND) :: size0, size1, g0, g1, r0, r1, tcycle
real(REAL_KIND) :: signal, g_rate, R_max, Hill_Km, Hill_n

real(REAL_KIND) :: tnow, C_O2, metab, dVdt, vol0	!, r_mean, c_rate
real(REAL_KIND) :: Vin_0, Vex_0, dV
real(REAL_KIND) :: Cin_0(MAX_CHEMO), Cex_0(MAX_CHEMO)
character*(20) :: msg
logical :: C_option = 1

ok = .true.
nlist0 = nlist
tnow = istep*DELTA_T
ndivide = 0

xmax = 0
xmin = 9999
do kcell = 1,nlist0
	xmax = max(xmax,cell_list(kcell)%site(1))
	xmin = min(xmin,cell_list(kcell)%site(1))
enddo

cellml_dt = 0.05
interval_dt = dt/3600.	! Cell_cycle uses time unit = hour
do kcell = 1,nlist0
!	write(*,*) 'kcell: ',kcell
	if (cell_list(kcell)%state == DEAD) cycle
	if (CellML_model == CELLML_CELL_CYCLE) then
		g0 = cell_list(kcell)%cellml_state(1)
		r0 = cell_list(kcell)%cellml_state(5)
		size0 = cell_list(kcell)%cellml_state(2)
		dx = xmax - cell_list(kcell)%site(1)
		signal = signal_max*exp(-signal_decay_coef*dx)	
		cell_list(kcell)%cellml_state(4) = signal
!		g_rate = 0.12	! 5
!		R_max = 0.1		! 6
!		Hill_Km = 0.33	! 7
!		Hill_n = 3		! 8
!		if (kcell == -5) then
!			write(nfout,'(a,2i6,9f6.2)') 'R: ',istep,kcell,cell_list(kcell)%cellml_state(0:8) 
!			write(*,'(a,2i6,9f6.2)') 'R: ',istep,kcell,cell_list(kcell)%cellml_state(0:8)
!		endif
	elseif (CellML_model == CELLML_WNT_CYCLIN) then
		tcycle = tnow - cell_list(kcell)%cell_cycle_entry_time
		cell_list(kcell)%cell_cycle_t = tcycle
		call UpdateCellCycle(kcell,tcycle)
		dx = xmax - cell_list(kcell)%site(1)
		signal = signal_max*exp(-signal_decay_coef*dx)	
		cell_list(kcell)%cellml_state(17) = signal	! Wnt
	endif
	call resetIntegrator
	call setState(cell_list(kcell)%cellml_state)
	call multiStep(0.0, interval_dt, cellml_dt);
	call getState(cell_list(kcell)%cellml_state)
	if (ReadyToDivide(kcell)) then
	    ndivide = ndivide + 1
	    divide_list(ndivide) = kcell
	endif
!	g1 = cell_list(kcell)%cellml_state(1)
!	r1 = cell_list(kcell)%cellml_state(5)
!	size1 = cell_list(kcell)%cellml_state(2)
enddo

do i = 1,ndivide
    kcell = divide_list(i)
	kcell_dividing = kcell
	call CellDivider(kcell, ok)
	if (.not.ok) return
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! The criterion for cell division currently needs to be hard-wired.
! cell.growth is variable(1)
!-----------------------------------------------------------------------------------------
logical function ReadyToDivide(kcell) result(ready)
integer :: kcell

if (CellML_model == CELLML_CELL_CYCLE) then
	if (cell_list(kcell)%cellml_state(1) > divide_growth) then
		ready = .true.
		! Here fix the divide_size overshoot.  This needs attention!
		cell_list(kcell)%cellml_state(2) = 2.0
	else
		ready = .false.
	endif
elseif (CellML_model == CELLML_WNT_CYCLIN) then
	ready = (cell_list(kcell)%cell_cycle_phase == CELL_CYCLE_D)
endif
end function

!-----------------------------------------------------------------------------------------
! Note: updating of concentrations is now done in extendODEdiff()
! The dividing cell, kcell0, is at site0.
! A neighbour site site01 is chosen randomly.  This becomes the site for the daughter cell.
!-----------------------------------------------------------------------------------------
subroutine CellDivider(kcell0, ok)
integer :: kcell0
logical :: ok
integer :: kpar=0
integer :: j, k, kcell1, site0(3), site1(3), site2(3), site01(3), site(3), ichemo, nfree, bestsite(3)
integer :: npath, path(3,200)
real(REAL_KIND) :: tnow, R, size ! v, vmax, V0, Cex(MAX_CHEMO), M0(MAX_CHEMO), M1(MAX_CHEMO), alpha(MAX_CHEMO)
logical :: freesite(27,3)
type (boundary_type), pointer :: bdry

if (kcell0 == -1) then
	dbug = .true.
	write(logmsg,*) 'CellDivider: ',kcell0
	call logger(logmsg)
endif
ok = .true.
tnow = istep*DELTA_T
cell_list(kcell0)%t_divide_last = tnow
!V0 = cell_list(kcell0)%volume
!cell_list(kcell0)%volume = V0/2
!R = par_uni(kpar)
!cell_list(kcell0)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
!cell_list(kcell0)%M = cell_list(kcell0)%M/2
!write(logmsg,'(a,f6.1)') 'Divide time: ',tnow/3600
!call logger(logmsg)

if (CellML_model == CELLML_CELL_CYCLE) then
	! Set kcell0 state to just-divided
	size = cell_list(kcell0)%cellml_state(2)/2
	cell_list(kcell0)%cellml_state(0) = 0		! time
	cell_list(kcell0)%cellml_state(1) = 1		! growth
	cell_list(kcell0)%cellml_state(2) = size	! size
	cell_list(kcell0)%volume = size
elseif (CellML_model == CELLML_WNT_CYCLIN) then
	cell_list(kcell0)%cell_cycle_phase = CELL_CYCLE_G0
	cell_list(kcell0)%volume = 1
	cell_list(kcell0)%cellml_state = state0
endif

site0 = cell_list(kcell0)%site
if (divide_option == DIVIDE_USE_CLEAR_SITE .or. &			! look for the best nearby clear site, if it exists use it
	divide_option == DIVIDE_USE_CLEAR_SITE_RANDOM) then		! make random choice from nearby clear sites
!	vmax = -1.0e10
	nfree = 0
	do j = 1,27
		if (j == 14) cycle
		site01 = site0 + jumpvec(:,j)
		if (occupancy(site01(1),site01(2),site01(3))%indx(1) == 0) then
			nfree = nfree + 1
			freesite(nfree,:) = site01
!			v = SiteValue(occupancy(site01(1),site01(2),site01(3))%C(:))
!			if (v > vmax) then
!				vmax = v
				bestsite = site01
!			endif
		endif
	enddo
	if (nfree > 0) then
		if (DIVIDE_USE_CLEAR_SITE) then	! use this site for the progeny cell
			site01 = bestsite
		else	! random choice
			j = random_int(1,nfree,0)
			site01 = freesite(j,:)
		endif
		call AddCell(kcell0,kcell1,site01,ok)
		if (.not.ok) then
			call logger('Error: AddCell: vacant site')
		endif
		Nreuse = Nreuse + 1
		return
	endif
endif

do
	j = random_int(1,6,kpar)
	site01 = site0 + neumann(:,j)
	if (site01(2) > ywall) exit
enddo
!if (dbug) write(*,*) 'CellDivider: ',kcell0,site0,occupancy(site0(1),site0(2),site0(3))%indx
if (occupancy(site01(1),site01(2),site01(3))%indx(1) == OUTSIDE_TAG) then	! site01 is outside, use it directly
	npath = 0
	site1 = site0
elseif (bdrylist_present(site01,bdrylist)) then	! site01 is on the boundary
	npath = 1
	site1 = site01
	path(:,1) = site01
else
	call ChooseBdrysite(site01,site1)
	if (occupancy(site1(1),site1(2),site1(3))%indx(1) == 0) then
		write(*,*) 'after ChooseBdrysite: site1 is VACANT: ',site1
		stop
	endif
	if (dbug) write(*,'(a,i6,9i4)') 'b4 ',kcell0,site0,site01,site1
	call SelectPath(site0,site01,site1,path,npath)
endif

if (npath > 0) then
	! Need to choose an outside or vacant site near site1
	call getOutsideSite(site1,site2)
	npath = npath+1
	path(:,npath) = site2
	! path(:,:) now goes from site01 to site2, which is an outside site next to site1
	if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
		Nsites = Nsites + 1
		call RemoveFromMedium
	endif
!	write(*,'(a,3i4,i6)') 'outside site: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
else
	! site01 is the destination site for the new cell
	site2 = site01
	Nsites = Nsites + 1
endif

!call GetPathMass(site0,site01,path,npath,M0)
if (npath > 0) then
	call PushPath(path,npath)
	if (dbug) write(*,*) 'did push_path'
endif
call AddCell(kcell0,kcell1,site01,ok)
if (.not.ok) then
	call logger('Error: AddCell: pushed site')
	return
endif
!if (vary_conc) then
!	call GetPathMass(site0,site01,path,npath,M1)
!	alpha = M0/M1	! scaling for concentrations on the path
!	call ScalePathConcentrations(site0,site01,path,npath,alpha)
!endif

call SetRadius(Nsites)
! Now need to fix the bdrylist.  
! site1 was on the boundary, but may no longer be.
! site2 may be now on the boundary
! First add site2
if (isbdry(site2)) then   ! add it to the bdrylist
	if (dbug) write(*,*) 'add site2 to bdrylist: ',site2
    allocate(bdry)
    bdry%site = site2
    nullify(bdry%next)
    call bdrylist_insert(bdry,bdrylist)
    occupancy(site2(1),site2(2),site2(3))%bdry => bdry
!    call SetBdryConcs(site2)
!else
!    write(logmsg,'(a,3i4,i6)') 'Added site is not bdry: ',site2,occupancy(site2(1),site2(2),site2(3))%indx(1)
!	call logger(logmsg)
!    call SetBdryConcs(site2)
!    stop
endif
if (dbug) write(*,*) 'Check for changed boundary status'
! Now check sites near site2 that may have changed their boundary status (including site1)
do j = 1,6
	site = site2 + neumann(:,j)
	if (dbug) write(*,*) j,site
	if (isbdry(site)) then
		if (dbug) write(*,*) 'isbdry'
		if (.not.bdrylist_present(site,bdrylist)) then	! add it
			if (dbug) write(*,*) 'not present, add it'
			allocate(bdry)
			bdry%site = site
			nullify(bdry%next)
			call bdrylist_insert(bdry,bdrylist)
			occupancy(site(1),site(2),site(3))%bdry => bdry
		endif
	else
		if (dbug) write(*,*) 'not isbdry'
		if (bdrylist_present(site,bdrylist)) then	! remove it
			if (dbug) write(*,*) 'present, remove it'
			call bdrylist_delete(site,bdrylist)
			nullify(occupancy(site(1),site(2),site(3))%bdry)
		endif
	endif
enddo
if (dbug) call logger('did CellDivider')
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine AdjustMedium

end subroutine

!-----------------------------------------------------------------------------------------
! Need to choose a site on the boundary in some sense near site0, and also to preserve
! the spherical shape
!-----------------------------------------------------------------------------------------
subroutine ChooseBdrysite(site0,site1)
integer :: site0(3), site1(3)
integer :: site(3), sitemin(3)
real(REAL_KIND) :: vc(3), v(3), r, rfrac, d, alpha_max, cosa_min, dmin, cosa
logical :: hit
type (boundary_type), pointer :: bdry

if (dbug) write(*,*) 'ChooseBdrysite: ',site0
vc = site0 - Centre
r = norm(vc)
vc = vc/r
if (ywall == 0) then
	rfrac = r/Radius
	alpha_max = getAlphaMax(rfrac)
	cosa_min = cos(alpha_max)
else
	cosa_min = 0	! for the half-plane
endif
dmin = 1.0e10
hit = .false.
bdry => bdrylist
do while ( associated ( bdry )) 
    site = bdry%site
    v = site - Centre
    d = norm(v)
    cosa = dot_product(v,vc)/d
    if (cosa > cosa_min) then
		hit = .true.
		if (d < dmin) then
			dmin = d
			sitemin = site
		endif
	endif
    bdry => bdry%next
enddo
if (.not.hit) then
	write(logmsg,*) 'Error: ChooseBdrysite: no candidate bdry site'
	call logger(logmsg)
	stop
endif
site1 = sitemin
if (dbug) write(*,*) 'site1: ',site1,occupancy(site1(1),site1(2),site1(3))%indx
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function getAlphaMax(r)
real(REAL_KIND) :: r
real(REAL_KIND) :: alf
real(REAL_KIND) :: r1 = 0.0, r2 = 0.8, alpha1 = PI/2, alpha2 = PI/6

if (r < r1) then
	getAlphaMax = alpha1
elseif (r > r2) then
	getAlphaMax = alpha2
else
	alf = (r-r1)/(r2-r1)
	getAlphaMax = (1-alf)*alpha1 + alf*alpha2
endif
end function


!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine test_get_path
integer :: site0(3), site1(3), site2(3), path(3,200), npath, kpar=0
integer :: x, y, z, it, k

site0 = 0
do it = 1,10
	do
		x = random_int(1,NX,kpar)
		y = random_int(1,NY,kpar)
		z = random_int(1,NZ,kpar)
		if (occupancy(x,y,z)%indx(1) > 0) exit
	enddo
	site1 = (/x,y,z/)
	do
		x = random_int(1,NX,kpar)
		y = random_int(1,NY,kpar)
		z = random_int(1,NZ,kpar)
		if (occupancy(x,y,z)%indx(1) > 0) exit
	enddo
	site2 = (/x,y,z/)
	call SelectPath(site0,site1,site2,path,npath)
	write(*,*) 'path: ',npath
	do k = 1,npath
		write(*,'(i3,2x,3i4)') path(:,k)
	enddo
enddo
end subroutine

!-----------------------------------------------------------------------------------------
! Find a path of sites leading from site1 to site2
! The first site in the list is site1, the last is site2
! Previously site2 was always a bdry site.  Now we stop the path when a vacant or outside site 
! is encountered.
! How do we avoid choosing a path through the site of the dividing cell? site0
!-----------------------------------------------------------------------------------------
subroutine SelectPath(site0,site1,site2,path,npath)
integer :: site0(3),site1(3), site2(3), path(3,200), npath
integer :: v(3), jump(3), site(3), k, j, jmin
real(REAL_KIND) :: r, d2, d2min
logical :: hit

if (occupancy(site2(1),site2(2),site2(3))%indx(1) == OUTSIDE_TAG) then
    call logger('site2 is OUTSIDE')
    stop
endif
if (occupancy(site2(1),site2(2),site2(3))%indx(1) == 0) then
    call logger('site2 is VACANT')
    stop
endif

k = 1
site = site1
path(:,k) = site
hit = .false.
do 
	d2min = 9999
	jmin = 0
	do j = 1,27
		if (j == 14) cycle
		jump = jumpvec(:,j)
		v = site + jump
		if (v(1)==site0(1) .and. v(2)==site0(2) .and. v(3)==site0(3)) cycle
!		if (occupancy(v(1),v(2),v(3))%indx(1) == OUTSIDE_TAG) then
		if (occupancy(v(1),v(2),v(3))%indx(1) <= 0) then	! outside or vacant - we'll use this!
			site2 = site
			hit = .true.
			exit
		endif
!		if (occupancy(v(1),v(2),v(3))%indx(1) < -100) cycle		! necrotic
		v = site2 - v
		d2 = v(1)*v(1) + v(2)*v(2) + v(3)*v(3)
		if (dbug) write(*,'(8i6,2f8.1)') k,j,site+jump,v,d2,d2min
		if (d2 < d2min) then
			d2min = d2
			jmin = j
		endif
	enddo
	if (hit) exit
	if (jmin == 0) then
	    call logger('get path: stuck')
	    stop
	endif
	site = site + jumpvec(:,jmin)
	k = k+1
	if (dbug) write(*,'(11i4,f8.1)') k,site1,site2,site,jmin,d2min
	if (k==30) then
		call logger('Error: SelectPath: too many iterations: k = 30')
		stop
	endif
	path(:,k) = site
	if (site(1) == site2(1) .and. site(2) == site2(2) .and. site(3) == site2(3)) exit
enddo
npath = k

end subroutine

!-----------------------------------------------------------------------------------------
! If there are any necrotic sites in the path, they are first shifted closer to the centre,
! then the path is adjusted.
!-----------------------------------------------------------------------------------------
subroutine ClearPath(path,npath)
integer :: npath, path(3,*)
logical :: clear
integer :: k, kcell, kvacant, site(3)

do
	clear = .true.
	do k = 1,npath
		site = path(:,k)
		kcell = occupancy(site(1),site(2),site(3))%indx(1)
		if (kcell < 0) then
			clear = .false.
			kvacant = k
			exit
		endif
	enddo
	if (clear) return
	
enddo			
end subroutine

!-----------------------------------------------------------------------------------------
! Need to modify the ODEdiff variables, at least %ivar(?)
!-----------------------------------------------------------------------------------------
subroutine PushPath(path,npath)
integer :: path(3,200),npath
integer :: k, site1(3), site2(3), kcell

do k = npath-1,1,-1
	site1 = path(:,k)
	kcell = occupancy(site1(1),site1(2),site1(3))%indx(1)
	if (dbug) write(*,*) k,' site1: ',site1,kcell
	site2 = path(:,k+1)
	if (dbug) write(*,*) 'site2: ',site2
	if (kcell > 0) then
    	cell_list(kcell)%site = site2
    endif
	occupancy(site2(1),site2(2),site2(3))%indx(1) = kcell
enddo
occupancy(site1(1),site1(2),site1(3))%indx = 0
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine getGrowthRateParameters(c_rate, r_mean)
real(REAL_KIND) :: c_rate, r_mean
real(REAL_KIND) :: divide_time, R
integer :: kpar = 0
real(REAL_KIND), parameter :: divide_time_variation = 0.2
if (randomise_divide_time) then
	R = par_uni(kpar)
	divide_time = ((1 - divide_time_variation) + 2*R*divide_time_variation)*divide_time_mean
else
	divide_time = divide_time_mean
endif
c_rate = log(2.0)/divide_time
r_mean = Vdivide0/(2*divide_time)
end subroutine

!-----------------------------------------------------------------------------------------
! The daughter cell kcell1 is given the same characteristics as kcell0 and placed at site1.
! Random variation is introduced into %divide_volume.
!-----------------------------------------------------------------------------------------
subroutine AddCell(kcell0,kcell1,site1,ok)
integer :: kcell0, kcell1, site1(3)
logical :: ok
integer :: kpar = 0
real(REAL_KIND) :: tnow, R, c_rate,r_mean

ok = .true.
tnow = istep*DELTA_T
nlist = nlist + 1
if (nlist > max_nlist) then
	call logger('Dimension of cell_list() has been exceeded: increase max_nlist and rebuild')
	ok = .false.
	return
endif
Ncells = Ncells + 1
kcell1 = nlist
if (site1(2) <= ywall) then
	write(*,*) 'AddCell: error: y < ywall: ',kcell0,kcell1,site1
	write(nflog,*) 'AddCell: error: y < ywall: ',kcell0,kcell1,site1
	stop
endif
allocate(cell_list(kcell1)%cellml_state(0:nvariables-1))
cell_list(kcell1)%cellml_state = cell_list(kcell0)%cellml_state
cell_list(kcell1)%celltype = cell_list(kcell0)%celltype
cell_list(kcell1)%state = cell_list(kcell0)%state
cell_list(kcell1)%site = site1
cell_list(kcell1)%ID = cell_list(kcell0)%ID
cell_list(kcell1)%exists = .true.
cell_list(kcell1)%active = .true.
cell_list(kcell1)%t_divide_last = tnow
cell_list(kcell1)%volume = cell_list(kcell0)%volume
if (CellML_model == CELLML_CELL_CYCLE) then
	R = par_uni(kpar)
	cell_list(kcell1)%divide_volume = Vdivide0 + dVdivide*(2*R-1)
endif
occupancy(site1(1),site1(2),site1(3))%indx(1) = kcell1
end subroutine

!-----------------------------------------------------------------------------------------
! Look at all bdry sites (those with a neumann neighbour outside)
!-----------------------------------------------------------------------------------------
subroutine check_bdry
integer :: kcell, i, site(3), site1(3), minv(3), maxv(3)
real(REAL_KIND) :: v(3), r, rmin, rmax
logical :: bdry

rmin = 1.0e10
rmax = 0
do kcell = 1,nlist
	site = cell_list(kcell)%site
	if (site(2) <= ywall) then
		write(*,*) 'check_bdry: Error: y <= ywall: ',kcell, site
		write(nflog,*) 'check_bdry: Error: y <= ywall: ',kcell, site
		stop
	endif
	bdry = .false.
	do i = 1,6
		site1 = site + neumann(:,i)
		if (occupancy(site1(1),site1(2),site1(3))%indx(1) == OUTSIDE_TAG) then
			bdry = .true.
			exit
		endif
	enddo
	if (bdry) then
		v = site - Centre
		r = norm(v)
		if (r < rmin) then
			rmin = r
			minv = site - Centre
		endif
		if (r > rmax) then
			rmax = r
			maxv = site - Centre
		endif
	endif
enddo
write(*,'(a,2(f7.1,3i4,2x))') 'rmin, rmax: ',rmin,minv,rmax,maxv
end subroutine


!-----------------------------------------------------------------------------------------
! The location site0 is just outside the blob.
! The aim is to find a cell near site0 that is further from the centre, and move it here.
!-----------------------------------------------------------------------------------------
subroutine Adjust(site0)
integer :: site0(3)
integer :: i, site(3), kcell, sitemax(3)
real(REAL_KIND) :: r0, r, rmax

if (occupancy(site0(1),site0(2),site0(3))%indx(1) /= OUTSIDE_TAG) then
	write(*,*) 'Error: adjust: site is not OUTSIDE: ',site0
	stop
endif
r0 = cdistance(site0)
!write(*,'(a,3i4,f6.2)') 'adjust: ',site0,r0
rmax = 0
do i = 1,6
	site = site0 + neumann(:,i)
	kcell = occupancy(site(1),site(2),site(3))%indx(1)
	if (kcell > 0) then
		r = cdistance(site)
		if (r > r0 .and. r > rmax) then	! move the cell here
			rmax = r
			sitemax = site
		endif
	endif
enddo
if (rmax > 0) then
	kcell = occupancy(sitemax(1),sitemax(2),sitemax(3))%indx(1)
!	write(*,'(i6,2x,3i4,f6.2)') kcell,sitemax,rmax
	cell_list(kcell)%site = site0
	occupancy(site0(1),site0(2),site0(3))%indx(1) = kcell
	occupancy(sitemax(1),sitemax(2),sitemax(3))%indx(1) = OUTSIDE_TAG
endif
end subroutine

end module
