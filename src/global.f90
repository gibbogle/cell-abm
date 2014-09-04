! Global definitions

module global

use omp_lib
use par_zig_mod
use winsock
use, intrinsic :: ISO_C_BINDING

implicit none

!INTEGER,  PARAMETER  ::  DP=SELECTED_REAL_KIND( 12, 60 )
INTEGER,  PARAMETER  ::  SP = kind(1.0), DP = kind(1.0d0)
integer, parameter :: REAL_KIND = DP
integer, parameter :: TCP_PORT_0 = 5000		! main communication port (logging)
integer, parameter :: TCP_PORT_1 = 5001		! data transfer port (plotting)
integer, parameter :: NORMAL_DIST      = 1
integer, parameter :: LOGNORMAL_DIST   = 2
integer, parameter :: EXPONENTIAL_DIST = 3
integer, parameter :: CONSTANT_DIST    = 4

integer, parameter :: DIVIDING  = 1
integer, parameter :: QUIESCENT = 2
integer, parameter :: DEAD      = 3

integer, parameter :: OUTSIDE_TAG  = -1
integer, parameter :: INACCESSIBLE_TAG  = -2

integer, parameter :: DIVIDE_ALWAYS_PUSH  = 1
integer, parameter :: DIVIDE_USE_CLEAR_SITE  = 2
integer, parameter :: DIVIDE_USE_CLEAR_SITE_RANDOM  = 3

integer, parameter :: nfin=10, nfout=11, nflog=12, nfres=13, nfrun=14, nfcell=15
integer, parameter :: neumann(3,6) = reshape((/ -1,0,0, 1,0,0, 0,-1,0, 0,1,0, 0,0,-1, 0,0,1 /), (/3,6/))

integer, parameter :: MAX_CHEMO = 3
integer, parameter :: OXYGEN = 1
integer, parameter :: GLUCOSE = 2
integer, parameter :: TRACER = 3
integer, parameter :: EXTRA = 1
integer, parameter :: INTRA = 2
integer, parameter :: MAX_CELLTYPES = 4
integer, parameter :: max_nlist = 500000

logical, parameter :: use_ODE_diffusion = .true.
logical, parameter :: compute_concentrations = .true.
logical, parameter :: use_division = .true.
logical, parameter :: use_death = .true.
logical, parameter :: use_react = .true.
logical, parameter :: use_migration = .false.	! causing an error with vacant site becoming bdry 
logical, parameter :: use_medium_flux = .true.	! flux of constituents between spheroid and medium is accounted for.
logical, parameter :: use_celltype_colour = .true.

logical, parameter :: use_Cex_Cin = .true.		! assume equilibrium to derive Cin from Cex
logical, parameter :: suppress_growth = .false.
logical, parameter :: randomise_divide_time = .true.
logical, parameter :: vary_conc = .false.

integer, parameter :: CELL_CYCLE_G0 = 0
integer, parameter :: CELL_CYCLE_G1 = 1
integer, parameter :: CELL_CYCLE_S  = 2
integer, parameter :: CELL_CYCLE_G2 = 3
integer, parameter :: CELL_CYCLE_M  = 4
integer, parameter :: CELL_CYCLE_D  = 5		! due to divide

real(REAL_KIND), parameter :: cell_cycle_fraction(4) = [0.4, 0.2, 0.3, 0.1]

integer, parameter :: CELLML_CELL_CYCLE = 0
integer, parameter :: CELLML_WNT_CYCLIN = 1

real(REAL_KIND), parameter :: PI = 4.0*atan(1.0)

type occupancy_type
	integer :: indx(2)
	real(REAL_KIND) :: C(MAX_CHEMO)
	type(boundary_type), pointer :: bdry
end type

type cell_type
	integer :: ID
	integer :: celltype
	integer :: site(3)
	integer :: iv
	logical :: active
	integer :: state
	integer :: cell_cycle_phase
	real(REAL_KIND) :: cell_cycle_entry_time
	real(REAL_KIND) :: cell_cycle_t
	logical :: divide_flag
	real(REAL_KIND) :: conc(MAX_CHEMO)
	real(REAL_KIND) :: dVdt
	real(REAL_KIND) :: c_rate, r_mean
	real(REAL_KIND) :: volume			! fractional volume (fraction of nominal cell volume Vcell_cm3)
	real(REAL_KIND) :: divide_volume
	real(REAL_KIND) :: t_divide_last
	real(REAL_KIND) :: t_hypoxic
	real(REAL_KIND) :: t_anoxia_die
	real(REAL_KIND) :: M
	logical :: anoxia_tag
	logical :: exists
	real(REAL_KIND), allocatable :: cellml_state(:)
end type

type boundary_type
    integer :: site(3)
    type (boundary_type), pointer :: next
end type

type dist_type
	integer :: class
	real(REAL_KIND) :: p1, p2, p3
end type

type, bind(C) :: field_data
	integer(c_int) :: site(3)
	integer(c_int) :: state
	real(c_double) :: dVdt
	real(c_double) :: volume
	real(c_double) :: conc(MAX_CHEMO+1)	! This needs to agree with the definition in field.h
end type

type(dist_type) :: divide_dist
type(occupancy_type), allocatable :: occupancy(:,:,:)
type(cell_type), allocatable :: cell_list(:)
real(REAL_KIND),allocatable :: state0(:)

integer :: NX, NY, NZ
integer :: initial_count
integer, allocatable :: zdomain(:),zoffset(:)
integer :: blobrange(3,2)
real(REAL_KIND) :: x0,y0,z0   ! centre in global coordinates (units = grids)
integer :: ywall
real(REAL_KIND) :: Radius, Centre(3)
real(REAL_KIND) :: adrop, bdrop, cdrop	! drop shape transformation parameters
real(REAL_KIND) :: z0drop				! drop centre is at z0drop = zmin + (bdrop-cdrop)*R
integer :: zmin							! drop lower bound at drop time
integer :: jumpvec(3,27)

integer :: nlist, Ncells, Ncells0, lastNcells, lastID, Ncelltypes
integer :: max_ngaps, ngaps, nadd_sites, Nsites, Nreuse
integer :: Nanoxia_tag, Nanoxia_dead
integer :: istep, nsteps, it_solve, NT_CONC, NT_GUI_OUT, show_progeny
integer :: Mnodes, ncpu_input
real(REAL_KIND) :: DELTA_T, DELTA_X, fluid_fraction, Vsite_cm3, Vextra_cm3, Vcell_cm3
real(REAL_KIND) :: medium_volume0, medium_volume, cell_radius, d_layer
real(REAL_KIND) :: celltype_fraction(MAX_CELLTYPES)
logical :: celltype_display(MAX_CELLTYPES)
real(REAL_KIND) :: MM_THRESHOLD, ANOXIA_THRESHOLD, t_anoxic_limit, anoxia_death_delay, Vdivide0, dVdivide
real(REAL_KIND) :: anoxia_tag_hours, anoxia_death_hours
real(REAL_KIND) :: divide_time_median, divide_time_shape, divide_time_mean
real(REAL_KIND) :: t_simulation, execute_t1
real(REAL_KIND) :: O2cutoff(3)
real(REAL_KIND) :: growthcutoff(3)
real(REAL_KIND) :: spcrad_value
real(REAL_KIND) :: signal_max, signal_decay_coef
logical :: bdry_changed
character*(128) :: inputfile
character*(128) :: outputfile
character*(128) :: cellmlfile
character*(2048) :: logmsg
logical :: test_case(4)

real(REAL_KIND) :: divide_growth = 2.2

real(REAL_KIND) :: cyclin_threshold
real(REAL_KIND) :: cell_cycle_duration
real(REAL_KIND) :: cell_cycle_rate
real(REAL_KIND) :: cell_cycle_endtime(4)

TYPE(winsockport) :: awp_0, awp_1
logical :: use_TCP = .true.         ! turned off in para_main()
logical :: use_CPORT1 = .false.
logical :: stopped, clear_to_send
logical :: simulation_start, par_zig_init, initialized
logical :: use_extracellular_O2
logical :: use_V_dependence
logical :: randomise_initial_volume
logical :: relax
logical :: use_parallel
logical :: dbug = .false.

real(REAL_KIND) :: ysave(100000),dCreactsave(100000)

integer :: divide_option = DIVIDE_USE_CLEAR_SITE
integer :: idbug = 0
integer :: Nbnd
integer :: seed(2)

real(REAL_KIND), allocatable :: omp_x(:), omp_y(:), omp_z(:)

integer, parameter :: CellML_model = CELLML_WNT_CYCLIN

!DEC$ ATTRIBUTES DLLEXPORT :: nsteps, DELTA_T

contains

!-----------------------------------------------------------------------------------------
! WTIME returns a reading of the wall clock time.
!-----------------------------------------------------------------------------------------
real(DP) function wtime()
!DEC$ ATTRIBUTES DLLEXPORT :: wtime
  integer :: clock_max, clock_rate, clock_reading

  call system_clock ( clock_reading, clock_rate, clock_max )
  wtime = real(clock_reading,kind=DP)/clock_rate
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine logger(msg)
character*(*) :: msg
integer :: error
logical :: logfile_isopen
character*(1) :: LF = char(94)

error = 0
inquire(unit=nflog,OPENED=logfile_isopen)
if (use_TCP) then
    if (awp_0%is_open) then
        call winsock_send(awp_0,trim(msg)//LF,len_trim(msg)+1,error)
    elseif (logfile_isopen) then
        write(nflog,*) trim(msg)
    else
        write(99,*) trim(msg)
    endif
else
	write(*,*) trim(msg)
endif
if (logfile_isopen) then
	write(nflog,*) 'msg: ',trim(msg)
	if (error /= 0) then
	    write(nflog,'(a,i4)') 'winsock_send error: ',error
	    close(nflog)
	endif
endif
if (error /= 0) stop
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine make_jumpvec
integer :: k,ix,iy,iz

k = 0
do ix = -1,1
	do iy = -1,1
		do iz = -1,1
			k = k+1
			jumpvec(:,k) = (/ix,iy,iz/)
		enddo
	enddo
enddo
end subroutine

!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine SetRadius(N)
integer :: N
Radius = (3.0*N/(4.0*PI))**(1./3.)
end subroutine

!---------------------------------------------------------------------
! Uniform randomly generates an integer I: n1 <= I <= n2
!---------------------------------------------------------------------
integer function random_int(n1,n2,kpar)
integer :: n1,n2,kpar
integer :: k,R

if (n1 == n2) then
    random_int = n1
elseif (n1 > n2) then
    write(logmsg,*) 'ERROR: random_int: n1 > n2: ',n1,n2
    call logger(logmsg)
    stop
endif
R = par_shr3(kpar)
k = abs(R)
random_int = n1 + mod(k,(n2-n1+1))
end function

!--------------------------------------------------------------------------------
! Make a random choice of an integer from 1 - N on the basis of probabilities in
! the array p(:) (assumed to be normalized).
!--------------------------------------------------------------------------------
integer function random_choice(p,N,kpar)
integer :: N,kpar
real(REAL_KIND) :: p(:)
integer :: k
real(REAL_KIND) :: R, psum

R = par_uni(kpar)
psum = 0
do k = 1,N
    psum = psum + p(k)
    if (R <= psum) then
        random_choice = k
        return
    endif
enddo
write(logmsg,*) 'ERROR: random_choice: ',N,p
call logger(logmsg)
stop
end function

!--------------------------------------------------------------------------------
! Returns a permutation of the elements of a()
!--------------------------------------------------------------------------------
subroutine permute(a,n,kpar)
integer :: a(*),n,kpar
integer :: i,k,tmp

do i = 1,n
    k = random_int(1,n,kpar)
	tmp = a(i)
	a(i) = a(k)
	a(k) = tmp
enddo
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine waste_time(n,dummy)
integer :: k, n
real(REAL_KIND) :: dummy
real(REAL_KIND) :: rsum,R
integer :: kpar=0

rsum = 0
do k = 1,n
    R = par_uni(kpar)
    rsum = rsum + R
enddo
dummy = rsum
end subroutine

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm(r)
real(REAL_KIND) :: r(3)

norm = sqrt(r(1)*r(1) + r(2)*r(2) + r(3)*r(3))
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function norm2(r)
real(REAL_KIND) :: r(3)

norm2 = r(1)*r(1) + r(2)*r(2) + r(3)*r(3)
end function

!-----------------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------------
subroutine normalize(r)
real(REAL_KIND) :: r(3)

r = r/norm(r)
end subroutine

!-----------------------------------------------------------------------------------------
! Distance from the blob centre (units = grids)
!-----------------------------------------------------------------------------------------
real(REAL_KIND) function cdistance(site)
integer :: site(3)
real(REAL_KIND) :: r(3)

r = site - Centre
cdistance = norm(r)
end function

!--------------------------------------------------------------------------------
! Crude check for a site inside, just looking at the allowable ranges of x, y and z.
!--------------------------------------------------------------------------------
logical function inside_xyz(site)
integer :: site(3)

if (site(1) < 1 .or. site(1) > NX .or. site(2) < 1 .or. site(2) > NY .or. site(3) < 1 .or. site(3) > NZ) then
    inside_xyz = .false.
else
    inside_xyz = .true.
endif
end function


!---------------------------------------------------------------------
!---------------------------------------------------------------------
subroutine get_vnorm(v,vnorm)
real(REAL_KIND) :: v(3), vnorm(3)
real(REAL_KIND) :: d

d = dot_product(v,v)
vnorm = v/sqrt(d)
end subroutine

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function DivideTime()
real(REAL_KIND) :: p1, p2
integer :: kpar = 0

dividetime = 0
p1 = divide_dist%p1
p2 = divide_dist%p2
select case (divide_dist%class)
case (NORMAL_DIST)
	DivideTime = rv_normal(p1,p2,kpar)
case (LOGNORMAL_DIST)
	DivideTime = rv_lognormal(p1,p2,kpar)
case (CONSTANT_DIST)
	DivideTime = p1
end select
!write(*,'(a,f8.2)') 'DivideTime: ',DivideTime/3600
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real function DivisionTime()
integer :: kpar = 0
real(REAL_KIND), parameter :: rndfraction = 0.2

DivisionTime = rv_lognormal(divide_dist%p1,divide_dist%p2,kpar)
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_normal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R

R = par_rnor(kpar)
rv_normal = p1+R*p2
end function

!--------------------------------------------------------------------------------------
! When Y is normal N(p1,p2) then X = exp(Y) is lognormal with
!   median = m = exp(p1)
!   shape  = s = exp(p2)
! Also median = m = mean/(s^2/2)
! kpar = parallel process number
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_lognormal(p1,p2,kpar)
integer :: kpar
real(REAL_KIND) :: p1,p2
real(REAL_KIND) :: R,z

R = par_rnor(kpar)
z = p1 + R*p2
rv_lognormal = exp(z)
end function

!--------------------------------------------------------------------------------------
! For testing.
!--------------------------------------------------------------------------------------
real(REAL_KIND) function my_rnor()
real(REAL_KIND) :: sum, R
integer :: k
integer :: kpar=0

sum = 0
do k = 1,12
    R = par_uni(kpar)
    sum = sum + R
enddo
my_rnor = sum - 6.0
end function

!--------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------------
real(REAL_KIND) function rv_exponential(p1)
real(REAL_KIND) :: p1
real(REAL_KIND) :: r
integer :: kpar = 0

r = par_rexp(kpar)
rv_exponential = p1*r
end function

!--------------------------------------------------------------------------------------
! Cumulative probability distribution for a lognormal variate with median m, shape s
! Computes Pr(X < a) where X = exp(Y) and Y is N(p1,p2), p1 = log(m), p2 = log(s)
! Pr(X < a) = Pr(Y < log(a)) = Pr(p1 + p2*R < log(a)) = Pr(R < (log(a)-p1)/p2)
! where R is N(0,1)
!--------------------------------------------------------------------------------------
real(REAL_KIND) function cum_prob_lognormal(a,p1,p2)
real(REAL_KIND) :: a, p1, p2
real(REAL_KIND) :: b, prob

b = (log(a) - p1)/p2
prob = 0.5 + 0.5*erf(b/sqrt(2.0))
cum_prob_lognormal = prob
end function

!--------------------------------------------------------------------------------------
! Determine real roots r(:) of the cubic equation:
! x^3 + a.x^2 + b.x + c = 0
! If there is one real root, n=1 and the root is r(1)
! If there are three distinct real roots, n=3 and the roots are r(1), r(2), r(3)
! If there is a repeated root, n=2 and the single root is r(1), the repeated root is r(2)
!--------------------------------------------------------------------------------------
subroutine cubic_roots(a, b, c, r, n)
real(REAL_KIND) :: a, b, c, r(3)
integer :: n
real(REAL_KIND) :: QQ, RR, theta, R2, Q3, AA, BB

QQ = (a*a - 3*b)/9
RR = (2*a*a*a - 9*a*b + 27*c)/54
Q3 = QQ*QQ*QQ
R2 = RR*RR
if (R2 < Q3) then
	n = 3
	theta = acos(RR/sqrt(Q3))
	r(1) = -2*sqrt(QQ)*cos(theta/3) - a/3
	r(2) = -2*sqrt(QQ)*cos((theta+2*PI)/3) - a/3
	r(3) = -2*sqrt(QQ)*cos((theta-2*PI)/3) - a/3
else
	n = 1
	AA = -sign(1.,RR)*(abs(RR) + sqrt(R2 - Q3))**(1.d0/3.d0)
	if (AA == 0) then
		BB = 0
	else
		BB = QQ/AA
	endif
	r(1) = AA + BB - a/3
endif
end subroutine

end module
