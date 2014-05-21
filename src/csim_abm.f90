module csim_abm

use, intrinsic :: iso_c_binding
implicit none

interface

    subroutine getNvariables(nv) BIND(C,NAME='getNvariables')
	use :: iso_c_binding
	integer(c_int) :: nv
	end subroutine
	
	subroutine setupCellml(URI) BIND(C,NAME='setupCellml')
	use :: iso_c_binding
	character(c_char),dimension(*) :: URI
	end subroutine
	
	subroutine setState(state) BIND(C,NAME='setState')
	use :: iso_c_binding
	real(c_double) :: state(*)
	end subroutine

	subroutine getState(state) BIND(C,NAME='getState')
	use :: iso_c_binding
	real(c_double) :: state(*)
	end subroutine

	subroutine setStateValue(k, val) BIND(C,NAME='setStateValue')
	use :: iso_c_binding
	integer(c_int),value :: k
	real(c_double),value :: val
	end subroutine

	subroutine getID(k,ID,clen) BIND(C,NAME='getID')
	use :: iso_c_binding
	integer(c_int), value :: k
	integer(c_int) :: clen
	character(c_char),dimension(*) :: ID
	end subroutine

	subroutine singleStep(stepSize) BIND(C,NAME='singleStep')
	use :: iso_c_binding
	real(c_double), value :: stepSize
	end subroutine

	subroutine multiStep(startTime, endTime, dt) BIND(C,NAME='multiStep')
	use :: iso_c_binding
	real(c_double), value :: startTime, endTime, dt
	end subroutine

	subroutine resetIntegrator() BIND(C,NAME='resetIntegrator')
	use :: iso_c_binding
	end subroutine

end interface

type ID_type
	character(64) :: comp_name
	character(64) :: var_name
	integer :: comp_index
end type

type(ID_type), allocatable :: IDlist(:)
integer :: nvariables, ncomponents
character(64), allocatable :: component(:)
real(c_double), allocatable :: csim_state(:), csim_savestate(:)
	
end module