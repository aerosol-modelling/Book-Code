module Mod_BAT_model    
    
use Mod_NumPrec_Types, only : wp

implicit none
private     !by default, make variables and procedures private

!public procedures (accessible by using this module):
public :: BAT_light, density_est !, ... 

!-----------------------------------------------------------
contains

pure subroutine BAT_light(x_org, M_org, OtoC, rho_org, ln_actcoeff, activity)
! ... (content not shown)
end subroutine BAT_light
!-----------------------------------------------------------

pure elemental function density_est(M, OtoC, HtoC) result(rho)
! ... (content not shown)
end function density_est
!-----------------------------------------------------------

! ... (additional procedures)
	
end module Mod_BAT_model