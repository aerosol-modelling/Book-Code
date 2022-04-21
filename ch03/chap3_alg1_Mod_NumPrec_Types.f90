!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module defining compiler-independent numerical precision parameters and derived    * 
!*   data types (for global use in a project).                                          * 
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
module Mod_NumPrec_Types

implicit none

!define a working precision (wp) level to be used with floating point (real) variables, e.g. 1.0 should be stated as 1.0_wp.
!number_of_digits = desired minimum level of precision in terms of number of floating point decimal digits.
integer,parameter,private :: number_of_digits = 12        
integer,parameter,public  :: wp = selected_real_kind(number_of_digits)

!definition of public derived types:
public
!type used for setting options of BAT_light wrapper functions
type :: foptions                    
    logical  :: determine_local_max
    integer  :: component_no
    real(wp) :: M_org
    real(wp) :: OtoC 
    real(wp) :: density_org
    real(wp) :: refval
end type foptions

end module Mod_NumPrec_Types