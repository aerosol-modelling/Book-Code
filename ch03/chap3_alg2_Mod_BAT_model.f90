!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing procedures for the Binary Activity Thermodynamics (BAT) model.   *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine  BAT_light                                                           *
!*   -  function    density_est                                                         *
!*   -  function    fnc_BAT_activity                                                    *
!*   -  function    fnc_BAT_act_dev                                                     *
!*   -  subroutine  BAT_miscibility_gap                                                 *
!*   -  subroutine  calc_BAT_LLE_dev                                                    *
!*                                                                                      *
!**************************************************************************************** 
module Mod_BAT_model    
    
use Mod_NumPrec_Types, only : wp

implicit none
private             !by default, make variables and procedures private to this module

!declaration of derived types and module variables:
!custom type used for setting additional BAT model options (for LLE computation)
type :: BAToptions
    real(wp),dimension(2) :: xtot
    real(wp) :: M_org
    real(wp) :: OtoC 
    real(wp) :: density_org
end type BAToptions

!declaration of module variables:
integer,private :: ncalls_LLE
type(BAToptions),private :: BATopt

!public procedures:
public :: BAT_light, BAT_miscibility_gap, density_est, fnc_BAT_activity, fnc_BAT_act_dev

!-----------------------------------------------------------
contains
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine for the evaluation of the core Duhem-Margules Binary Activity           *
!*   Thermodynamics (BAT) model expressions as outlined by Gorkowski et al., (2019;     *
!*   https://doi.org/10.5194/acp-19-13383-2019).                                        *
!*   This procedure calculates the natural log of the mole-fraction based activity      *
!*   coefficients and the associated activities of water (component 1) and the organic  *
!*   component (2) of binary liquid mixture.                                            *
!*                                                                                      *
!*   Notes: equation and table numbers refer to those from the book chapter             *
!*          (section 3.1.2);                                                            *
!*          a single set of model parameters is used in this simplified "light"         *
!*          version of the BAT model.                                                   *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!****************************************************************************************
pure subroutine BAT_light(x_org, M_org, OtoC, rho_org, ln_actcoeff, activity)
    
implicit none

!interface arguments:
!x_org       [-], mole fraction of the organic component (input)  
!M_org       [g/mol], molar mass of the organic (input)   
!OtoC        [-], elemental oxygen-to-carbon ratio of organic (input)
!rho_org     [g/cm^3], liquid-state density of organic (input)
!ln_actcoeff [-], ln[activity coeff.] (output); array elements: 1 water, 2 org
!activity    [-], mole-fraction-based activity (output) 
real(wp),intent(in) :: x_org                         
real(wp),intent(in) :: M_org                  
real(wp),intent(in) :: OtoC                         
real(wp),intent(in) :: rho_org                      
real(wp),dimension(2),intent(out) :: ln_actcoeff    
real(wp),dimension(2),intent(out) :: activity       
    
!local constants:
![g/mol], molar mass of water
real(wp),parameter :: M_w = 18.015_wp
![g/cm^3], liquid-state density of water (here assumed 
! temperature-independent)
real(wp),parameter :: rho_w = 0.997_wp               
!fitted BAT model parameters [-] (as listed in Table 3.1):
real(wp),parameter :: s1 = 4.06991_wp
real(wp),parameter :: s2 = -1.23723_wp
!2-D array of a_n,i coefficients (column-major order)
real(wp),dimension(4,2),parameter :: apar = &                   
 & reshape([ 5.88511_wp, -4.73125_wp, -5.20165_wp, -30.8230_wp, &
 &          -0.98490_wp, -6.22721_wp, 2.32029_wp, -25.8404_wp], &
 &         shape = [4,2]) 
!local variables:
real(wp) :: dGEbyRT_dxorg, GEbyRT, M_ratio, one_minus_2phi, &
            & phi_param, phi_org, phi_by_x
real(wp),dimension(2) :: cpar   !array for c1, c2 parameters defined by Eq. (3.15)
!...............................................
    
!step 1) calculate recurring equation terms and model coefficients
M_ratio = M_w/M_org
phi_param = (rho_org/rho_w) * M_ratio * s1*(1.0_wp + OtoC)**s2
!Eq. (3.13):
phi_org = x_org / ( x_org + (1.0_wp - x_org)*phi_param )    
!Eq. (3.15):
cpar(:) = apar(1,:)*exp(apar(2,:)*OtoC) + &
            & apar(3,:)*exp(apar(4,:)*M_ratio)                 
    
!step 2) calculate normalized Gibbs excess energy term and its derivative
one_minus_2phi = 1.0_wp - 2.0_wp*phi_org
GEbyRT = phi_org*(1.0_wp - phi_org)*( cpar(1) + cpar(2)*one_minus_2phi )
    
! express (phi_org/x_org) in form avoiding division by zero:
phi_by_x = 1.0_wp/( x_org + (1.0_wp - x_org)*phi_param )  
dGEbyRT_dxorg = ( one_minus_2phi*(cpar(1) + cpar(2)*one_minus_2phi) &
    & -2.0_wp*cpar(2)*phi_org*(1.0_wp - phi_org) )* phi_param*phi_by_x**2
    
!step 3) calculate natural log of activity coefficients and activities
!Eqs. (3.11, 3.12):
ln_actcoeff(1) = GEbyRT - x_org*dGEbyRT_dxorg          
ln_actcoeff(2) = GEbyRT + (1.0_wp - x_org)*dGEbyRT_dxorg        
activity(1) = (1.0_wp - x_org)*exp(ln_actcoeff(1))  !water
activity(2) = x_org*exp(ln_actcoeff(2))             !organic
    
end subroutine BAT_light
!-----------------------------------------------------------
    
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Function to estimate the pure-component liquid-state density of organic compounds  *
!*   based on the simple method by Girolami (J. Chem. Educ., 71, 962–964, 1994).        *
!*   This implementation is suitable for applications with limited input information.   *
!*   Inputs: molar mass, O:C and H:C ratios. If the H:C ratio is unknown, enter a       *
!*   negative value. In that case the actual H:C will be estimated based on an initial  *
!*   assumption of H:C = 2 and futher corrected by oxygen content.                      *
!*                                                                                      *
!*   Notes: this implementation is for compounds containing only the elements O, C, H;  *
!*          the weak temperature dependence of density is ignored (a simplification);   *
!*          no correction applied for rings and/or aromatic compounds (limited info);   *
!*          the applied scaling implies the assumption that most of the oxygen atoms    *
!*          are able to make H-bonds (as donor or acceptor).                            *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
pure elemental function density_est(M, OtoC, HtoC) result(rho)
    
implicit none

!interface arguments:
real(wp),intent(in) :: M            ![g/mol] compound molar mass in 
real(wp),intent(in) :: OtoC         ![-] elemental oxygen-to-carbon ratio of molecule
real(wp),intent(in) :: HtoC         ![-] elemental hydrogen-to-carbon ratio of molecule
real(wp)            :: rho          ![g/cm^3] function return value
!local variables:
real(wp),parameter :: MC = 12.01_wp, MO = 16.0_wp, MH = 1.008_wp    !the molar masses of the carbon, oxygen and hydrogen atoms in [g/mol]
real(wp) :: n_carb, HtoC_est, rho_init
!...............................................

!estimate the HtoC value if not provided as input:
if (HtoC < 0.0_wp) then
    !estimate HtoC assuming an aliphatic compound with HtoC = 2 in absence of oxygen-bearing groups, 
    !then correct for oxygen content assuming a -1 slope (based on Van Krevelen Diagram of typical SOA).
    HtoC_est = 2.0_wp -OtoC
else
    HtoC_est = HtoC
endif

!determine approximate number of carbon atoms per organic molecule:
n_carb = M / (MC + HtoC_est*MH + OtoC*MO)

!compute density estimate based on method by Girolami (1994):
rho_init = M / (5.0_wp*n_carb*(2.0_wp + HtoC_est + OtoC*2.0_wp))
rho = rho_init*( 1.0_wp + min(n_carb*OtoC*0.1_wp, 0.3_wp) )     !density in [g/cm^3];   

end function density_est
!-----------------------------------------------------------
    
    
!-----------------------------------------------------------
!** wrapper function for BAT_light with interface as required by local_min function. **
pure function fnc_BAT_activity(x_org, fopt)
    
use Mod_NumPrec_Types, only : wp, foptions
implicit none
    
!interface arguments:
real(wp),intent(in) :: x_org
type(foptions),optional,intent(in) :: fopt
real(wp) :: fnc_BAT_activity
!local variables:
real(wp),dimension(2) :: ln_actcoeff
real(wp),dimension(2) :: activity
!....................................
    
call BAT_light(x_org, fopt%M_org, fopt%OtoC, fopt%density_org, ln_actcoeff, activity)
fnc_BAT_activity = activity(fopt%component_no)  !transfer water activity (component_no = 1) or organic activity (component_no = 2)
    
if (fopt%determine_local_max) then
    fnc_BAT_activity = -fnc_BAT_activity        !flip sign (to find a local maximum with a minimum search)
endif
    
end function fnc_BAT_activity
!-----------------------------------------------------------
    
    
!-----------------------------------------------------------
!** wrapper function for determining an activity deviation from a target value (using BAT_light). **
pure function fnc_BAT_act_dev(x_org, fopt)
    
use Mod_NumPrec_Types, only : wp, foptions
implicit none
    
!interface arguments:
real(wp),intent(in) :: x_org
type(foptions),optional,intent(in) :: fopt
real(wp) :: fnc_BAT_act_dev
!local variables:
real(wp),dimension(2) :: ln_actcoeff
real(wp),dimension(2) :: activity
!....................................
    
call BAT_light(x_org, fopt%M_org, fopt%OtoC, fopt%density_org, ln_actcoeff, activity)
    
!return relative deviation between computed activity and a given reference value:
fnc_BAT_act_dev = (activity(fopt%component_no) / fopt%refval) -1.0_wp   
    
end function fnc_BAT_act_dev
!-----------------------------------------------------------
    
    
!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Subroutine to compute the extent of a miscibility gap for a binary (!) system in   *
!*   terms of x_org using a nonlinear equation solver in two unknowns.                  *
!*                                                                                      *
!*   Approach: Use the p1, p2 points as initial guess for solving the phase separation  * 
!*   problem. Here we use an initial molar phase size ratio, omega, of 1.0 for finding  *
!*   the two x_org compositions denoting the binodal points of the LLE. This omega      *
!*   value is equivalent to picking an overall xtot_org input value in the middle of    *
!*   the phase separation range:  xtot_org = 0.5_wp*(xA_org + xB_org). A sigmoidal map  *
!*   is used to convert between q_alpha and apparent ln(affinity) values used by the    *
!*   solver.                                                                            *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend,                                                                        *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
subroutine BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, density_org, x_org_LLE_limits)
    
use Mod_MINPACK, only : hybrd1
    
implicit none
!interface arguments:
real(wp),dimension(2),intent(in) :: x_org_p1p2
real(wp),intent(in) :: M_org, OtoC, density_org 
real(wp),dimension(2),intent(out) :: x_org_LLE_limits
!local variables:
integer :: info
real(wp),parameter :: tol_LLE = max(1.0E-5_wp, sqrt(epsilon(1.0_wp)))
real(wp) :: sum_qx, xtot_org
real(wp),dimension(2) :: aff, lnaff, fvec, qA
real(wp),dimension(2,2) :: activities, ln_actcoeff
!...............................................

!omega = 1.0 initially (and for sigmoidal mapping)
xtot_org = 0.5_wp*(x_org_p1p2(1) + x_org_p1p2(2))   

!set BATopt values to be passed via private module variable 
!in Mod_BAT_model to the 'calc_BAT_LLE_dev' subroutine:
BATopt%xtot = [1.0_wp - xtot_org, xtot_org]     !x_w, x_org
BATopt%M_org = M_org
BATopt%OtoC = OtoC
BATopt%density_org = density_org

!evaluate BAT at initial guess points for phases alpha 'A' 
!and beta 'B'; 
!activities(:,1) = val. for water and org. of phase A (= 1):
write(*,'(A,2(ES13.6,1X),/)') 'initial guess for xA_org, xB_org = ', x_org_p1p2(2), x_org_p1p2(1)
call BAT_light(x_org_p1p2(1), M_org, OtoC, density_org, ln_actcoeff(:,1), activities(:,1))
call BAT_light(x_org_p1p2(2), M_org, OtoC, density_org, ln_actcoeff(:,2), activities(:,2))
ncalls_LLE = 2      !count BAT_light calls

!use vector lnaff(:) as unknown solver variables:
lnaff = ln_actcoeff(:,1) - ln_actcoeff(:,2)

!use Powell's hybrid method from module Mod_MINPACK to solve  
!the system of equations of the LLE isoactivity conditions 
!(implemented in subroutine calc_BAT_LLE_dev below), starting 
!with the guess for lnaff:
call hybrd1(calc_BAT_LLE_dev, 2, lnaff, fvec, tol_LLE, info)

!for information, report diagnostics and relative deviations
!from isoactivity condition:
write(*,'(A,I0)') 'number of BAT_light calls during LLE calc.: ', ncalls_LLE
write(*,'(A,I0)') 'after hybrd1 solver; info = ', info
write(*,'(A,2(ES13.6,1X))') 'after hybrd1 solver; fvec = ', fvec

if (info < 0) then 
    !terminated due to convergence to trivial solution
    !based on test within calc_BAT_LLE_dev
    qA = 0.5_wp
else
    aff(:) = exp(lnaff) 
    qA = aff/(1.0_wp + aff)
endif

!mole fractions in phases alpha and beta from determined qA:
sum_qx = sum(qA*BATopt%xtot)
x_org_LLE_limits(1) = qA(2)*BATopt%xtot(2)/sum_qx
sum_qx = sum((1.0_wp - qA)*BATopt%xtot)
x_org_LLE_limits(2) = (1.0_wp - qA(2))*BATopt%xtot(2)/sum_qx
!make x_org_LLE_limits(1) = xA_org the water-rich phase (by choice):
if (x_org_LLE_limits(1) > x_org_LLE_limits(2)) then     !switch
    sum_qx = x_org_LLE_limits(1)
    x_org_LLE_limits(1) = x_org_LLE_limits(2)
    x_org_LLE_limits(2) = sum_qx
endif
write(*,'(A,2(ES13.6,1X),/)') 'after hybrd1 solver; xA_org, xB_org = ', x_org_LLE_limits(1:2)

end subroutine BAT_miscibility_gap
!-----------------------------------------------------------

!-----------------------------------------------------------
!** subroutine for computing the deviation from isoactivity condition for LLE computation. **
subroutine calc_BAT_LLE_dev(n, lnaff, fvec, iflag)
    
use Mod_NumPrec_Types, only : wp, foptions
implicit none
    
!interface arguments:
integer,intent(in) :: n                         !input : number of variables & equations
real(wp),dimension(n),intent(in) :: lnaff       !input : set solver variables
real(wp),dimension(n),intent(out) :: fvec       !output: evaluated activity deviations
integer,intent(out) :: iflag                    !output: indicator flag to report issues
!local variables:
real(wp) :: sum_qx, xA_org, xB_org
real(wp),dimension(2) :: aff, qA, ln_actcoeff_A, ln_actcoeff_B
real(wp),dimension(2) :: activity_A, activity_B
!....................................
    
iflag = 0 
!sigmoidal map from lnaff solver vector to q_alpha vector:
aff = exp(lnaff)
qA = aff/(1.0_wp + aff)
    
!mole fractions in phases alpha and beta based on qA values:
!BATopt%xtot(2) is the overall mole frac. of component 2 (= org).
sum_qx = sum(qA*BATopt%xtot)
xA_org = qA(2)*BATopt%xtot(2)/sum_qx 
sum_qx = sum((1.0_wp - qA)*BATopt%xtot)
xB_org = (1.0_wp - qA(2))*BATopt%xtot(2)/sum_qx 
    
!compute activities in each phase:
!(additional BAT model component properties are provided via module variable BATopt)
call BAT_light(xA_org, BATopt%M_org, BATopt%OtoC, BATopt%density_org, ln_actcoeff_A, activity_A)
call BAT_light(xB_org, BATopt%M_org, BATopt%OtoC, BATopt%density_org, ln_actcoeff_B, activity_B)
ncalls_LLE = ncalls_LLE + 2 
    
!return vector of relative deviations between component activities in both phases:
fvec(1:2) = (activity_A - activity_B) / (0.5_wp*(activity_A + activity_B))

!check for signs of convergence to trivial solution:
if (all(qA < 0.5_wp) .OR. all(qA > 0.5_wp)) then
    if (sum(abs(fvec)) < 1.0E-2_wp) then
        iflag = -1      !signal to hybrd1 to stop and return;
    endif
endif

end subroutine calc_BAT_LLE_dev
!-----------------------------------------------------------
    
end module Mod_BAT_model