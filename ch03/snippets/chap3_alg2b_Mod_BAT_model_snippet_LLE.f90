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
subroutine BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, &
							& density_org, x_org_LLE_limits)
    
use Mod_MINPACK, only : hybrd1
    
implicit none
!interface arguments:
real(wp),dimension(2),intent(in) :: x_org_p1p2
real(wp),intent(in) :: M_org, OtoC, density_org 
real(wp),dimension(2),intent(out) :: x_org_LLE_limits
!local variables:
integer :: info
real(wp),parameter :: tol_LLE = max(1.0E-5_wp, &
	& sqrt(epsilon(1.0_wp)))
real(wp) :: sum_qx, xtot_org
real(wp),dimension(2) :: aff, lnaff, fvec, qA
real(wp),dimension(2,2) :: activities, ln_actcoeff
!...............................................

!omega = 1.0 initially (and fixed for sigmoidal mapping)
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
call BAT_light(x_org_p1p2(1), M_org, OtoC, density_org, &
    & ln_actcoeff(:,1), activities(:,1))
call BAT_light(x_org_p1p2(2), M_org, OtoC, density_org, &
    & ln_actcoeff(:,2), activities(:,2))
ncalls_LLE = 2      !count BAT_light calls

!initialize lnaff(:) as unknown (solver variable):
lnaff = ln_actcoeff(1,:) - ln_actcoeff(2,:)

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
!make x_org_LLE_limits(1) (= xA_org) the water-rich phase 
!(by choice):
if (x_org_LLE_limits(1) > x_org_LLE_limits(2)) then
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
integer,intent(in) :: n                     	!input : number of variables & equations
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