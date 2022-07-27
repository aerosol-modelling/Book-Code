!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Chapter 3, exercise 2 solution to determining the O:C ratio of LLPS onset.         *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
subroutine BAT_exercise_2()

use Mod_NumPrec_Types, only : wp, foptions
use Mod_BAT_model, only : BAT_light, BAT_miscibility_gap, density_est, &
    & fnc_BAT_activity, fnc_BAT_act_dev, fnc_LLPS_det
use Mod_NumMethods, only : local_min, Ridders_zero

implicit none
!local variables: 
integer :: i, k, np, nfcalls, iloc1, iloc2
integer,dimension(2) :: imin, imax
logical :: phase_sep_detected
real(wp) :: density_org, M_org, OtoC, OtoC_crit, xinc
real(wp) :: OtoC_low, OtoC_high, fnc_low, fnc_high, fnc_crit, tol
real(wp) :: ax, bx, fnc_loc1, fnc_loc2, flocmin, tol_init
real(wp),dimension(2) :: activity_locmin, x_org_min, x_org_LLE_limits, &
    & x_org_p1p2
real(wp),dimension(:),allocatable :: x_org_vec
real(wp),dimension(:,:),allocatable :: activities, ln_actcoeff
type(foptions) :: fopt, fopt_LLPS
!...............................................

!## Chapter 3, exercise 2 ##

! For a given binary mixture (here  water + 1,6-hexanediol), 
! determine the presence/absence of a liquid--liquid phase separation.

! set general BAT function options for this organic compound 
! via elements of custom derived type fopt_LLPS:
M_org = 118.18_wp                       ![g/mol], molar mass
fopt_LLPS%M_org = M_org
fopt_LLPS%HtoC = 14.0_wp/6.0_wp
! set number of data points to be evaluated during each LLPS 
! detection function call:
np = 15
fopt_LLPS%np = np

! set and compute function values at the limits of the O:C interval 
! in which we aim to find OtoC_crit; we need to provide bounds that 
! bracket a zero crossing of fnc_LLPS_det:
OtoC_low = 0.0_wp
OtoC_high = 0.9_wp          !set to a reasonably high value such that no LLPS occurs
fnc_low = fnc_LLPS_det(OtoC_low, fopt_LLPS) 
fnc_high = fnc_LLPS_det(OtoC_high, fopt_LLPS)
tol = 1.0E-5_wp             !or use max( 1.0E-5_wp, sqrt(epsilon(1.0_wp)) )

OtoC_crit = Ridders_zero(fnc_LLPS_det, fopt_LLPS, OtoC_low, OtoC_high, fnc_low, &
            & fnc_high, tol, nfcalls)

! the computed OtoC_crit value should be within numerical tolerance of the onset of LLPS;
! check whether it is refers to an LLPS case or an O:C just outside that range:
fnc_crit = fnc_LLPS_det(OtoC_crit, fopt_LLPS)
! determine critical O:C value for which LLPS just occurs:
do while (fnc_crit < 0.0_wp)            !should take 0 or only 1 iteration
    OtoC_crit = OtoC_crit -tol
    fnc_crit = fnc_LLPS_det(OtoC_crit, fopt_LLPS)
enddo

write(*,'(A,I0,/)') '## Exercise 2 solution -- critical O:C for LLPS ##'
write(*,*) ''
write(*,'(A,I0)') 'Number of points for activity curve resolution: ', np
write(*,'(A,ES12.5)') 'The critical O:C ratio for LLPS onset is: ', OtoC_crit

! set general properties of the organic component with the determined limiting case O:C ratio:
OtoC = OtoC_crit
density_org = density_est(M_org, OtoC, fopt_LLPS%HtoC)    ![g/cm^3]

! with the determined O:C ratio, call the BAT_light model again to compute 
! activities, looping over x_org points:
xinc = 1.0_wp/real(np-1, kind=wp)               !the x_org increment
allocate( x_org_vec(np), ln_actcoeff(np,2), activities(np,2) )
x_org_vec = [(i*xinc, i = 0,np-1)]              !population of array values via implied loop
do i = 1,np
    call BAT_light(x_org_vec(i), M_org, OtoC, density_org, & 			
        & ln_actcoeff(i,:), activities(i,:))
enddo

!(2.4, optional): generate a plot of the activities vs x_org using the DISLIN  
!       graphics library (https://www.dislin.de/index.html);
!       this requires an installed and linked DISLIN library.
!       Otherwise comment-out the following [block ... end block] code section.
block
    use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
    character(len=75) :: xlabel, ylabel
    integer,dimension(3),parameter :: rgb_blue = [40, 40, 255], rgb_green = [0, 178, 0] 
    !....................................
    
    !first, add ideal mixing activity data for comparison:
    call add_plot_xydata(xv=x_org_vec, yv=(1.0_wp - x_org_vec), ltext='ideal, $a_w = x_w$', &
            & pen_wid=2.0_wp, lstyle='dotted', plot_symb='curve')
    call add_plot_xydata(xv=x_org_vec, yv=x_org_vec, ltext='ideal, $a_{\rm org} = x_{\rm org}$', &
            & pen_wid=2.0_wp, lstyle='dashed', plot_symb='curve')
    
    !second, add non-ideal water and organic activity curves:
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,1), ltext='water activity, $a_w$', &
            & pen_wid=4.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='both', symb_id=15)
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,2), ltext='organic activity, $a_{\rm org}$', &
            & pen_wid=4.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='both', symb_id=4) 
    
    !set overall plot properties and generate plot:
    xlabel = 'mole fraction of organic, $x_{\rm org}$'
    ylabel = 'activity'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.75_wp, legend_position=3, &
            & metafile='pdf', out_file_name='exerc2_activity_curves')
end block

!(2.7) find a local minimum and a local maximum of the activity 
!   curves -- if bracketing intervals are present.
!   first determine approximate index locations of extrema 
!   for both water and organic activity curves:
phase_sep_detected = .false.
imin = 0
imax = 0
do k = 1,2   !loop over the two mixture components
    do i = 2,np-1
        if ( activities(i,k) < activities(i-1,k) ) then
            if ( activities(i,k) < activities(i+1,k) ) then 
                !found a local minimum
                imin(k) = i
                phase_sep_detected = .true.
            endif
        else
            if ( activities(i,k) > activities(i+1,k) ) then 
                !found a local maximum
                imax(k) = i
                phase_sep_detected = .true.
            endif
        endif
        if (imin(k) > 0 .AND. imax(k) > 0) then
            exit
        endif
    enddo
enddo !k

!report on whether LLPS was detected for the given input:
if (phase_sep_detected) then
    write(*,'(A)') 'LLPS was detected. Composition limits of LLPS will need to be determined.'
else
    write(*,'(A)') 'LLPS was not detected. A single liquid phase is likely stable over the whole composition range.'
endif

!Note: if any local extremum is found for any curve, then both 
!a local max and min must be present.  Due to the functional 
!relation between aw and aorg, we know that the x-coordinate of 
!a local min of one curve is the same as the local max x_org 
!for the other curve and analogous for a local maximum.  
!Because limited curve resolution was used, make sure that all 
!indices are correctly set when local extrema are present:
if (phase_sep_detected) then
    !possibly adjust based on local min/max from other curve:
    if (imin(1) /= imax(2)) then
        imin(1) = max(imin(1), imax(2))
        imax(2) = imin(1)
    endif
    if (imin(2) /= imax(1)) then
        imin(2) = max(imin(2), imax(1))
        imax(1) = imin(2)
    endif
    !also consider a local min/max being near the activity = 
    !1.0 end point (one index off) due to low resolution of 
    !curve (missing extrema) and knowing that x_org_vec values 
    !increase with increasing index number.
    if (imin(1) == 0 .AND. imax(2) == 0) then  
        imin(1) = 2
        imax(2) = 2
    endif
    if (imin(2) == 0 .AND. imax(1) == 0) then  
        imin(2) = np-1
        imax(1) = np-1
    endif
endif

if (phase_sep_detected) then
    tol_init = max( 1.0E-4_wp, sqrt(epsilon(1.0_wp)) )
    write(*,'(A,ES12.5)') 'set tolerance for local_min: ', tol_init
    !set general BAT wrapper function options via elements 
    !of custom derived type fopt:
    fopt%M_org = fopt_LLPS%M_org 
    fopt%OtoC = OtoC 
    fopt%density_org = density_org
    do k = 1,2
        !determine local minimum within interval [ax, bx]:
        ax = x_org_vec(imin(k) -1)
        bx = x_org_vec(imin(k) +1)
        !set fopt type components for finding a local minimum 
        !in water (k=1) or organic (k=2) activity:
        fopt%determine_local_max = .false.
        fopt%component_no = k
        !use local_min method:
        flocmin = local_min(ax, bx, tol_init, fnc_BAT_activity, &
            & fopt, x_org_min(k), nfcalls)
        !save activity value at determined local min:
        activity_locmin(k) = flocmin 
        write(*,'(A,I0,A,ES13.6,A,ES13.6)') 'found local min &
            & of activity curve at x_org_min(',k,') = ', &
            & x_org_min(k), ', flocmin = ', flocmin
        write(*,'(A,I0,/)') 'number of BAT_light calls during &
            & local min search = ', nfcalls
    enddo  !k        

    
    !(2.8) Determine approximate x_org coordinates of the points p1 and p2.
    !      We are interested in finding the approximate intersection of the activity = activity_locmin 
    !      value with the components' activity curve located on the "other" side of the local maximum 
    !      of the curve. This means we can constrain the search range to the curve data by using the 
    !      determined x_org array index values corresponding to the local activity min and max values.
    do k = 1,2
        !determine approx. array index value (iloc) of the p1 
        !or p2 point:
        if (imin(k) <= imax(k)) then
            iloc1 = imax(k)-1 &
                & + minloc( abs(activities(imax(k):np,k) &
                & -activity_locmin(k)), DIM = 1 ) 
        else
            iloc1 = minloc( abs(activities(1:imax(k),k) &
                & -activity_locmin(k)), DIM = 1 ) 
        endif
    
        !determine curve data interval indices [iloc1, iloc2] 
        !to bracket the point:
        if (iloc1 > 1 .AND. iloc1 < np) then
            if (activities(iloc1,k) > activity_locmin(k)) then
                if (activities(iloc1+1, k) < activity_locmin(k)) then
                    iloc2 = iloc1 +1
                else 
                    iloc2 = iloc1 -1
                endif
            else
                if (activities(iloc1+1, k) > activity_locmin(k)) then
                    iloc2 = iloc1 +1
                else 
                    iloc2 = iloc1 -1
                endif
            endif 
        else 
            !iloc1 is first or last index, so iloc2 must be  
            !the only neighbor
            if (iloc1 == 1) then
                iloc2 = 2
            else
                iloc2 = np-1
            endif
        endif
        write(*,'(A,I0,A,ES12.5,/)') 'approx x_org of p',k,' &
            & point is: ', 0.5_wp*( x_org_vec(iloc1) + &
            & x_org_vec(iloc2) )
        
        !use Ridders' method to determine the x_org value of the
        !p1 or p2 point (within set tolerance):
        tol_init = max( 1.0E-3_wp, sqrt(epsilon(1.0_wp)) )
        fopt%component_no = k
        fopt%refval = activity_locmin(k)
        fnc_loc1 = (activities(iloc1,k) / activity_locmin(k)) -1.0_wp
        fnc_loc2 = (activities(iloc2,k) / activity_locmin(k)) -1.0_wp
        
        x_org_p1p2(k) = Ridders_zero(fnc_BAT_act_dev, fopt, &
            & x_org_vec(iloc1), x_org_vec(iloc2), fnc_loc1, &
            & fnc_loc2, tol_init, nfcalls)
        
        !for debugging output to screen:
        write(*,'(A,I0,A,ES12.5)') 'determined x_org of p',k,' &
            & point is: ', x_org_p1p2(k)   
        write(*,'(A,ES12.5)') 'relative deviation from targeted &
            & activity: ', fnc_BAT_act_dev(x_org_p1p2(k), fopt)
        write(*,'(A,I0,/)') 'number of BAT_light calls during &
            & Ridders_zero = ', nfcalls
    enddo !k
    
    !(2.9)  Use the p1, p2 points as initial guess to 
    !compute the composition limits of phase separation 
    !(i.e. determine binodal points denoting the 
    !miscibility gap).
    call BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, density_org, x_org_LLE_limits)   
else
    write(*,'(A)') 'no LLPS detected'
    x_org_LLE_limits = -1.0_wp
endif

end subroutine BAT_exercise_2  
!----------------------------------------------------------- 