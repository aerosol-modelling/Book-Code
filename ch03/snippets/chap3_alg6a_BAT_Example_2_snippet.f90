!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Example 2 demonstrating finding local minima/maxima of BAT activity curves and     * 
!*   determining potential phase separation and the related x_org composition range.    *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
subroutine BAT_example_2()

use Mod_NumPrec_Types, only : wp, foptions
use Mod_BAT_model, only : BAT_light, BAT_miscibility_gap, density_est, &
    & fnc_BAT_activity, fnc_BAT_act_dev
use Mod_NumMethods, only : local_min, Ridders_zero

implicit none
!local variables: 
integer :: i, k, np, nfcalls, iloc1, iloc2
integer,dimension(2) :: imin, imax
logical :: phase_sep_detected
real(wp) :: density_org, HtoC, M_org, OtoC, xinc
real(wp) :: ax, bx, fnc_loc1, fnc_loc2, flocmin, tol_init
real(wp),dimension(2) :: activity_locmin, x_org_min, x_org_LLE_limits, &
    & x_org_p1p2
real(wp),dimension(:),allocatable :: x_org_vec, delta_Gmix_by_nRT, &
    & delta_Gmix_ideal_by_nRT
real(wp),dimension(:,:),allocatable :: activities, ln_actcoeff
type(foptions) :: fopt
!...............................................

!## Example 2 -- consideration of LLPS ##

!For a given binary mixture (here  water + 1-pentanol), determine the presence/absence
!of a liquid--liquid phase separation.
!First we compute a reasonable number of points of the activity curves, then we
!explore whether local extrema are present on those curves; if present, we determine the
!approximate x_org coordinates of the local min/max.

!(2.1): set input parameters for the organic component, 
!       1-pentanol, C5H12O, SMILES: CCCCCO
!       in the units required by the BAT model (see the interface declaration of BAT_light);
!       the liquid-state density will be estimated using a simple, limited-information method 
!       based on that by Girolami (1994).
M_org       = 88.15_wp                       ![g/mol], molar mass
OtoC        = 1.0_wp/5.0_wp                  ![-], O:C ratio
HtoC        = 12.0_wp/5.0_wp
density_org = density_est(M_org, OtoC, HtoC) ![g/cm^3] 

!(2.2): allocate and initialize input and output arrays to cover composition range with a few points, as in example 1.
np = 1001                                    !no. of points
xinc = 1.0_wp/real(np-1, kind=wp)            !the x_org increment
allocate( x_org_vec(np), ln_actcoeff(np,2), activities(np,2) )
x_org_vec = [(i*xinc, i = 0,np-1)]              

!(2.3): call the BAT model to compute activities, looping over x_org points:
do i = 1,np
    call BAT_light(x_org_vec(i), M_org, OtoC, density_org, &
		& ln_actcoeff(i,:), activities(i,:))
enddo

write(*,'(A,I0,/)') '## Example 2 -- consideration of LLPS ##'
write(*,'(A,I0,/)') 'number of BAT_light calls in step (2.3): ', np

!(2.4, optional): generate a plot of the activities vs x_org using the DISLIN  
!       graphics library (https://www.dislin.de/index.html);
!       this requires an installed and linked DISLIN library.
!       Otherwise comment-out the following [block ... end block] code section.
block
    use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
    character(len=75) :: xlabel, ylabel
    integer,dimension(3),parameter :: rgb_blue =  [40, 40, 255], rgb_green = [0, 178, 0] 
    !....................................
    
    !first, add ideal mixing activity data for comparison:
    call add_plot_xydata(xv=x_org_vec, yv=(1.0_wp - x_org_vec), ltext='ideal, $a_w = x_w$', &
            & pen_wid=2.0_wp, lstyle='dotted', plot_symb='curve')
    call add_plot_xydata(xv=x_org_vec, yv=x_org_vec, ltext='ideal, $a_{\rm org} = x_{\rm org}$', &
            & pen_wid=2.0_wp, lstyle='dashed', plot_symb='curve')
    
    !second, add non-ideal water and organic activity curves:
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,1), ltext='water activity, $a_w$', &
            & pen_wid=8.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve', symb_id=15)
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,2), ltext='organic activity, $a_{\rm org}$', &
            & pen_wid=8.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='curve') 
    
    !set overall plot properties and generate plot:
    xlabel = 'mole fraction of organic, $x_{\rm org}$'
    ylabel = 'activity'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.75_wp, legend_position=3, &
            & metafile='pdf', out_file_name='example2_activity_curves')
end block

!(2.5)  Compute normalized molar Gibbs energy of mixing from 
!		the computed activity data.
!       Since we are interested in the non-ideal mixing effects 
!		and the phase separation region, all we need to know about
!		is the change in the Gibbs energy due to mixing, 
!       normalized by RT and total molar amount nt.
allocate( delta_Gmix_by_nRT(np), delta_Gmix_ideal_by_nRT(np) )
where (x_org_vec(:) > 0.0_wp .AND. (1.0_wp - x_org_vec(:)) &
& > 0.0_wp)
    delta_Gmix_by_nRT = &
        & (1.0_wp - x_org_vec)*log(activities(:,1)) &
        & + x_org_vec*log(activities(:,2))
    !ideal mixing case for comparison:
    delta_Gmix_ideal_by_nRT = &
        & (1.0_wp - x_org_vec)*log((1.0_wp - x_org_vec)) &
        & + x_org_vec*log(x_org_vec)
elsewhere
    delta_Gmix_by_nRT = 0.0_wp
    delta_Gmix_ideal_by_nRT = 0.0_wp
endwhere

!(2.6)  plot normalized molar Gibbs energy of mixing:
block
    use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
    character(len=75) :: xlabel, ylabel
    !....................................
    
    !add curves/data sets to be plotted
    call add_plot_xydata(xv=x_org_vec, yv=delta_Gmix_ideal_by_nRT, ltext='ideal mixing', pen_wid=4.0_wp, rgb_col=[50, 50, 50], lstyle='dotted', plot_symb='curve')
    call add_plot_xydata(xv=x_org_vec, yv=delta_Gmix_by_nRT, ltext='non-ideal mixing', pen_wid=8.0_wp) 
            !use default values for the other (optional) data set properties;
    
    !set overall plot properties and generate plot:
    xlabel = 'mole fraction of organic, $x_{\rm org}$'
    ylabel = '$\Delta G_{\rm mix} / (n R T)$'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, legend_position=3, &
            & metafile='pdf', out_file_name='example2_Delta_Gmix')
end block


!(2.7) find a local minimum and a local maximum of the activity 
!   curves -- if bracketing intervals are present.
!   first determine approximatactive index location of extrema 
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
    fopt%M_org = M_org 
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
        write(*,'(A,I0,/)') 'number of BAT_light function calls &
            & during local min search = ', nfcalls
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
		&point is: ', 0.5_wp*( x_org_vec(iloc1) + &
		& x_org_vec(iloc2) )
	
	!use Ridders' method to determine the x_org value of the
	!p1 or p2 points (within set tolerance):
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
		&point is: ', x_org_p1p2(k)   
	write(*,'(A,ES12.5)') 'relative deviation from targeted &
		& activity: ', fnc_BAT_act_dev(x_org_p1p2(k), fopt)
	write(*,'(A,I0,/)') 'number of BAT_light function calls &
		&during Ridders_zero = ', nfcalls
enddo !k
    
!!x_org_p1p2(1:2) = [0.1, 0.9]  !for test calculations only

!(2.9)  Use the p1, p2 points to generate an initial guess for solving the phase separation problem.
call BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, density_org, x_org_LLE_limits) 


!(2.10) Generate activity curve data with consideration of LLE, when within the miscibility gap.

!example to be added...