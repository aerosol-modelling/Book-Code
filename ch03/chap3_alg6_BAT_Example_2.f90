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
    & x_org_p1p2, act_LLE_A, act_LLE_B, ln_ac_LLE_A, ln_ac_LLE_B
real(wp),dimension(:),allocatable :: x_org_vec, delta_Gmix_by_nRT, &
    & delta_Gmix_ideal_by_nRT, fracA, qA_org, qA_w, x_org_A, x_org_B
real(wp),dimension(:,:),allocatable :: activities, activities_B, &
    & ln_actcoeff, ln_actcoeff_B
type(foptions) :: fopt
!...............................................

!## Example 2 -- consideration of LLPS ##

!For a given binary mixture (here  water + 1-pentanol), determine the presence/absence
!of a liquid--liquid phase separation.
!First we compute a reasonable number of points of the activity curves, then we
!explore whether local extrema are present on those curves; if present, we determine the
!approximate x_org coordinates of the local min/max.

!(2.1): set input parameters for the organic component, 1-pentanol, C5H12O, SMILES: CCCCCO,
!       in the units required by the BAT model (see the interface declaration of BAT_light);
!       the liquid-state density will be estimated using a simple, limited-information method 
!       based on that by Girolami (1994).
M_org       = 88.15_wp                          ![g/mol], molar mass
OtoC        = 1.0_wp/5.0_wp                     ![-], O:C ratio
HtoC        = 12.0_wp/5.0_wp
density_org = density_est(M_org, OtoC, HtoC)    ![g/cm^3] 

!(2.2): allocate and initialize input and output arrays to cover composition range
!       with a few points, as in example 1.
np = 15 !1001                                   !set the number of BAT evaluation points (> 5) for x_org within [0.0, 1.0]
xinc = 1.0_wp/real(np-1, kind=wp)               !the x_org increment
allocate( x_org_vec(np), ln_actcoeff(np,2), activities(np,2) )
x_org_vec = [(i*xinc, i = 0,np-1)]              !population of array values via implied loop

!(2.3): call the BAT_light model to compute activities, looping over x_org points:
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
    integer,dimension(3),parameter :: rgb_blue = [40, 40, 255], rgb_green = [0, 178, 0] 
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
            & pen_wid=8.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='curve', symb_id=4) 
    
    !set overall plot properties and generate plot:
    xlabel = 'mole fraction of organic, $x_{\rm org}$'
    ylabel = 'activity'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.75_wp, legend_position=3, &
            & metafile='pdf', out_file_name='example2_activity_curves')
end block

!(2.5)  Compute normalized molar Gibbs energy of mixing from the computed activity data.
!       Since we are interested in the non-ideal mixing effects and the phase separation region,
!       all we need to know about is the change in the Gibbs energy due to mixing, 
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
    call add_plot_xydata(xv=x_org_vec, yv=delta_Gmix_ideal_by_nRT, ltext='ideal mixing', pen_wid=4.0_wp, &
            & rgb_col=[50, 50, 50], lstyle='dotted', plot_symb='curve')
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
    
    !!x_org_p1p2(1:2) = [0.1, 0.9]  !set initial guess for test calculations only
    
    !(2.9)  Use the p1, p2 points as initial guess to 
    !compute the composition limits of phase separation 
    !(i.e. determine binodal points denoting the 
    !miscibility gap).
    call BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, density_org, x_org_LLE_limits)   
else
    write(*,'(A)') 'no LLPS detected'
    x_org_LLE_limits = -1.0_wp
endif


!(2.10) Generate equilibrium-state phase fractions and activity 
!       data with consideration of potential presence of LLPS.
deallocate ( x_org_vec, ln_actcoeff, activities, delta_Gmix_by_nRT )

np = 1001
allocate( x_org_vec(np), ln_actcoeff(np,2), activities(np,2), &
    & ln_actcoeff_B(np,2), activities_B(np,2), &
    & delta_Gmix_by_nRT(np), fracA(np), qA_org(np), qA_w(np), &
    & x_org_A(np), x_org_B(np) )

xinc = 1.0_wp/real(np-1, kind=wp)
x_org_vec = [(i*xinc, i = 0, np-1)]     !the input (total) x_org

if (phase_sep_detected) then
    !compute activity coeff. and activities for LLE compositions;
    !phase alpha:
    call BAT_light(x_org_LLE_limits(1), M_org, OtoC, &
        & density_org, ln_ac_LLE_A, act_LLE_A)  
    !phase beta:
    call BAT_light(x_org_LLE_limits(2), M_org, OtoC, &
        & density_org, ln_ac_LLE_B, act_LLE_B)
endif

!For equilibrium computations, call the BAT_light model to 
!compute activities *only* when outside the miscibility gap:
do i = 1,np
    if ( x_org_vec(i) < x_org_LLE_limits(1) .OR. &
        & x_org_vec(i) > x_org_LLE_limits(2) ) then
        ! --> outside miscibility gap:
        call BAT_light(x_org_vec(i), M_org, OtoC, &
            & density_org, ln_actcoeff(i,:), activities(i,:))
        if (x_org_vec(i) < x_org_LLE_limits(1)) then
            fracA(i) = 1.0_wp
            qA_org(i) = 1.0_wp
            qA_w(i) = 1.0_wp
        else
            fracA(i) = 0.0_wp
            qA_org(i) = 0.0_wp
            qA_w(i) = 0.0_wp
        endif
        x_org_A(i) = x_org_vec(i)
        x_org_B(i) = x_org_vec(i)
        ln_actcoeff_B(i,:) = ln_actcoeff(i,:)
        activities_B(i,:) = activities(i,:)
    else
        ! --> initial input x_org is within miscibility gap:
        x_org_A(i) = x_org_LLE_limits(1)
        x_org_B(i) = x_org_LLE_limits(2)
        activities(i,:) = act_LLE_A
        activities_B(i,:) = act_LLE_B
        ln_actcoeff(i,:) = ln_ac_LLE_A
        ln_actcoeff_B(i,:) = ln_ac_LLE_B
        
        !use lever rule to determine total fraction in phase A:
        fracA(i) = (x_org_vec(i) - x_org_B(i)) / &
            & (x_org_A(i) - x_org_B(i))
        
        !calculate fractions qA_org and qA_w:
        qA_org(i) = fracA(i) * x_org_A(i) / x_org_vec(i)
        qA_w(i)   = fracA(i) * (1.0_wp - x_org_A(i)) &
            & / (1.0_wp - x_org_vec(i))
    endif
enddo

!compute equilibrium-state (normalized) Gibbs energy of mixing:
where (activities(:,1) > 0.0_wp .AND. activities(:,2) > 0.0_wp)
    !the generally valid version based on Eq. (3.22):
    delta_Gmix_by_nRT = fracA*( (1.0_wp - x_org_A) * &
      & log(activities(:,1)) + x_org_A*log(activities(:,2)) ) &
      & + (1.0_wp - fracA)*( (1.0_wp - x_org_B) * &
      & log(activities_B(:,1)) + x_org_B * &
      & log(activities_B(:,2)) )
elsewhere
    delta_Gmix_by_nRT = 0.0_wp
endwhere


!(2.11) (optional) generate several plots of equilibrium calculation outputs:
block
    use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
    character(len=75) :: xlabel, ylabel
    integer,dimension(3),parameter :: rgb_blue = [40, 40, 255], &
        & rgb_green = [0, 178, 0], rgb_purple = [112, 51, 173]
    real(wp),dimension(2) :: DGmixbyRT_binodal
    !....................................
 
    !(a) Plot of equilibrium-state delta Gibbs energy of mixing (normalized)
    call add_plot_xydata(xv=x_org_vec, yv=delta_Gmix_by_nRT, &
        & ltext='non-ideal mixing, LLE considered', pen_wid=9.0_wp)
    if (phase_sep_detected) then
        !add open diamond symbols denoting binodal points:
        DGmixbyRT_binodal(1:2) = (1.0_wp - x_org_LLE_limits(1:2))*log(act_LLE_A(1)) &
            & + x_org_LLE_limits(1:2)*log(act_LLE_A(2))
        call add_plot_xydata(xv=x_org_LLE_limits, yv=DGmixbyRT_binodal, ltext='LLPS limits', &
            & pen_wid=3.0_wp, plot_symb='symbols', symb_id=5)
    endif
    !set overall plot properties and generate plot:
    xlabel = 'input mole fraction of organic, $x_{\rm org}^t$'
    ylabel = 'equil. $\Delta G_{\rm mix} / (n R T)$'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.3_wp, legend_position=3, &
        & metafile='pdf', out_file_name='ex2_equil_Delta_Gmix')
    !------------------------------------------------------------
    
    !(b) Plot of equilibrium-state fraction of cumulative molar amounts in phase alpha:
    !   plot is only useful in LLPS case;
    if (phase_sep_detected) then
        call add_plot_xydata(xv=x_org_vec, yv=fracA, ltext='$f^{\;\alpha}$', &
            & pen_wid=9.0_wp, rgb_col=rgb_purple, lstyle='solid', plot_symb='curve')
        !add open diamond symbols denoting binodal points:
        call add_plot_xydata(xv=x_org_LLE_limits, yv=[1.0_wp, 0.0_wp], ltext='LLPS limits', &
            & pen_wid=3.0_wp, rgb_col=rgb_purple, plot_symb='symbols', symb_id=5)
        !set overall plot properties and generate plot:
        xlabel = 'input mole fraction of organic, $x_{\rm org}^t$'
        ylabel = '$f^{\;\alpha}$'   !'fractional amount in $\alpha$'
        call dislin_plot(xlabel, ylabel, yaxis_mod=0.25_wp, legend_position=3, &
            & metafile='pdf', out_file_name='ex2_equil_fracA')
    endif
    !------------------------------------------------------------
    
    !(c) Plot of equilibrium-state activity curves:
    !first, add ideal mixing activity data for comparison:
    call add_plot_xydata(xv=x_org_vec, yv=(1.0_wp - x_org_vec), ltext='ideal, $a_w = x_w$', &
            & pen_wid=2.0_wp, lstyle='dotted', plot_symb='curve')
    call add_plot_xydata(xv=x_org_vec, yv=x_org_vec, ltext='ideal, $a_{\rm org} = x_{\rm org}$', &
            & pen_wid=2.0_wp, lstyle='dashed', plot_symb='curve')
    !water:
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,1), &
        & ltext='water activity, $a_w$, LLE considered', &
        & pen_wid=9.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve', symb_id=15)
    if (phase_sep_detected) then
        !add open diamond symbols denoting binodal points:
        call add_plot_xydata(xv=x_org_LLE_limits, yv=[act_LLE_A(1), act_LLE_B(1)], &
            & ltext='LLPS limits', pen_wid=3.0_wp, rgb_col=rgb_blue, plot_symb='symbols', symb_id=5) 
    endif
    !organic:
    call add_plot_xydata(xv=x_org_vec, yv=activities(:,2), &
        & ltext='organic activity, $a_{\rm org}$, LLE considered', &
        & pen_wid=8.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='curve', symb_id=4) 
    if (phase_sep_detected) then
        !add open diamond symbols denoting binodal points:
        call add_plot_xydata(xv=x_org_LLE_limits, yv=[act_LLE_A(2), act_LLE_B(2)], &
            & ltext='LLPS limits', pen_wid=3.0_wp, rgb_col=rgb_green, plot_symb='symbols', symb_id=5) 
    endif
    !set overall plot properties and generate plot:
    xlabel = 'input mole fraction of organic, $x_{\rm org}^t$'
    ylabel = 'equil. activities'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.3_wp, legend_position=3, &
        & metafile='pdf', out_file_name='ex2_equil_activity_curves')
    !------------------------------------------------------------
    
    !(d) Plot fractions of org. or water in phase alpha (qA_j):
    !   plot is only useful in LLPS case;
    if (phase_sep_detected) then
        !water (qA_w):
        call add_plot_xydata(xv=x_org_vec, yv=qA_w, ltext='$q^\alpha_{\rm w}$', &
            & pen_wid=9.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve')
        !organic (qA_org):
        call add_plot_xydata(xv=x_org_vec, yv=qA_org, ltext='$q^\alpha_{\rm org}$', &
            & pen_wid=7.0_wp, rgb_col=rgb_green, lstyle='dashed_medium', plot_symb='curve')
        !add open diamond symbols denoting binodal points:
        call add_plot_xydata(xv=x_org_LLE_limits, yv=[1.0_wp, 0.0_wp], ltext='LLPS limits', &
            & pen_wid=3.0_wp, plot_symb='symbols', symb_id=5)
        !set overall plot properties and generate plot:
        xlabel = 'input mole fraction of organic, $x_{\rm org}^t$'
        ylabel = 'fraction of $j$ in $\alpha$, $q^\alpha_{\rm j}$'
        call dislin_plot(xlabel, ylabel, yaxis_mod=0.25_wp, legend_position=3, &
            & metafile='pdf', out_file_name='ex2_equil_qA')
    endif
    !------------------------------------------------------------ 
        
    if (phase_sep_detected) then
        !determine the array index values of points close to and 
        !just within the LLPS limits:
        iloc1 = minloc(abs(x_org_vec(:) - x_org_LLE_limits(1)), &
            & mask = x_org_vec(:) > x_org_LLE_limits(1), dim=1)
        iloc2 = minloc(abs(x_org_vec(:) - x_org_LLE_limits(2)), &
            & mask = x_org_vec(:) < x_org_LLE_limits(2), dim=1)
    endif
    
    !(e) Plot equilibrium-state total mole fraction of water vs. water activity:
    if (phase_sep_detected) then
        call add_plot_xydata(xv=activities(1:iloc1,1), yv=(1.0_wp - x_org_vec(1:iloc1)), &
            & ltext='equil.-state $x_w$', pen_wid=9.0_wp, rgb_col=rgb_blue, &
            & lstyle='solid', plot_symb='curve')
        
        call add_plot_xydata(xv=activities(iloc2:,1), yv=(1.0_wp - x_org_vec(iloc2:)), &
            & ltext='equil.-state $x_w$', pen_wid=9.0_wp, rgb_col=rgb_blue, &
            & lstyle='solid', plot_symb='curve')
        
        call add_plot_xydata(xv=activities(iloc1:iloc2,1), yv=(1.0_wp - x_org_vec(iloc1:iloc2)), &
            & ltext='equil.-state $x_w$ within LLPS', pen_wid=5.0_wp, rgb_col=rgb_blue, &
            & lstyle='dotted', plot_symb='curve')
        
        !add open diamond symbols denoting binodal points:
        call add_plot_xydata(xv=[act_LLE_A(1), act_LLE_B(1)], yv=(1.0_wp - x_org_LLE_limits), &
            & ltext='LLPS limits', pen_wid=3.0_wp, rgb_col=rgb_blue, plot_symb='symbols', symb_id=5) 
    else    !a single liquid phase only
        call add_plot_xydata(xv=activities(:,1), yv=(1.0_wp - x_org_vec), &
            & ltext='equil.-state $x_w$', pen_wid=9.0_wp, rgb_col=rgb_blue, &
            & lstyle='solid', plot_symb='curve')
    endif
    !set overall plot properties and generate plot:
    xlabel = 'water activity'
    ylabel = 'equil. $x_w^t$'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.6_wp, legend_position=3, &
        & metafile='pdf', out_file_name='ex2_equil_xwtot')
    !------------------------------------------------------------ 
    
    !(f) Plot of equilibrium-state ln(activity coefficient) curves:
    if (phase_sep_detected) then
        !water:
        call add_plot_xydata(xv=x_org_vec(1:iloc2), yv=ln_actcoeff(1:iloc2,1), & 
            & ltext='water act. coeff., $\ln{(\gamma_w^\alpha)}$', &
            & pen_wid=6.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve')
        call add_plot_xydata(xv=x_org_vec(iloc1:), yv=ln_actcoeff_B(iloc1:,1), &
            & ltext='water act. coeff., $\ln{(\gamma_w^\beta)}$', &
            & pen_wid=4.0_wp, rgb_col=rgb_blue, lstyle='dashed_medium', plot_symb='curve')
        !organic:
        call add_plot_xydata(xv=x_org_vec(1:iloc2), yv=ln_actcoeff(1:iloc2,2), &
            & ltext='org. act. coeff., $\ln{(\gamma_{\rm org}^\alpha)}$', &
            & pen_wid=6.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='curve')
        call add_plot_xydata(xv=x_org_vec(iloc1:), yv=ln_actcoeff_B(iloc1:,2), &
            & ltext='org. act. coeff., $\ln{(\gamma_{\rm org}^\beta)}$', &
            & pen_wid=4.0_wp, rgb_col=rgb_green, lstyle='dashed_medium', plot_symb='curve')
    else
        !water:
        call add_plot_xydata(xv=x_org_vec, yv=ln_actcoeff(:,1), & 
            & ltext='water act. coeff., $\ln{(\gamma_w^\alpha)}$', &
            & pen_wid=6.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve')
        !organic:
        call add_plot_xydata(xv=x_org_vec, yv=ln_actcoeff(:,2), &
            & ltext='org. act. coeff., $\ln{(\gamma_{\rm org}^\alpha)}$', &
            & pen_wid=6.0_wp, rgb_col=rgb_green, lstyle='solid', plot_symb='curve')
    endif
    
    !set overall plot properties and generate plot:
    xlabel = 'input mole fraction of organic, $x_{\rm org}^t$'
    ylabel = 'equil. ln(activity coefficient)'
    call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, legend_position=3, &
        & metafile='pdf', out_file_name='ex2_equil_lnActCoeff')
    !------------------------------------------------------------
end block

!## end of example 2 ##

end subroutine BAT_example_2