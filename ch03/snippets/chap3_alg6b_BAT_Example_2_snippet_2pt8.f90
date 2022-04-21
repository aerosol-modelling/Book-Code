!(2.8) Determine approximate x_org coordinates of the points p1 and p2.
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
        & point is: ', x_org_p1p2(k)   
    write(*,'(A,ES12.5)') 'relative deviation from targeted &
        & activity: ', fnc_BAT_act_dev(x_org_p1p2(k), fopt)
    write(*,'(A,I0,/)') 'number of BAT_light function calls &
        & during Ridders_zero = ', nfcalls
enddo !k

!(2.9)  Use the p1, p2 points as initial guess to 
!compute the composition limits of phase separation 
!(i.e. determine binodal points denoting the 
!miscibility gap).
call BAT_miscibility_gap(x_org_p1p2, M_org, OtoC, density_org, &
    & x_org_LLE_limits) 
    
    