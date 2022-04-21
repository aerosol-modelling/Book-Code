!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Example 1 demonstrating running the BAT model for set organic properties.          *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
subroutine BAT_example_1()

use Mod_NumPrec_Types, only : wp
use Mod_BAT_model, only : BAT_light, density_est

implicit none
!local variables: 
integer :: i, np, unt
real(wp) :: density_org, HtoC, M_org, OtoC, xinc
real(wp),dimension(:),allocatable :: x_org_vec
real(wp),dimension(:,:),allocatable :: activities, ln_actcoeff
!...............................................

!## Example 1 ## 

!We will run the BAT calculation for a series of organic mole 
!fractions, x_org, covering the composition range from pure water 
!to pure organic.

!(1.1): set input parameters for the organic component, pinic  
!       acid, C9H14O4, SMILES: CC1(C)C(C(O)=O)CC1CC(O)=O, in  
!       the units required by the BAT_light model 
!       (see the interface of BAT_light); the liquid-state density will be estimated 
!       will be estimated using a simple, limited-information method based on that  
!       by Girolami (1994).
M_org       = 186.207_wp                        ![g/mol], molar mass
OtoC        = 4.0_wp/9.0_wp                     ![-] O:C ratio
HtoC        = 14.0_wp/9.0_wp
density_org = density_est(M_org, OtoC, HtoC)    ![g/cm^3] 

!(1.2): allocate and initialize input and output arrays to cover
!       composition range; 
!       e.g., for a mole fraction increment of 0.01, we need 101 points
!       for the range from 0.0 to 1.0. The activities and ln_actcoeff 
!       arrays are 2-dimensional, with the first dimension denoting the
!       data point index while the second dimension saves the two entries
!       for water (1) and organic (2) aligned with an x_org value.
np = 101
allocate( x_org_vec(np), ln_actcoeff(np,2), activities(np,2) )
xinc = 1.0_wp/real(np-1, kind=wp)           !the x_org increment
!population of array values via implied loop:
x_org_vec = [(i*xinc, i = 0,np-1)]          

!(1.3): use BAT_light to compute activities, looping over x_org:
do i = 1,np
    call BAT_light(x_org_vec(i), M_org, OtoC, density_org, ln_actcoeff(i,:), activities(i,:))
enddo

!(1.4): write the data to a text file (comma-separated values):
open(newunit = unt, file = './example1_BAT_output.csv', status = 'unknown')
write(unt,'(A)') 'x_org,  act_coeff_water,  act_coeff_org,  a_w,  a_org,'
do i = 1,np
    write(unt,'(4(ES12.5,","),ES12.5)') x_org_vec(i), exp(ln_actcoeff(i,:)), activities(i,:) 
enddo
close(unt)

write(*,'(A)') 'Output saved to file example1_BAT_output.csv'
write(*,*)

!!(1.5, optional): generate a quick, simple plot of the activities 
!!       vs x_org using the DISLIN graphics library 
!!       (https://www.dislin.de/index.html);
!!       this requires an installed and linked DISLIN library.
!block
!   integer,parameter :: dp = kind(1.0D0)
!   real(dp),dimension(:),allocatable :: xval, yval
!   external :: metafl, qplot
!   !................
!   xval = real(x_org_vec, kind=dp)          !for correct floting point number conversion to double precision kind used here with Dislin;
!   yval = real(activities(:,1), kind=dp)    !water activity
!   !simple x--y scatter plot:
!   call metafl('xwin')                      !'xwin' or 'pdf'
!   call qplot(xval, yval, np)
!   yval = real(activities(:,2), kind=dp)    !organic activity
!   call qplot(xval, yval, np)
!end block


!!(1.6, optional): Generate a nicer x--y curve plot, showing
!!   several curves in the same plot (using Dislin but without use
!!   of quick plots). 
!!   Here we use a set of module procedures from Mod_Dislin_plots 
!!   to make such plots with relatively few statements, while 
!!   being able to set various curve and plot properties.
!!   The following [block ... end block] code section requires an 
!!   installed and linked DISLIN library (double precision version). 
!!   Otherwise comment-out this block of code.
!block
!   use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
!   character(len=75) :: xlabel, ylabel
!   integer,dimension(3),parameter :: rgb_blue =  [40, 40, 255], &
!                                 &   rgb_green = [0, 178, 0] 
!   !....................................
!   
!   !first, add ideal mixing activity curve data for comparison:
!   !all arguments, including optional ones, and their meaning 
!   !are outlined near the top of subroutine 'add_plot_xydata' 
!   !in module Mod_Dislin_plots.
!   call add_plot_xydata(xv=x_org_vec, yv=(1.0_wp - x_org_vec), &
!     & ltext='ideal, $a_w = x_w$', pen_wid=2.0_wp, &
!     & lstyle='dotted', plot_symb='curve')
!   call add_plot_xydata(xv=x_org_vec, yv=x_org_vec, &
!     & ltext='ideal, $a_{\rm org} = x_{\rm org}$', &
!     & pen_wid=2.0_wp, lstyle='dashed', plot_symb='curve')
!   
!   !second, add non-ideal water and organic activity curves:
!   call add_plot_xydata(xv=x_org_vec, yv=activities(:,1), &
!     & ltext='water activity, $a_w$', pen_wid=8.0_wp, &
!     & rgb_col=rgb_blue, lstyle='solid', plot_symb='curve', &
!     & symb_id=15)
!   call add_plot_xydata(xv=x_org_vec, yv=activities(:,2), &
!     & ltext='organic activity, $a_{\rm org}$', &
!     & pen_wid=8.0_wp, rgb_col=rgb_green, lstyle='solid', &
!    & plot_symb='curve') 
!   
!   !set overall plot properties and generate Dislin plot:
!   xlabel = 'mole fraction of organic, $x_{\rm org}$'
!   ylabel = 'activity'
!   call dislin_plot(xlabel, ylabel, yaxis_mod=0.67_wp, &
!     & legend_position=3, metafile='pdf', &
!     & out_file_name='example1_activity_curves')
!end block

end subroutine BAT_example_1