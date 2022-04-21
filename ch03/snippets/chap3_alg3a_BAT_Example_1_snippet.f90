!(1.5, optional): generate a quick, simple plot of the activities 
!       vs x_org using the DISLIN graphics library 
!       (https://www.dislin.de/index.html);
!       this requires an installed and linked DISLIN library.
block
   integer,parameter :: dp = kind(1.0D0)
   real(dp),dimension(:),allocatable :: xval, yval
   external :: metafl, qplot
   !................
   !apply floating point number conversion to double 
   !precision kind used here with Dislin;
   xval = real(x_org_vec, kind=dp)          
   yval = real(activities(:,1), kind=dp)        !water activity
   !simple x--y scatter plot:
   call metafl('xwin')                          !'xwin' or 'pdf'
   call qplot(xval, yval, np)
   yval = real(activities(:,2), kind=dp)        !organic activity
   call qplot(xval, yval, np)
end block


!(1.6, optional): Generate a nicer x--y curve plot, showing
!   several curves in the same plot (using Dislin but without use
!   of quick plots). 
!   Here we use a set of module procedures from Mod_Dislin_plots 
!   to make such plots with relatively few statements, while 
!   being able to set various curve and plot properties.
!   The following [block ... end block] code section requires an 
!   installed and linked DISLIN library (double precision version). 
!   Otherwise comment-out this block of code.
block
   use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot 
   character(len=75) :: xlabel, ylabel
   integer,dimension(3),parameter :: rgb_blue =  [40, 40, 255], &
                                 &   rgb_green = [0, 178, 0] 
   !....................................
   
   !first, add ideal mixing activity curve data for comparison:
   !all arguments, including optional ones, and their meaning 
   !are outlined near the top of subroutine 'add_plot_xydata' 
   !in module Mod_Dislin_plots.
   call add_plot_xydata(xv=x_org_vec, yv=(1.0_wp - x_org_vec), &
     & ltext='ideal, $a_w = x_w$', pen_wid=2.0_wp, &
     & lstyle='dotted', plot_symb='curve')
   call add_plot_xydata(xv=x_org_vec, yv=x_org_vec, &
     & ltext='ideal, $a_{\rm org} = x_{\rm org}$', &
     & pen_wid=2.0_wp, lstyle='dashed', plot_symb='curve')
   
   !second, add non-ideal water and organic activity curves:
   call add_plot_xydata(xv=x_org_vec, yv=activities(:,1), &
     & ltext='water activity, $a_w$', pen_wid=8.0_wp, &
     & rgb_col=rgb_blue, lstyle='solid', plot_symb='curve', &
     & symb_id=15)
   call add_plot_xydata(xv=x_org_vec, yv=activities(:,2), &
     & ltext='organic activity, $a_{\rm org}$', &
     & pen_wid=8.0_wp, rgb_col=rgb_green, lstyle='solid', &
     & plot_symb='curve') 
   
   !set overall plot properties and generate Dislin plot:
   xlabel = 'mole fraction of organic, $x_{\rm org}$'
   ylabel = 'activity'
   call dislin_plot(xlabel, ylabel, yaxis_mod=0.67_wp, &
     & legend_position=3, metafile='pdf', &
     & out_file_name='example1_activity_curves')
end block
