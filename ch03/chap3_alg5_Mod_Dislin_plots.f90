!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module offering a simplified procedure for Dislin-based x--y scatter (curve) plots *
!*   using a few preset common plot settings. Several curves/data sets can be plotted   *
!*   on the same axis system with curve-specific properties provided.                   *
!*                                                                                      *
!*  This module requires an installed and linked DISLIN graphics library                *
!*  (https://www.dislin.de/index.html). See website for details about Dislin and        *
!*  installation instructions on different platforms.                                   *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   :: License ::                                                                      *
!*   This program is free software: you can redistribute it and/or modify it under the  *
!*   terms of the GNU General Public License as published by the Free Software          *
!*   Foundation, either version 3 of the License, or (at your option) any later         *
!*   version.                                                                           *
!*   The AIOMFAC model code is distributed in the hope that it will be useful, but      *
!*   WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or      *
!*   FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more      *
!*   details.                                                                           *
!*   You should have received a copy of the GNU General Public License along with this  *
!*   program. If not, see <http://www.gnu.org/licenses/>.                               *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  subroutine  add_plot_xydata                                                     *
!*   -  subroutine  dislin_plot                                                         *
!*                                                                                      *
!**************************************************************************************** 
module Mod_Dislin_plots

use Mod_NumPrec_Types, only : wp

implicit none

!public module variables and types:
type, private :: plot_xydata                                    !the (second) array dimension denotes the data set (curve) no.
    integer,dimension(:),allocatable        :: npoints          !stores number of data points of a data set to be plotted
    real(wp),dimension(:,:),allocatable     :: xval             !x coordinates of x--y data
    real(wp),dimension(:,:),allocatable     :: yval             !y coordinates of x--y data
    character(len=75),dimension(:),allocatable :: legtext       !legend entry text for a curve
    real(wp),dimension(:),allocatable :: pen_width              !pen width; typical values are between 1.0 and 9.0
    integer,dimension(:,:),allocatable :: rgb_color             !first dimension is red, green, blue values 
                                                                !(array of 3 values each in range [0, 255])
    character(len=13),dimension(:),allocatable :: line_style    !options: 'solid', 'dotted', 'dashed', 'dashed_medium'
    character(len=7),dimension(:),allocatable :: plot_symb      !plot curves w/ or w/o symbols or symbols only 
                                                                !(options: 'symbols', 'curve', 'both')
    integer,dimension(:),allocatable :: symbol_id               !ID (number) of symbol type to be plotted (see Dislin manual)
end type plot_xydata

type(plot_xydata),allocatable,private :: xy_data

!public procedures:
public :: add_plot_xydata, dislin_plot
!-----------------------------------------------------------

    contains
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Public subroutine to define the x--y data arrays and a selection of curve/symbol   *
    !*   properties for a specific data set (curve) to be added to a Dislin plot page.      *
    !*   The data will be stored in xy_data and later used within 'dislin_plot' to generate *
    !*   the plot. xy_data will be increased in the number of curves dimension with each    *
    !*   call of this subroutine (until a reset at the end of the subroutine dislin_plot)   *
    !*                                                                                      *
    !*   Note: Dislin offers many additional plot types and controls of plot properties.    *
    !*         Here we use a typically useful set of options for x--y scatter/curve plots.  *
    !*                                                                                      *
    !*   :: Authors ::                                                                      *
    !*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2021-06-29                                                      *
    !*   -> latest changes: 2022-03-23                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine add_plot_xydata(xv, yv, ltext, pen_wid, rgb_col, lstyle, plot_symb, symb_id)
    
    implicit none
    !interface arguments:
    real(wp),dimension(:),intent(in) :: xv, yv                  !arrays of x and y values of data set to be added
    character(len=*),intent(in) :: ltext                        !legend entry text for data set
    real(wp),optional,intent(in) :: pen_wid                     !pen width; typical values are between 1.0 and 9.0
    integer,dimension(3),optional,intent(in) :: rgb_col         !color array of the 3 r,g,b values, each in range [0, 255])
    character(len=*),optional,intent(in) :: lstyle              !line style; options: 'solid', 'dotted', 'dashed', 'dashed_medium'
    character(len=*),optional,intent(in) :: plot_symb           !plot curves w/ or w/o symbols or symbols only 
                                                                !(options: 'symbols', 'curve', 'both')
    integer,optional,intent(in) :: symb_id                      !the symbol ID for dislin (15 = open circle, 5 = open diamond, 
                                                                !3 = +, 4 = X, 16 = filled square, 21 = filled circle, etc.)
    
    !local variables:
    integer :: k, ndset, np, np_old, istat
    logical :: is_allocated
    type(plot_xydata),allocatable :: temp_data
    !...................................
    
    k = size(xv)
    
    !check whether a xy_data variable has been allocated already: 
    if ( allocated(xy_data) ) then  !xy_data%xval
        is_allocated = .true.
        ndset = size(xy_data%xval, dim=2)
        np_old = size(xy_data%xval, dim=1)
        !update ndset and np (with room for data set to be added to xy_data):
        ndset = ndset +1
        np = max(k, np_old)
    else
        is_allocated = .false.
        ndset = 1
        np = k
    endif
    
    allocate( temp_data, stat=istat)
    allocate( temp_data%npoints(ndset), temp_data%xval(np,ndset), temp_data%yval(np,ndset),  &
        & temp_data%legtext(ndset), temp_data%pen_width(ndset), temp_data%rgb_color(3,ndset), &
        & temp_data%line_style(ndset), temp_data%plot_symb(ndset), temp_data%symbol_id(ndset), stat=istat )
    
    if (is_allocated) then  
        !copy the data for the first ndset-1 entries to the temporary array
        temp_data%npoints(1:ndset-1)       = xy_data%npoints(1:ndset-1)
        temp_data%xval(1:np_old,1:ndset-1) = xy_data%xval(1:np_old,1:ndset-1)
        temp_data%yval(1:np_old,1:ndset-1) = xy_data%yval(1:np_old,1:ndset-1)
        temp_data%legtext(1:ndset-1)       = xy_data%legtext(1:ndset-1)
        temp_data%pen_width(1:ndset-1)     = xy_data%pen_width(1:ndset-1)
        temp_data%rgb_color(:,1:ndset-1)   = xy_data%rgb_color(:,1:ndset-1)
        temp_data%line_style(1:ndset-1)    = xy_data%line_style(1:ndset-1)
        temp_data%plot_symb(1:ndset-1)     = xy_data%plot_symb(1:ndset-1) 
        temp_data%symbol_id(1:ndset-1)     = xy_data%symbol_id(1:ndset-1) 
    endif
    !move the data into enlarged xy_data using memory-friendly move_alloc:
    call move_alloc(temp_data, xy_data)     !temp_data will be deallocated
    
    !finally add the new data set to xy_data, including accouting for optional arguments:
    xy_data%npoints(ndset) = k
    xy_data%xval(1:k,ndset) = xv
    xy_data%yval(1:k,ndset) = yv
    xy_data%legtext(ndset) = ''                                 !initialize whole string
    xy_data%legtext(ndset)(1:min(75, len_trim(ltext))) = ltext(1:min(75, len_trim(ltext)))
    if (present(pen_wid)) then
        xy_data%pen_width(ndset) = pen_wid
    else    
        xy_data%pen_width(ndset) = 1.0_wp                       !set a default value
    endif
    if (present(rgb_col)) then
        xy_data%rgb_color(1:3,ndset) = rgb_col(1:3)
    else
        xy_data%rgb_color(1:3,ndset) = [0, 0, 0]                !default (black)
    endif
    if (present(lstyle)) then
        xy_data%line_style(ndset) = ''
        xy_data%line_style(ndset)(1:min(75, len_trim(lstyle))) = lstyle(1:min(75, len_trim(lstyle)))
    else
        xy_data%line_style(ndset) = 'solid'
    endif
    if (present(plot_symb)) then
        xy_data%plot_symb(ndset) = ''
        xy_data%plot_symb(ndset)(1:min(75, len_trim(plot_symb))) = plot_symb(1:min(75, len_trim(plot_symb)))
    else
        xy_data%plot_symb(ndset) = 'curve'
    endif
    if (present(symb_id)) then
        xy_data%symbol_id(ndset) = symb_id
    else
        xy_data%symbol_id(ndset) = 15
    endif
    
    end subroutine add_plot_xydata
    !-----------------------------------------------------------
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Public subroutine to generate a x--y scatter plot with the previously loaded data  *
    !*   (via calls to subroutine 'add_plot_xydata'                                         *
    !*                                                                                      *
    !*   Note: Dislin offers many additional plot types and controls of plot properties.    *
    !*         Here we use a typically useful set of options for x--y scatter plots.        *
    !*                                                                                      *
    !*   :: Authors ::                                                                      *
    !*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        2021-06-29                                                      *
    !*   -> latest changes: 2022-03-23                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine dislin_plot(xlabel, ylabel, yaxis_mod, xaxis_limits, yaxis_limits, &
        & legend_position, metafile, out_file_name)

    implicit none
    !interface arguments:
    character(len=*),intent(in) :: xlabel, ylabel       !text for axis labels;
    real(wp),intent(in) :: yaxis_mod                    !set scaling of default length of y-axis on graph 
                                                        !(e.g. 0.5 indicates a plot height of half the x-axis length);
    real(wp),dimension(2),optional,intent(in) :: &
        & xaxis_limits, yaxis_limits                    !(optional) targeted x-axis and y-axis limits for plot;
    integer,intent(in) :: legend_position               !legend position value options:
                                                        != 0 for no legend
                                                        != 1 is the lower left corner of the page.
                                                        != 2 is the lower right corner of the page.
                                                        != 3 is the upper right corner of the page.
                                                        != 4 is the upper left corner of the page.
                                                        != 5 is the lower left corner of the axis system.
                                                        != 6 is the lower right corner of the axis system.
                                                        != 7 is the upper right corner of the axis system.
                                                        != 8 is the upper left corner of the axis system.
    character(len=*),intent(in) :: metafile             !typically set as 'xwin', 'cons' or 'pdf';
    character(len=*),intent(in) :: out_file_name        !relative path/file name of output file (e.g. for a .pdf file, 
                                                        !without stating file extension);
    
    !local variables:
    character(len=:),allocatable :: cbuff
    integer,parameter :: dp = kind(1.0D0)               !use this real kind with the double precision version of Dislin.
    integer :: i, istat, ndsets, npts, nxl, nyl, nzl
    real(dp) :: axis_min, axis_max, axis_range, xa, xe, xor, xstep, ya, ye, yor, ystep
    ! the external Dislin procedures (available via linked static library):
    external :: metafl, scrmod, disini, pagfll, setclr, psfont, chaspc, height, hwfont, &
        & texmod, name, setscl, graf, incmrk, legini, leglin, legend, lncap, dot, dash, &
        & dashm, solid, penwid, color, curve, setrgb, endgrf, disfin, setfil, filmod, &
        & legtit, frame, linesp, getlen, axslen, psmode, marker, nochek   
    !...............................................
    
    !The object xy_data contains stored data for a curve or for several curves, 
    !including curve color, thickness, legend text entry, etc.
    
    ndsets = size(xy_data%pen_width)            !number of data sets
    allocate(character(len=75*min(ndsets, 30)) :: cbuff)

    call metafl(metafile)       
    call setfil(trim(out_file_name)//'.'//trim(metafile))
    call filmod('delete')                       !overwrite if file already exists                
    call scrmod('norev')
    call disini
    call pagfll(255)                            !set page background color to white
    if (trim(metafile) == 'pdf') then
        call psfont('Helvetica')                !for 'pdf' use a postscript font; e.g. 'Times-Roman' or 'Helvetica'
        call psmode('both')                     !allow both Greek and italic modes
    else
        call hwfont()                           !use a hardware font; choice depends on operating system
    endif
    call chaspc(-0.06_dp)                       !slightly adjust character spacing
    call texmod('on')                           !allow for TeX-style statments in figure text 
    call color('black')                         !set text/axis/curve color by name
    
    !potentially modify axis system properties
    call getlen(nxl, nyl, nzl)
    nxl = ceiling(0.8_dp*nxl)
    call axslen(nxl, ceiling(nxl*yaxis_mod))    !modify the aspect ratio of the axis system via y-axis scaling
    
    !set axis system properties
    call name(trim(xlabel), 'X')                !set x-axis label text
    call name(trim(ylabel), 'Y')
    
    !set automatic scaling for x-axis and y-axis based on slighly scaled input ranges:
    axis_min = minval( [( minval(xy_data%xval(1:xy_data%npoints(i), i)), i=1,ndsets )] ) 
    axis_max = maxval( [( maxval(xy_data%xval(1:xy_data%npoints(i), i)), i=1,ndsets )] ) 
    if (present(xaxis_limits)) then
        axis_min = xaxis_limits(1)  
        axis_max = xaxis_limits(2)   
    endif
    axis_range = axis_max - axis_min
    axis_min = axis_min -0.008_dp*axis_range
    axis_max = axis_max +0.008_dp*axis_range
    call setscl([axis_min, axis_max], 2, 'X')
    
    axis_min = minval( [( minval( xy_data%yval(1:xy_data%npoints(i), i)), i=1,ndsets )] ) 
    axis_max = maxval( [( maxval( xy_data%yval(1:xy_data%npoints(i), i)), i=1,ndsets )] )     
    if (present(yaxis_limits)) then
        axis_min = yaxis_limits(1)  
        axis_max = yaxis_limits(2)   
    endif
    axis_range = axis_max - axis_min
    axis_min = axis_min -(0.015_dp/yaxis_mod)*axis_range
    axis_max = axis_max +(0.015_dp/yaxis_mod)*axis_range
    call setscl([axis_min, axis_max], 2, 'Y')   !set automatic scaling for y-axis based on targeted range
    
    !initialize graph axis system (using automatic scaling):
    call nochek()                               !suppress warning about points outside of the plotting area
    call graf(xa, xe, xor, xstep, ya, ye, yor, ystep)  
    call legini(cbuff, min(ndsets, 30), 75)     !initialize legend
    
    !plot x--y data for all curves/symbols with the curve-specific properties:
    do i = 1,ndsets
        npts = xy_data%npoints(i)
        !select plotting only symbols or curves or both:
        select case( trim(xy_data%plot_symb(i)) )
        case('curve')
            call incmrk(0)
        case('symbols')
            call incmrk(-1)    
            call marker(xy_data%symbol_id(i)) 
        case('both')
            call incmrk(1)        
            call marker(xy_data%symbol_id(i))                       
        case default
            call incmrk(0)              
        end select
        !select line style:
        select case( trim(xy_data%line_style(i)) )
        case('solid')
            call lncap('long')
            call solid()
        case('dotted')
            call lncap('round')     !rounded line caps
            call dot()                             
        case('dashed')  
            call lncap('round')
            call dash()                          
        case('dashed_medium')
            call lncap('cut')
            call dashm()                              
        case default
            call incmrk(0)              
        end select
        call penwid(real(xy_data%pen_width(i), kind=dp))    !set pen / curve width (especially for pdf output) 
        call setrgb(xy_data%rgb_color(1,i)/255.0_dp, xy_data%rgb_color(2,i)/255.0_dp, &
            & xy_data%rgb_color(3,i)/255.0_dp)              !set color by RGB value
        call leglin(cbuff, xy_data%legtext(i), i)           !call leglin(cbuff, trim(xy_data%legtext(i)), i)
        call curve(real(xy_data%xval(1:npts,i), kind=dp), real(xy_data%yval(1:npts,i), kind=dp), npts)
    enddo
    
    !plot legend and finalize plot:
    if (legend_position /= 0) then
        call height(26)                             !set font size (height) for legend
        call color('black')
        call penwid(1.0_dp)  
        call legtit('')                             
        call frame(1)                               !set legend frame thickness
        call linesp(2.0_dp)                         !modify line spacing of legend entries
        call legend(cbuff, legend_position)         !plot legend; e.g. 3 = position at upper right corner of page
    endif
    call endgrf
    call disfin
    deallocate(cbuff)
    
    !deallocate xy_data content so the next plot is not including current data:
    deallocate( xy_data, stat=istat )    

    end subroutine dislin_plot
    !------------------------------------------------------------ 
    
end module Mod_Dislin_plots