block !for chapter 3, exercise 4 solution;
	integer :: iH
	!(f) Plot equil-state phase pH vs. water activity (for
	!    each phase in LLPS case), where applicable:
	iH = findloc(int(out_data_A(7,1,:)), VALUE=205, DIM=1)
	!iH = index location of H+ ion
	if (iH > 0) then
		call add_plot_xydata(xv=out_data_A(6,:,1), &
			& yv=-log10( out_data_A(6,:,iH) ), &
			& ltext='pH of phase $\alpha$', &
			& pen_wid=8.0_wp, lstyle='solid', plot_symb='curve')

		if (LLPSpresent) then
			call add_plot_xydata(xv=out_data_B(6,iloc1:iloc2,1), &
				& yv=-log10( out_data_B(6,iloc1:iloc2,iH) ), &
				& ltext='pH of phase $\beta$', &
				& pen_wid=8.0_wp, lstyle='dashed_medium', &
				& plot_symb='curve')
		endif

		!set overall plot properties and generate plot:
		xlabel = 'water activity, $a_w$'
		ylabel = 'pH'
		yax_lim = [-1.0_wp, 3.0_wp]
		call dislin_plot(xlabel, ylabel, yaxis_mod=0.6_wp, &
			& yaxis_limits=yax_lim, legend_position=3, &
			& metafile='pdf', out_file_name=trim(fname)//'_pH')
	endif
end block