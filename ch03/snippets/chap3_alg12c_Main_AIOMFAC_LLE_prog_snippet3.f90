!==== AIOMFAC-LLE calculation section ==========================
!
if (errorind == 0) then !perform AIOMFAC-LLE calculations
                        !; else jump to termination section;
allocate(inputconc(nindcomp), outnames(NKNpNGS), &
& out_data_A(8,npoints,NKNpNGS), out_data_B(8,npoints,NKNpNGS), &
& out_LLEprop(7,npoints), LLEprop(8), LLEoutvars(8,NKNpNGS,2), &
& name_species_TeX(NKNpNGS), STAT=allocstat)
inputconc = 0.0_wp
out_data_A = 0.0_wp
out_data_B = 0.0_wp
!--
if (verbose) then
    write(unito,*) ""
    write(unito,'(A)') "MESSAGE from AIOMFAC: mixture defined, &
    & calculating composition points... "
    write(unito,*) ""
endif

!set AIOMFAC input and call the main AIOMFAC-LLE subroutine for 
!all composition points;
!note: xinputtype is set .true. for mole fraction input (on 
!component-basis, i.e. electrolytes undissociated) or .false. to 
!indicate mass fraction input (setting from input file);
ignore_LLPS = .false.   !(set .true. to force a single-phase 
                        !calculation)

do pointi = 1,npoints   !loop over mixture points, changing 
                        !composition and/or temperature       
    inputconc(1:ncp) = compos2(pointi,1:ncp)
    !--
    call AIOMFAC_LLE_inout(inputconc, xinputtype, T_K(pointi), &
    & ignore_LLPS, nspecies, LLEprop, LLEoutvars, outnames, &
    & errorflag, warningflag)
    !--
    if (warningflag > 0 .OR. errorflag > 0) then
        call RepErrorWarning(unito, errorflagmix, warningflag, &
        & errorflag, pointi, errorind, warningind)
    endif
    !--
    !save properties of this input point:
    out_LLEprop(:,pointi) = LLEprop
    !"nspecies" denotes the maximum number of different species 
    !in the mixture (accounting for ions and for the possibility
    !of HSO4- dissociation after input)
    !out_data_A array structure: | LLEoutvars data columns <1:7>
    ! | data point no. | species no.|
    do nc = 1,nspecies
        out_data_A(1:7,pointi,nc) = LLEoutvars(1:7,nc,1)   
        out_data_B(1:7,pointi,nc) = LLEoutvars(1:7,nc,2)       
        out_data_A(8,pointi,nc) = real(errorflag, kind=wp)      
        out_data_B(8,pointi,nc) = real(errorflag, kind=wp)
        if (errorflag == 0 .AND. warningflag > 0) then  
            !do not overwrite an errorflag if present
            if (warningflag == 16) then                     
                !a warning that only affects viscosity calc. 
                !(here not used)
            else
                out_data_A(8,pointi,nc) = real(warningflag, &
                & kind=wp)
                out_data_B(8,pointi,nc) = real(warningflag, &
                & kind=wp)
            endif
        endif
    enddo !nc
    !--
enddo !pointi