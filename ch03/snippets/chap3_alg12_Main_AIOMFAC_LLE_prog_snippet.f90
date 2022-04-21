!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main program section to interface the AIOMFAC model with input parameters received * 
!*   from a webpage or command line.  The program gets the name of an input data file   *
!*   via a command line argument and calls AIOMFAC-LLE using these input parameters.    *
!*   Output is produced in the form of an output txt-file and/or a related .html file   *
!*   for webpage display.                                                               *
!*                                                                                      *
!*   Relative folder structure expected for input/output; directories 'Inputfiles' and  *
!*   'Outputfiles' need to be present at run time with read & write permissions:        *
!*   [location of executable *.out] -> [./Inputfiles/input_0???.txt]                    *                       
!*   [location of executable *.out] -> [./Outputfiles/output_0???.txt]                  *
!*                                                                                      *
!*   The AIOMFAC model expressions and parameters are described in Zuend et al. (2008,  * 
!*   Atmos. Chem. Phys.) and Zuend et al. (2011, Atmos. Chem. Phys.). Interaction       *
!*   parameters of Zuend et al. (2011) are used where they differ from the previous     *
!*   version. Additional parameters from Zuend and Seinfeld (2012), e.g. for peroxides, *
!*   are included as well. The LLE method is an improved and extended version of the    *
!*   method by Zuend and Seinfeld (2013).                                               *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
!*   Div. Chemistry and Chemical Engineering, Caltech, Pasadena, CA, USA (2009 - 2012)  *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
!*                                                                                      *
!*   -> created:        2011  (original web version)                                    *
!*   -> latest changes: 2021/08/23                                                      *
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
!****************************************************************************************
    
program Main_AIOMFAC_LLE_prog

use Mod_NumPrec, only : wp
use ModSystemProp, only : errorflagmix, nindcomp, NKNpNGS, SetSystem, topsubno, waterpresent
use ModSubgroupProp, only : SubgroupAtoms, SubgroupNames
use ModMRpart, only : MRdata
use ModSRunifac, only : SRdata
use Mod_InputOutput, only : OutputLLE_TXT, OutputLLE_plots, ReadInputFile, RepErrorWarning

implicit none
!set preliminary input-related parameters:
integer,parameter :: maxpoints = 101     
integer,parameter :: ninpmax = 51        !set the maximum number of mixture components allowed (arbitrary parameter)
!local variables:
character(len=4) :: VersionNo
character(len=200) :: filename
character(len=3000) :: filepath, folderpathout, fname, txtfilein  
character(len=200),dimension(:),allocatable :: cpnameinp
character(len=200),dimension(:),allocatable :: outnames, &
    & name_species_TeX
integer :: allocstat, errorflag, errorind, i, nc, ncp, npoints, &
    & nspecies, nspecmax, pointi, unito, warningflag, &
    & warningind, watercompno
integer,dimension(:,:),allocatable :: cpsubg   !list of input component subgroups and corresponding subgroup quantities
real(wp),dimension(:),allocatable :: T_K
real(wp),dimension(:),allocatable :: inputconc
real(wp),dimension(:,:),allocatable :: composition, compos2, &
    & out_LLEprop
real(wp),dimension(:,:,:),allocatable :: out_data_A, out_data_B
real(wp),dimension(:),allocatable :: LLEprop         !selected properties of predicted LLE; structure ( #-phases, omega, phi, Gibbs-E_diff, LstarA, LstarB, reladiffbest )
real(wp),dimension(:,:,:),allocatable :: LLEoutvars  !3-D array
    !with computed compositions and activities for each species
    !and phase; structure is:  (| input_mass-frac, mass-frac, 
    !mole-frac, molality, act.coeff., activity, ion-indicator |
    !species-no |phase-no |)
logical :: ignore_LLPS, filevalid, verbose, xinputtype
!..............................................................

!
!==== initialization section ===================================
!
VersionNo = "3.0" 
verbose = .true.    !if true, some debugging information will be 
                    !printed to the unit "unito" (errorlog file)
nspecmax = 0
errorind = 0
warningind = 0 
!
!==== input data section =======================================
!
!read command line for text-file name (which contains the input 
!parameters to run the AIOMFAC program):
call get_command_argument(1, txtfilein) 
if (len_trim(txtfilein) < 4) then
    !no command line argument stated; use specific input file 
    !for tests:
    txtfilein = './Inputfiles/input_0404.txt'
endif
filepath = adjustl(trim(txtfilein))
write(*,*) ""
write(*,'(A,A)') "MESSAGE from AIOMFAC-LLE: program started, command line argument 1 = ", trim(filepath)
write(*,*) ""
allocate(cpsubg(ninpmax,topsubno), cpnameinp(ninpmax), &
    & composition(maxpoints,ninpmax), T_K(maxpoints), &
    & STAT=allocstat)
!--
call ReadInputFile(filepath, folderpathout, filename, ninpmax, &
    & maxpoints, unito, verbose, ncp, npoints, warningind, &
    & errorind, filevalid, cpnameinp, cpsubg, T_K, composition, &
    & xinputtype)
!--
if (filevalid) then
    !
    !==== AIOMFAC model parameters and system initialization ======
    !
    if (verbose) then
        write(unito,*) ""
        write(unito,'(A)') "MESSAGE from AIOMFAC: input file &
        & read, starting AIOMFAC mixture definitions and &
        & initialization... "
        write(unito,*) ""
    endif
    
    !load the MR and SR interaction parameter data:
    call MRdata()
    call SRdata()
    call SubgroupNames()
    call SubgroupAtoms()

    !set system properties based on the data from the input file:
    call SetSystem(1, .true., ncp, cpnameinp(1:ncp), &
        & cpsubg(1:ncp,1:topsubno) )

    !check whether water is present in the mixture and as which 
    !component number:
    watercompno = 0
    if (waterpresent) then
        watercompno = findloc(cpsubg(1:ncp,16), value=1, dim=1)
    endif
    !transfer composition data to adequate array size:
    allocate(compos2(npoints,ncp), STAT=allocstat)
    do nc = 1,ncp
        compos2(1:npoints,nc) = composition(1:npoints,nc)
    enddo
    deallocate(cpsubg, composition, STAT=allocstat)
    
    if (errorflagmix /= 0) then  !a mixture-related error occurred:
        call RepErrorWarning(unito, errorflagmix, warningflag, &
            & errorflag, i, errorind, warningind)
    endif
    !
    !==== AIOMFAC-LLE calculation section ==========================
    !
    if (errorind == 0) then !perform AIOMFAC-LLE calculations
                            !else jump to termination section
        allocate(inputconc(nindcomp), outnames(NKNpNGS), &
        & out_data_A(8,npoints,NKNpNGS), out_data_B(8,npoints,NKNpNGS), out_LLEprop(7,npoints), LLEprop(7), LLEoutvars(7,NKNpNGS,2), name_species_TeX(NKNpNGS), STAT=allocstat)
        inputconc = 0.0_wp
        out_data_A = 0.0_wp
        out_data_B = 0.0_wp
        !--
        if (verbose) then
            write(unito,*) ""
            write(unito,'(A)') "MESSAGE from AIOMFAC: mixture defined, calculating composition points... "
            write(unito,*) ""
        endif
        
        !set AIOMFAC input and call the main AIOMFAC-LLE subroutine for all composition points:
        !note: xinputtype is .true. for mole fraction input (on component-basis, i.e. electrolytes undissociated) 
        !or .false. to indicate mass fraction input (setting from reading the input file);
        ignore_LLPS = .false.   !.false. is the default setting (set .true. to force a single-phase calculation)
        
        do pointi = 1,npoints   !loop over points, changing composition
                                !and/or temperature       
            inputconc(1:ncp) = compos2(pointi,1:ncp)
            !--
            call AIOMFAC_LLE_inout(inputconc, xinputtype, T_K(pointi), &
                & ignore_LLPS, nspecies, LLEprop, LLEoutvars, &
                & outnames, errorflag, warningflag)
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
                        out_data_A(8,pointi,nc) = real(warningflag, kind=wp)
                        out_data_B(8,pointi,nc) = real(warningflag, kind=wp)
                    endif
                endif
            enddo !nc
            !--
        enddo !pointi
        
        !
        !==== output data to file section ==================================================
        !
        !use the name of the input file to create a corresponding output file name from a string like "inputfile_0004.txt"
        i = index(filename, ".txt")
        filename = "AIOMFAC-LLE_output_"//filename(i-4:)
        !create an output text file with an overall mixture header and individual tables 
        !for all components / species (in case of ions) and separate phases:
        fname = trim(folderpathout)//trim(filename)
        call OutputLLE_TXT(fname, VersionNo, nspecies, npoints, watercompno, cpnameinp(1:nspecies), &
            & T_K(1:npoints), out_LLEprop, out_data_A, out_data_B, name_species_TeX)
        
        !-- for debugging
        write(unito,'(A)') "................................................................................"
        write(unito,'(A)') "MESSAGE from AIOMFAC: computations successfully performed."
        write(unito,'(4(A))') "Output file, ", trim(filename), " created at path: ", trim(folderpathout)
        write(unito,'(A)') "................................................................................"
        
        !
        !==== (optional) generate a few plots using Dislin ========================================
        !
        !note: this requires a linked Dislin library        
        i = index(filename, ".txt")
        filename = "AIOMFAC-LLE_"//filename(i-4:i-1)
        fname = trim(folderpathout)//trim(filename)
        call OutputLLE_plots(fname, npoints, name_species_TeX, out_LLEprop, out_data_A, out_data_B)
        
        !
        !==== Program termination section ==========================================================
        !
        deallocate(inputconc, outnames, out_LLEprop, out_data_A, out_data_B, T_K, cpnameinp, &
            & name_species_TeX, STAT=allocstat)
        if (allocated(compos2)) then
            deallocate(compos2, STAT=allocstat)
        endif
    endif !errorind
endif !file valid

write(unito,*) "+-+-+-+-+"
write(unito,'(A)') "Final warning indicator (an entry '00' means no warnings found):"
write(unito,'(I2.2)') warningind
write(unito,*) "+-+-+-+-+"
write(unito,*) ""
write(unito,*) "########"
write(unito,'(A)') "Final error indicator (an entry '00' means no errors found):"
write(unito,'(I2.2)') errorind
write(unito,*) "########"
close(unito)    !close the error log-file

write(*,*) ""
write(*,'(A,I2.2)') "MESSAGE from AIOMFAC: end of program; final error indicator: ", errorind
write(*,*) ""
!read(*,*)  !Pause; just for debugging and testing.
!
!==== the end ======================================================================
!
end program Main_AIOMFAC_LLE_prog