!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module containing input (file) data reading subroutines and output processing      *
!*   subroutines, including HTML output for the web model.                              *
!*                                                                                      *
!*   :: Author & Copyright ::                                                           *
!*   Andreas Zuend,                                                                     *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!*   -> created:        2021/07/26                                                      *
!*   -> latest changes: 2021/08/26                                                      *
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
!*   -  SUBROUTINE ReadInputFile                                                        *
!*   -  SUBROUTINE OutputLLE_TXT                                                        *
!*   -  SUBROUTINE RepErrorWarning                                                      *
!*   -  SUBROUTINE OutputLLE_plots                                                      *
!*                                                                                      *
!****************************************************************************************
module Mod_InputOutput

implicit none
private

public :: OutputLLE_TXT, OutputLLE_plots, ReadInputFile, RepErrorWarning

!============================================================================================
    contains
!============================================================================================

    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Read an input text file with information about the subgroups of each system        *
    !*   component and the compositions and temperatures for specific mixture calculations. *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2019/10/17                                                      *
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
    SUBROUTINE ReadInputFile(filepath, folderpathout, filename, ninpmax, maxpoints, unito, verbose, &
        & ncp, npoints, warningind, errorind, filevalid, cpnameinp, cpsubg, T_K, composition, xinputtype)

    USE ModSystemProp, ONLY : topsubno

    IMPLICIT NONE

    !interface variables:
    CHARACTER(LEN=3000),INTENT(INOUT) :: filepath, folderpathout
    CHARACTER(LEN=200),INTENT(INOUT) :: filename
    INTEGER(4),INTENT(IN) :: ninpmax, maxpoints
    INTEGER(4),INTENT(INOUT) :: unito
    LOGICAL(4),INTENT(IN) :: verbose
    INTEGER(4),INTENT(OUT) :: ncp, npoints
    INTEGER(4),INTENT(INOUT) :: warningind, errorind
    LOGICAL(4),INTENT(OUT) :: filevalid
    CHARACTER(LEN=200),DIMENSION(:),INTENT(OUT) :: cpnameinp  !list of assigned component names (from input file)
    INTEGER(4),DIMENSION(:,:),INTENT(OUT) :: cpsubg          !list of input component subgroups and corresponding subgroup quantities
    REAL(8),DIMENSION(:),INTENT(OUT) :: T_K                  !temperature of data points in Kelvin
    REAL(8),DIMENSION(:,:),INTENT(OUT) :: composition        !array of mixture composition points for which calculations should be run
    LOGICAL(4),INTENT(OUT) :: xinputtype
    !--
    !local variables:
    CHARACTER(LEN=:),ALLOCATABLE :: cn   !this assumes a maximum four-digit component number in the system (max. 9999); to be adjusted otherwise.
    CHARACTER(LEN=4) :: dashes, equalsigns, pluses
    CHARACTER(LEN=20) :: dummy, cnformat
    CHARACTER(LEN=50) :: txtcheck, inpfolder, outpfolder
    CHARACTER(LEN=3000) :: errlogfile, fname
    CHARACTER(LEN=20),DIMENSION(ninpmax) :: txtarray
    INTEGER(4) :: cpno, i, inpfilesize, istat, k, kinpf, qty, subg, unitx
    LOGICAL(4) :: fileopened, fileexists
    !...................................................................................

    !initialize variables
    ncp = 0
    cpsubg = 0
    cpnameinp = "none"
    composition = 0.0D0
    T_K = 298.15D0      !default temperature in [K]
    dashes = "----"
    equalsigns = "===="
    pluses = "++++"
    fileexists = .false.
    filevalid = .false.

    !==== READ INPUT data ===========================================================

    !extract filepath from the input filepath (which could include a path to the directory):
    k = LEN_TRIM(filepath)
    IF (k < 1) THEN !no file was provided
        WRITE(*,*) ""
        WRITE(*,*) "ERROR from AIOMFAC-web: no input file was provided!"
        WRITE(*,*) "An input file path needs to be provided via command line argument 2"
        WRITE(*,*) ""
        READ(*,*)
        STOP
    ENDIF

    i = INDEX(filepath, "/input_", BACK = .true.)  !find index position of input file, assuming that all input files start with "input_"
    IF (i < k .AND. i > 0) THEN !found a direcory path in UNIX file path style
        filename = TRIM(filepath(i+1:k))
        filepath = TRIM(filepath(1:k))
        kinpf = INDEX(filepath(1:i-1), "/", BACK = .true.)  !find folder length of folder containing the input file
    ELSE
        kinpf = 0
        i = INDEX(filepath, "\input_", BACK = .true.) 
        IF (i > 0 .AND. i < k) THEN !found a direcory path in WINDOWS file path style
            filename = TRIM(filepath(i+1:k))
            filepath = TRIM(filepath(1:k))
            kinpf = INDEX(filepath(1:i-1), "\", BACK = .true.)  !find folder length of folder containing the input file
        ENDIF
    ENDIF
    IF (kinpf > i-1) THEN
        kinpf = 0
    ENDIF
    inpfolder = TRIM(filepath(kinpf+1:i-1))
    kinpf = LEN_TRIM(inpfolder)

    !check / create for associated output directory:
    outpfolder = "Outputfiles/"  !initialize
    istat = INDEX(inpfolder, "inp")
    IF (istat > 0 .AND. istat < kinpf) THEN
        outpfolder = "outp"//TRIM(inpfolder(4:))//"/"  
    ELSE
        istat = INDEX(inpfolder, "Inp")
        IF (istat > 0) THEN
            outpfolder = "Outp"//TRIM(inpfolder(4:))//"/"   
        ENDIF
    ENDIF
    folderpathout = TRIM(filepath(1:i-(kinpf+1)))//TRIM(outpfolder)

    !use the name of the input file to create a corresponding output file name:
    i = INDEX(filename, ".txt") !returns starting position of string input within string filename
    !create an error-logfile associated with the input file name:
    errlogfile = "Errorlog_"//filename(i-4:)
    fname = TRIM(folderpathout)//TRIM(errlogfile)
    OPEN (NEWUNIT = unito, FILE = fname, STATUS ='UNKNOWN') !unito is the error / logfile unit
    !-----
    !check if file exists and read its content if true:
    fname = TRIM(filepath)
    INQUIRE(FILE = fname, EXIST = fileexists, SIZE = inpfilesize) !inpfilesize is the file size in [bytes]
    !Delete very large files that can only mean uploaded spam content and not actual input:
    IF (fileexists) THEN
        IF (REAL(inpfilesize, KIND=8) > 50.0D0*(ninpmax +ninpmax*maxpoints) ) THEN !likely not a valid input file
            fname = TRIM(filepath)
            OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='READ', STATUS='OLD')
            READ(unitx,*) dummy, dummy, dummy, txtcheck
            CLOSE(unitx)
            IF (.NOT. (txtcheck(1:11) == "AIOMFAC-web")) THEN !invalid file (likely spam)
                OPEN (NEWUNIT = unitx, FILE = fname, STATUS='OLD')
                CLOSE(unitx, STATUS='DELETE')  !close and delete the file
                fileexists = .false.
            ENDIF
        ENDIF
    ENDIF

    !attempt to read a supposedly valid input file:
    IF (fileexists) THEN
        IF (verbose) THEN
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: input file exists."
        ENDIF
        fname = TRIM(filepath)
        OPEN (NEWUNIT = unitx, FILE = fname, IOSTAT=istat, ACTION='READ', STATUS='OLD')
        IF (istat /= 0) THEN ! an error occurred
            WRITE(unito,*) ""
            WRITE(unito,*) "MESSAGE from AIOMFAC: an error occurred while trying to open the file! IOSTAT: ", istat
        ENDIF
        !validate file as a correct input text-file (no spam):
        READ(unitx,*) dummy, dummy, dummy, txtcheck
        IF (txtcheck(1:11) == "AIOMFAC-web") THEN
            filevalid = .true.
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: input file has passed the first line text validation and will be read."
            ENDIF
            !read input data from file:
            BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
            READ(unitx,*) dummy, dummy, dummy, dummy, dummy !read first line
            READ(unitx,*)               !read empty line 2
            READ(unitx,*) dummy, dummy  !read line 3
            READ(unitx,*) dummy         !read line 4
        ELSE !file not valid. It will be closed and deleted below
            filevalid = .false.
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: input file does not pass first line text validation and will be deleted."
            ENDIF
        ENDIF
        !loop over mixture components with variable numbers of subgroups to read (using inner loop):
        IF (filevalid) THEN
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: reading component data from input file."
            ENDIF
            k = MAX(2, CEILING(LOG10(REAL(ninpmax))) ) !determine order of magnitude digits
            WRITE(dummy,'(I0)') k
            cnformat = "(I"//TRIM(dummy)//"."//TRIM(dummy)//")" !variable format specifier
            ALLOCATE(CHARACTER(LEN=k) :: cn)
            DO ncp = 1,ninpmax
                IF (verbose .AND. istat /= 0) THEN
                    WRITE(unito,*) ""
                    WRITE(unito,*) "MESSAGE from AIOMFAC: file end found in input file at ncp = ", ncp
                ENDIF
                READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
                IF (txtcheck(1:4) == pluses .OR. istat /= 0) THEN !"++++" indicates end of this components definition part
                    EXIT !exit ncp DO-loop
                ELSE !in this case, argument 3 of txtcheck is the component no.:
                    BACKSPACE unitx !jump back to beginning of record (to the beginning of the line)
                    READ(unitx,*) dummy, dummy, cpno !read line with component's number
                ENDIF
                READ(unitx,*) dummy, dummy, cpnameinp(cpno) !read line with component's name
                IF (LEN_TRIM(cpnameinp(cpno)) < 1) THEN !no component name was assigned, generate a default name
                    WRITE(cn,cnformat) cpno
                    cpnameinp(cpno) = "cp_"//cn
                ENDIF
                IF (verbose) THEN
                    WRITE(unito,*) ""
                    WRITE(unito,*) "MESSAGE from AIOMFAC: found a component cpno =", cpno
                ENDIF
                DO !until exit
                    READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line
                    !check whether another subgroup is present or not
                    IF (txtcheck(1:4) == dashes .OR. istat /= 0) THEN !"----" indicates no more subgroups of this component
                        EXIT !leave the inner DO-loop
                    ELSE !else continue reading the next subgroup
                        BACKSPACE unitx
                        READ(unitx,*) dummy, dummy, dummy, subg, qty  !continue reading line with subgroup no. and corresp. quantity
                        cpsubg(cpno,subg) = cpsubg(cpno,subg)+qty
                    ENDIF
                ENDDO
            ENDDO
            ncp = ncp-1 !ncp-1 is the number of different components in this mixture
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: component data read."
                WRITE(unito,*) "MESSAGE from AIOMFAC: number of components: ", ncp
            ENDIF
            IF (ncp == ninpmax) THEN
                READ(unitx,*,IOSTAT=istat) txtcheck !read only first argument on this line for subsequent check
                IF (txtcheck(1:4) /= pluses) THEN
                    errorind = 34
                    filevalid = .false.
                    WRITE(*,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
                    WRITE(*,*) "AIOMFAC ERROR 34: check whether ninpmax value is too small."
                    WRITE(unito,*) ""
                    WRITE(unito,*) "AIOMFAC ERROR 34: maximum number of input components reached while reading input file."
                ENDIF
            ENDIF
            IF (filevalid) THEN
                READ(unitx,*) dummy, dummy     !read mixture composition line
                READ(unitx,*) dummy, dummy, i  !read mass fraction? line
                IF (i == 1) THEN !composition in mass fractions
                    xinputtype = .false.
                ELSE !composition in mole fractions
                    xinputtype = .true.
                ENDIF
                READ(unitx,*) dummy, dummy, i !read mole fraction? line
                READ(unitx,*) dummy !read dashes line
                !now read the lines with the mixture composition points:
                READ(unitx,*) txtarray(1:ncp) !read header line of composition table
                DO npoints = 1,maxpoints !or until exit
                    READ(unitx,*,IOSTAT=istat) txtcheck !read only the first argument on this line
                    IF (txtcheck(1:4) == equalsigns .OR. istat /= 0) THEN !"====" indicates no more composition points (and last line of input file)
                        EXIT !leave the DO-loop (normal exit point)
                    ELSE IF (IACHAR(txtcheck(1:1)) > 47 .AND. IACHAR(txtcheck(1:1)) < 58) THEN !validate that the data is actual intended input and not some sort of text field spam).
                        BACKSPACE unitx
                        READ(unitx,*) txtcheck, dummy !read only the first two arguments on this line
                        IF (IACHAR(dummy(1:1)) > 47 .AND. IACHAR(dummy(1:1)) < 58) THEN !it is a number
                            BACKSPACE unitx
                            READ(unitx,*) i, T_K(npoints), composition(npoints,2:ncp) !read the temperature in [K] and composition values of the components [2:ncp] into the array
                        ELSE
                            IF (npoints == 1) THEN
                                filevalid = .false. !file is not valid due to incorrect input in text field
                                EXIT
                            ELSE
                                warningind = 31
                                EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                            ENDIF
                        ENDIF
                    ELSE
                        IF (npoints == 1) THEN
                            filevalid = .false. !file is not valid due to incorrect input in text field
                            EXIT
                        ELSE
                            warningind = 31
                            EXIT !abnormal exit point of loop due to non-number entries at a certain line in the text field
                        ENDIF
                    ENDIF
                ENDDO
                npoints = npoints-1 !the number of composition points
                IF (verbose) THEN
                    WRITE(unito,*) ""
                    WRITE(unito,*) "MESSAGE from AIOMFAC: composition points read."
                    WRITE(unito,*) "MESSAGE from AIOMFAC: number of points: ", npoints
                ENDIF
                IF (.NOT. filevalid) THEN
                    !close and delete file from server:
                    CLOSE(unitx, STATUS = 'DELETE')
                ELSE
                    CLOSE(unitx)
                ENDIF
                !assign component 01 its composition values based on the rest of the component's data:
                DO i = 1,npoints
                    composition(i,1) = 1.0D0-SUM(composition(i,2:ncp))
                ENDDO
            ENDIF !filevalid2
        ENDIF !filevalid1
    ELSE
        WRITE(unito,*) ""
        WRITE(unito,*) "ERROR in AIOMFAC: Input file does not exist at expected location: ", TRIM(filename)
        WRITE(unito,*) ""
    ENDIF !fileexists

    IF (fileexists) THEN
        IF (.NOT. filevalid) THEN !input file contains errors or is completely invalid (submitted spam etc.)
            IF (verbose) THEN
                WRITE(unito,*) ""
                WRITE(unito,*) "MESSAGE from AIOMFAC: input file did not pass full validation and may be a spam file."
                WRITE(unito,*) "MESSAGE from AIOMFAC: the input file will be deleted to prevent spam files and malicious code from occupying the server."
                WRITE(unito,*) ""
            ENDIF
            INQUIRE(FILE = fname, OPENED = fileopened)
            IF (errorind == 0) THEN
                errorind = 32
                !close and delete file from server:
                IF (fileopened) THEN
                    CLOSE(unitx, STATUS = 'DELETE')
                ENDIF
            ELSE
                IF (fileopened) THEN
                    CLOSE(unitx)
                ENDIF
            ENDIF
        ELSE
            !check for any valid points for calculations before initializing AIOMFAC:
            IF (npoints < 1) THEN !no points input for calculations
                errorind = 33 !no input in text field..
                IF (verbose) THEN
                    WRITE(unito,*) ""
                    WRITE(unito,*) "MESSAGE from AIOMFAC: no composition points have been entered in the text field. There is nothing to calculate."
                    WRITE(unito,*) ""
                ENDIF
                filevalid = .false.
            ENDIF
        ENDIF
    ENDIF

    END SUBROUTINE ReadInputFile
    !========================================================================================
    
        
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Create an output ASCII text file with an overall mixture data header and           *
    !*   individual tables for all components / species (in case of ions).                  *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2020/07/18                                                      *
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
    SUBROUTINE OutputLLE_TXT(fname, VersionNo, nspecmax, npoints, watercompno, cpnameinp, T_K, &
        & out_LLEprop, out_data_A, out_data_B, name_species_TeX)

    USE ModSystemProp, ONLY : compname, compsubgroups, compsubgroupsTeX, NGS, NKNpNGS, ninput, nneutral
    USE ModSubgroupProp, ONLY : subgrname, subgrnameTeX

    IMPLICIT NONE
    !interface arguments:
    CHARACTER(LEN=*),INTENT(IN) :: fname 
    CHARACTER(LEN=*),INTENT(IN) :: VersionNo
    INTEGER(4),INTENT(IN) ::  nspecmax, npoints, watercompno
    CHARACTER(LEN=*),DIMENSION(:),INTENT(IN) :: cpnameinp           !component names (from input file)
    REAL(8),DIMENSION(npoints),INTENT(IN) :: T_K
    REAL(8),DIMENSION(8,npoints),INTENT(IN) :: out_LLEprop
    REAL(8),DIMENSION(8,npoints,NKNpNGS),INTENT(IN) :: out_data_A, out_data_B
    CHARACTER(LEN=*),DIMENSION(:),INTENT(OUT) :: name_species_TeX   !species names (for output)
    !local variables:
    CHARACTER(LEN=:),ALLOCATABLE :: cn
    CHARACTER(LEN=5) :: tlen
    CHARACTER(LEN=1) :: phase_char
    CHARACTER(LEN=50) :: subntxt, Iformat
    CHARACTER(LEN=150) :: cnformat, horizline, txtn, tablehead
    CHARACTER(LEN=3000) :: txtsubs, mixturestring
    INTEGER(4) ::  i, k, kms, phase_no, pointi, qty, unitx
    REAL(8) :: RH
    REAL(8),DIMENSION(8,npoints,NKNpNGS) :: out_data
    !...................................................................................

    k = MAX(2, CEILING(LOG10(REAL(nspecmax))) )  
    WRITE(tlen,'(I0)') k
    Iformat = "I"//TRIM(tlen)//"."//TRIM(tlen)  !variable integer specifier
    cnformat = "("//Iformat//")"                !constructed format specifier
    ALLOCATE(CHARACTER(LEN=k) :: cn)
    
    out_data = out_data_A           !first use the data for phase alpha (A);

    !create a character string of the mixture as a series of its components:
    kms = LEN(mixturestring)
    mixturestring = TRIM(cpnameinp(1)) !first component
    !loop over all further components / species:
    DO i = 2,nspecmax
        IF (INT(out_data(7,1,i)) == 0) THEN   !neutral component
            txtn = TRIM(ADJUSTL(cpnameinp(i)))
        ELSE    !ion
            IF (INT(out_data(7,npoints,i)) > 1) THEN
                txtn = TRIM( ADJUSTL(subgrname(INT(out_data(7,npoints,i) ))) )
            ELSE
                txtn = "unknown_sub"
            ENDIF
        ENDIF
        k = LEN_TRIM(mixturestring) +LEN_TRIM(' + '//TRIM(txtn) )
        IF (k < kms - 50) THEN
            mixturestring = TRIM(mixturestring)//' + '//TRIM(ADJUSTL(txtn))
        ELSE
            qty = nspecmax -i
            WRITE(subntxt,'(I0)') qty
            mixturestring = TRIM(mixturestring)//' + '//TRIM(subntxt)//' additional components ...'
            EXIT
        ENDIF
    ENDDO
    mixturestring = ADJUSTL(mixturestring)

    OPEN (NEWUNIT = unitx, FILE = fname, STATUS = "UNKNOWN")
    WRITE(unitx,'(A)') "==========================================================================================================="
    WRITE(unitx,'(A)') "AIOMFAC-LLE (web), version "//VersionNo
    WRITE(unitx,'(A)') "==========================================================================================================="
    WRITE(unitx,*) ""
    mixturestring = "'"//TRIM(mixturestring)//"'"
    WRITE(unitx, '(A,A)') "Mixture name:  ", TRIM(mixturestring)
    WRITE(unitx, '(A,I0.2)') "Number of independent input components: ", ninput
    WRITE(unitx, '(A,I0.2)') "Number of different neutral components: ", nneutral
    WRITE(unitx, '(A,I0.2)') "Number of different inorganic ions    : ", NGS
    WRITE(unitx,*) ""
    WRITE(unitx,'(A)') "The AIOMFAC-LLE phase composition data are tabulated for each component/species individually."
    WRITE(unitx,'(A)') "Selected system properties (LLE metrics) are listed first, followed by data for phases A, B."
    WRITE(unitx,*) ""
    WRITE(unitx,'(A)') '---  Table key  -------------------------------------------------------------------------------------------'
    WRITE(unitx,'(A)') 'no.               : composition point number;                                                        '
    WRITE(unitx,'(A)') 'T [K]             : absolute temperature;                                                            '
    WRITE(unitx,'(A)') 'RH [%]            : relative humidity in equilibrium with the liquid mixture (bulk conditions);      '
    WRITE(unitx,'(A)') '#ph [-]           : number of coexisting liquid phases (1 or 2);                                     '
    WRITE(unitx,'(A)') 'fracMassA [-]     : fraction of total mass residing in phase alpha (< 1 in LLE case);                '
    WRITE(unitx,'(A)') 'omega [-]         : (in LLE case) molar phase amount ratio (sum[n_beta] : sum[n_alpha]);             '
    WRITE(unitx,'(A)') 'phi [-]           : (in LLE case) solvent mass ratio (phase beta : alpha);                           '
    WRITE(unitx,'(A)') 'Gmix_diff [J/mol] : (in LLE case) molar Gibbs energy difference (1-phase - 2-phase state);           '
    WRITE(unitx,'(A)') 'LLE_res [-]       : weighted cumulative relative activity deviations between phases A and B serving  '
    WRITE(unitx,'(A)') '                    as check for the numerical resolution of an LLE state (tiny abs. value is good); '
    WRITE(unitx,'(A)') '                                                                                                     '
    WRITE(unitx,'(A)') 'w_inp(j) [-]      : input weight fraction (mass fraction) of species "j" in the overall mixture;     '
    WRITE(unitx,'(A)') 'w(j) [-]          : mass fraction of species "j" in indicated phase;                                 '
    WRITE(unitx,'(A)') 'x_i(j) [-]        : mole fraction of species "j" in phase, calculated on the basis of completely     '
    WRITE(unitx,'(A)') '                    dissociated inorganic ions; exception: the partial dissociation of bisulfate     '
    WRITE(unitx,'(A)') '                    (HSO4- <--> H+ + SO4--) is explicitly considered when present in the mixture;    '
    WRITE(unitx,'(A)') 'm_i(j) [mol/kg]   : molality of species "j" [mol(j)/(kg solvent mixture)] in phase, where            '
    WRITE(unitx,'(A)') '                    "solvent mixture" refers to the electrolyte-free mixture (water + organics);     '
    WRITE(unitx,'(A)') 'a_coeff_x(j) [-]  : activity coefficient of "j" in phase, defined on mole fraction basis (used for   '
    WRITE(unitx,'(A)') '                    non-ionic components) with pure (liquid) component "j" reference state;          '
    WRITE(unitx,'(A)') 'a_coeff_m(j) [-]  : activity coefficient of "j" in phase, defined on molality basis (used for inorg. '
    WRITE(unitx,'(A)') '                    ions) with reference state of infinite dilution of "j" in pure water;            '
    WRITE(unitx,'(A)') 'a_x(j) [-]        : activity of "j", defined on mole fraction basis (pure component reference state);'
    WRITE(unitx,'(A)') 'a_m(j) [-]        : activity of "j", defined on molality basis (used for inorg. ions) with reference '
    WRITE(unitx,'(A)') '                    state of infinite dilution of "j" in pure water;                                 '
    WRITE(unitx,'(A)') 'flag              : error/warning flag, a non-zero value (error/warning number) indicates that a     '
    WRITE(unitx,'(A)') '                    numerical issue or a warning occurred at this data point, suggesting evaluation  '
    WRITE(unitx,'(A)') '                    with caution (warnings) or exclusion (errors) of this data point.                '
    WRITE(unitx,'(A)') '-----------------------------------------------------------------------------------------------------------'
    WRITE(unitx,*) ''
    WRITE(unitx,*) ''
    !--
    horizline = "-----------------------------------------------------------------------------------------------------------"
    !--
    
    !output selected overall LLE properties data:
    WRITE(unitx,'(A)') ""
    WRITE(unitx,'(A)') "Selected system properties in equilibrium state (LLE or single phase) "
    !write column headers:
    WRITE(unitx,'(A)') TRIM(horizline)
    tablehead = "no.  T_[K]  RH_[%]  #ph   fracMassA    omega        phi      Gmix_diff_[J/mol]  LLE_res  flag"
    WRITE(unitx,'(1X,A)') TRIM(tablehead)
    !--
    WRITE(unitx,'(A)') TRIM(horizline)
    !write data to table:
    DO pointi = 1,npoints  !loop over composition points
        IF (watercompno > 0) THEN
            RH = out_data(6,pointi,watercompno)*100.0D0 !RH in %
            IF (RH > 1000.0D0 .OR. RH < 0.0D0) THEN
                RH = -99.99D0
            ENDIF
        ELSE
            RH = 0.0D0
        ENDIF
        WRITE(unitx,'(1X,I0.3,1X,F7.2,1X,F6.2,3X,I0,2X,5(ES12.5,1X),1X,I2)') pointi, T_K(pointi), RH, &
            & INT(out_LLEprop(1,pointi)), out_LLEprop(2:5,pointi), out_LLEprop(8,pointi), INT(out_data(8,pointi,1))
    ENDDO !pointi
    WRITE(unitx,'(A)') TRIM(horizline)
    WRITE(unitx,*) ""
    WRITE(unitx,'(A)') ""
    
    DO phase_no = 1,2       !loop over phases alpha and beta for output in same file;
        
        SELECT CASE(phase_no)
        CASE(1)
            phase_char = "A"
            WRITE(unitx,'(A)') "================="
            WRITE(unitx,'(A)') " Phase alpha (A) "
            WRITE(unitx,'(A)') "================="
        CASE(2)
            IF (ANY(out_LLEprop(1,:) > 1.0D0)) THEN
                WRITE(unitx,*) ""
                phase_char = "B"
                out_data = out_data_B
                WRITE(unitx,'(A)') "================="
                WRITE(unitx,'(A)') " Phase beta (B)  "
                WRITE(unitx,'(A)') "================="
            ELSE
                WRITE(unitx,*) ""
                WRITE(unitx,'(A)') "============================================================================="
                WRITE(unitx,'(A)') " No LLE predicted for any of the input compositions (==> no phase beta data) "
                WRITE(unitx,'(A)') "============================================================================="
                EXIT        !leave the phase_no loop
            ENDIF
        END SELECT
        
        !write individual data tables for each component / ionic species:
        DO i = 1,nspecmax
            WRITE(cn,cnformat) i !component / species number as character string
            WRITE(unitx,*) ""
            !distinguish between neutral components and ionic species:
            IF (INT(out_data(7,1,i)) == 0) THEN !neutral component
                WRITE(unitx,'(A,A,",",1X,I0.2)') "Phase, component #    : ", phase_char, i
                name_species_TeX(i) = TRIM(compname(i))     !save processed component name
                txtn = "'"//TRIM(ADJUSTL(compname(i)))//"'"
                WRITE(unitx,'(A,A)') "Component's name      : ", TRIM(txtn)
                txtsubs = "'"//TRIM(compsubgroups(i))//"'"
                WRITE(unitx,'(A,A)') "Component's subgroups : ", TRIM(txtsubs)
                txtsubs = "'"//TRIM(compsubgroupsTeX(i))//"'"
                WRITE(unitx,'(A,A)') "Subgroups, TeX format : ", TRIM(txtsubs)
                !write table column headers:
                WRITE(unitx,'(A)') TRIM(horizline)
                tablehead = "no.   T_[K]   RH_[%]  w_inp("//cn//")   w("//cn//")        x_i("//cn//")      m_i("//cn//")     a_coeff_x("//cn//")  a_x("//cn//")    flag"
                WRITE(unitx,'(1X,A)') TRIM(tablehead)
                !--
            ELSE IF (INT(out_data(7,1,i)) < 240) THEN !cation
                WRITE(unitx,'(A,A,",",1X,I0.2)') "Phase, species #      : ", phase_char, i
                subntxt = TRIM(ADJUSTL(subgrname(INT(out_data(7,1,i)))))
                qty = LEN_TRIM(subntxt)
                txtn = ADJUSTL(subntxt(2:qty-1)) !to print the ion name without the enclosing parenthesis ()
                txtn = "'"//TRIM(txtn)//"'"
                WRITE(unitx,'(A,A)')  "Cation's name         : ", TRIM(txtn)
                txtsubs = "'"//TRIM(subntxt)//"'"
                WRITE(unitx,'(A,A)') "Cation's subgroups    : ", TRIM(txtsubs)
                subntxt = TRIM(ADJUSTL(subgrnameTeX(INT(out_data(7,1,i)))))
                name_species_TeX(i) = subntxt(2:LEN_TRIM(subntxt)-1)    !save ion species name TeX string
                txtsubs = "'"//TRIM(subntxt)//"'"
                WRITE(unitx,'(A,A)') "Subgroups, TeX format : ", TRIM(txtsubs)
                !write table column headers:
                WRITE(unitx,'(A)') TRIM(horizline)
                tablehead = "no.   T_[K]   RH_[%]  w_inp("//cn//")   w("//cn//")        x_i("//cn//")      m_i("//cn//")     a_coeff_m("//cn//")  a_m("//cn//")    flag"
                WRITE(unitx,'(1X,A)') TRIM(tablehead)
                !--
            ELSE IF (INT(out_data(7,1,i)) > 240) THEN !anion
                WRITE(unitx,'(A,A,",",1X,I0.2)') "Phase, species #      : ", phase_char, i
                subntxt = TRIM( ADJUSTL( subgrname(INT(out_data(7,1,i))) ) )
                qty = LEN_TRIM(subntxt)
                txtn = ADJUSTL(subntxt(2:qty-1))
                txtn = "'"//TRIM(txtn)//"'"
                WRITE(unitx,'(A,A)')  "Anion's name          : ", TRIM(txtn)
                txtsubs = "'"//TRIM(subntxt)//"'"
                WRITE(unitx,'(A,A)') "Anion's subgroups     : ", TRIM(txtsubs)
                subntxt = TRIM( ADJUSTL( subgrnameTeX(INT(out_data(7,1,i))) ) )
                name_species_TeX(i) = subntxt(2:LEN_TRIM(subntxt)-1)    !save ion species name TeX string
                txtsubs = "'"//TRIM(subntxt)//"'"
                WRITE(unitx,'(A,A)') "Subgroups, TeX format : ", TRIM(txtsubs)
                !write table column headers:
                WRITE(unitx,'(A)') TRIM(horizline)
                tablehead = "no.   T_[K]   RH_[%]  w_inp("//cn//")   w("//cn//")        x_i("//cn//")      m_i("//cn//")     a_coeff_m("//cn//")  a_m("//cn//")    flag"
                WRITE(unitx,'(1X, A)') TRIM(tablehead)
                !--
            ELSE
                !error
            ENDIF
            !--
            WRITE(unitx,'(A)') TRIM(horizline)
            !write data to table:
            DO pointi = 1,npoints  !loop over composition points
                IF (watercompno > 0) THEN
                    RH = out_data(6,pointi,watercompno)*100.0D0 !RH in %
                    IF (RH > 1000.0D0 .OR. RH < 0.0D0) THEN
                        RH = -99.99D0
                    ENDIF
                ELSE
                    RH = 0.0D0
                ENDIF
                WRITE(unitx,'(1X,I0.3,2X,F7.2,1X,F6.2,1X,6(ES12.5,1X),1X,I2)') pointi, T_K(pointi), RH, out_data(1:6,pointi,i), INT(out_data(8,pointi,i))
            ENDDO !pointi
            WRITE(unitx,'(A)') TRIM(horizline)
            WRITE(unitx,*) ""
        ENDDO
        WRITE(unitx,*) ""
    ENDDO !phase_no
    !--
    WRITE(unitx,*) ""
    WRITE(unitx,'(A)') "==========================================================================================================="
    CLOSE(unitx)
        
    END SUBROUTINE OutputLLE_TXT
    !========================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Report error and warning messages to the error logfile from a list of defined      *
    !*   cases.                                                                             *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2011                                                            *
    !*   -> latest changes: 2018/08/08                                                      *
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
    SUBROUTINE RepErrorWarning(unito, errorflagmix, warningflag, errorflagcalc, pointi, errorind, warningind)

    IMPLICIT NONE
    !interface variables:
    INTEGER(4),INTENT(IN) :: unito, errorflagmix, warningflag, errorflagcalc, pointi
    INTEGER(4),INTENT(OUT) :: errorind, warningind
    !...................................................................................

    IF (errorflagmix /= 0) THEN !some mixture related error occurred:
        SELECT CASE(errorflagmix)
        CASE(1)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 1: mixture related."
            WRITE(unito,*) "An organic main group <-> cation interaction parameter"
            WRITE(unito,*) "is not defined for the requested mixture. "
            WRITE(unito,*) "Please check your mixture components for available parameters"
            WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(2)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 2: mixture related."
            WRITE(unito,*) "An organic main group <-> anion interaction parameter"
            WRITE(unito,*) "is not defined for the requested mixture. "
            WRITE(unito,*) "Please check your mixture components for available parameters"
            WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(9)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 9: mixture related."
            WRITE(unito,*) "At least one cation <-> anion interaction parameter is "
            WRITE(unito,*) "not defined for the requested mixture. "
            WRITE(unito,*) "Please check all ion combinations for available parameters"
            WRITE(unito,*) "stated in the AIOMFAC interaction matrix."
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(13)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 13: Incorrect hydroxyl group assignment."
            WRITE(unito,*) "At least one component containing (CH_n[(OH)]) groups"
            WRITE(unito,*) "has been assigned an incorrect number of (OH) groups."
            WRITE(unito,*) "Note that the notation of a CHn group bonded to an OH"
            WRITE(unito,*) "group does not include the OH group; rather the hydroxyl"
            WRITE(unito,*) "groups have to be defined separately.                   "
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(14)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 14: Missing short-range ARR parameter."
            WRITE(unito,*) "A neutral main group <-> main group interaction coeff."
            WRITE(unito,*) "of this particular mixture is not available in the SR."
            WRITE(unito,*) "part of the model."
            WRITE(unito,*) "Check your organic components and their subgroups in"
            WRITE(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE(15)
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR 15: Missing short-range BRR parameter."
            WRITE(unito,*) "A neutral main group <-> main group interaction coeff."
            WRITE(unito,*) "of this particular mixture is not available in the SR."
            WRITE(unito,*) "part of the model for 3-parameter temperature dependence."
            WRITE(unito,*) "Check your organic components and their subgroups in"
            WRITE(unito,*) "comparison to available subgroups in the AIOMFAC matrix."
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        CASE DEFAULT
            WRITE(unito,*) ""
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) "AIOMFAC ERROR XX: an undefined mixture error occurred!"
            WRITE(unito,*) "errorflagmix = ", errorflagmix
            WRITE(unito,*) "======================================================="
            WRITE(unito,*) ""
        END SELECT
        errorind = errorflagmix
    ELSE 
        !check warnings and errors related to occurences during specific data point calculations:
        IF (warningflag > 0) THEN
            SELECT CASE(warningflag)
            CASE(10)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC WARNING 10: Temperature range related."
                WRITE(unito,*) "At least one data point has a set temperature outside of"
                WRITE(unito,*) "the recommended range for model calculations of "
                WRITE(unito,*) "electrolyte-containing mixtures. This may be intended, "
                WRITE(unito,*) "but caution is advised as AIOMFAC is not designed to  "
                WRITE(unito,*) "perform well at this temperature."
                WRITE(unito,*) "Data point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(11)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC WARNING 11: Temperature range related."
                WRITE(unito,*) "At least one data point has a set temperature outside of"
                WRITE(unito,*) "the recommended range for model calculations of "
                WRITE(unito,*) "electrolyte-free organic mixtures. This may be intended,"
                WRITE(unito,*) "but caution is advised as AIOMFAC is not designed to "
                WRITE(unito,*) "perform well at this temperature."
                WRITE(unito,*) "Data point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(16)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC-VISC WARNING 16: Mixture viscosity issue.      "
                WRITE(unito,*) "Note that mixture viscosity is currently not computed  " 
                WRITE(unito,*) "for electrolyte-containing mixtures. Therefore, an     "
                WRITE(unito,*) "unrealistic mixture viscosity of                       "
                WRITE(unito,*) "log_10(eta/[Pa.s]) = -999.999 is output.               "
                WRITE(unito,*) "Data point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE DEFAULT
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC WARNING XX: an undefined WARNING occurred!"
                WRITE(unito,*) "warningflag = ", warningflag
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            END SELECT
            warningind = warningflag
        ENDIF !warningflag
        IF (errorflagcalc > 0) THEN
            SELECT CASE(errorflagcalc)
            CASE(3)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR 3: Mixture composition related."
                WRITE(unito,*) "Composition data for this point is missing or incorrect!"
                WRITE(unito,*) "The sum of the mole or mass fractions of all components"
                WRITE(unito,*) "has to be equal to 1.0 and individual mole or mass "
                WRITE(unito,*) "fractions have to be positive values <= 1.0!"
                WRITE(unito,*) "Composition point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(4,5)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR 4: Mixture composition related."
                WRITE(unito,*) "Composition data for this point is incorrect!"
                WRITE(unito,*) "The sum of the mole or mass fractions of all components"
                WRITE(unito,*) "has to be equal to 1.0 and individual mole or mass "
                WRITE(unito,*) "fractions have to be positive values <= 1.0!"
                WRITE(unito,*) "Composition point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(6,7)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR 6: Numerical issue."
                WRITE(unito,*) "A numerical issue occurred during computation of the"
                WRITE(unito,*) "data points flagged in the output tables."
                WRITE(unito,*) "This error was possibly caused due to input of very"
                WRITE(unito,*) "high electrolyte concentrations."
                WRITE(unito,*) "Composition point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(8)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR 8: Mixture composition related."
                WRITE(unito,*) "At least one neutral component must be present in the  "
                WRITE(unito,*) "system! "
                WRITE(unito,*) "Composition point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE(12)
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR 12: Charge neutrality violated."
                WRITE(unito,*) "The mixture violates the electrical charge neutrality  "
                WRITE(unito,*) "condition (moles cation*[cation charge] =              "
                WRITE(unito,*) "                          moles anion*[anion charge]). "
                WRITE(unito,*) "Make sure that selected integer amounts of cation and  "
                WRITE(unito,*) "anion 'subgroups' fulfill the charge balance (in the   "
                WRITE(unito,*) "inorganic component definition of the input file).     "
                WRITE(unito,*) "Composition point no.: ", pointi
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            CASE DEFAULT
                WRITE(unito,*) ""
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) "AIOMFAC ERROR XX: an undefined calculation error occurred!"
                WRITE(unito,*) "errorflagcalc = ", errorflagcalc
                WRITE(unito,*) "======================================================="
                WRITE(unito,*) ""
            END SELECT
            errorind = errorflagcalc
        ENDIF !errorflagcalc
    ENDIF !errorflagmix
        
    END SUBROUTINE RepErrorWarning
    !========================================================================================
    
    
    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Generate several plot examples using equilibrium (LLE or single-phase) calculation *
    !*   outputs. Note: this requires a linked Dislin library as part of the program        *
    !*   building steps.                                                                    *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Andi Zuend, (andi.zuend@gmail.com)                                                 *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University (2013 - present)         *
    !*                                                                                      *
    !*   -> created:        2021/08/01                                                      *
    !*   -> latest changes: 2021/08/26                                                      *
    !*                                                                                      *
    !****************************************************************************************
    subroutine OutputLLE_plots(fname, npoints, name_species_TeX, out_LLEprop, out_data_A, out_data_B)

        use Mod_NumPrec, only : wp
        use Mod_Dislin_plots, only : add_plot_xydata, dislin_plot
    
        implicit none
        !interface arguments:
        character(len=*),intent(in) :: fname 
        integer,intent(in) :: npoints
        character(len=200),dimension(:),intent(in) :: name_species_TeX
        real(wp),dimension(:,:),intent(in) :: out_LLEprop
        real(wp),dimension(:,:,:),intent(in) :: out_data_A, out_data_B
    
        !local variables:    
        character(len=75) :: xlabel, ylabel, legend_text
        character(len=1) :: char_phase
        character(len=3) :: act_basis
        integer :: i, iloc1, iloc2, istat, nphases, nc, ncolors, nspecies, nph, u1
        integer,dimension(3),parameter :: rgb_blue = [40, 40, 255], &
            & rgb_lightblue = [0, 170, 255], rgb_purple = [112, 51, 173]
        integer,dimension(3) :: rgb_set
        real(wp),dimension(npoints) :: mfracA
        real(wp),dimension(2) :: xax_lim, yax_lim
        real(wp),dimension(:,:),allocatable :: col_palette
        logical :: LLPSpresent
        !....................................
    
        !determine presence of phase separation in output data:
        if (any(out_LLEprop(1,:) > 1.0_wp)) then
            LLPSpresent = .true.
            nphases = 2
            !compute fraction of total mass in phase alpha:
            mfracA = out_LLEprop(2,:)
            !determine index limits of LLPS (assuming contiguous range):
            iloc1 = minloc(out_LLEprop(1,:), MASK = out_LLEprop(1,:) > 1.0_wp, DIM=1)
            iloc2 = minloc(out_LLEprop(1,:), MASK = out_LLEprop(1,:) > 1.0_wp, BACK=.true., DIM=1)
        else
            LLPSpresent = .false.
            nphases = 1
        endif
        
        !load RGB color table for use with certain plots:
        allocate(col_palette(3,20))
        !./Color_palettes/RGBColTabBlackPurpleRedYellow.dat !256 colors
        !./Color_palettes/RGBColTabViridisPurpleBlueGreenYellow.dat !256 colors
        !./Color_palettes/RGBTwentyDistinctCol.dat !22 colors
        ncolors = 20
        open(newunit=u1, file='./Color_palettes/RGBTwentyDistinctCol.dat', iostat=istat, action='read', status='old')
        do i = 1,ncolors
            read(u1,*,iostat=istat) col_palette(1:3,i)
            if (istat /= 0) exit
        enddo
        close(u1)
        if (.NOT. any(col_palette > 1.1_wp)) then
            !scale to 0 to 255 range:
            col_palette = col_palette*255.0_wp
        endif

        !(a) Plot equilibrium-state fraction of total mass in phase alpha vs. input mass frac. of water:
        if (LLPSpresent) then  !plot only of interest in case of LLPS being present
            call add_plot_xydata(xv = out_data_A(1,:,1), yv = out_LLEprop(2,:), &
                & ltext='mass frac. in phase $\alpha$', pen_wid=8.0_wp, rgb_col=rgb_purple, &
                & lstyle='solid', plot_symb='curve')

            !set overall plot properties and generate plot:
            xlabel = 'input mass fraction of water, $w_w^t$'
            ylabel = 'fraction of mass in $\alpha$'
            !set/limit x-and y-axis ranges for plot:
            xax_lim = [0.0_wp, 1.0_wp]
            yax_lim = [0.0_wp, 1.0_wp]
            call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, xaxis_limits=xax_lim, yaxis_limits=yax_lim, &
                & legend_position=3, metafile='pdf', out_file_name=trim(fname)//'_massfracA_vs_ww')
        endif
        !------------------------------------------------------------

        !(b) Plot equilibrium-state fraction of total mass in phase alpha vs. water activity:
        if (LLPSpresent) then  !plot only of interest in LLPS case
            call add_plot_xydata(xv = out_data_A(6,:,1), yv = out_LLEprop(2,:), &
                & ltext='mass frac. in phase $\alpha$', pen_wid=8.0_wp, rgb_col=rgb_purple, &
                & lstyle='solid', plot_symb='curve')

            !set overall plot properties and generate plot:
            xlabel = 'water activity, $a_w$'
            ylabel = 'fraction of mass in $\alpha$'
            !set/limit x-and y-axis ranges for plot:
            xax_lim = [0.0_wp, 1.0_wp]
            yax_lim = [0.0_wp, 1.0_wp]
            call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, xaxis_limits=xax_lim, yaxis_limits=yax_lim, &
                & legend_position=3, metafile='pdf', out_file_name=trim(fname)//'_massfracA_vs_aw')
        endif
        !------------------------------------------------------------

        !(c) Plot mass fractions of water vs. water activity
        !    (total and for each phase in LLPS case):
        call add_plot_xydata(xv=out_data_A(6,:,1), yv=out_data_A(1,:,1), &
            & ltext='input mass fraction $w_w^t$', pen_wid=4.0_wp, &
            & lstyle='dotted', plot_symb='curve')

        call add_plot_xydata(xv=out_data_A(6,:,1), yv=out_data_A(2,:,1), &
            & ltext='$w_w^\alpha$ (phase $\alpha$)', &
            & pen_wid=8.0_wp, rgb_col=rgb_blue, lstyle='solid', plot_symb='curve')

        if (LLPSpresent) then
            call add_plot_xydata(xv=out_data_B(6,iloc1:iloc2,1), yv=out_data_B(2,iloc1:iloc2,1), &
                & ltext='$w_w^\beta$ (phase $\beta$)', &
                & pen_wid=8.0_wp, rgb_col=rgb_lightblue, lstyle='dashed_medium', plot_symb='curve')
        endif

        !set overall plot properties and generate plot:
        xlabel = 'water activity, $a_w$'
        ylabel = 'mass fraction of water'
        yax_lim = [0.0_wp, 1.0_wp]
        call dislin_plot(xlabel, ylabel, yaxis_mod=0.6_wp, yaxis_limits=yax_lim, legend_position=3, &
            & metafile='pdf', out_file_name=trim(fname)//'_wwater_vs_aw')
        !------------------------------------------------------------
        
        !(d) Plot activity coefficients of species, separately per phase:
        nspecies = size(name_species_TeX)
        do nph = 1,nphases
            select case(nph)
            case(1) !phase alpha
                char_phase = 'A'
                do nc = 1,nspecies
                    if (int(out_data_A(7,1,nc)) == 0) then  !neutral (non-ion) component
                        act_basis = '(x)'
                    else
                        act_basis = '(m)'
                    endif
                    legend_text = 'log$_{10}$[$\gamma^{\,\alpha , '//act_basis//'}$('//trim(name_species_TeX(nc))//')]'
                    rgb_set = int(col_palette(1:3, min(nc,ncolors)))
                    call add_plot_xydata(xv = out_data_A(1,:,1), yv = log10(out_data_A(5,:,nc)), &
                        & ltext=legend_text, pen_wid=6.0_wp, rgb_col=rgb_set, &
                        & lstyle='solid', plot_symb='curve')
                enddo
                !..
            case(2) !phase beta
                char_phase = 'B'
                do nc = 1,size(name_species_TeX)
                    if (int(out_data_B(7,1,nc)) == 0) then  !neutral (non-ion) component
                        act_basis = '(x)'
                    else
                        act_basis = '(m)'
                    endif
                    legend_text = 'log$_{10}$[$\gamma^{\,\beta , '//act_basis//'}$('//trim(name_species_TeX(nc))//')]'
                    rgb_set = int(col_palette(1:3, min(nc,ncolors)))
                    call add_plot_xydata(xv = out_data_B(1,iloc1:iloc2,1), yv = log10(out_data_B(5,iloc1:iloc2,nc)), &
                        & ltext=legend_text, pen_wid=6.0_wp, rgb_col=rgb_set, &
                        & lstyle='solid', plot_symb='curve')
                enddo
            end select

            !set overall plot properties and generate plot:
            xlabel = 'input mass fraction of water, $w_w^t$'
            ylabel = 'log$_{10}$(act. coeff.)'
            !set/limit x-and y-axis ranges for plot:
            xax_lim = [0.0_wp, 1.0_wp]
            yax_lim = [-2.0_wp, 2.0_wp]
            call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, xaxis_limits=xax_lim, yaxis_limits=yax_lim, &
                & legend_position=3, metafile='pdf', out_file_name=trim(fname)//'_actcoeff_'//char_phase ) 
        enddo !nphases
        !------------------------------------------------------------
        
        !(e) Plot activities of species, separately per phase:
        nspecies = size(name_species_TeX)
        do nph = 1,nphases
            select case(nph)
            case(1) !phase alpha
                char_phase = 'A'
                do nc = 1,nspecies
                    if (int(out_data_A(7,1,nc)) == 0) then  !neutral (non-ion) component
                        act_basis = '(x)'
                    else
                        act_basis = '(m)'
                    endif
                    legend_text = '$a^{\,\alpha ,'//act_basis//'}$('//trim(name_species_TeX(nc))//')'
                    rgb_set = int(col_palette(1:3, min(nc,ncolors)))
                    call add_plot_xydata(xv = out_data_A(1,:,1), yv = out_data_A(6,:,nc), &
                        & ltext=legend_text, pen_wid=6.0_wp, rgb_col=rgb_set, &
                        & lstyle='solid', plot_symb='curve')
                enddo
                !..
            case(2) !phase beta
                char_phase = 'B'
                do nc = 1,size(name_species_TeX)
                    if (int(out_data_B(7,1,nc)) == 0) then  !neutral (non-ion) component
                        act_basis = '(x)'
                    else
                        act_basis = '(m)'
                    endif
                    legend_text = '$a^{\,\beta ,'//act_basis//'}$('//trim(name_species_TeX(nc))//')'
                    rgb_set = int(col_palette(1:3, min(nc,ncolors)))
                    call add_plot_xydata(xv = out_data_B(1,iloc1:iloc2,1), yv = out_data_B(6,iloc1:iloc2,nc), &
                        & ltext=legend_text, pen_wid=6.0_wp, rgb_col=rgb_set, &
                        & lstyle='solid', plot_symb='curve')
                enddo
            end select

            !set overall plot properties and generate plot:
            xlabel = 'input mass fraction of water, $w_w^t$'
            ylabel = 'activitiy'
            !set/limit x-and y-axis ranges for plot:
            xax_lim = [0.0_wp, 1.0_wp]
            yax_lim = [0.0_wp, 5.0_wp]
            call dislin_plot(xlabel, ylabel, yaxis_mod=0.5_wp, xaxis_limits=xax_lim, yaxis_limits=yax_lim, &
                & legend_position=3, metafile='pdf', out_file_name=trim(fname)//'_activities_'//char_phase ) 
        enddo !nphases
        !------------------------------------------------------------
        
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
        !------------------------------------------------------------
        
        block !for chapter 3, exercise 4b solution;
            integer :: iH, iHSO4, iSO4
            real(wp),dimension(size(out_data_A(4,:,1))) :: alpha_HSO4, mmax_HSO4
            
            !Compute the degree of HSO4- dissociation
            !first determine the array index of relevant ions:
            iHSO4 = findloc(int(out_data_A(7,1,:)), VALUE=248, DIM=1)
            iSO4 = findloc(int(out_data_A(7,1,:)), VALUE=261, DIM=1)
            iH = findloc(int(out_data_A(7,1,:)), VALUE=205, DIM=1)
            
            if (iHSO4 > 0) then !bisulfate is present
                mmax_HSO4(:) = out_data_A(4,:,iHSO4) + &
                    & min(out_data_A(4,:,iSO4), out_data_A(4,:,iH))
                alpha_HSO4(:) = 1.0_wp - out_data_A(4,:,iHSO4)/mmax_HSO4(:)
                
                call add_plot_xydata(xv=out_data_A(6,:,1), yv=alpha_HSO4, &
                    & ltext='${\rm HSO_4^-}$ dissoc. deg., phase $\alpha$', &
                    & pen_wid=8.0_wp, lstyle='solid', plot_symb='curve')

                if (LLPSpresent) then
                    mmax_HSO4(iloc1:iloc2) = out_data_B(4,iloc1:iloc2,iHSO4) &
                        & + min(out_data_B(4,iloc1:iloc2,iSO4), out_data_B(4,iloc1:iloc2,iH))
                    alpha_HSO4(iloc1:iloc2) = 1.0_wp - out_data_B(4,iloc1:iloc2,iHSO4)/mmax_HSO4(iloc1:iloc2)
                    call add_plot_xydata(xv=out_data_B(6,iloc1:iloc2,1), &
                        & yv=alpha_HSO4(iloc1:iloc2) , &
                        & ltext='${\rm HSO_4^-}$ dissoc. deg., phase $\beta', &
                        & pen_wid=8.0_wp, lstyle='dashed_medium', &
                        & plot_symb='curve')
                endif

                !set overall plot properties and generate plot:
                xlabel = 'water activity, $a_w$'
                ylabel = '${\rm HSO_4^-}$ dissociation degree'
                xax_lim = [0.65_wp, 1.0_wp]
                yax_lim = [0.0_wp, 1.0_wp]
                call dislin_plot(xlabel, ylabel, yaxis_mod=0.4_wp, &
                    & xaxis_limits=xax_lim, yaxis_limits=yax_lim, legend_position=3, &
                    & metafile='pdf', out_file_name=trim(fname)//'_HSO4_dissoc')
            endif
        end block
        !------------------------------------------------------------

    end subroutine OutputLLE_plots
    !========================================================================================

end module Mod_InputOutput