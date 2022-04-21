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