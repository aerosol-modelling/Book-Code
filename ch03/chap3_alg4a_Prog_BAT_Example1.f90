program Prog_BAT_Example_1

implicit none
external :: BAT_example_1
!.........................................
call BAT_example_1()
write(*,*) 'Done; press Enter to end the program'
read(*,*)    !pause to wait for user action.

end program Prog_BAT_Example_1