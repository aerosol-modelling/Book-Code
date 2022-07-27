!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main program for launching the example solutions to the BAT-related exercises from *
!*   Chapter 3 of the book.                                                             *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
program Prog_BAT_exercises

implicit none
!local variables: 
integer :: exercise_no
external :: BAT_exercise_1, BAT_exercise_2, BAT_exercise_3
!...............................................

!### parameters: set example number to be run ##
exercise_no = 1
!###

select case(exercise_no)
case(1)
    call BAT_exercise_1()
case(2)
    call BAT_exercise_2()
case(3)
    call BAT_exercise_3()
case default
    write(*,*) 'No case exists for the provided input.'
end select

write(*,*) ''
write(*,*) 'Done; press Enter to end the program Prog_BAT_Examples'
read(*,*)   !pause to wait for user action.

end program Prog_BAT_exercises