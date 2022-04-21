!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Main program providing several examples of BAT model applications via subroutine   *
!*   calls (see settings below).                                                        *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*                                                                                      *
!**************************************************************************************** 
program Prog_BAT_Examples

implicit none
!local variables: 
integer :: example_no
external :: BAT_example_1, BAT_example_2
!...............................................

!### parameters: set example number to be run ##
example_no = 2
!###

select case(example_no)
case(1)
    call BAT_example_1()
case(2)
    call BAT_example_2()
case default
    write(*,*) 'No example program exists for the provided input.'
end select

write(*,*) ''
write(*,*) 'Done; press Enter to end the program Prog_BAT_Examples'
read(*,*)   !pause to wait for user action.

end program Prog_BAT_Examples