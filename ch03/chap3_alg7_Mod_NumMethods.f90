!****************************************************************************************
!*   :: Purpose ::                                                                      *
!*   Module providing numerical equation solver and optimization methods.               *
!*                                                                                      *
!*   :: Authors ::                                                                      *
!*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
!*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
!*   For original authors of methods, see information at top of specific procedures.    *
!*                                                                                      *
!*   :: List of subroutines and functions contained in this module:                     *
!*   --------------------------------------------------------------                     *
!*   -  function    local_min                                                           *
!*   -  function    linear_interpolation                                                *
!*   -  function    Ridders_zero                                                         *
!*                                                                                      *
!****************************************************************************************
module Mod_NumMethods

use Mod_NumPrec_Types, only : wp, foptions

implicit none
private

!public procedures:
public :: linear_interpolation, local_min, Ridders_zero

contains


    !*****************************************************************************
    !  ** function local_min(a, b, tolinit, fnc, fopt, x_min, calls) **
    !
    !   local_min() seeks a local minimum of a function fnc(x, fopt) in an interval [a, b].
    !
    !  Discussion:
    !
    !    If the function fnc is defined on the interval [a, b], then local_min
    !    finds an approximation x_min to the point at which fnc(x, fopt) attains its minimum
    !    (or the appropriate limit point), and returns the value of fnc at x_min.
    !    fopt is a customized derived type providing function parameters (if needed);
    !
    !    tolinit and eps define a tolerance tol = eps*abs(x) + tolinit,
    !    where eps is here 5 times the rel. machine epsilon (for set real kind).
    !    fnc is never evaluated at two points closer than tol.
    !
    !    If fnc is delta-unimodal for some delta less than tol, the x_min approximates
    !    the global minimum of fnc with an error less than 3*tol.
    !
    !    If fnc is not delta-unimodal, then x_min may approximate a local, but
    !    perhaps non-global, minimum.
    !
    !    The method used is a combination of golden section search and
    !    successive (inverse) parabolic interpolation.  Convergence is never much
    !    slower than that for a Fibonacci search.  If 'f' has a continuous second
    !    derivative which is positive at the minimum (which is not at limit a or
    !    b), then, ignoring rounding errors, convergence is superlinear, and
    !    usually of the order of about 1.3247.
    !
    !    Thanks to Jonathan Eggleston for pointing out a correction to the
    !    golden section step, 01 July 2013.
    !
    !  Licensing:
    !
    !    This code is distributed under the GNU LGPL license.
    !
    !  Modified:
    !
    !    30 May 2021 (by John Burkardt)
    !
    !    26 June 2021: a few modifications by Andi Zuend, including:
    !           use of Mod_NumPrec_Types for real kinds (wp) instead of
    !           real( kind = 8 ) and use of literal floating point constants of
    !           wp kind. Renaming of some variables;
    !           eps has been made a local parameter rather than an input;
    !           addition of fopt to pass additional function parameters to fnc.
    !
    !
    !  Author:
    !
    !    Original FORTRAN 77 version by Richard Brent.
    !    Fortran 90 version by John Burkardt.
    !
    !  Reference:
    !
    !    Richard Brent,
    !    Algorithms for Minimization Without Derivatives,
    !    Dover, 2002,
    !    ISBN: 0-486-41998-3,
    !    LC: QA402.5.B74.
    !
    !*****************************************************************************

    function local_min(a, b, tolinit, fnc, fopt, x_min, calls)

    implicit none

    !interface; input/output arguments:
    real(kind=wp),intent(in) :: a, b            !input: the endpoints of the interval.
    real(kind=wp),intent(in) :: tolinit         !input: a positive absolute error tolerance; typical value is sqrt(epsilon(1.0_wp))
    interface
        function fnc(x, fopt)                   !fnc, the name of a user-supplied function to be evalulated
        use Mod_NumPrec_Types, only : wp, foptions
        implicit none
        real(wp),intent(in) :: x
        type(foptions),optional,intent(in) :: fopt
        real(wp) :: fnc
        end function fnc
    end interface
    type(foptions),optional,intent(in) :: fopt  !input: function option arguments to be passed to fnc.
    real(kind=wp),intent(inout) :: x_min        !output: the estimated value of the abscissa of the minimum in fnc within [a, b]
    integer,intent(out) :: calls                !output: the number of calls of function fnc.
    real(kind=wp) :: local_min                  !output: the value of fnc(x=x_min, fopt).
    !
    !   eps is a positive relative error tolerance. eps should be no smaller than twice the relative machine precision,
    !   and preferably not much less than the square root of the relative machine precision.
    real(kind=wp),parameter :: eps = 5.0_wp*epsilon(1.0_wp)
    !local variables:
    real(kind=wp) :: c, d, e, fu, fv, fw, fx, m, p, q, r, sa, sb, t2, tol, u, v, w
    !.............................................

    calls = 0
    !  c is the square of the inverse of the golden ratio:
    c = 0.5_wp * (3.0_wp - sqrt(5.0_wp))

    sa = a
    sb = b
    x_min = sa + c * (b - a)
    w = x_min
    v = w
    e = 0.0_wp
    fx = fnc(x_min, fopt)
    calls = calls + 1
    fw = fx
    fv = fw

    do  !loop until "exit"
        m = 0.5_wp *(sa + sb)
        tol = eps *abs(x_min) + tolinit
        t2 = 2.0_wp * tol
        !
        !  Check the stopping criterion.
        !
        if ( abs(x_min - m) <= t2 - 0.5_wp * (sb - sa) ) then
            exit
        endif
        !
        !  Fit a parabola.
        !
        r = 0.0_wp
        q = r
        p = q

        if ( tol < abs(e) ) then
            r = (x_min - w) * (fx - fv)
            q = (x_min - v) * (fx - fw)
            p = (x_min - v) * q - (x_min - w) * r
            q = 2.0_wp * (q - r)
            if ( 0.0_wp < q) then
                p = - p
            endif
            q = abs(q)
            r = e
            e = d
        endif

        if ( abs(p) < abs(0.5_wp *q*r) .AND. &
            & q*(sa - x_min) < p .AND. p < q*(sb - x_min) ) then
            !
            !  Take the inverse parabolic interpolation step:
            !
            d = p / q
            u = x_min + d
            !
            !  fnc must not be evaluated too close to a or b.
            !
            if ( (u - sa) < t2 .OR. (sb - u) < t2 ) then
                if (x_min < m) then
                    d = tol
                else
                    d = - tol
                endif
            endif
            !
            !  A golden section step.
            !
        else
            if (x_min < m) then
                e = sb - x_min
            else
                e = sa - x_min
            endif
            d = c * e
        endif
        !
        !  fnc must not be evaluated too close to x_min.
        !
        if ( tol <= abs(d) ) then
            u = x_min + d
        elseif ( 0.0_wp < d ) then
            u = x_min + tol
        else
            u = x_min - tol
        endif

        fu = fnc(u, fopt)
        calls = calls + 1
        !
        !  Update a, b, v, w, and x_min.
        !
        if ( fu <= fx ) then
            if ( u < x_min ) then
                sb = x_min
            else
                sa = x_min
            endif
            v = w
            fv = fw
            w = x_min
            fw = fx
            x_min = u
            fx = fu
        else
            if ( u < x_min ) then
                sa = u
            else
                sb = u
            endif
            if ( fu <= fw .OR. w == x_min ) then
                v = w
                fv = fw
                w = u
                fw = fu
            elseif ( fu <= fv .OR. v == x_min .OR. v == w ) then
                v = u
                fv = fu
            endif
        endif
    enddo

    local_min = fx

    end function local_min
    !-----------------------------------------------------------


    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Function for computing the linear interpolation value of an x-axis value at a      *
    !*   targeted y-axis value based on only two given (x1, y1), (x2, y2) points that       *
    !*   enclose the interval of interest (assumed to be fulfilled at input).               *
    !*                                                                                      *
    !*   :: Authors ::                                                                      *
    !*   Andi Zuend (andreas.zuend@mcgill.ca)                                               *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !****************************************************************************************
    pure function linear_interpolation(xval, yval, target_yval)  result(intp_xval)

    implicit none

    real(wp),dimension(2),intent(in) :: xval, yval
    real(wp),intent(in) :: target_yval
    real(wp) :: intp_xval, y1m2
    !...............................................

    y1m2 = yval(1) - yval(2)
    if (abs(y1m2) > 0.0_wp) then
        intp_xval = ( xval(1)*(yval(1) - target_yval) + xval(2)*(target_yval - yval(2)) ) / y1m2
    else
        intp_xval = sqrt(xval(1)*xval(2))           !use geometric mean in this exceptional case
    endif

    end function linear_interpolation
    !-----------------------------------------------------------


    !****************************************************************************************
    !*   :: Purpose ::                                                                      *
    !*   Function to locate the root (zero) bounded within a given interval [x1, x2]        *
    !*   where the function changes sign -- by means of Ridders' method as outlined in      *
    !*   Numerical Recipes for Fortran 90 (2nd edition by Press et al., 1996).              *
    !*   The original Fortran 90 code from Numerical Recipes has been slightly adapted      *
    !*   here (e.g., using real(wp) instead of single precision and using alternative exit  *
    !*   points). Also, additional function parameters can be passed via fopt.              *
    !*                                                                                      *
    !*   Original short description: "Using Ridders’ method, return the root of a function  *
    !*   fnc known to lie between x1 and x2. The root, returned as Ridder_zero, will be     *
    !*   refined to an approximate relative accuracy set by xacc".                          *
    !*                                                                                      *
    !*   :: Author & Copyright ::                                                           *
    !*   Numerical Recipes in Fortran 90 (Press et al., 1996)                               *
    !*   Some modifications by Andi Zuend,                                                  *
    !*   Dept. Atmospheric and Oceanic Sciences, McGill University                          *
    !*                                                                                      *
    !*   -> created:        1996 (?)                                                        *
    !*   -> latest changes: 2021/06/30                                                      *
    !*                                                                                      *
    !****************************************************************************************
    function Ridders_zero(fnc, fopt, x1, x2, f1, f2, xacc, calls)

    implicit none
    !interface arguments:
    interface
        function fnc(x, fopt)                           !input: function 'fnc' to be evaluated.
        use Mod_NumPrec_Types, only : wp, foptions
        implicit none
        real(wp),intent(in) :: x
        type(foptions),optional,intent(in) :: fopt
        real(wp) :: fnc
        end function fnc
    end interface
    type(foptions),optional,intent(in) :: fopt          !input: (optional) function option arguments to be passed to fnc
    real(wp),intent(in) :: x1, x2, f1, f2               !input: interval limits [x1, x2] and corresponding fnc values (f1, f2)
    real(wp),intent(in) :: xacc                         !input: relative accuracy of x-value determination
    integer,intent(out) :: calls                        !output: the number of calls of function fnc.
    real(wp) :: Ridders_zero                            !output: function return value (estimate of x for which fnc is ~ zero)
    !local variables:
    integer,parameter :: maxit = 50
    real(wp),parameter :: xdev = 0.1_wp*sqrt(epsilon(1.0_wp))
    integer :: j
    real(wp) :: fl, fh, fm, fnew, s, xh, xl, xm, xnew
    !...............................................

    calls = 0
    fl = f1
    fh = f2

    if ((fl > 0.0_wp .and. fh < 0.0_wp) .or. (fl < 0.0_wp .and. fh > 0.0_wp)) then
        xl = x1
        xh = x2
        Ridders_zero = 0.5_wp*(xl + xh)             !initial guess
        do j = 1,maxit
            xm = 0.5_wp*(xl + xh)
            fm = fnc(xm, fopt)
            calls = calls +1
            s = sqrt(fm**2 - fl*fh)
            if (s == 0.0_wp) return
            xnew = xm + (xm - xl)*(sign(1.0_wp, fl - fh)*fm/s)
            if (abs(xnew - Ridders_zero)/sqrt(abs(xm)) < xacc) then
                Ridders_zero = xnew
                exit
            endif
            Ridders_zero = xnew
            fnew = fnc(Ridders_zero, fopt)
            calls = calls +1
            if (sign(fm, fnew) /= fm) then
                xl = xm
                fl = fm
                xh = Ridders_zero
                fh = fnew
            else if (sign(fl, fnew) /= fl) then
                xh = Ridders_zero
                fh = fnew
            else if (sign(fh, fnew) /= fh) then
                xl = Ridders_zero
                fl = fnew
            else
                !$OMP critical
                write(*,*) 'WARNING from Ridder_zero: never get here'
                !read(*,*)
                !$OMP end critical
            endif
            if (abs(xh - xl) < xdev) then   !alternative exit point for absolute bracket width smaller than a tolerance value
                exit
            endif
        enddo
        if (j > maxit) then
            !$OMP critical
            write(*,*) 'WARNING from Ridder_zero: exceeded maximum number of iterations'
            !read(*,*)
            !$OMP end critical
        endif
    else if (fl == 0.0_wp) then
        Ridders_zero = x1
    else if (fh == 0.0_wp) then
        Ridders_zero = x2
    else
        !$OMP critical
        write(*,*) 'WARNING from Ridder_zero: root must be bracketed'
        read(*,*)
        !$OMP end critical
        Ridders_zero = x1
    endif

    end function Ridders_zero
    !-----------------------------------------------------------


end module Mod_NumMethods