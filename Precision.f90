module Precision_Module
implicit none
private
public precision

! 6, 15 or 33 decimal digits.
integer, parameter :: singlePrecision    = selected_real_kind(6)
integer, parameter :: doublePrecision    = selected_real_kind(15)
integer, parameter :: quadruplePrecision = selected_real_kind(33)

! Global precision of real and complex numbers.
integer, parameter :: precision = doublePrecision  

real(precision), parameter :: pi = 3.14159265358979323846264338327950288_precision
complex(precision), parameter :: im = cmplx(0,1,precision)   ! Imaginary unit.

end module Precision_Module