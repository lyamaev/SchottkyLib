program RamificationPts
! Construct Schottky group associated with M-curve given in algebraic model and check that projections
! of ramification points in Schottky model coincide with ramification points in algebraic model.

    use Precision_Module, only: p => precision
    use DiffsAndFuncsOnCurve_Module, only: SchottkyNumerics_Type
    use Uniformization_Module, only: PoincareUniformization

    type(SchottkyNumerics_Type) :: schottkyGroup
    real(p), allocatable :: x(:), ramInSchM(:), ramInAlgM(:)

    ! Pick vector x(:) of even size. 
    x = [0.1_p, 0.2_p, 1.0_p, 1.3_p, 1.5_p, 2.9_p] 

    ! Construct Schottky group of genus size(x)/2 that corresponds to the M-curve with affine part given by the equation 
    ! w^2 = z(z-x(1))(z-x(2))...(z-x(size(x))), (z,w)\in\C. We use Poincare opening process with 5 generations of slots.
    schottkyGroup = PoincareUniformization(ramInAlgM = x, slG = 5)

    ! Ramification points in Schottky model.
    allocate(ramInSchM(size(x)))
    forall(j = 1:size(x)/2)
        ramInSchM(2*j-1) = schottkyGroup%c(j) - schottkyGroup%r(j)
        ramInSchM(2*j)   = schottkyGroup%c(j) + schottkyGroup%r(j)
    end forall

    ! Chi is the projection from Schottky model to algebraic model normalized by 0->0, infty->infty, b->1.
    ramInAlgM = schottkyGroup%Chi(cmplx(ramInSchM,0,p), b = 1.0_p)
    ramInAlgM = ramInAlgM * x(1)/ramInAlgM(1)   ! Stretch to the same scale as x(:).

    ! Check that ramification points in algebraic model for constructed Schottky group coincide with x(:).
    print *, 'Ramification points in algebraic model:', real(ramInAlgM, kind = 4)

! Output:
! Ramification points in algebraic model:  0.1000000   0.2000000   1.000000   1.300000   1.500000   2.900000    
end program RamificationPts