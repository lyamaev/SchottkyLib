module Uniformization_Module
! Poincare algorithm for numerical uniformization of hyperelliptic M-curves. Implementation is based on the 
! article M. Seppala, "Myrberg's numerical uniformization of hyperelliptic curves". 
!
! Numerical uniformization is converting algebraic model of a given curve into Schottky model. Namely, given 
! a set of ramification points $0 < x_1^- < x_1^+ < x_2^- < x_2^+ < ... < x_g^- < x_g^+ < \infty$, one needs 
! to compute Schottky moduli c(1:g), r(1:g), sigma(1:g) = 1, which correspond to hyperelliptic M-curve with 
! affine part $w^2 = x\prod_{j=1}^g(x-x_j^-)(x-x_j^+)$, $(x,w)\in\C^2$.

! Note that there is alternative approach to numerical Schottky uniformization: solve system of equations 
! [Chi(*ramification pts in Schottky model*)](*moduli*) = *ramification pts in algebraic model*, using 
! variational formulae for Abelian integrals and their periods. 

use Precision_Module, only: precision
use DiffsAndFuncsOnCurve_Module, only: SchottkyNumerics_Type
implicit none
public PoincareUniformization
private

integer, parameter :: slN = 100000   ! Maximum number of slots.
real(precision), parameter :: uEps = Epsilon(1.0_precision)

contains

function PoincareUniformization(ramInAlgM, slG) result(schottkyGroup)   ! Opening process with 'slG' gens of slots. 
! ramInAlgM := [x_1^-, x_1^+, x_2^-, x_2^+, ..., x_g^-, x_g^+] -- flat array of ramification points in algebraic 
! model sotred in ascending order.
    real(precision), intent(in) :: ramInAlgM(:)
    integer, intent(in) :: slG                    ! Number of slot generations we want to process.
    real(precision) :: p(size(ramInAlgM)/2+1,2)   ! Current approximation of intersection points $c_j \pm r_j$.
    real(precision) :: q(2)                       ! q := [c_0,r_0] in thesis notations.
    real(precision) :: slot(slG,slN,2)   ! Array of slots: sl(a,b,c) -- a-th end of b-th slot of c-th generation.
    type(SchottkyNumerics_Type) :: schottkyGroup    
    integer :: n(slG), j, k, i, g

    g = size(ramInAlgM)/2   ! Genus of curve.

    forall(j = 1:g)
        p(j+1,1) = ramInAlgM(2*j-1) 
        p(j+1,2) = ramInAlgM(2*j)
    end forall
    p = p / (2*p(g+1,2) - p)   ! Normalize: 0->0, $x_g^+$->1, $\infty$->-1.
    p(1,1) = -1
    p(1,2) = 0

    n(1) = g+1
    n(2:slG) = 0

    slot(1,1:n(1),:) = p(1:n(1),:)

    ! Notations:
    ! k    -- (number of) current generation.
    ! i    -- (number of) current slot inside k-th generation.
    ! n(a) -- number of slots of a-th generation.
    ! slG  -- max number of generations we want to process.
    ! j    -- index to access generations with numbers > k.
    OpeningProcess: do k = 1,slG
        do i = 1,n(k)
            if (abs(slot(k,i,1)-slot(k,i,2)) <= uEps .or. maxval(n(k+1:slG)) > slN) EXIT OpeningProcess
            q = [(slot(k,i,1) + slot(k,i,2))/2, (slot(k,i,2) - slot(k,i,1))/4]
            call ApplyPhi(p)   ! Update approximation of intersection points.

            ! Apply opening transform Phi(*,q) to closed slots. 
            call ApplyPhi(slot(k,i+1:n(k),:))
            do j = k+1,slG 
                call ApplyPhi(slot(j,1:n(j),:))
            end do

            ! Add new slots -- images of closed slots by involution G.
            do j = slG,k+2,-1
                slot(j,n(j)+1:n(j)+n(j-1),:) = Invl(slot(j-1,1:n(j-1),:))
                n(j) = n(j) + n(j-1)
            end do    
            if (k < slG) then
                slot(k+1,n(k+1)+1:n(k+1)+n(k)-i,:) = Invl(slot(k,i+1:n(k),:))
                n(k+1) = n(k+1) + n(k) - i
            end if
        end do
    end do OpeningProcess

    ! Normalize: p(1,1)->\infty, p(1,2)->0, c(g)+r(g):=p(g+1,2)->1.
    p = (p - p(1,2))/(p - p(1,1)) / ((p(g+1,2) - p(1,2))/(p(g+1,2) - p(1,1)))
    
    ! Construct M-curve with class costructor. sigma=1 by default.
    schottkyGroup = SchottkyNumerics_Type(c=(p(2:g+1,2)+p(2:g+1,1))/2, r=(p(2:g+1,2)-p(2:g+1,1))/2)   

    contains

        subroutine ApplyPhi(z)   ! Apply opening slot transform $\Phi(z)$. 
            real(precision) :: z(:,:)   ! size(z,2)=2. 
            z = z + ((z-q(1))/2) * (sqrt(abs(1 - 4*(q(2)/(z-q(1)))**2)) - 1)
        end subroutine ApplyPhi

        function Invl(z)   ! Real elliptic involution $G(z)$.
            real(precision) :: z(:,:), Invl(size(z,1),2)   ! size(z,2)=2. 
            Invl = q(1) + q(2)**2/(z-q(1))
        end function Invl

end function PoincareUniformization

end module Uniformization_Module