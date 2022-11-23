module SchottkyGroup_Module

! Real hyperelliptic Schottky group construction: 
! 1. Fix $g,k\in\N$ such that 1<=k<=g+1. Let $\{1,2,...,g\} = $\Sigma^+\sqcup\Sigma^-$ with $|\Sigma^+| = k-1$.
! 2. Draw $g$ circles $C_1,\dots,C_g$ with disjoint diameters on positive half-line $(0;+\infty)$.
! 3. On each circle $C_j$ mark a pair of points: 
!       if $j\in\Sigma^+$ then mark intersection points $c_j\pm r_j$ of $C_j$ with real line,
!       if $j\in\Sigma^-$ then mark any conjugate pair $c_j \pm i r_j \in C_j$.
! 4. Let $\sigma_j := \pm 1$ for $j \in \Sigma^\pm$ and $S_j(u) := c_j - \sigma_j r_j^2/(u+c_j)$ for 1<=j<=g. 
!    Hyperbolic transformations $S_1,\dots,S_g$ freely generate Schottky group $\mathfrak{S}$. 
!
! Theorem: Orbit manifold of the Schottky group $\mathfrak{S}$ is a real hyperelliptic curve of 
!          genus $g$ with $k$ real and $k$ coreal ovals and all such curves can be obtained in this way.      
! 
! Let $S_{-j} := S_j^{-1}$, $C_{-j} := -C_j$, $c_{-j} := -c_j$, $r_{-j} = r_j$, $\sigma_{-j} = \sigma_j$ for 1<=j<=g.

use Precision_Module, only: precision
implicit none
private 

type, public :: SchottkyGroup_Type                
    integer :: g                                 ! Genus of curve associated with Schottky group.
    real(precision), allocatable :: c(:), r(:)   ! (-g:g)  Schottky moduli. c(-j) = -c(j), r(-j) = r(j), 1<=j<=g.
    integer, allocatable :: sigma(:)             ! (-g:g)  sigma(j) = sigma(-j) := $\pm1$, $j \in \Sigma^\pm$. 
    real(precision), allocatable :: rr(:)        ! (-g:g)  rr(j) = rr(-j) := $\sigma_j r_j^2$, 1<=j<=g.
contains 
    procedure :: S_Real, S_Complex           
    generic   :: S => S_Real, S_Complex          ! Generators and their inverses overloaded for real/complex inputs.
    procedure :: IsInsideModuliSpace             ! Is moduli+shift inside moduli space?
end type SchottkyGroup_Type

contains

real(precision) elemental function S_Real(self, j, u)
! Generators of Schottky group and their inverses: S(j,u) := $S_j(u)$, S(-j,u) := $S_j^{-1}(u)$, 1<=j<=g.
    class(SchottkyGroup_Type), intent(in) :: self
    real(precision), intent(in) :: u
    integer, intent(in) :: j
    S_Real = self%c(j) - self%rr(j)/(u + self%c(j))     
end function S_Real


complex(precision) elemental function S_Complex(self, j, u)
! Generators of Schottky group and their inverses: S(j,u) := $S_j(u)$, S(-j,u) := $S_j^{-1}(u)$, 1<=j<=g.
    class(SchottkyGroup_Type), intent(in) :: self
    complex(precision), intent(in) :: u
    integer, intent(in) :: j
    S_Complex = self%c(j) - self%rr(j)/(u + self%c(j))      
end function S_Complex


logical pure function IsInsideModuliSpace(self, shift)  
! Check if shifted moduli define a point inside moduli space. Shift is optional, zero by default.
! Shift formula: c(1:g) --> c(1:g) + shift(1:g), r(1:g) --> r(1:g) + shift(g+1:2*g). 
    class(SchottkyGroup_Type), intent(in) :: self
    real(precision), optional, intent(in) :: shift(2*self%g)
    real(precision) :: curr, rightPrev, c_new(self%g), r_new(self%g)
    integer :: j
    
    associate(g => self%g, c => self%c, r => self%r, sigma => self%sigma)
        c_new = c(1:g)
        r_new = r(1:g)
        if (present(shift)) then
            c_new = c_new + shift(1:g)
            r_new = r_new + shift(g+1:2*g)
        end if

        IsInsideModuliSpace = .false.
        if (any(r_new <= 0)) RETURN
        rightPrev = 0
        do j = 1,g
            if (sigma(j) == 1) then
                curr = c_new(j) - r_new(j)
                if (curr <= rightPrev) RETURN
                rightPrev = c_new(j) + r_new(j)
            else
                curr = c_new(j) + sigma(j)*r_new(j)**2/(rightPrev - c_new(j)) 
                if (curr <= rightPrev) RETURN
                rightPrev = curr
            end if
        end do
        IsInsideModuliSpace = .true.   ! Here if all tests passed.
    end associate
end function IsInsideModuliSpace

end module SchottkyGroup_Module