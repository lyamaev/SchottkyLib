module SchottkyGroup_Module

! Real hyperelliptic Schottky group construction: 
! 1. Fix $g,k\in\N$ such that 1<=k<=g+1. Let $\{1,2,...,g\} = $\Sigma^+\sqcup\Sigma^-$ 
!    with $|\Sigma^+| = k-1$.
! 2. Draw $g$ circles $C_1,\dots,C_g$ with disjoint diameters on positive half-line $(0;+\infty)$.
! 3. On each circle $C_j$ mark a pair of points: 
!       if $j\in\Sigma^+$ then mark intersection points $c_j\pm r_j$ of $C_j$ with real line,
!       if $j\in\Sigma^-$ then mark any conjugate pair $c_j \pm i r_j \in C_j$.
! 4. Let $\sigma_j := \pm 1$ for $j \in \Sigma^\pm$ and $S_j(u) := c_j - \sigma_j r_j^2/(u+c_j)$ 
!    for 1<=j<=g. Hyperbolic transformations $S_1,\dots,S_g$ freely generate Schottky group 
!    $\mathfrak{S}$. 
!
! Theorem: Orbit manifold of the Schottky group $\mathfrak{S}$ is a real hyperelliptic curve of 
!    genus $g$ with $k$ real and $k$ coreal ovals and all such curves can be obtained in this way.      
! 
! Let $S_{-j} := S_j^{-1}$, $C_{-j} := -C_j$, $c_{-j} := -c_j$, $r_{-j} = r_j$, 
!    $\sigma_{-j} = \sigma_j$ for 1<=j<=g.

use Precision_Module, only: precision
implicit none
private 
public SchottkyGroup_Type, SchottkyGroupFromFile


type SchottkyGroup_Type(g)
    integer, len :: g                     ! Genus of curve associated with Schottky group.
    real(precision) :: c(-g:g), r(-g:g)   ! Schottky moduli. c(-l) = -c(l), r(-l) = r(l), 1<=l<=g.
    integer :: sigma(-g:g)                ! sigma(l) = sigma(-l) := $\pm1$, $l \in \Sigma^\pm$. 
    real(precision) :: rr(-g:g)           ! rr(l) = rr(-l) := $\sigma_l r_l^2$, 1<=l<=g.
contains 
    procedure :: Init => InitSchottkyGroup   ! Set new moduli and initialize all fields of type.
    procedure :: IsInsideModuliSpace         ! Is moduli+shift inside moduli space?
    procedure :: S_real, S_complex           ! Overloading for real/complex inputs.
    generic   :: S => S_real, S_complex      ! Generators $S_j$, 1<=j<=g, and their inverses. 
end type SchottkyGroup_Type

contains


subroutine InitSchottkyGroup(self, c_new, r_new, sigma_new)   
! Set new moduli and initialize all fields of type.
    class(SchottkyGroup_Type(*)), intent(inout) :: self
    real(precision), intent(in) :: c_new(self%g), r_new(self%g)
    integer, intent(in) :: sigma_new(self%g)
    associate(g => self%g, c => self%c, r => self%r, sigma => self%sigma, rr => self%rr)
        c(1:g) = c_new(1:g)
        r(1:g) = r_new(1:g)
        sigma(1:g) = sigma_new(1:g)
    
        if (.not. self%IsInsideModuliSpace()) then
            error stop 'SchottkyGroup_Module/InitSchottkyGroup: outside of moduli space.'
        end if

        c(-g:-1) = -c(g:1:-1)
        r(-g:-1) = r(g:1:-1)
        sigma(-g:-1) = sigma(g:1:-1)
        rr(-g:g) = sigma(-g:g) * r(-g:g)**2
    end associate
end subroutine InitSchottkyGroup


function SchottkyGroupFromFile(pathToInputModuli) result(schottkyGroup)
! Read Schottky moduli from file and return fully initialized Schottky group object.
! Input file must contain g+1 lines and be organized by the following template:
!       g
!       c_1  r_1
!       c_2  r_2
!       ...
!       c_g  r_g
    character(*), intent(in) :: pathToInputModuli   ! String with full path to input file.
    class(SchottkyGroup_Type(:)), allocatable :: schottkyGroup
    integer :: j, g, inputModuli_File
    real(precision), allocatable :: c(:), r(:) 
    integer, allocatable :: sigma(:)

    open(file = pathToInputModuli, newunit = inputModuli_File, action = 'read')
    read(inputModuli_file, *) g
    allocate(SchottkyGroup_Type(g) :: schottkyGroup)
    allocate(c(g), r(g), sigma(g))
    do j = 1, g
        read(inputModuli_file, *) c(j), r(j), sigma(j)
    end do
    call schottkyGroup%Init(c, r, sigma)
    close(inputModuli_File)
end function SchottkyGroupFromFile


real(precision) elemental function S_real(self, j, u)
! Generators of Schottky group and their inverses: 
! S(j,u) := $S_j(u)$, S(-j,u) := $S_j^{-1}(u)$ for 1<=l<=g.
    class(SchottkyGroup_Type(*)), intent(in) :: self
    real(precision), intent(in) :: u
    integer, intent(in) :: j
    S_real = self%c(j) - self%rr(j)/(u + self%c(j))     
end function S_real


complex(precision) elemental function S_complex(self, j, u)
! Generators of Schottky group and their inverses: 
! S(j,u) := $S_j(u)$, S(-j,u) := $S_j^{-1}(u)$ for 1<=l<=g.
    class(SchottkyGroup_Type(*)), intent(in) :: self
    complex(precision), intent(in) :: u
    integer, intent(in) :: j
    S_complex = self%c(j) - self%rr(j)/(u + self%c(j))      
end function S_complex


logical pure function IsInsideModuliSpace(self, shift)  
! Check if shifted moduli define a point inside moduli space. 
! Shift is optional, zero by default.
! Shift: c(1:g) --> c(1:g) + shift(1:g), r(1:g) --> r(1:g) + shift(g+1:2*g). 
    class(SchottkyGroup_Type(*)), intent(in) :: self
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