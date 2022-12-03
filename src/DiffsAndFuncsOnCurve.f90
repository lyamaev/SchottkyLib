module DiffsAndFuncsOnCurve_Module   ! Differentials and functions on a curve in Schottky model.

! This module provides 'SchottkyNumerics_Type' -- extension of 'SchottkyGroup_Type' with methods for computing
! a number of differentials and functions on a curve (all needed in Chebyshev ansatz). For almost all of them 
! reduced sum/product is used (reduction is possible due to symmetry $S^{-}(-u) = -S(u)$, $S\in\mathfrak{S}$).

! If some diff/func you need is missing in the current implementation, it's not difficult to add it using
! universal Cayley tree traversal interface provided by CayleyTreeTraveral_Module.

! Notations:
! $\eta_{zw}$   -- normalized linear Abelian differential of the 3rd kind. 
! $\zeta_j$     -- normalized linear holomorphic differential.
! $\theta_j$    -- even holomorphic quadratic differential.
! $\chi(u)$     -- projection from Schottky model to algebraic model.
! $b_j$         -- b-cycles of symplectic basis on a curve.
! $\zeta := (i\pi)^{-1}(zetaCoeff(1)*\zeta_1 + ... + zetaCoeff(g)*\zeta_g$).

! Define alpha(j) as follows:
! 1<=j<=g       --> alpha(j) := attractive fixed point (AFP) of $S_j$; 
! g+1<=j<=2*g-1 --> alpha(j) := attractive fixed point of $S_{j-g+1}S_1^{-1}S_{j-g+1}$ (these transforms are 
! used to form theta-basis of even holomorphic quadratic differentials). 

use Precision_Module, only: precision, pi
use SchottkyGroup_Module, only: SchottkyGroup_Type
use ChildrenOrder_Module, only: ChildrenOrder, Child_Type
use CayleyTreeTraversal_Module, only: Traversal_Interface, TraversalParameters_Type, & 
                                      Traversal_BogatyrevAlg, Traversal_NewAlg

implicit none 
private
public SchottkyNumerics_Type

procedure(Traversal_Interface), pointer :: Traversal   ! Module-scale traversal function.

type, extends(SchottkyGroup_Type) :: SchottkyNumerics_Type
    real(precision)  :: eps = 1000*Epsilon(1.0_precision)    ! Parameter to control accuracy of (a posteriori) algs.
    logical :: useNewAlg = .true.                            ! Use NewAlg (or BogatyrevAlg)?
    real(precision), allocatable :: alpha(:)                 ! (2*g-1)  AFP of $S_j$ and $S_jS_1^{-1}S_j$. 
    real(precision), allocatable :: lambda(:)                ! (g)  lambda(j) := $S_j'\alpha_j$ -- multiplier of $S_j$.
    type(Child_Type), allocatable :: children(:,:)           ! (-g:g,1:2*g-1)  childEps for linear differentials.
    type(Child_Type), allocatable :: children_sqrtEps(:,:)   ! (-g:g,1:2*g-1)  childEps for quadratic differentials.
contains
    procedure :: Eta            ! $\eta_{zw}(u)/du$ -- coeff of Abelian diff of 3rd kind.
    procedure :: EtaSymm        ! $(\eta_{zw}(u)+\eta_{-z-w}(u)/)du$.
    procedure :: EtaZeroInfty   ! $\eta_{0\infty}(u)/du$.
    procedure :: ZetaJ          ! $\zeta_j(u)/du$ -- coeff of Abelian diff of 1st kind.  
    procedure :: Zeta           ! $\zeta(u)/du$.
    procedure :: Theta          ! $\theta_j(u)/(du)^2$ -- coeff of even holomorphic quadratic diff.
    procedure :: ExpIntEta      ! $\exp\int_\infty^u\eta_{zw}$.
    procedure :: ExpIntZetaJ    ! $\exp\int_\infty^u\zeta_j$.
    procedure :: ExpPeriod      ! $\exp\int_{b_l}\zeta_j$. 
    procedure :: PeriodMatrix   ! $\{\int_{b_j}\zeta_l, 1<=l,j<=g\}$.
    procedure :: IntZeta        ! $\int_\infty^u\zeta$.
    procedure :: Chi            ! $\chi(u) := (\int_b^u\eta_{0\infty)^2}$. 
    procedure :: SetParams
    final :: Destruct
end type SchottkyNumerics_Type

interface SchottkyNumerics_Type
    procedure :: Constructor
end interface SchottkyNumerics_Type

contains

function Constructor(c, r, sigma, eps, useNewAlg) result(a)
    real(precision), intent(in) :: c(:), r(size(c))
    integer, intent(in), optional :: sigma(size(c))
    real(precision), intent(in), optional :: eps
    logical, intent(in), optional :: useNewAlg
    type(SchottkyNumerics_Type) :: a

    a%g = size(c)
    allocate(a%c(-a%g:a%g), a%r(-a%g:a%g), a%rr(-a%g:a%g), a%sigma(-a%g:a%g), a%alpha(2*a%g-1), &
             a%lambda(a%g), a%children(-a%g:a%g,1:2*a%g-1), a%children_sqrtEps(-a%g:a%g,1:2*a%g-1))
    call a%SetParams(c, r, sigma, eps, useNewAlg)
end function Constructor


impure elemental subroutine Destruct(self)
    type(SchottkyNumerics_Type), intent(inout) :: self
    if (allocated(self%c)) deallocate(self%c)
    if (allocated(self%r)) deallocate(self%r)
    if (allocated(self%rr)) deallocate(self%rr)
    if (allocated(self%sigma)) deallocate(self%sigma)
    if (allocated(self%alpha)) deallocate(self%alpha)
    if (allocated(self%lambda)) deallocate(self%lambda)
    if (allocated(self%children)) deallocate(self%children)
    if (allocated(self%children_sqrtEps)) deallocate(self%children_sqrtEps)
end subroutine Destruct


subroutine SetParams(self, c_new, r_new, sigma_new, eps_new, useNewAlg_new)   ! Set parameters of SchottkyNumerics.
    class(SchottkyNumerics_Type), intent(inout) :: self
    real(precision), intent(in), optional :: c_new(self%g), r_new(self%g)
    integer, intent(in), optional :: sigma_new(self%g)
    real(precision), intent(in), optional :: eps_new
    logical, intent(in), optional :: useNewAlg_new

    if (present(c_new) .and. present(r_new)) call self%SetNewModuli(c_new, r_new, sigma_new)
    if (present(eps_new)) self%eps = eps_new
    if (present(useNewAlg_new)) self%useNewAlg = useNewAlg_new

	if (self%useNewAlg) then
        Traversal => Traversal_NewAlg
        self%children = ChildrenOrder(self, self%eps)
        self%children_sqrtEps%childIndex = self%children%childIndex
        self%children_sqrtEps%childEps = sqrt(self%children%childEps)
    else 
        Traversal => Traversal_BogatyrevAlg
    end if

    associate(g => self%g, c => self%c, rr => self%rr, alpha => self%alpha, lambda => self%lambda)
        alpha(1:g) = sqrt(c(1:g)**2 - rr(1:g))
        alpha(g+1:2*g-1) = sqrt( ((alpha(2:g)**2-c(2:g)*c(1))**2 - rr(1)*c(2:g)**2) / &
                                 ((c(2:g)-c(1))**2 - rr(1)) )
        lambda(1:g) = rr(1:g)/(c(1:g) + alpha(1:g))**2 
    end associate
end subroutine SetParams


pure function Eta(self, u, z)   ! $\eta_{z_1z_2}(u)/du$.
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    real(precision), intent(in) :: z(2)
    type(TraversalParameters_Type) :: params
    complex(precision) :: Eta(size(u))

    ! Sum is over the whole Schottky group, so we allocate arrays leftCoset and rightCoset of zero size. 
    ! Note: ifort works fine without these dummy allocations, but gfortran occasionally throws segfault.
    allocate(params%leftCoset(0), params%rightCoset(0))

    params%Operation => Sum
    params%isReduced = .false.
    params%u = u;  params%z = z
    params%Term => Eta_Term;  params%idTransformTerm = Eta_Term(u,z)
    params%eps = self%eps;  params%children = self%children
    
    Eta = Traversal(self, params)

    contains 
        pure function Eta_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: Eta_Term(size(u))
            Eta_Term = (Sz(1)-Sz(2)) / ((u-Sz(1))*(u-Sz(2)))
        end function Eta_Term
end function Eta


pure function EtaSymm(self, u, z)   ! $ (\eta_{z_1z_2}+\eta_{-z_1-z_2})(u)/du$.
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    real(precision), intent(in) :: z(2)
    type(TraversalParameters_Type) :: params
    complex(precision) :: EtaSymm(size(u))

    params%Operation => Sum
    params%isReduced = .true.;  allocate(params%leftCoset(0), params%rightCoset(0))
    params%u = u**2   ! Compute u**2 straightaway and once.
    params%z = [z(1),z(2),-z(1),-z(2)]
    params%Term => EtaSymm_Term;  params%idTransformTerm = (z(1)**2-z(2)**2)/((u**2-z(1)**2)*(u**2-z(2)**2))
    params%eps = self%eps;  params%children = self%children

    EtaSymm = 2*u * Traversal(self, params)

    contains 
        pure function EtaSymm_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)    ! dim(Sz) = 4.
            complex(precision) :: EtaSymm_Term(size(u)), v(size(Sz)) 
            v = Sz**2
            EtaSymm_Term = 1/(u-v(1)) - 1/(u-v(2)) + 1/(u-v(3)) - 1/(u-v(4))
        end function EtaSymm_Term
end function EtaSymm


pure function EtaZeroInfty(self, u)   ! $\eta_{0\infty}(u)/du$.
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    real(precision) :: zero
    type(TraversalParameters_Type) :: params
    complex(precision) :: EtaZeroInfty(size(u))

    zero = 0
    params%u = u**2                 ! Compute u**2 straightaway and once.
    params%z = [zero, -log(zero)]   ! -log(zero) = infty.
    params%Operation => Sum
    params%isReduced = .true.;  allocate(params%leftCoset(0), params%rightCoset(0))
    params%Term => EtaZeroInfty_Term;  params%idTransformTerm = 1/(2*u**2)
    params%eps = self%eps;  params%children = self%children

    EtaZeroInfty = 2*u * Traversal(self, params)

    contains 
        pure function EtaZeroInfty_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: EtaZeroInfty_Term(size(u)), v(size(Sz)) 
            v = Sz**2
            EtaZeroInfty_Term = (v(1)-v(2)) / ((u-v(1))*(u-v(2)))
        end function EtaZeroInfty_Term
end function EtaZeroInfty


pure function ZetaJ(self, j, u)   ! $\zeta_j(u)/du$.
    class(SchottkyNumerics_Type), intent(in) :: self
    integer, intent(in) :: j
    complex(precision), intent(in) :: u(:)
    type(TraversalParameters_Type) :: params
    complex(precision) :: ZetaJ(size(u))

    associate(alpha => self%alpha)
        params%Operation => Sum
        params%rightCoset = [j];  allocate(params%leftCoset(0));  params%isReduced = .true.
        params%u = u**2
        params%z = [alpha(j), -alpha(j)]
        params%Term => ZetaJ_Term;  params%idTransformTerm = alpha(j)/(u**2-alpha(j)**2)
        params%eps = self%eps;  params%children = self%children
    end associate

    ZetaJ = 2*Traversal(self, params)

    contains 
        pure function ZetaJ_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: ZetaJ_Term(size(u))
                ZetaJ_term = (u+Sz(1)*Sz(2))*(Sz(1)-Sz(2)) / ((u-Sz(1)**2)*(u-Sz(2)**2))
        end function ZetaJ_Term
end function ZetaJ


pure function Zeta(self, u, zetaCoeff)
! $\zeta/du$, $\zeta = (i\pi)^{-1}(zetaCoeff(1)*\zeta_1 + ... + zetaCoeff(g)*\zeta_g$).
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    integer, intent(in) :: zetaCoeff(self%g)
    complex(precision) :: Zeta(size(u)), zetaJVal(size(u),self%g)
    integer :: j
    forall(j = 1:self%g) zetaJVal(:,j) = self%ZetaJ(j,u)
    Zeta = cmplx(0,-1/pi,precision) * matmul(zetaJVal,zetaCoeff)
end function Zeta


pure function Theta(self, j, u)
! \theta_T(u)/(du)^2 -- coefficient of even holomorphic quadratic differential: 
! 1<=j<=g -> T = S_j, 1+g<=j<=2*g-1 -> T = S_{j-g+1}S_1^{-1}S_{j-g+1}.
    class(SchottkyNumerics_Type), intent(in) :: self
    integer, intent(in) :: j
    complex(precision), intent(in) :: u(:)
    type(TraversalParameters_Type) :: params
    complex(precision) :: Theta(size(u))

    associate(alpha => self%alpha, g => self%g)
        params%Operation => Sum
        if (j <= g) then 
            params%rightCoset = [j]
        else 
            params%rightCoset = [j-g+1,-1,j-g+1]
        end if
        allocate(params%leftCoset(0))
        params%isReduced = .true.
        params%u = u
        params%z = [alpha(j), -alpha(j)]
        params%Term => Theta_Term;  params%idTransformTerm = (2*alpha(j)/(u**2-alpha(j)**2))**2
        params%eps = sqrt(self%eps);  params%children = self%children_sqrtEps

        Theta = (alpha(j)/2) * Traversal(self, params)
    end associate

    contains 
        pure function Theta_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision), dimension(size(u)) :: v1, v2, Theta_Term
            v1 = (((u-Sz(1))*(u-Sz(2))))**2
            v2 = (((u+Sz(1))*(u+Sz(2))))**2
            Theta_term = (v1+v2)*(Sz(1)-Sz(2))**2/(v1*v2)
        end function Theta_Term
end function Theta


pure function ExpIntEta(self, u, z, w)   ! $\exp\int_\infty^u\eta_{zw}$.
    class(SchottkyNumerics_Type), intent(in) :: self
    real(precision), intent(in) :: z, w
    complex(precision), intent(in) :: u(:)
    complex(precision) :: ExpIntEta(size(u))
    type(TraversalParameters_Type) :: params

    params%Operation => Product
    params%isReduced = .false.;  allocate(params%leftCoset(0), params%rightCoset(0))
    params%u = u;  params%z = [z,w]
    params%Term => ExpIntEta_Term;  params%idTransformTerm = (u-z)/(u-w)
    params%eps = self%eps;  params%children = self%children
    ExpIntEta = Traversal(self, params)

    contains 
        pure function ExpIntEta_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: ExpIntEta_Term(size(u))
            ExpIntEta_term = (u-Sz(1))/(u-Sz(2))
        end function ExpIntEta_Term
end function ExpIntEta


pure function ExpIntZetaJ(self, j, u)   ! $\exp\int_\infty^u\zeta_j$.
    class(SchottkyNumerics_Type), intent(in) :: self
    integer, intent(in) :: j
    complex(precision), intent(in) :: u(:)
    complex(precision) :: ExpIntZetaJ(size(u))
    type(TraversalParameters_Type) :: params

    associate(alpha => self%alpha)
        params%Operation => Product
        params%rightCoset = [j];  allocate(params%leftCoset(0))
        params%isReduced = .true.
        params%u = u;  params%z = [alpha(j),-alpha(j)]
        params%Term => ExpIntZetaJ_Term;  params%idTransformTerm = (u-alpha(j))/(u+alpha(j))
        params%eps = self%eps;  params%children = self%children
    end associate
    ExpIntZetaJ = Traversal(self, params)

    contains 
        pure function ExpIntZetaJ_Term(u, Sz)
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: ExpIntZetaJ_Term(size(u))
            ExpIntZetaJ_term = (u-Sz(1))*(u+Sz(2)) / ((u-Sz(2))*(u+Sz(1)))
        end function ExpIntZetaJ_Term
end function ExpIntZetaJ


real(precision) pure function ExpPeriod(self, l, j)   ! $\exp\int_{b_j}\zeta_l$, 1<=l,j<=g.
    class(SchottkyNumerics_Type), intent(in) :: self
    integer, intent(in) :: j, l
    complex(precision) :: v, ExpPeriodVec(1)
    type(TraversalParameters_Type) :: params

    associate(alpha => self%alpha, lambda => self%lambda)
        params%Operation => Product
        params%leftCoset = [l];  params%rightCoset = [j]
        params%isReduced = .true.
        params%u = [complex(precision) :: alpha(l)];  params%z = [alpha(j), -alpha(j)]
        params%Term => ExpPeriod_Term
        if (l /= j) then 
            v = ((alpha(j)-alpha(l))/(alpha(j)+alpha(l)))**2
        else 
            v = lambda(l)
        end if
        params%idTransformTerm = [v]
        params%eps = self%eps;  params%children = self%children
    end associate

    ExpPeriodVec = Traversal(self, params)
    ExpPeriod = ExpPeriodVec(1)

    contains 
        pure function ExpPeriod_term(u, Sz)   ! Squared cross ratio.
            complex(precision), intent(in) :: u(:)
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: ExpPeriod_term(size(u))
            ExpPeriod_term = (((u-Sz(1))*(u+Sz(2))) / ((u-Sz(2))*(u+Sz(1))))**2
        end function ExpPeriod_term
end function ExpPeriod


pure function PeriodMatrix(self)   ! ${\int_{b_j}\zeta_l, 1<=l,j<=g}$.
    class(SchottkyNumerics_Type), intent(in) :: self
    integer :: j, l
    real(precision) :: PeriodMatrix(self%g,self%g)
    associate(g => self%g)
        forall(j = 1:g)
            forall(l = j:g) PeriodMatrix(l,j) = log(self%ExpPeriod(j,l))
            PeriodMatrix(j,j+1:g) = PeriodMatrix(j+1:g,j)
        end forall
    end associate
end function PeriodMatrix


pure function IntZeta(self, u, zetaCoeff)   ! $\int_\infty^u\zeta$.
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    integer, intent(in) :: zetaCoeff(self%g)
    complex(precision) :: IntZeta(size(u)), ExpIntZetaJVec(size(u),self%g)
    integer :: j
    forall(j = 1:self%g) ExpIntZetaJVec(:,j) = self%ExpIntZetaJ(j,u)
    IntZeta = cmplx(0,-1/pi,precision) * matmul(log(ExpIntZetaJVec),zetaCoeff)
end function IntZeta


pure function Chi(self, u, b)   ! Projection $\chi(u): (b,0,\infty)\mapsto(1,0,\infty)$.
    class(SchottkyNumerics_Type), intent(in) :: self
    complex(precision), intent(in) :: u(:)
    real(precision), intent(in) :: b   ! basePoint.
    complex(precision) :: Chi(size(u)), uu(size(u))
    real(precision) :: zero, bb
    type(TraversalParameters_Type) :: params

    zero = 0
    uu = u**2
    bb = b**2
    params%Operation => Product
    params%isReduced = .true.;  allocate(params%leftCoset(0), params%rightCoset(0))
    params%u = uu;  params%z = [zero,-log(zero)]   ! -log(zero) = infty.
    params%Term => Chi_Term;  params%idTransformTerm = u/b
    params%eps = self%eps;  params%children = self%children

    Chi = Traversal(self, params)**2

    contains 
        pure function Chi_Term(u, Sz)
            complex(precision), intent(in) :: u(:) 
            real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
            complex(precision) :: Chi_Term(size(u))
            real(precision) :: v(size(Sz)) 
            v = Sz**2
            Chi_Term = (u-v(1))*(bb-v(2)) / ((bb-v(1))*(u-v(2)))
        end function Chi_Term
end function Chi


pure function Sum(arg1, arg2)
    complex(precision), intent(in) :: arg1(:), arg2(size(arg1))
    complex(precision) :: Sum(size(arg1))
    Sum = arg1 + arg2
end function Sum


pure function Product(arg1, arg2)
    complex(precision), intent(in) :: arg1(:), arg2(size(arg1))
    complex(precision) :: Product(size(arg1))
    Product = arg1 * arg2
end function Product

end module DiffsAndFuncsOnCurve_Module