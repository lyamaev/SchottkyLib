module DiffsAndFuncsOnCurve_Module   ! Differentials and functions on a curve in Schottky model.

! This module provides 'SchottkyWithNumerics_Type' -- extension of 'SchottkyGroup_Type' with methods for computing
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

use Precision_Module, only: precision, invPi
use SchottkyGroup_Module, only: SchottkyGroup_Type
use ChildrenOrder_Module, only: ChildrenOrder, Child_Type
use CayleyTreeTraversal_Module, only: Traversal_Interface, TraversalParameters_Type, & 
                                      Traversal_BogatyrevAlg, Traversal_NewAlg

implicit none 
private
public SchottkyWithNumerics_Type


procedure(Traversal_Interface), pointer :: Traversal   ! Module-scale traversal function.


type, extends(SchottkyGroup_Type) :: SchottkyWithNumerics_Type
	real(precision)  :: alpha(2*g-1)            ! AFP of $S_j$ and $S_jS_1^{-1}S_j$. 
	real(precision)  :: lambda(g)               ! lambda(j) := $S_j'\alpha_j$ -- multiplier of $S_j$.
	real(precision)  :: eps = 1e-12_precision   ! Parameter to control accuracy of (a posteriori) algs.	
	logical :: useNewAlg = .false.              ! Use Bogatyrev or NewAlg?

	! Children order (with per child eps) to use in NewAlg: 'children' -- for linear differentials and 
	! functions associated with them, 'children_WithSqrtEps' -- for quadratic differentials.
	type(Child_Type), dimension(-g:g,1:2*g-1), private :: children, children_WithSqrtEps   
contains
	procedure :: Eta            ! $\eta_{zw}(u)/du$ -- coeff of Abelian diff of 3rd kind.
	procedure :: EtaSymm        ! $(\eta_{zw}(u)+\eta_{-z-w}(u)/)du$.
	procedure :: EtaZeroInfty   ! $\eta_{0\infty}(u)/du$.
	procedure :: ZetaJ          ! $\zeta_j(u)/du$ -- coeff of Abelian diff of 1st kind.  
	procedure :: Zeta           ! $\zeta(u)/du$.
	procedure :: Theta          ! $\theta_j(u)/(du)^2$ -- coeff of even holomorphic quadratic diff.
	procedure :: ExpIntZetaJ    ! $\exp\int_\infty^u\zeta_j$.
	procedure :: ExpPeriod      ! $\exp\int_{b_l}\zeta_j$. 
	procedure :: PeriodMatrix   ! $\{\int_{b_j}\zeta_l, 1<=l,j<=g\}$.
	procedure :: IntZeta        ! $\int_\infty^u\zeta$.
	procedure :: Chi            ! $\chi(u) := (\int_b^u\eta_{0\infty)^2}$. 
	procedure :: Init => InitSchottkyWithNumerics
end type SchottkyWithNumerics_Type


contains


subroutine InitSchottkyWithNumerics(self, c_new, r_new, sigma_new)
	class(SchottkyWithNumerics_Type(*)), intent(inout) :: self
	real(precision), intent(in) :: c_new(self%g), r_new(self%g)
	integer, intent(in) :: sigma_new(self%g)

	call self%SchottkyGroup_Type%Init(c_new, r_new, sigma_new)   ! Call parent initializer.

	if (self%useNewAlg) then
		Traversal => Traversal_NewAlg
		self%children = ChildrenOrder(self, self%eps)
		self%children_WithSqrtEps(:,:)%childIndex = self%children(:,:)%childIndex
		self%children_WithSqrtEps(:,:)%childEps = sqrt(self%children(:,:)%childEps)
	else 
		Traversal => Traversal_BogatyrevAlg
	end if

	associate(g => self%g, c => self%c, rr => self%rr, alpha => self%alpha, lambda => self%lambda)
		alpha(1:g) = sqrt(c(1:g)**2 - rr(1:g))
		alpha(g+1:2*g-1) = sqrt( ((alpha(2:g)**2-c(2:g)*c(1))**2 - rr(1)*c(2:g)**2) / &
			                     ((c(2:g)-c(1))**2 - rr(1)) )
		lambda(1:g) = rr(1:g)/(c(1:g) + alpha(1:g))**2 
	end associate
end subroutine InitSchottkyWithNumerics


function Eta(self, u, z)   ! $\eta_{z_1z_2}(u)/du$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	real(precision), intent(in) :: z(2)
	type(TraversalParameters_Type) :: params
	complex(precision) :: Eta(size(u))

	params%Operation => Sum
	params%isReduced = .false.
	allocate(params%u, source=u)
	allocate(params%z, source=z)
	params%Term => Eta_Term
	allocate(params%idTransformTerm, source=Eta_Term(u,z))
	params%eps = self%eps; 	allocate(params%children, source=self%children) 

	Eta = Traversal(self, params)

	contains 
		pure function Eta_Term(u, Sz)
			complex(precision), intent(in) :: u(:)
			real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
			complex(precision) :: Eta_Term(size(u))
			Eta_Term = (Sz(1)-Sz(2)) / ((u-Sz(1))*(u-Sz(2)))
		end function Eta_Term
end function Eta


function EtaSymm(self, u, z)   ! $ (\eta_{z_1z_2}+\eta_{-z_1-z_2})(u)/du$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	real(precision), intent(in) :: z(2)	
	type(TraversalParameters_Type) :: params
	complex(precision) :: EtaSymm(size(u))

	params%Operation => Sum
	params%isReduced = .true.
	allocate(params%u, source = u**2)   ! Compute u**2 straightaway and once.
	allocate(params%z, source = [z(1),z(2),-z(1),-z(2)])
	params%Term => EtaSymm_Term
	allocate(params%idTransformTerm, source = (z(1)**2-z(2)**2)/((u**2-z(1)**2)*(u**2-z(2)**2)))
	params%eps = self%eps; 	allocate(params%children, source=self%children) 

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


function EtaZeroInfty(self, u)   ! $\eta_{0\infty}(u)/du$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	real(precision) :: zero = 0	
	type(TraversalParameters_Type) :: params
	complex(precision) :: EtaZeroInfty(size(u))

	params%Operation => Sum
	params%isReduced = .true.
	allocate(params%u, source = u**2)                 ! Compute u**2 straightaway and once.
	allocate(params%z, source = [zero, -log(zero)])   ! -log(zero) = infty.
	params%Term => EtaZeroInfty_Term
	allocate(params%idTransformTerm, source = 1/(2*u**2))
	params%eps = self%eps; 	allocate(params%children, source=self%children) 

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


function ZetaJ(self, j, u)   ! $\zeta_j(u)/du$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	integer, intent(in) :: j
	complex(precision), intent(in) :: u(:)
	type(TraversalParameters_Type) :: params
	complex(precision) :: ZetaJ(size(u))

	associate(alpha => self%alpha)
		params%Operation => Sum
		allocate(params%rightCoset, source = [j])
		params%isReduced = .true.
		allocate(params%u, source = u**2)
		allocate(params%z, source = [alpha(j), -alpha(j)])
		params%Term => ZetaJ_Term
		allocate(params%idTransformTerm, source = alpha(j)/(u**2-alpha(j)**2))
		params%eps = self%eps; 	allocate(params%children, source=self%children) 
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


function Zeta(self, u, zetaCoeff)
! $\zeta/du$, $\zeta = (i\pi)^{-1}(zetaCoeff(1)*\zeta_1 + ... + zetaCoeff(g)*\zeta_g$).
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	integer, intent(in) :: zetaCoeff(self%g)
	complex(precision) :: Zeta(size(u)), zetaJVal(size(u),self%g)
	integer :: j
	forall(j = 1:self%g) zetaJVal(:,j) = self%ZetaJ(j,u)
	Zeta = cmplx(0,-invPi,precision) * matmul(zetaJVal,zetaCoeff)
end function Zeta


function Theta(self, j, u)
! \theta_T(u)/(du)^2 -- coefficient of even holomorphic quadratic differential: 
! 1<=j<=g -> T = S_j, 1+g<=j<=2*g-1 -> T = S_{j-g+1}S_1^{-1}S_{j-g+1}.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	integer, intent(in) :: j
	complex(precision), intent(in) :: u(:)
	type(TraversalParameters_Type) :: params
	complex(precision) :: Theta(size(u))

	associate(alpha => self%alpha, g => self%g)
		params%Operation => Sum
		if (j <= g) then 
			allocate(params%rightCoset, source = [j])
		else 
			allocate(params%rightCoset, source = [j-g+1,-1,j-g+1])
		end if
		params%isReduced = .true.
		allocate(params%u, source = u)
		allocate(params%z, source = [alpha(j), -alpha(j)])
		params%Term => Theta_Term
		allocate(params%idTransformTerm, source = (2*alpha(j)/(u**2-alpha(j)**2))**2)
		params%eps = sqrt(self%eps);  allocate(params%children, source=self%children_WithSqrtEps) 

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


function ExpIntZetaJ(self, j, u)   ! $\exp\int_\infty^u\zeta_j$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	integer, intent(in) :: j
	complex(precision), intent(in) :: u(:)
	complex(precision) :: ExpIntZetaJ(size(u))
	type(TraversalParameters_Type) :: params

	associate(alpha => self%alpha)
		params%Operation => Product
		allocate(params%rightCoset, source=[j])
		params%isReduced = .true.
		allocate(params%u, source=u)
		allocate(params%z, source=[alpha(j),-alpha(j)])
		params%Term => ExpIntZetaJ_Term
		allocate(params%idTransformTerm, source=(u-alpha(j))/(u+alpha(j)))
		params%eps = self%eps; 	allocate(params%children, source=self%children) 
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


real(precision) function ExpPeriod(self, l, j)   ! $\exp\int_{b_j}\zeta_l$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	integer, intent(in) :: j, l
	complex(precision) :: v, ExpPeriodVec(1)
	type(TraversalParameters_Type) :: params

	associate(alpha => self%alpha, lambda => self%lambda)
		params%Operation => Product
		allocate(params%leftCoset, source=[l])
		allocate(params%rightCoset, source=[j])
		params%isReduced = .true.
		allocate(params%u, source = [complex(precision) :: alpha(l)])
		allocate(params%z, source = [alpha(j), -alpha(j)])
		params%Term => ExpPeriod_Term
		if (l /= j) then 
			v = ((alpha(j)-alpha(l))/(alpha(j)+alpha(l)))**2
		else 
			v = lambda(l)
		end if
		allocate(params%idTransformTerm, source=[v])
		params%eps = self%eps; 	allocate(params%children, source=self%children) 
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


function PeriodMatrix(self)   ! ${\int_{b_j}\zeta_l, 1<=l,j<=g}$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	integer :: j, l
	real(precision) :: PeriodMatrix(self%g,self%g)
	associate(g => self%g)
		forall(j = 1:g)
			forall(l = j:g) PeriodMatrix(l,j) = log(self%ExpPeriod(j,l))
			PeriodMatrix(j,j+1:g) = PeriodMatrix(j+1:g,j)
		end forall
	end associate
end function PeriodMatrix


function IntZeta(self, u, zetaCoeff)   ! $\int_\infty^u\zeta$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	integer, intent(in) :: zetaCoeff(self%g)
	complex(precision) :: IntZeta(size(u)), ExpIntZetaJVec(size(u),self%g)
	integer :: j
	forall(j = 1:self%g) ExpIntZetaJVec(:,j) = self%ExpIntZetaJ(j,u)
	IntZeta = cmplx(0,-invPi,precision) * matmul(log(ExpIntZetaJVec),zetaCoeff)
end function IntZeta


function Chi(self, u, b)   ! Projection $\chi(u): (b,0,\infty)\mapsto(1,0,\infty)$.
	class(SchottkyWithNumerics_Type(*)), intent(in) :: self	
	complex(precision), intent(in) :: u(:)
	real(precision), intent(in) :: b   ! basePoint.
	complex(precision) :: Chi(size(u)), uu(size(u))
	real(precision) :: zero = 0, bb
	type(TraversalParameters_Type) :: params

	uu = u**2
	bb = b**2
	params%Operation => Product
	params%isReduced = .true.
	allocate(params%u, source=uu)
	allocate(params%z, source=[zero,-log(zero)])   ! -log(zero) = infty.
	params%Term => Chi_Term
	allocate(params%idTransformTerm, source = uu/bb)
	params%eps = sqrt(self%eps);  allocate(params%children, source=self%children_WithSqrtEps) 


	Chi = Traversal(self, params)**2

	contains 
		function Chi_Term(u, Sz)
			complex(precision), intent(in) :: u(:) 
			real(precision), intent(in) :: Sz(:)   ! dim(Sz) = 2.
			complex(precision) :: Chi_Term(size(u))
			real(precision) :: v(size(Sz)) 
			v = Sz**2
			Chi_Term = (u-v(1))*(bb-v(2)) / ((bb-v(1))*(u-v(2)))
		end function Chi_Term
end function Chi


function Sum(arg1, arg2)
	complex(precision), intent(in) :: arg1(:), arg2(size(arg1))
	complex(precision) :: Sum(size(arg1))
	Sum = arg1 + arg2
end function Sum


function Product(arg1, arg2)
	complex(precision), intent(in) :: arg1(:), arg2(size(arg1))
	complex(precision) :: Product(size(arg1))
	Product = arg1 * arg2
end function Product


end module DiffsAndFuncsOnCurve_Module