module CayleyTreeTraversal_Module

! Provides universal Cayley tree traversal interface to compute functions and differentials on
! real hyperelliptic curve in Schottky model. 
!
! Parameters 'isReduced', 'leftCoset', 'rightCoset' specify subset of Schottky group to traverse over: 
! isReduced == .false. --> $\mathfrak{S} / (S \sim T_1S \sim S T_2)$, 
! isReduced == .true.  --> $\mathfrak{S} / (S \sim S^{-} \sim T_1 S \sim S T_2)$,
! where $T_1$ and $T_2$ are specified by parameters 'leftCoset' and 'rightCoset' correspondingly.
! For each $S = S_{j_1}^{i_1}\dots S_{j_k}^{i_k} \in \mathfrak{S}$ we define $S^-$ by the formula
! $S^- := S_{j_1}^{-i_1}\dots S_{j_k}^{-i_k}$.
! 
! Arrays 'u' and 'z' are arguments of 'Term' function and have slightly different meanings.
! In numerical applications (for example, in Chebyshev ansatz) one often wants to compute the same 
! object (function or differential) at several points: 
! --  'u' is an array of points at which we want to compute the same object;
! --  'z' contains all term parameters of this (single) object. 
! Sz is computed for each visited node S (and is stored in 'Sz' array), Su is not. 
! REQUIREMENTS: dim(u)>=1, dim(z)>=2, z(1)/=z(2).
! Example: if one wants to compute $\eta_{z_1z_2}(u)/du (that is coefficient of normalized Abelian 
! differential of the 3rd kind with poles $z_1$ and $z_2$) at points $u_1$, $u_2$ and $u_3$, then 
! z = [z_1,z_2], u = [u_1,u_2,u_3], Term(u,z) = 1/(u(1:3)-z(1)) - 1/(u(1:3)-z(2)).
! We use real type for 'z'-array and complex type for 'u'-array because it's the most convenient
! choice for Chebyshev anzats implementation.

! 'Operation_Interface' gives template for operation committed on each iteration: that's Sum in case 
! of series and Product in case of infinite product.

use Precision_Module, only: precision
use SchottkyGroup_Module, only: SchottkyGroup_Type
use ChildrenOrder_Module, only: Child_Type
implicit none
private
public Traversal_Interface, TraversalParameters_Type, Traversal_BogatyrevAlg, Traversal_NewAlg

integer, parameter :: maxLen = 10000   ! Max length of branch. All nodes S: |S|>maxLen are excluded from traversal.

interface 
    pure function Term_Interface(u, z)
        import :: precision
        complex(precision), intent(in) :: u(:)
        real(precision), intent(in) :: z(:)
        complex(precision) :: Term_Interface(size(u))
    end function Term_Interface
end interface

interface 
    pure function Operation_Interface(arg1, arg2)
        import :: precision
        complex(precision), intent(in) :: arg1(:), arg2(size(arg1))
        complex(precision) :: Operation_Interface(size(arg1))
    end function Operation_Interface
end interface

type :: TraversalParameters_Type
    procedure(Operation_Interface), pointer, nopass :: Operation   ! Sum or Product.
    real(precision), allocatable :: z(:)                    ! $Sz$ is computed for each visited node $S$.
    complex(precision), allocatable :: u(:)                 ! $Su$ is not computed.
    integer, allocatable :: leftCoset(:), rightCoset(:)     ! Cosets for relative series/products.
    logical :: isReduced                                    ! Is sum/prod reduced with equity $S^{-}(-u) = -Su$?
    procedure(Term_Interface), pointer, nopass :: Term      ! Term of sum/product.
    complex(precision), allocatable :: idTransformTerm(:)   ! (size(u))  Computed term for root.
    real(precision) :: eps                                  ! Parameter to control accuracy (for BogatyrevAlg).
    type(Child_Type), allocatable :: children(:,:)          ! (-g:g,1:2*g+1)  Children order and childEps (for NewAlg).
end type TraversalParameters_Type

interface 
    pure function Traversal_Interface(schottkyGroup, params)
        import :: precision, SchottkyGroup_Type, TraversalParameters_Type
        class(SchottkyGroup_Type), intent(in) :: schottkyGroup
        type(TraversalParameters_Type), intent(in) :: params
        complex(precision) :: Traversal_Interface(size(params%u))
    end function Traversal_Interface
end interface

contains

pure function Traversal_BogatyrevAlg(schottkyGroup, params) result(res)
! Bogatyrev's algorithm. Implementation is based on depth-first search and lexicographical order.
    class(SchottkyGroup_Type), intent(in) :: schottkyGroup
    type(TraversalParameters_Type), intent(in) :: params
    complex(precision) :: res(size(params%idTransformTerm))
    integer :: len                                    ! Length of current branch: len := $|T|$.
    integer :: letter(0:maxLen)                       ! Word representation of current node T.
    real(precision) :: Sz(size(params%z), 0:maxLen)   ! Sz for all nodes S in current branch.
    logical :: isNotFound, shouldContribute 

    associate(g => schottkyGroup%g, lc => params%leftCoset, rc => params%rightCoset)
        len = 0
        letter(0) = 0
        Sz(:,0) = params%z 
        res = params%idTransformTerm

        Traversal_Loop: do 
            ! Find next node.
            if (abs(Sz(1,len)-Sz(2,len)) >= params%eps .and. len < maxLen) then
                len = len + 1          ! Grow current branch. 
                letter(len) = -g - 1   ! We will add +1 in FindNextNode_Loop.
            end if
            isNotFound = .true.       ! Next node is not found yet.
            FindNextNode_Loop: do while (isNotFound)
                do while (letter(len) == g)
                    len = len - 1          ! Layer is full: go one layer backward.
                    if (len == 0) RETURN   ! Finish traversal if we are at root again.
                end do
                letter(len) = letter(len) + 1
                if (letter(len) == 0) letter(len) = 1   ! Jump over 0 (no such letter).
                isNotFound = (letter(len) == -letter(len-1)) .or.                        &  
                             (len == 1 .and. params%isReduced .and. letter(1) < 0) .or. &   
                             (len == size(rc) .and. AreEqualOrInverse(letter(1:len),rc))   
            end do FindNextNode_Loop

            ! Contribute node to sum/prod.  
            Sz(:,len) = schottkyGroup%S(letter(len), Sz(:,len-1))
            shouldContribute = (size(lc) == 0) .or. (len < size(lc)) .or. &
                               (.not. AreEqualOrInverse(letter(len-size(lc)+1:len),lc))
            if (shouldContribute) res = params%Operation(res, params%Term(params%u, Sz(:,len))) 
        end do Traversal_Loop
    end associate
end function Traversal_BogatyrevAlg


pure function Traversal_NewAlg(schottkyGroup, params) result(res)   
    class(SchottkyGroup_Type), intent(in) :: schottkyGroup
    type(TraversalParameters_Type), intent(in) :: params
    complex(precision) :: res(size(params%idTransformTerm))
    integer :: len                                    ! Length of current branch: len := $|T|$.
    integer :: letter(0:maxLen)                       ! Word representation of current node T.
    integer :: child_num(0:maxLen)                    ! Representation of current node T by child order numbers.
    real(precision) :: Sz(size(params%z), 0:maxLen)   ! Sz for all nodes S in current branch.
    real(precision) :: dist(0:maxLen)                 ! |Sz_1-Sz_2| for all nodes S in current branch.
    logical :: shouldContribute 

    associate(g => schottkyGroup%g, lc => params%leftCoset, rc => params%rightCoset)
        len = 0
        letter(0) = 0
        Sz(:,0) = params%z 
        dist(0) = abs(Sz(1,0)-Sz(2,0))
        res = params%idTransformTerm

        Traversal_Loop: do 
            len = len + 1   ! Grow current branch.
            child_num(len) = 0

            FindNextNode_Loop: do
                do while (len > 1 .and. child_num(len) == 2*g-1)
                    len = len - 1   ! Сurrent layer is full, go one layer backward.
                end do
                if (len == 1 .and. child_num(len) == 2*g) RETURN   ! Finish if the first layer is full.
                
                child_num(len) = child_num(len) + 1   ! Check next child of letter(1:len-1). 

                if (len == 1) then
                    if (child_num(len) <= g) then  
                        if (params%isReduced) CYCLE FindNextNode_Loop
                        letter(len) = -child_num(len)
                    else
                        letter(len) = child_num(len) - g
                    end if
                else 
                    if (dist(len-1) < params%children(letter(len-1),child_num(len))%childEps) then 
                        len = len - 1   ! Throw away current node together with all its younger brothers.
                        CYCLE FindNextNode_Loop
                    end if
                    letter(len) = params%children(letter(len-1),child_num(len))%childIndex
                end if

                if (len /= size(rc) .or. .not. AreEqualOrInverse(letter(1:len),rc)) EXIT FindNextNode_Loop
            end do FindNextNode_Loop

            ! Contribute node to sum/prod.  
            Sz(:,len) = schottkyGroup%S(letter(len), Sz(:,len-1))
            dist(len) = abs(Sz(1,len) - Sz(2,len))
            shouldContribute = (size(lc) == 0) .or. (len < size(lc)) .or. &
                               (.not. AreEqualOrInverse(letter(len-size(lc)+1:len),lc))
            if (shouldContribute) res = params%Operation(res, params%Term(params%u, Sz(:,len)))
        end do Traversal_Loop
    end associate
end function Traversal_NewAlg


logical pure function AreEqualOrInverse(word1, word2)   
! Returns true if $S = T$ or $S = T^{-1}$ where $S,T\in\mathfrak{S}$ are given by word1 and word2.
    integer, intent(in) :: word2(:), word1(size(word2))
    AreEqualOrInverse = all(word1 == word2) .or. all(word1 == -word2(size(word2):1:-1)) 
end function AreEqualOrInverse

end module CayleyTreeTraversal_Module