module ChildrenOrder_Module
! Provides child(-g:g,1:2*g-1) -- table of children ordered "by age". Namely, any node $T\ne\id$ has 2g-1 
! children -- they are sorted in array child(t,1:2*g-1) by ascending of childEps(child) := eps/M(child), 
! where t is an index of the leftmost letter of T. 
! M(child S_jT) := M(j,t) := L(-j,t)(K(t)+1), where L(j,t) -- upper estimate of $f(S_{-j}T)/f(T)$, 
! K(t) -- upper estimate of $\sum_{S: S>T} f(S)/f(T)$. Here $f(S) := |Sz_1-Sz_2|$ (or $f(S) := |S'u|$). 
! For child S_jT define parentIndex := t, childIndex := j.

use Precision_Module, only: precision
use SchottkyGroup_Module, only: SchottkyGroup_Type
implicit none
private 
public ChildrenOrder, Child_Type

type :: Child_Type   
    integer :: childIndex         ! Equals j for child $S_jT$, $j\ne -t$.
    real(precision) :: childEps   ! eps/M(j,t).
end type Child_Type

contains

function ChildrenOrder(a, eps) result(children)
    class(SchottkyGroup_Type), intent(in) :: a   ! Schottky group.
    type(Child_Type) :: children(-a%g:a%g,1:2*a%g-1)   
    real(precision), intent(in) :: eps
    integer :: j, t, cnt
    real(precision) :: L(-a%g:a%g,-a%g:a%g)   ! L(j,t) -- estimate of $f(S_{-j}T)/f(T)$.
    real(precision) :: K(-a%g:a%g)            ! K(t) -- estimate of $\sum_{S: S>T} f(S)/f(T)$.
    real(precision) :: M(-a%g:a%g,-a%g:a%g)   ! M(j,t) := L(-j,t)(K(t)+1).
    real(precision) :: ip(-1:1,-a%g:a%g)      ! ip(\pm1,j) -- intersection points $C_j\cap\R$.
    real(precision) :: sumOfDiamsEstimate     ! Estimate of $\sum_{S\ne\id} diam(S\cF)$.
    real(precision) :: leftmostRight(a%g), rightPrev, lambda(a%g), lambdaMax, gamma(a%g), h(-a%g:a%g,1:a%g)

    associate(g => a%g, c => a%c, r => a%r, sigma => a%sigma)	
        ! Generate ip(\pm1,j) := $p^\pm_j$ -- intersection points of $C_j$ with real line. If $\sigma_j = -1$ then 
        ! circle $C_j$ is not uniqly determined by moduli, and so intersection points $p^\pm_j$.
        forall(j = 1:g, sigma(j) == 1) 
            ip(-1,j) = c(j) - r(j);  ip(+1,j) = c(j) + r(j)
        end forall
        if (any(sigma(1:g) == -1)) then
            rightPrev = 0
            do j = 1,g
                if (sigma(j) == -1) then
                    leftmostRight(j) = a%S(j,-rightPrev)  ! Leftmost position of $p^+_j$.
                    rightPrev = leftmostRight(j)
                else 
                    rightPrev = ip(+1,j)
                end if
            end do
            do j = g,1,-1
                if (sigma(j) == -1) then
                    if (j == g) then 
                        ip(+1,j) = leftmostRight(g) + 1
                    else
                        ip(+1,j) = (leftmostRight(j) + ip(-1,j+1))/2
                    end if
                    ip(-1,j) = a%S(j, -ip(+1,j))
                end if
            end do
        end if
        ip(+1,-g:-1) = -ip(-1,g:1:-1);  ip(-1,-g:-1) = -ip(+1,g:1:-1)

        ! Compute L(j,t) and lambda(t) := sum(L(:,t)). Schmies' estimate is only applicable, if max(lambda(:))<1.
        forall(j = -g:g, t = -g:g, j /= t .and. j /= 0 .and. t /= 0) 
            L(j,t) = (r(j) / (ip(sign(1,j-t),t) - c(j)))**2
        end forall
        forall(j = -g:g) 
            L(j,j) = 0;  L(j,0) = 0;  L(0,j) = 0
        end forall
        forall(t = 1:g) lambda(t) = sum(L(-g:g,t))
        lambdaMax = maxval(lambda(1:g))

        if (lambdaMax < 1) then
            K(1:g) = lambda(1:g)/(1-lambdaMax)   ! Schimes' estimate.
            K(-g:-1) = K(g:1:-1)   
        else
            ! $\gamma_k$ is an upper estimate of $\sum_{l: SS_l>S_l}\diam(SS_l\cF)/diam(S\cF)$, where 
            ! $k$ is an index of the rightmost letter of $S\ne\id$, $\cF$ -- fundamental domain.
            forall(j = 2:g-1) gamma(j) = ((ip(-1,j+1) - ip(-1,j)) / (ip(-1,j+1) - ip(+1,j))) / &
                                         ((ip(+1,j-1) - ip(-1,j)) / (ip(+1,j-1) - ip(+1,j)))
            gamma(1) = ((ip(-1,2) - ip(-1,1)) / (ip(-1,2) - ip(+1,1))) / &
                       ((ip(+1,-1) - ip(-1,1)) / (ip(+1,-1) - ip(+1,1)))    
            gamma(g) = ((ip(-1,-g) - ip(-1,g)) / (ip(-1,-g) - ip(+1,g))) / &
                       ((ip(+1,g-1) - ip(-1,g)) / (ip(+1,g-1) - ip(+1,g)))	

            ! Upper estimate of $\sum_{S\ne\id} diam(S\cF)$.
            sumOfDiamsEstimate = (sqrt(maxval(gamma)) + 1) * sum(ip(+1,1:g) - ip(-1,1:g))

            forall(j = -g:g, t = 1:g, j /= t .and. j /= 0) 
                h(j,t) = (ip(+1,j) - ip(-1,j)) / (4*(ip(sign(1,j-t),t) - ip(-1,j))*(ip(sign(1,j-t),t) - ip(+1,j)))
            end forall
            forall(t = -g:g) h(t,t) = 0
            h(0,1:g) = 0
        
            forall(t = 1:g) K(t) = sumOfDiamsEstimate * maxval(h(-g:g,t))   ! New estimate (Schmies' isn't applicable).
            K(-g:-1) = K(g:1:-1)
        end if

        M = 0
        forall(t = -g:g, t /= 0) M(-g:g,t) = (1 + K(-g:g)) * L(g:-g:-1,t)
        M = M/maxval(M)   ! Normalize.

        do t = 1,g 
            cnt = 1
            do j = -g,g 
                if (j /= 0 .and. j /= -t) then 
                    children(t,cnt)%childIndex = j
                    children(t,cnt)%childEps = eps/M(j,t)
                    cnt = cnt + 1
                end if
            end do
            call QSort(children(t,1:2*g-1))
            children(-t,1:2*g-1)%childIndex = -children(t,1:2*g-1)%childIndex
            children(-t,1:2*g-1)%childEps = children(t,1:2*g-1)%childEps
        end do
    end associate
end function ChildrenOrder


recursive subroutine QSort(array)   ! Quick sort by ascending of array(:)%childEps.
    type(Child_Type), intent(inout) :: array(:)
    type(Child_Type) :: temp
    integer :: n, left, right, marker
    real(precision) :: random, pivot

    n = size(array)
    if (n > 1) then
        call random_seed()
        call random_number(random)
        pivot = array(int(random*(n-1))+1)%childEps
        left = 1;  right = n

        Main_Loop: do
            if (left >= right) exit
            FindRight_Loop: do
                if (array(right)%childEps <= pivot) exit
                right = right - 1
            end do FindRight_Loop
            FindLeft_Loop: do
                if (array(left)%childEps >= pivot) exit
                left = left + 1
            end do FindLeft_Loop
            if (left < right) then
                temp = array(left);  array(left) = array(right);  array(right) = temp
            end if
        end do Main_Loop

        if (left == right) then
            marker = left + 1
       	else
            marker = left
        end if

        call QSort(array(:marker-1))
        call QSort(array(marker:))
    end if
end subroutine QSort

end module ChildrenOrder_Module