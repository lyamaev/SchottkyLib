program Example1   
! Check numerically that reciprocity law for Abelian differentials holds: $\int_{b_j}\eta_{zw}=\int_w^z\zeta_j$, 
! where $\zeta_j$ and $\eta_{zw}$ are the normalized Abelian differentials of the 1st and the 3rd kinds, respectively.

use Precision_Module, only: p => precision
use DiffsAndFuncsOnCurve_Module, only: SchottkyNumerics_Type

type(SchottkyNumerics_Type) :: schottkyGroup
complex(p) :: lhs(2), rhs(2), a
integer :: j
real(p) :: z, w

! Construct Schottky group of genus 3 from the given Schottky moduli using class constructor.
schottkyGroup = SchottkyNumerics_Type(c = [0.1_p, 0.55_p, 1.0_p], r = [0.01_p, 0.1_p, 0.05_p], sigma = [1, -1, 1])

j = 2
z = -2.1_p
w = 4.3_p

! Pick any point $a$ in fundamental domain of Schottky group, result doesn't depend on it.
a = cmplx(1.7_p, -2.8_p, p)   

! lhs(1) = $\exp(\int_\infty^{S_j(a)} \eta_{zw})$,  lhs(2) = $\exp(\int_\infty^{a} \eta_{zw})$.
! $b_j$ is the cycle that connects points $S_j(a)$ and $a$.
lhs = schottkyGroup%ExpIntEta([schottkyGroup%S(j,a), a], z, w)

! rhs(1) = $\exp(\int_\infty^{z} \zeta_j)$,  rhs(2) = $\exp(\int_\infty^{w} \zeta_j)$.
rhs = schottkyGroup%ExpIntZetaJ(j, [cmplx(z,0,p), cmplx(w,0,p)])

print *, 'LHS and RHS of reciprocity law:', log(real([lhs(1)/lhs(2), rhs(1)/rhs(2)], kind=4))

! Output:
! LHS and RHS of reciprocity law:  0.8132047      0.8132047 
end program Example1