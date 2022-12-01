# SchottkyLib

SchottkyLib is a tool for computation of analytical objects (meromorphic functions, differentials of different orders, spinors, etc) on real hyperelliptic curves using Schottky uniformization. 

*Real hyperelliptic curve* (over the complex numbers) is a compact curve with affine part given by the equation $y^2 = P(x),~(x,y)\in\mathbb{C}^2$, where $P(x)$ is a real polynomial without multiple roots. Computations with such curves and their moduli space arise in numerous applications, such as Chebyshev optimization (computation of multiband filters using Chebyshev ansatz, computation of Chebyshev polynomials for several intervals), finding algebro-geometric solutions of nonlinear partial differential equations, description of magnetic states in planar magnetic nanoelements, simulation of the water flow under a stepped dam, and others. 

In *Schottky model*, a real hyperelliptic curve is represented as an orbit manifold of action of a suitable Schottky group, analytical objects on the curve can then be explicitly represented in terms of Poincaré $\theta$-series. Schottky model is free of the complicated modular dependencies inherent in other approaches to computations with such curves (Riemann $\theta$-functions, Weierstrass $\sigma$-functions) and allows one to deal with curves of large genera (up to 100 and larger).

Modules of SchottkyLib:

- `Precision` — Set precision of floating-point numbers (real and complex): single, double or quadruple.
- `SchottkyGroup` — A class `SchottkyGroup_Type` is defined which contains main parameters of Schottky group: genus of associated real hyperelliptic curve and its Schottky moduli. Methods of this class are generators of Schottky group and their inverses, as well as indicator function of moduli space.
- `CayleyTreeTraversal` — Provides a universal interface for computation of sums and products over Schottky group and its quotients. Two traversal algorithms are present: Bogatyrev’s algorithm and a new faster algorithm[^2] proposed by the author.
- `ChildrenOrder` — Auxiliary computations for the new algorithm: estimate needed in the new traversal algorithm and order of children nodes by seniority for each parent.
- `DiffAndFuncsOnCurve` — Contains `SchottkyNumerics_Type`, the main class of SchottkyLib, which extends `SchottkyGroup_Type`. Added methods are linear and quadratic differentials (all needed in Chebyshev ansatz), Schottky functions, period matrix, and projection from Schottky model to algebraic model. If some analytical object you need is missing in current implementation it is not difficult to add it using the universal traversal interface provided by `CayleyTreeTraversal` module.
- `Uniformization` — Implementation of the classical algorithm (A. Poincaré, 1884) of numerical Schottky uniformization for M-curves: translation of an M-curve given in algebraic model to Schottky model.

SchottkyLib is a part of a bigger package developed by the author during his PhD research. The bigger library is an implementation of Chebyshev ansatz[^4][^5] for practical synthesis of quasi-optimal multiband filters (of different types: digital, analog and microwave) and also contains evaluation of gradients of Abelian integrals and their periods by Schottky moduli using Hejhal’s approach. 

SchottkyLib is written in Fortran 2018 and is tested to compile and run successfully with gfortran 12.2.0 and ifort 2022.2.1. 

Here are two short examples (wrapped in blocks inside one program) to try it out:

```Fortran
program Examples

use Precision_Module, only: p => precision
use DiffsAndFuncsOnCurve_Module, only: SchottkyNumerics_Type
use Uniformization_Module, only: PoincareUniformization


Example1: block
    type(SchottkyNumerics_Type) :: schottkyGroup
    complex(p) :: lhs(2), rhs(2), a
    integer :: j
    real(p) :: z, w

    ! Construct Schottky group of genus 3 from the given Schottky moduli using 
    ! SchottkyNumerics_Type constructor.
    schottkyGroup = SchottkyNumerics_Type(c=[0.1_p, 0.55_p, 1.0_p], &
                                          r=[0.01_p, 0.1_p, 0.05_p], &
                                          sigma=[1, -1, 1])

    ! Check numerically that reciprocity law for Abelian differentials of 1st and 
    ! 3rd kinds holds. 
    ! a is any point in fundamental domain of Schottky group,
    ! result doesn't depend on it.
    a = cmplx(1.7_p, -2.8_p, p) 
    j = 2;  z = -2.1_p;  w = 4.3_p
    lhs = schottkyGroup%ExpIntEta([schottkyGroup%S(j,a), a], z, w)
    rhs = schottkyGroup%ExpIntZetaJ(j, [cmplx(z,0,p), cmplx(w,0,p)])
    print *, 'Example1. LHS and RHS of reciprocity law:', &
             log(real([lhs(1)/lhs(2), rhs(1)/rhs(2)], kind=4))
end block Example1


Example2: block
    type(SchottkyNumerics_Type) :: schottkyGroup
    complex(p), allocatable :: ramInSchM(:), ramInAlgM(:)

    ! Construct Schottky group of genus 2 that corresponds to M-curve given by the 
    ! equation w^2 = x(x-0.1)(x-0.2)(x-1)(x-1.3) using Poincare opening process with 
    ! 5 generations of slots.
    schottkyGroup = PoincareUniformization(ramInAlgM=[0.1_p, 0.2_p, 1.0_p, 1.3_p], &
                                           slotsGen=4)

    associate(c => schottkyGroup%c, r => schottkyGroup%r)
        ! Ramification points in Schottky model.
        ramInSchM = cmplx([c(1)-r(1), c(1)+r(1), c(2)-r(2), c(2)+r(2)], 0, p)
    end associate

    ! Chi is the projection from Schottky model to algebraic model, normalized by 0->0, 
    ! infty->infty, b->1.
    ramInAlgM = schottkyGroup%Chi(ramInSchM, b=1.0_p)

    ! Check if ramification points in algebraic model for obtained Schottky group coincide
    ! with array [0.1,0.2,1.0,1.3] passed to PoincareUniformization.
    print *, 'Example2. Ramification pts (in algebraic model):', real(ramInAlgM, kind=4)
end block Example2

end program Examples
```

Output:

```Fortran
Example1. LHS and RHS of reciprocity law:          0.8132047   0.8132047
Example2. Ramification pts (in algebraic model):   0.1000000   0.2000000   1.000000   1.300000
```

[^2]: S.Yu. Lyamaev, [Summation of Poincaré theta-series in Schottky model](https://link.springer.com/article/10.1134/S0965542522070053), *Computational Mathematics and Mathematical Physics,* **62**:7 (2022), 1059–1073.

[^4]: A. B. Bogatyrev, [Chebyshev representation for rational functions](https://doi.org/10.1070/SM2010v201n11ABEH004123), *Sb. Math.*, **201**:11 (2010), 1579–1598.

[^5]: A.B. Bogatyrev, S.A. Goreinov, S.Yu. Lyamaev, [Analytical approach to multiband filter synthesis and comparison to other approaches](https://link.springer.com/article/10.1134/S0032946017030073), *Problems Inform. Transmission*, **53**:3 (2017), 260–273.