# SchottkyLib

SchottkyLib is a tool for computation of analytical objects (meromorphic functions, differentials of different orders, spinors, etc) on real hyperelliptic curves using Schottky uniformization. 

*Real hyperelliptic curve* (over the complex numbers) is a compact curve with affine part given by the equation $y^2 = P(x),~(x,y)\in\mathbb{C}^2$, where $P(x)$ is a real polynomial without multiple roots. Computations with such curves and their moduli space arise in numerous applications, such as Chebyshev optimization (computation of multiband filters using Chebyshev ansatz, computation of Chebyshev polynomials for several intervals), finding algebro-geometric solutions of nonlinear partial differential equations, description of magnetic states in planar magnetic nanoelements, simulation of the water flow under a stepped dam, and others. 

In *Schottky model*, a real hyperelliptic curve is represented as an orbit manifold of action of a suitable Schottky group, analytical objects on the curve can then be explicitly represented in terms of Poincaré $\theta$-series. Schottky model is free of the complicated modular dependencies inherent in other approaches to computations with such curves (Riemann $\theta$-functions, Weierstrass $\sigma$-functions) and allows one to deal with curves of large genera (up to 100 and larger).

Modules of SchottkyLib (folder `./src`):

- `Precision` — Set precision of floating-point numbers (real and complex): single, double or quadruple.
- `SchottkyGroup` — A class `SchottkyGroup_Type` is defined which contains main parameters of Schottky group: genus of associated real hyperelliptic curve and its Schottky moduli. Methods of this class are generators of Schottky group and their inverses, as well as indicator function of moduli space.
- `CayleyTreeTraversal` — Provides a universal interface for computation of sums and products over Schottky group and its quotients. Two traversal algorithms are present: Bogatyrev’s algorithm and a new faster algorithm[^2] proposed by the author.
- `ChildrenOrder` — Auxiliary computations for the new algorithm: estimate needed in the new traversal algorithm and order of children nodes by seniority for each parent.
- `DiffAndFuncsOnCurve` — Contains `SchottkyNumerics_Type`, the main class of SchottkyLib, which extends `SchottkyGroup_Type`. Added methods are linear and quadratic differentials (all needed in Chebyshev ansatz), Schottky functions, period matrix, and projection from Schottky model to algebraic model. If some analytical object you need is missing in current implementation it is not difficult to add it using the universal traversal interface provided by `CayleyTreeTraversal` module.
- `Uniformization` — Implementation of the classical algorithm (A. Poincaré, 1884) of numerical Schottky uniformization for M-curves: translation of an M-curve given in algebraic model to Schottky model.

Folder `./examples` contains two short examples to try the library out:

- `Example1` — Check numerically that reciprocity law for Abelian differentials holds: $\int_{b_j}\eta_{zw}=\int_w^z\zeta_j$, where $\zeta_j$ and $\eta_{zw}$ are the normalized Abelian differentials of the 1st and the 3rd kinds, respectively.
- `Example2` — Construct Schottky group associated with an M-curve given in algebraic model and check that projections of ramification points in Schottky model coincide with ramification points in algebraic model.

SchottkyLib is a part of a bigger package developed by the author during his PhD research. The bigger library is an implementation of Chebyshev ansatz[^4][^5] for practical synthesis of quasi-optimal multiband filters (of different types: digital, analog and microwave) and also contains evaluation of gradients of Abelian integrals and their periods by Schottky moduli using Hejhal’s approach. 

SchottkyLib is written in Fortran 2018 and is tested to compile and run successfully with gfortran 12.2.0 and ifort 2022.2.1. 

[^2]: S.Yu. Lyamaev, [Summation of Poincaré theta-series in Schottky model](https://link.springer.com/article/10.1134/S0965542522070053), *Computational Mathematics and Mathematical Physics,* **62**:7 (2022), 1059–1073.

[^4]: A. B. Bogatyrev, [Chebyshev representation for rational functions](https://doi.org/10.1070/SM2010v201n11ABEH004123), *Sb. Math.*, **201**:11 (2010), 1579–1598.

[^5]: A.B. Bogatyrev, S.A. Goreinov, S.Yu. Lyamaev, [Analytical approach to multiband filter synthesis and comparison to other approaches](https://link.springer.com/article/10.1134/S0032946017030073), *Problems Inform. Transmission*, **53**:3 (2017), 260–273.