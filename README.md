# SchottkyLib

SchottkyLib is a tool for computation of analytical objects (meromorphic functions, differentials of different orders, spinors, etc) on real hyperelliptic curves using Schottky model. 

*Real hyperelliptic curve* (over the complex numbers) is a compact curve with affine part given by the equation $y^2 = P(x),~(x,y)\in\mathbb{C}^2$, where $P(x)$ is a real polynomial without multiple roots. Computations with such curves and their moduli space arise in numerous applications, such as Chebyshev optimization (computation of multiband filters using Chebyshev ansatz, computation of Chebyshev polynomials for several intervals), finding algebro-geometric solutions of nonlinear partial differential equations, description of magnetic states in planar magnetic nanoelements, simulation of the water flow under a stepped dam, and others. 

In *Schottky model*, a real hyperelliptic curve is represented as an orbit manifold of action of a suitable Schottky group. Analytical objects on the curve can then be explicitly represented in terms of Poincaré $\theta$-series. Schottky model is free of the complicated modular dependencies inherent in other approaches to numerical computations on curves (Riemann $\theta$-functions, Weierstrass $\sigma$-functions) and allows one to deal with curves of large genera (up to 100 and even larger).
