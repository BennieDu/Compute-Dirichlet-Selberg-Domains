# Compute Dirichlet-Selberg domains in the symmetric space $SL(3,\mathbb{R})/SO(3)$

## Background

This program is motivated by the Poincaré's Algorithm for the Lie group $SL(n,\mathbb{R})$ of non-compact type, aiming to determine if a given subgroup is discrete and obtain a group presentation.

Specifically, the corresponding symmetric space $\mathcal{X}_n = SL(n,\mathbb{R})/SO(n)$ is realized through a projective model:

$$\mathcal{X}_n = \{X\in \mathbf{P}(Sym_n(\mathbb{R}))\mid X>0\},$$

and the entries of $X$ compose a projective coordinate of the point. Furthermore, the Selberg's two-point invariant,

$$ s(X,Y) = \mathrm{tr}(X^{-1}Y),$$

is a linear function in the coordinates of $Y$. Hence, for a discrete subgroup $\Gamma<SL(n,\mathbb{R})$, the Dirichlet-Selberg domain

$$ DS(X,\Gamma) = \{Y\in \mathcal{X}_n\mid s(X,Y)\leq s(g.X,Y),\ \forall g\in\Gamma\},$$

has a projective polytope structure.

The Poincaré's Algorithm decides the discreteness and obtains the presentation by contructing the Dirichlet-Selberg domain. More accurately, the algorithm tries to compute the polytope structure of the Dirichlet-Selberg domain $$DS(X,\Gamma_0)$$ for a finite subset $\Gamma_0\subset\Gamma$, then determine if it equals to the actual Dirichlet-Selberg domain by checking certain conditions and applying Poincaré's Fundamental Polyhedron Theorem:
