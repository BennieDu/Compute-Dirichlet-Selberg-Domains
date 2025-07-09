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

- Assume that a subgroup $\Gamma<SL(n,\mathbb{R})$ is given by generators $g_1,\dots,g_m$, with relators initially unknown. We begin by selecting a point $X\in\mathcal{X}_n$, setting $l = 1$, and computing the finite subset $\Gamma_l\subset \Gamma$, which consists of elements represented by words of length $\leq l$ in the letters $g_i$ and $g_i^{-1}$.
- Compute the face poset of the Dirichlet-Selberg domain $DS(X,\Gamma_l)$, which forms a finitely-sided polytope in $\mathcal{X}_n$.
- Utilizing this face poset data, check if $DS(X,\Gamma_l)$ satisfies the following conditions:
  - Verify that $DS(X,\Gamma_l)$ is an exact convex polytope. For each $w\in \Gamma_l$, confirm that the isometry $w$ pairs the two facets contained in $\mathrm{Bis}(X,w.X)$ and $\mathrm{Bis}(X,w^{-1}.X)$, provided these facets exist.
  - Verify that $D(X,\Gamma_l)$ satisfies the angle sum condition for each ridge cycle.
  - Verify that each element $g_i$ can be expressed as a product of the facet pairings of $DS(X,\Gamma_l)$.
- If any of these conditions are not met, increment $l$ by $1$ and repeat the initialization, computation, and verification processes.
- If all conditions are satisfied, by Poincar\'e's Fundamental Polyhedron Theorem, $DS(X,\Gamma_l)$ is a fundamental domain for $\Gamma$, and $\Gamma$ is geometrically finite. Specifically, $\Gamma$ is discrete and has a finite presentation derived from the ridge cycles of $DS(X,\Gamma_l)$.

We implement all steps of the algorithm in this program for the $SL(3,\mathbb{R})$ case. Specifically:

- Given a center and generators in $SL(3,R)$ (both as numpy.array matrices) and maximum length of words, compute the polytope structure of the Dirichlet-Selberg domain.
- Given polytope data for a Dirichlet-Selberg domain, check if it is exact.
- Given polytope data for a Dirichlet-Selberg domain and knowing exactness, compute the ridge cycles and corresponding angle sums.
- Given polytope data for a Dirichlet-Selberg domain, check if a given word (especially a generator) can be recovered by the Dirichlet-Selberg facet pairings.

## Quick Start

## Usage/API

## Examples

## License

## Contact/Acknowledgments
