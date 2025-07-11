# Compute Dirichlet-Selberg domains in the symmetric space $SL(3,\mathbb{R})/SO(3)$

## Background

This program is motivated by the generalized Poincaré's Algorithm, aiming to determine if a given subgroup of the Lie group $SL(n,\mathbb{R})$ is discrete and obtain a group presentation[^Kap23][^Du24].

[^Kap23]: Michael Kapovich. Geometric algorithms for discreteness and faithfulness. In *Computational Aspects of Discrete Subgroups of Lie Groups*, Contemporary Mathematics, pages 87–112. AMS, 2023.
[^Du24]: Yukun Du. Geometry of Selberg’s bisectors in the symmetric space $SL(n,R)/SO(n,R)$. *J. Lond. Math. Soc.*, 110(4), Oct 2024.

Specifically, the corresponding symmetric space $\mathcal{X}_n = SL(n,\mathbb{R})/SO(n)$ is realized through a projective model[^Ebe96]:

[^Ebe96]: Patrick Eberlein. *Geometry of nonpositively curved manifolds*. University of Chicago Press, 1996.

$$\mathcal{X}_n = \lbrace X\in \mathbf{P}(Sym_n(\mathbb{R}))\mid X>0\rbrace ,$$

and the entries of $X$ compose a projective coordinate of the point. Furthermore, the Selberg's two-point invariant[^Sel62],

[^Sel62]: Atle Selberg. On discontinuous groups in higher-dimensional symmetric spaces. *Matematika*, 6(3):3–16, 1962.

$$ s(X,Y) = \mathrm{tr}(X^{-1}Y),$$

is a linear function in the coordinates of $Y$. Hence, for a discrete subgroup $\Gamma<SL(n,\mathbb{R})$, the Dirichlet-Selberg domain

$$ DS(X,\Gamma) = \lbrace Y\in \mathcal{X}_n\mid s(X,Y)\leq s(g.X,Y),\ \forall g\in\Gamma\rbrace ,$$

has a projective polytope structure.

Our algorithm generalizes the original Poincaré's Algorithm[^Ril83][^EP94], deciding the discreteness and obtains the presentation of a $SL(n,\mathbb{R})$ subgroup by contructing a Dirichlet-Selberg domain. More accurately, the algorithm tries to compute the polytope structure of the Dirichlet-Selberg domain $$DS(X,\Gamma_0)$$ for a finite subset $\Gamma_0\subset\Gamma$, then determine if it equals to the actual Dirichlet-Selberg domain by checking certain conditions and applying Poincaré's Fundamental Polyhedron Theorem[^Du24]:

- Assume that a subgroup $\Gamma<SL(n,\mathbb{R})$ is given by generators $g_1,\dots,g_m$, with relators initially unknown. We begin by selecting a point $X\in\mathcal{X}_n$, setting $l = 1$, and computing the finite subset $\Gamma_l\subset \Gamma$, which consists of elements represented by words of length $\leq l$ in the letters $g_i$ and $g_i^{-1}$.
- Compute the face poset of the Dirichlet-Selberg domain $DS(X,\Gamma_l)$, which forms a finitely-sided polytope in $\mathcal{X}_n$.
- Utilizing this face poset data, check if $DS(X,\Gamma_l)$ satisfies the following conditions:
  - Verify that $DS(X,\Gamma_l)$ is an exact convex polytope. For each $w\in \Gamma_l$, confirm that the isometry $w$ pairs the two facets contained in $\mathrm{Bis}(X,w.X)$ and $\mathrm{Bis}(X,w^{-1}.X)$, provided these facets exist.
  - Verify that $D(X,\Gamma_l)$ satisfies the angle sum condition for each ridge cycle[^Rat94].
  - Verify that each element $g_i$ can be expressed as a product of the facet pairings of $DS(X,\Gamma_l)$[^Ril83] .
- If any of these conditions are not met, increment $l$ by $1$ and repeat the initialization, computation, and verification processes.
- If all conditions are satisfied, by Poincar\'e's Fundamental Polyhedron Theorem, $DS(X,\Gamma_l)$ is a fundamental domain for $\Gamma$, and $\Gamma$ is geometrically finite. Specifically, $\Gamma$ is discrete and has a finite presentation derived from the ridge cycles of $DS(X,\Gamma_l)$.

[^Rat94]: John G. Ratcliffe. *Foundations of hyperbolic manifolds*, volume 149. Springer, 1994.
[^Ril83]: Robert Riley. Applications of a computer implementation of Poincar´e’s theorem on fundamental polyhedra. *Math. Comput.*, 40(162):607–632, 1983.
[^EP94]: David B. A. Epstein and Carlo Petronio. An exposition of Poincar´e’s polyhedron theorem. *Enseign. Math.*, 40(1-2):113–170, 1994.

We implement all steps of the algorithm in this program for the $SL(3,\mathbb{R})$ case. Specifically:

- Given a center and generators in $SL(3,R)$ (both as numpy.array matrices) and maximum length of words, compute the polytope structure of the Dirichlet-Selberg domain.
- Given polytope data for a Dirichlet-Selberg domain, check if it is exact.
- Given polytope data for a Dirichlet-Selberg domain and knowing exactness, compute the ridge cycles and corresponding angle sums.
- Given polytope data for a Dirichlet-Selberg domain, check if a given word (especially a generator) can be recovered by the Dirichlet-Selberg facet pairings.

## Usage

### Install dependencies

Make sure you have a Python 3 environment, then:

<pre markdown>
  pip install -r requirements.txt
</pre>

### Launch the notebook

<pre markdown>
  jupyter notebook Dirichlet_Selberg.ipynb
</pre>

### Run cells in sequence

In the notebook, execute the cells from top to bottom. In particular:

- Import the libraries
- Configuration of the tolerance
- Data classes
- Generic helper functions
- Question-specific helper functions
- Special Utilities
- Core solver functions

### Now you are ready to explore the program!

## Quick Start
After excuted the prerequisite cells, copy the following to a new cell and run:
<pre markdown>
  generators = [np.diag([3,2,1/6])]
  center = np.eye(3)
  my_wbs, my_face_list = compute_selberg_domain_short(generators, 5, center)
  for i in range(len(my_face_list)):
      print("Codimension of the", i, "th face:", my_face_list[i].codim)
      print("Indices of equations defining the", i, "th face's spanned plane:", my_face_list[i].equs)
      print("Indices of the", i, "th face's subfaces:", my_face_list[i].subfaces)
</pre>

You will expect output in the following form:

<pre markdown>
  Codimension of the 0 th face: 1
  Indices of equations defining the 0 th face's spanned plane: [0]
  Indices of the 0 th face's subfaces: []
  Codimension of the 1 th face: 1
  Indices of equations defining the 1 th face's spanned plane: [1]
  Indices of the 1 th face's subfaces: []
  Codimension of the 2 th face: 0
  Indices of equations defining the 2 th face's spanned plane: []
  Indices of the 2 th face's subfaces: [0, 1]
</pre>

This implies that the Dirichlet-Selberg domain for center $X = diag(1,1,1)$ and cyclic group generated by $g = diag(3,2,\frac{1}{6})$ is just a dihedron bounded by $Bis(X,g.X)$ and $Bis(X,g^{-1}.X)$.

## API Reference

### `Word_Bis`

<pre markdown>
  @dataclass
  class Word_Bis:
      word: np.ndarray
      bis: np.ndarray
</pre>

- Description: A bisector $Bis(X,g^{-1}.X)\subset \mathcal{X}_3$ along with the isometry $g\in SL(3,\mathbb{R})$.
- Attributes:
  - `word`: the 3*3 numpy matrix of the isometry $g\in SL(3,\mathbb{R})$.
  - `bis`: the 3*3 numpy matrix of the normal matrix for the bisector: $Bis(X,g^{-1}.X) = A^\perp = \lbrace Y\mid \mathrm{tr}(AY) = 0\rbrace$.

### `Poly_Face`

<pre markdown>
  @dataclass
  class Poly_Face:
      equs: list[int]
      codim: int
      subfaces: list[int]
      sample_point: np.ndarray
</pre>

- Description:
- Attributes:
  - `equs`:
  - `codim`:
  - `subfaces`:
  - `sample_point`:

### `Ridge_Cycles`

<pre markdown>
  @dataclass
  class Ridge_Cycles:
      ridge: list[int]
      pairing: list[int]
</pre>

- Description:
- Attributes:
  - `ridge`:
  - `pairing`:

## Examples

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

## Contact/Acknowledgments

- Yukun Du, Postdoc Researcher, University of Georgia
  - Email: [yukun.du@uga.edu](mailto:yukun.du@uga.edu)
  - Website (departmental): [Yukun Du](https://math.franklin.uga.edu/directory/people/yukun-du)
- I would express my deepest gratitude to my Ph.D. advisor, Michael Kapovich, for introducing me to the area of Geometric Group Theory and for his invaluable mentoring throughout my graduate studies.
