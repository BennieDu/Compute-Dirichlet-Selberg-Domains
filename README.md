# Compute Dirichlet-Selberg domains in the symmetric space $SL(3,\mathbb{R})/SO(3)$

## Table of Contents

- [Background](#background)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
  - [Important Data Classes](#important-data-classes)
  - [Core Solver Functions](#core-solver-functions)
  - [Other Important Utility Functions](#other-important-utility-functions)
- [Examples](#examples)
  - [Elementary Subgroups](#elementary-subgroups)
  - [Lattices](#lattices)
  - [Subgroups of $SL(2,\mathbb{R})$](#subgroups-of-sl2mathbbr)
- [License](#license)
- [Contact/Acknowledgments](#contactacknowledgments)

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

### Important Data Classes

#### `Word_Bis`

<pre markdown>
  @dataclass
  class Word_Bis:
      word: np.ndarray
      bis: np.ndarray
</pre>

- Description: A bisector $Bis(X,g^{-1}.X)\subset \mathcal{X}_3$ along with the isometry $g\in SL(3,\mathbb{R})$. A list `my_wbs` of `Word_Bis` elements describes the hyperplanes defining a convex polytope in $\mathcal{X}_3$.
- Attributes:
  - `word`: The 3*3 numpy matrix of the isometry $g\in SL(3,\mathbb{R})$.
  - `bis`: The 3*3 numpy matrix of the normal matrix for the bisector: $Bis(X,g^{-1}.X) = A^\perp = \lbrace Y\mid \mathrm{tr}(AY) = 0\rbrace$.

#### `Poly_Face`

<pre markdown>
  @dataclass
  class Poly_Face:
      equs: list[int]
      codim: int
      subfaces: list[int]
      sample_point: np.ndarray
</pre>

- Description: A face of a convex polytope in $\mathcal{X}_3$. A list `my_face_list` along with the list `my_wbs` describes the polytope structure of a convex polytope in $\mathcal{X}_3$.
- Attributes:
  - `equs`: A list of indices of equations defining the minimal plane in $\mathcal{X}_3$ that contains the face. Equivalently, the indices `i` such that `my_wbs[i].bis` is normal to the face.
  - `codim`: The codimension of the face.
  - `subfaces`: The indices of proper subfaces of the current face. Equivalently, the indices `i` such that `my_face_list[i]` represents a proper subface.
  - `sample_point`: A 3*3 numpy matrix representing a point lying in the face interior.

#### `Ridge_Cycle`

<pre markdown>
  @dataclass
  class Ridge_Cycles:
      ridge: list[int]
      pairing: list[int]
</pre>

- Description: A ridge cycle $[r] = \lbrace r_0,r_1,\dots,r_{m-1}\rbrace$ of the convex polytope (described by lists such as `my_face_list` and `my_wbs`).
- Attributes:
  - `ridge`: the indices of the ridges (as in `my_face_list`) in $[r]$, following the order. That is, the $i$-th index is for the ridge $r_i$.
  - `pairing`: the indices of the words taking a ridge to the next one. That is to say, `my_wbs[i].word` takes $r_i$ to $r_{i+1}$ (indices realized in modulo $m$).

### Core Solver Functions

#### `compute_selberg_domain`

- Description: Compute the polytope structure of the (pre-)Dirichlet-Selberg domain $DS(X,\Gamma_0)$ from generators and center.
- Parameters:
  - `generators`: Isometries $g_1,\dots, g_k\in SL(3,\mathbb{R})$ (as `numpy.array` matrices) that generates the subgroup $\Gamma$.
  - `length_1`: Maximal length of words in generators $g_1,\dots, g_k$ that will be added to the subset $\Gamma_0\subset \Gamma$.
  - `length_2`: Maximal length of words, typically greater than `length_1` that may be added to $\Gamma_0$ to eliminate unpaired ridges.
  - `loop_times`: Number of loops to pick a word up to length `length_2` and add it to $\Gamma_0$.
  - `center`: The center $X\in\mathcal{X}_3$ (as a numpy matrix) of the Dirichlet-Selberg domain.
- Returns:
  - `my_wbs`: A list of data class `Word_Bis` elements, each describing a bisector that bounds the polytope $DS(X,\Gamma_0)$.
  - `my_face_list`: A list of data class `Poly_Face` elements, each describing a face of the polytope $DS(X,\Gamma_0)$. The two lists together provide a thorough description of the polyhedral structure of the Dirichlet-Selberg domain.

#### `polytope_is_exact`

- Description: Decide if a given (pre-)Dirichlet-Selberg domain is exact with respect to the canonical facet pairings.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
- Returns:
  - `is_exact`: A boolean variable deciding if the given domain is exact.
  - `facet_indices`: A list of indices (in `my_face_list`, in order) of the facets of the domain.
  - `paired_indices`: A list of indices of the facets so the corresponding indices in `facet_indices` and `paired_indices` correspond to paired facets.

#### `compute_ridge_cycle`

- Description: Compute all ridge cycles of a given exact (pre-)Dirichlet-Selberg domain.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
- Returns:
  - `ridge_cycle_list`: A list of `Ridge_Cycle` data class elements, each is a ridge cycle of the given domain.

#### `angle_sum`

- Description: Compute the quotient of $2\pi$ by the angle sum of a given ridge cycle in a given exact (pre-)Dirichlet-Selberg domain. The angle-sum condition is satisfied if this quotient is an integer.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `ridge_cycle`: A `Ridge_Cycle` element describing a ridge cycle of the given domain.
- Returns:
  - `angle_sum_quotient`: A natural number $k$ such that the angle sum of the given ridge cycle equals $2\pi/k$. Return `None` if no such natural number is satisfied.

#### `word_is_recovered`

- Description: Decide if an isometry in $SL(3,\mathbb{R})$ is generated by the facet pairings of a given Dirichlet-Selberg domain.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `word`: An isometry in $SL(3,\mathbb{R})$ as a `numpy.array` matrix.
- Returns:
  - `is_recovered`: Boolean variable determining if the given isometry is generated by the facet pairings of the given domain.

### Other Important Utility Functions

#### `word_bisectors`

- Description: Compute all pairs of bisector $\mathrm{Bis}(X,g^{-1}.X)$ and isometry $g$ (as a `Word.Bis` data class element) generated by given generators and up to a given word length.
- Parameters
  - `generators`: A list of isometries $g_1,\dots, g_k\in SL(3,\mathbb{R})$ as `numpy.array` matrices.
  - `length`: The maximum word length of output words.
  - `center`: The center $X\in \mathcal{X}_3$ as a `numpy.array` matrix.
- Returns:
  - `my_wbs`: A list of `Word.Bis` elements consisting of bisectors $\mathrm{Bis}(X,g^{-1}.X)$ and isometry $g$ for all words $g$ in given generators up to the given length.

#### `find_positive_definite_intersection`

- Description: Decide if the intersection $\cap A_i^\perp$ of given hyperplanes is non-empty; if so, obtain a sample point in this intersection plane.
- Parameters:
  - `words`: A list of indefinite `numpy.array` matrices $A_i$, each represents a hyperplane $A_i^\perp\subset \mathcal{X}_3$.
- Returns: a data class `Find_Intersection` consisting of the following attributes:
  - `sample_point`: A `numpy.array` matrix representing a sample point in the intersection $\cap A_i^\perp$.
  - `is_intersection`: A boolean deciding if the intersection is non-empty.

#### `selberg_domain_add_facet`

- Description: Intersect the given (pre-)Dirichlet-Selberg domain with a new half space and obtain the polyhedral structure of the resulting domain.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `new_wb`: A `numpy.array` matrix $A$ representing the normal matrix of the new half space $\lbraceX\mid \mathrm{tr}(AX)\geq 0\rbrace$.
- Returns:
  - `my_wbs`: A list of bisectors with facet pairings bounding the new domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the new domain.

#### `face_is_paired`

- Description: Decide if two certain faces in a given (pre-)Dirichlet-Selberg domain are isometrically paired by the given isometry in $SL(3,\mathbb{R})$.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `old_face_ind`: The index of the first face in `my_face_list`.
  - `new_face_ind`: The index of the second face in `my_face_list`.
  - `word`: An isometry in $SL(3,\mathbb{R})$ as a `numpy.array` matrix.
- Returns:
  - `is_paired`: Boolean variable showing if thw two faces are paired by the isometry.

#### `path_word`

- Description: Find the word in letters of the facet pairings that takes a given point in $\mathcal{X}_3$ into the given Dirichlet-Selberg domain.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `dest_point`: A point $X\in \mathcal{X}_3$ as a `numpy.array` matrix.
- Returns:
  - `my_path`: A list of indices $i_k$, such that the product of the isometries `my_wbs[i_k].word` takes $X$ into the given Dirichlet-Selberg domain.

## Examples

### Elementary Subgroups

### Lattices

### Subgroups of $SL(2,\mathbb{R})$

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

## Contact/Acknowledgments

- Yukun Du, Postdoc Researcher, University of Georgia
  - Email: [yukun.du@uga.edu](mailto:yukun.du@uga.edu)
  - Website (departmental): [Yukun Du](https://math.franklin.uga.edu/directory/people/yukun-du)
- I would express my deepest gratitude to my Ph.D. advisor, Michael Kapovich, for introducing me to the area of Geometric Group Theory and for his invaluable mentoring throughout my graduate studies.
