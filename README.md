# Compute Dirichlet-Selberg domains in $SL(3,\mathbb{R})/SO(3)$

## Table of Contents

- [Background](#background)
- [Usage](#usage)
- [Quick Start](#quick-start)
- [API Reference](#api-reference)
  - [Data Classes](#data-classes)
  - [Core Solvers](#core-solvers)
  - [Utility Functions](#utility-functions)
- [Examples](#examples)
  - [Elementary Subgroups](#elementary-subgroups)
  - [Lattices](#lattices)
  - [$SL(2,\mathbb{R})$ Subgroups](#sl2mathbbr-subgroups)
- [License](#license)
- [Contact & Acknowledgments](#contactacknowledgments)

## Background

We implemented an generalized version of Poincaré's Algorithm[^Kap23][^Du24] to decide the discreteness and compute a presentation of subgroups $\Gamma < SL(n,\mathbb{R})$ by constructing a Dirichlet–Selberg domain in the symmetric space

$$\mathcal{X}_n = SL(n,\mathbb{R})/SO(n),$$

realized projectively as a region of positive definite symmetric matrices[^Ebe96]

$$\lbrace X\in \mathbf{P}(Sym_n(\mathbb{R}))\mid X>0\rbrace.$$

The key is the Selberg's two‑point invariant[^Sel62],

$$ s(X,Y) = \mathrm{tr}(X^{-1}Y), $$

that turns bisectors (and half-spaces) into linear equations (inequalities, respectively) in the entries of $Y$, so the Dirichlet-Selberg domain

$$ DS(X,\Gamma_0) = \lbrace Y\in \mathcal{X}_n\mid s(X,Y)\leq s(g.X,Y),\ \forall g\in\Gamma_0\rbrace $$

is a convex projective polytope for any finite $\Gamma_0\subset\Gamma$.

Our algorithm, for $n = 3$, generalizes the original Poincaré's Algorithm[^Ril83][^EP94] and proceeds by:

- **Incremental word search**: build

$$ \Gamma_l = \lbrace g\in \Gamma\mid \lvert g\rvert \leq l\rbrace. $$
  
- **Polytope construction**: compute the face poset of $DS(X,\Gamma_l)$.
- **Conditions checks**: verify exactness, angle sum condition[^Rat94] and generator recovery[^Ril83].
- **Iteration**: if any test fails, increment $l$ and repeat.

When all checks pass, $DS(X,\Gamma_l)$ is a fundamental domain and $\Gamma$ is discrete with a finite presentation.

[^Kap23]: Michael Kapovich. Geometric algorithms for discreteness and faithfulness. In *Computational Aspects of Discrete Subgroups of Lie Groups*, Contemporary Mathematics, pages 87–112. AMS, 2023.
[^Du24]: Yukun Du. Geometry of Selberg’s bisectors in the symmetric space $SL(n,R)/SO(n,R)$. *J. Lond. Math. Soc.*, 110(4), Oct 2024.
[^Ebe96]: Patrick Eberlein. *Geometry of nonpositively curved manifolds*. University of Chicago Press, 1996.
[^Sel62]: Atle Selberg. On discontinuous groups in higher-dimensional symmetric spaces. *Matematika*, 6(3):3–16, 1962.
[^Rat94]: John G. Ratcliffe. *Foundations of hyperbolic manifolds*, volume 149. Springer, 1994.
[^Ril83]: Robert Riley. Applications of a computer implementation of Poincar´e’s theorem on fundamental polyhedra. *Math. Comput.*, 40(162):607–632, 1983.
[^EP94]: David B. A. Epstein and Carlo Petronio. An exposition of Poincar´e’s polyhedron theorem. *Enseign. Math.*, 40(1-2):113–170, 1994.

## Usage

### Install Dependencies

Ensure you have a Python 3 Jupyter environment, then install all required packages:

<pre markdown>
  pip install numpy sympy scipy cvxpy tqdm ipywidgets widgetsnbextension
</pre>

### Launch the Notebook

<pre markdown>
  jupyter notebook Dirichlet_Selberg.ipynb
</pre>

### Execute Cells Sequentially

In the notebook, run each cell in order to ensure all dependencies and definitions are loaded:

- Import the libraries
- Configuration of the tolerance
- Data classes
- Generic helpers
- Question-specific helpers
- Special utilities
- Core solvers

### Now you are ready to explore the program!

## Quick Start
Once you’ve run all prerequisite cells, copy & paste the following into a new cell and execute:
<pre markdown>
  generators = [np.diag([3,2,1/6])]
  center = np.eye(3)
  wbs, faces = compute_selberg_domain(generators, center, 5)
  for ind, face in enumerate(faces):
      print(f"Face #{ind}:")
      print(f"  Codimension:              {face.codim}")
      print(f"  Defining equations (ids): {face.equs}")
      print(f"  Subfaces (ids):           {face.subfaces}")
      print()
</pre>

Expected output (for this cyclic example):

<pre markdown>
  Face #0:
    Codimension:              1
    Defining equations (ids): [0]
    Subfaces (ids):           []
  
  Face #1:
    Codimension:              1
    Defining equations (ids): [1]
    Subfaces (ids):           []
  
  Face #2:
    Codimension:              0
    Defining equations (ids): []
    Subfaces (ids):           [0, 1]
</pre>

This confirms that, for

$$ X = \mathrm{diag}(1,1,1),\ g = \mathrm{diag}(3,2,\frac{1}{6}), $$

the Dirichlet-Selberg domain is the **dihedron** bounded by the two bisectors $Bis(X,g.X)$ and $Bis(X,g^{-1}.X)$.

## API Reference

### Data Classes

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

### Core Solvers

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
  - `pairing_dict`: A dictionary whose keys are indices (in `my_face_list`) of the facets, and the values are the incides of the corresponding paired facets. It is `None` if the given domain is not exact.

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

### Utility Functions

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
- Returns:
  - `sample_point`: A `numpy.array` matrix representing a sample point in the intersection $\cap A_i^\perp$, or `None` if the intersection is empty.

#### `selberg_domain_add_facet`

- Description: Intersect the given (pre-)Dirichlet-Selberg domain with a new half space and obtain the polyhedral structure of the resulting domain.
- Parameters:
  - `my_wbs`: A list of bisectors with facet pairings bounding the Dirichlet-Selberg domain.
  - `my_face_list`: A list of faces describing the polyhedral structure of the Dirichlet-Selberg domain.
  - `new_wb`: A `numpy.array` matrix $A$ representing the normal matrix of the new half space $\lbrace X\mid \mathrm{tr}(AX)\geq 0\rbrace$.
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

Consider the center

<pre markdown>
  [2 3 -2]
  [3 5 -3]
  [-2 -3 3]
</pre>

and the group $\Gamma$ generated by two elements

<pre markdown>
  [2 1 1]  [2 0 -1]
  [1 2 0]  [0 1 1]
  [1 0 1]  [-1 1 2]
</pre>

This is a two-generated Abelian subgroup of $SL(3,\mathbb{Z})$. To compute its Dirichlet-Selberg domain, run the following script in a new cell:

<pre markdown>
  generators = [np.array([[2, 1, 1],
              [1, 2, 0],
              [1, 0, 1]]),
              np.array([[2, 0, -1],
              [0, 1, 1],
              [-1, 1, 2]])]
  center = np.array([[2, 3, -2],
              [3, 5, -3],
              [-2, -3, 3]])
  my_wbs, my_face_list = compute_selberg_domain(generators, 1, 3, 20, center)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in my_face_list if face.codim == 1))
  print("Number of ridges:", sum(1 for face in my_face_list if face.codim == 2))
  print("Number of peaks:", sum(1 for face in my_face_list if face.codim == 3))
  print("Number of edges:", sum(1 for face in my_face_list if face.codim == 4))
  print("Number of vertices:", sum(1 for face in my_face_list if face.codim == 5))
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  Number of facets: 10
  Number of ridges: 18
  Number of peaks: 8
  Number of edges: 0
  Number of vertices: 0
  ```
</details>

Verify the exactness and angle-sum condition:

<pre markdown>
  is_exact, original_facets, paired_facets = polytope_is_exact(my_wbs, my_face_list)
  if not is_exact:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycle(my_wbs, my_face_list)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          pairings = [my_wbs[wb].word for wb in ridge_cycle.pairing]
          my_angle_sum = angle_sum(my_wbs, my_face_list, ridge_cycle)
          if my_angle_sum == None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", my_angle_sum)
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  The indices of ridges in the 0 th cycle: [8, 23, 21]
  The angle sum divisor for the 0 th ridge cycle equals 1
  The indices of ridges in the 1 th cycle: [10, 13, 16]
  The angle sum divisor for the 1 th ridge cycle equals 1
  The indices of ridges in the 2 th cycle: [12, 18, 25]
  The angle sum divisor for the 2 th ridge cycle equals 1
  The indices of ridges in the 3 th cycle: [15, 9, 14]
  The angle sum divisor for the 3 th ridge cycle equals 1
  The indices of ridges in the 4 th cycle: [20, 11, 22]
  The angle sum divisor for the 4 th ridge cycle equals 1
  The indices of ridges in the 5 th cycle: [24, 17, 19]
  The angle sum divisor for the 5 th ridge cycle equals 1
  ```
</details>

Verifying that all generators are recovered by the facet pairings:

<pre markdown>
  for ind, matrix in enumerate(generators):
      if word_is_recovered(my_wbs, my_face_list, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", path_word(my_wbs, my_face_list, matrix.T @ center @ matrix))
      else:
          print("the", ind, "th generator is not recovered by the facet pairings." )
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  the 0 th generator is recovered by the product of facet pairings with indices: [2]
  the 1 th generator is recovered by the product of facet pairings with indices: [3]
  ```
</details>

### Lattices

Consider the congruence subgroup of level $2$,

$$ \Gamma = \langle I + 2\mathbf{e}_i\otimes 2\mathbf{e}_j,\ i\neq j\rangle. $$

To compute its Dirichlet-Selberg domain, run the following script in a new cell:

<pre markdown>
  generators = []
  for i, j in itertools.permutations(range(3), 2):
      generator = np.eye(3)
      generator[i, j] = 2
      generators.append(generator)
  center = np.eye(3)
  my_wbs, my_face_list = compute_selberg_domain(generators, 1, 2, 20, center)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in my_face_list if face.codim == 1))
  print("Number of ridges:", sum(1 for face in my_face_list if face.codim == 2))
  print("Number of peaks:", sum(1 for face in my_face_list if face.codim == 3))
  print("Number of edges:", sum(1 for face in my_face_list if face.codim == 4))
  print("Number of vertices:", sum(1 for face in my_face_list if face.codim == 5))
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  Number of facets: 24
  Number of ridges: 84
  Number of peaks: 96
  Number of edges: 0
  Number of vertices: 0
  ```
</details>

Indeed, all the thirteen vertices of the domain are Satake point of rank $1$, hence no vertices or edges exist in the interior of $\mathcal{X}_3$.

Verify the exactness and angle-sum condition:

<pre markdown>
  is_exact, original_facets, paired_facets = polytope_is_exact(my_wbs, my_face_list)
  if not is_exact:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycle(my_wbs, my_face_list)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          pairings = [my_wbs[wb].word for wb in ridge_cycle.pairing]
          my_angle_sum = angle_sum(my_wbs, my_face_list, ridge_cycle)
          if my_angle_sum == None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", my_angle_sum)
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  The indices of ridges in the 0 th cycle: [96, 102, 106, 104]
  The angle sum divisor for the 0 th ridge cycle equals 1
  The indices of ridges in the 1 th cycle: [100, 154, 136]
  The angle sum divisor for the 1 th ridge cycle equals 1
  The indices of ridges in the 2 th cycle: [116, 122, 168]
  The angle sum divisor for the 2 th ridge cycle equals 1
  The indices of ridges in the 3 th cycle: [124, 129, 143, 169]
  The angle sum divisor for the 3 th ridge cycle equals 1
  The indices of ridges in the 4 th cycle: [133, 118, 153]
  The angle sum divisor for the 4 th ridge cycle equals 1
  The indices of ridges in the 5 th cycle: [165, 103, 123]
  The angle sum divisor for the 5 th ridge cycle equals 1
  The indices of ridges in the 6 th cycle: [97, 144, 131]
  The angle sum divisor for the 6 th ridge cycle equals 1
  The indices of ridges in the 7 th cycle: [109, 126, 172]
  The angle sum divisor for the 7 th ridge cycle equals 1
  The indices of ridges in the 8 th cycle: [120, 134, 155, 166]
  The angle sum divisor for the 8 th ridge cycle equals 1
  The indices of ridges in the 9 th cycle: [128, 111, 142]
  The angle sum divisor for the 9 th ridge cycle equals 1
  The indices of ridges in the 10 th cycle: [170, 105, 127]
  The angle sum divisor for the 10 th ridge cycle equals 1
  The indices of ridges in the 11 th cycle: [98, 107, 112, 110]
  The angle sum divisor for the 11 th ridge cycle equals 1
  The indices of ridges in the 12 th cycle: [99, 160, 149]
  The angle sum divisor for the 12 th ridge cycle equals 1
  The indices of ridges in the 13 th cycle: [113, 140, 176]
  The angle sum divisor for the 13 th ridge cycle equals 1
  The indices of ridges in the 14 th cycle: [125, 145, 130, 171]
  The angle sum divisor for the 14 th ridge cycle equals 1
  The indices of ridges in the 15 th cycle: [147, 115, 159]
  The angle sum divisor for the 15 th ridge cycle equals 1
  The indices of ridges in the 16 th cycle: [174, 108, 141]
  The angle sum divisor for the 16 th ridge cycle equals 1
  The indices of ridges in the 17 th cycle: [138, 148, 161, 175]
  The angle sum divisor for the 17 th ridge cycle equals 1
  The indices of ridges in the 18 th cycle: [101, 114, 119, 117]
  The angle sum divisor for the 18 th ridge cycle equals 1
  The indices of ridges in the 19 th cycle: [121, 156, 135, 167]
  The angle sum divisor for the 19 th ridge cycle equals 1
  The indices of ridges in the 20 th cycle: [139, 162, 150, 177]
  The angle sum divisor for the 20 th ridge cycle equals 1
  The indices of ridges in the 21 th cycle: [132, 178, 146]
  The angle sum divisor for the 21 th ridge cycle equals 2
  The indices of ridges in the 22 th cycle: [163, 173, 151]
  The angle sum divisor for the 22 th ridge cycle equals 2
  The indices of ridges in the 23 th cycle: [137, 179, 157]
  The angle sum divisor for the 23 th ridge cycle equals 2
  The indices of ridges in the 24 th cycle: [164, 158, 152]
  The angle sum divisor for the 24 th ridge cycle equals 2
  ```
</details>

The output implying that $21$ of the ridge cycles satisfy the angle sum condition with angle sum $2\pi$, and $4$ of them satisfy the condition with angle sum $\pi$.

Verifying that all generators are recovered by the facet pairings:

<pre markdown>
  for ind, matrix in enumerate(generators):
      if word_is_recovered(my_wbs, my_face_list, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", path_word(my_wbs, my_face_list, matrix.T @ center @ matrix))
      else:
          print("the", ind, "th generator is not recovered by the facet pairings." )
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  the 0 th generator is recovered by the product of facet pairings with indices: [6]
  the 1 th generator is recovered by the product of facet pairings with indices: [7]
  the 2 th generator is recovered by the product of facet pairings with indices: [8]
  the 3 th generator is recovered by the product of facet pairings with indices: [9]
  the 4 th generator is recovered by the product of facet pairings with indices: [10]
  the 5 th generator is recovered by the product of facet pairings with indices: [11]
  ```
</details>

### $SL(2,\mathbb{R})$ Subgroups

We can recognize discrete subgroups of $SL(2,\mathbb{R})$ as $SL(3,\mathbb{R})$-subgroups following the block‐diagonal embedding,

$$ SL(2,\mathbb{R})\hookrightarrow SL(3,\mathbb{R}),\ g\mapsto \mathrm{diag}(g,1). $$

The Dirichlet-Selberg domains (with block-diagonal centers) of the images of such embeddings are equivalent to the product of $\mathbb{R}^3$ with the Dirichlet domain of the corresponding $SL(2,\mathbb{R})$-subgroups. Their polyhedral structures are also computable using our program. For example the fundamental group of the double-torus:

<pre markdown>
  generators = []
  for i, j in itertools.product([-1,1], repeat=2):
      generator = np.array([[(np.sqrt(2) + 2 + i*(2 ** (3/4)))/2, (j*(np.sqrt(2) + 2) - i*j*(2 ** (3/4)  + 2 ** (5/4)))/2, 0],
                      [(-j*(np.sqrt(2) + 2) - i*j*(2 ** (3/4)  + 2 ** (5/4)))/2, (np.sqrt(2) + 2 - i*(2 ** (3/4)))/2, 0],
                      [0, 0, 1]])
      generators.append(generator)
  center = np.eye(3)
  my_wbs, my_face_list = compute_selberg_domain(generators, 1, 1, 20, center)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in my_face_list if face.codim == 1))
  print("Number of ridges:", sum(1 for face in my_face_list if face.codim == 2))
  print("Number of peaks:", sum(1 for face in my_face_list if face.codim == 3))
  print("Number of edges:", sum(1 for face in my_face_list if face.codim == 4))
  print("Number of vertices:", sum(1 for face in my_face_list if face.codim == 5))
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  Number of facets: 8
  Number of ridges: 8
  Number of peaks: 0
  Number of edges: 0
  Number of vertices: 0
  ```
</details>

This corresponds to an octagon, the standard fundamental domain of the double-torus group.

Verify the exactness and angle-sum condition:

<pre markdown>
  is_exact, original_facets, paired_facets = polytope_is_exact(my_wbs, my_face_list)
  if not is_exact:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycle(my_wbs, my_face_list)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          pairings = [my_wbs[wb].word for wb in ridge_cycle.pairing]
          my_angle_sum = angle_sum(my_wbs, my_face_list, ridge_cycle)
          if my_angle_sum == None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", my_angle_sum)
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  The indices of ridges in the 0 th cycle: [1, 2, 4, 3, 0, 5, 7, 6]
  The angle sum divisor for the 0 th ridge cycle equals 1
  ```
</details>

This implies that the angle sum of the unique ridge-cycle is $2\pi$, also agrees with the well-known fact for the hyperbolic double-torus.

Verifying that all generators are recovered by the facet pairings:

<pre markdown>
  for ind, matrix in enumerate(generators):
      if word_is_recovered(my_wbs, my_face_list, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", path_word(my_wbs, my_face_list, matrix.T @ center @ matrix))
      else:
          print("the", ind, "th generator is not recovered by the facet pairings." )
</pre>

<details>
  <summary><strong>Expected Output (click to expand)</strong></summary>

  ```text
  the 0 th generator is recovered by the product of facet pairings with indices: [4]
  the 1 th generator is recovered by the product of facet pairings with indices: [5]
  the 2 th generator is recovered by the product of facet pairings with indices: [6]
  the 3 th generator is recovered by the product of facet pairings with indices: [7]
  ```
</details>

## License

This project is licensed under the [MIT License](https://opensource.org/licenses/MIT).

## Contact & Acknowledgments

- Yukun Du, Postdoc Researcher, University of Georgia
  - Email: [yukun.du@uga.edu](mailto:yukun.du@uga.edu)
  - Website (departmental): [Yukun Du](https://math.franklin.uga.edu/directory/people/yukun-du)
- I would express my deepest gratitude to my Ph.D. advisor, Michael Kapovich, for introducing me to the area of Geometric Group Theory and for his invaluable mentoring throughout my graduate studies.
