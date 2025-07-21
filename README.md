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
- [Contact & Acknowledgments](#contact--acknowledgments)

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
      word: np.ndarray   # g ∈ SL(3,ℝ)
      bis:  np.ndarray   # normal matrix A for Bis(X, g⁻¹·X)
</pre>

- What it does:
  Capture an isometry $g$ and its bisector.
#### `Poly_Face`

<pre markdown>
  @dataclass
  class Poly_Face:
      equs:         list[int]   # indices into `wbs` whose `bis` vanish on this face
      codim:        int         # codimension of the face (0 = interior, 1 = facet, 2 = ridge)
      subfaces:     list[int]   # indices into `faces` of strictly smaller faces contained in this face
      sample_point: np.ndarray  # a point in the interior of this face
</pre>

- What it does:
  Records one face of the polytope (facet, ridge, etc.)

#### `Ridge_Cycle`

<pre markdown>
  @dataclass
  class Ridge_Cycles:
      ridge:   list[int]    # ordered indices into `wbs` of ridges in the cycle
      pairing: list[int]    # for each ridge i, the index into `wbs` pairing it to the next one
</pre>

- What it does:
  Encodes one cycle of ridges and their facet‑pairings

### Core Solvers

#### `compute_selberg_domain`

<pre markdown>
  compute_selberg_domain(
      generators: list[np.ndarray],
      center:     np.ndarray,
      length_1:   int,
      length_2:   Optional[int] = None,
      loop_times: int           = 0
  ) -> tuple[list[Word_Bis], list[Poly_Face]]
</pre>

- What it does:
  Builds the (pre‑)Dirichlet–Selberg domain $DS(X,\Gamma_0)$ by
  - Enumerating words up to `length_1`,
  - Adding words up to `length_2` to pair unpaired ridges (`loop_time` rounds).
- Returns:
  - `wbs` - A list of `Word_Bis` bisectors
  - `faces` - A list of `Poly_Face` describing the resulting polytope

#### `polytope_exactness`

<pre markdown>
  polytope_exactness(
      wbs:   list[Word_Bis],
      faces: list[Poly_Face]
  ) -> Optional[dict[int, int]]
</pre>

- What it does:
  Checks that every facet in `faces` is paired canonically by a word in `wbs`.
- Returns:
  - `pairing_dict` - A dict `{"facets": [...], "paired_facets": [...]}` if exact, or `None` if any facet is unpaired

#### `compute_ridge_cycles`

<pre markdown>
  compute_ridge_cycles(
      wbs:   list[Word_Bis],
      faces: list[Poly_Face]
  ) -> Optional[list[Ridge_Cycle]]
</pre>

- What it does:
  For an exact domain, traces each cycle of ridges via facet pairings.
- Returns:
  - `ridge_cycles` - A list of `Ridge_Cycle` objects, or `None` if the domain is not exact.

#### `compute_angle_sum`

<pre markdown>
  compute_angle_sum(
      wbs:         list[Word_Bis],
      faces:       list[Poly_Face],
      ridge_cycle: Ridge_Cycle
  ) -> Optional[int]
</pre>

- What it does: Computes the integer $k$ such that the sum of dihedral angles around `ridge_cycle` is $2\pi/k$.

- Returns:
  - `angle_sum_quotient` - The natural number $k$ if the condition holds, or `None` otherwise.

#### `is_word_recovered`

<pre markdown>
  is_word_recovered(
      wbs:   list[Word_Bis],
      faces: list[Poly_Face],
      isom:  np.ndarray
  ) -> bool
</pre>

- What it does:
  Tests whether the isometry `isom` is realized by composing facet pairings of the domain.
- Returns:
  - `is_recovered` - `True` if recovered, `False` otherwise.

### Utility Functions

#### `compute_word_bisectors`

<pre markdown>
  compute_word_bisectors(
    gens:       list[np.ndarray],
    max_length: int,
    center:     np.ndarray
  ) -> list[WordBis]
</pre>

- What it does:
  Enumerates all words in `gens` (and inverses) up to `max_length`, computes associated bisectors at `center`.
- Returns:
  - `wbs` - A list of `Word_Bis` bisectors
#### `find_feasible_point`

<pre markdown>
  find_feasible_point(
      plane_eqs: list[np.ndarray]
  ) -> Optional[np.ndarray]
</pre>

- What it does:
  Solves a small CVXPY feasibility problem to find any point in $\cap A_i^\perp$ for $A_i$ given in `plane_eqs`.
- Returns:
  - `feasible_point` - A 3×3 `np.ndarray` if feasible, or `None` if the hyperplane intersection is empty

#### `add_facet_to_domain`

<pre markdown>
  add_facet_to_domain(
      wbs:      list[Word_Bis],
      faces:    list[Poly_Face],
      new_wb:   Word_Bis
  ) -> tuple[list[Word_Bis], list[Poly_Face]]
</pre>

- What it does:
  Intersects the current domain with the half‑space given by `new_wb`, updating both `wbs` and `faces`.

#### `are_faces_paired`

<pre markdown>
  are_faces_paired(
      wbs:    list[Word_Bis],
      faces:  list[Poly_Face],
      i_old:  int,
      i_new:  int,
      isom:   np.ndarray
  ) -> bool
</pre>

- What it does:
  Checks if the facets `faces[i_old]` and `faces[i_new]` are paired by the isometry `isom`.

#### `find_path_word`

<pre markdown>
  find_path_word(
      wbs:   list[Word_Bis],
      faces: list[Poly_Face],
      dest:  np.ndarray
  ) -> list[int]
</pre>

- What it does:
  Finds a sequence of facet‑pairing indices that carry the point `dest` into the domain.
- Returns:
  - `path_word` - A list of indices into `wbs` whose composed word maps `dest` inside

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
  wbs, faces = compute_selberg_domain(generators, center, 1, 3, 20)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in faces if face.codim == 1))
  print("Number of ridges:", sum(1 for face in faces if face.codim == 2))
  print("Number of peaks:", sum(1 for face in faces if face.codim == 3))
  print("Number of edges:", sum(1 for face in faces if face.codim == 4))
  print("Number of vertices:", sum(1 for face in faces if face.codim == 5))
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
  pairing_dict = polytope_exactness(wbs, faces)
  if pairing_dict is None:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycles(wbs, faces)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          angle_sum_quotient = angle_sum(wbs, faces, ridge_cycle)
          if angle_sum_quotient is None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", angle_sum_quotient)
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
      if is_word_recovered(wbs, faces, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", find_path_word(wbs, faces, matrix.T @ center @ matrix))
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
  wbs, faces = compute_selberg_domain(generators, center, 1, 2, 20)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in faces if face.codim == 1))
  print("Number of ridges:", sum(1 for face in faces if face.codim == 2))
  print("Number of peaks:", sum(1 for face in faces if face.codim == 3))
  print("Number of edges:", sum(1 for face in faces if face.codim == 4))
  print("Number of vertices:", sum(1 for face in faces if face.codim == 5))
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
  pairing_dict = polytope_exactness(wbs, faces)
  if pairing_dict is None:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycles(wbs, faces)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          angle_sum_quotient = angle_sum(wbs, faces, ridge_cycle)
          if angle_sum_quotient is None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", angle_sum_quotient)
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
      if is_word_recovered(wbs, faces, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", find_path_word(wbs, faces, matrix.T @ center @ matrix))
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
  wbs, faces = compute_selberg_domain(generators, center, 1, 1, 20)
</pre>

Counting the number of faces:

<pre markdown>
  print("Number of facets:", sum(1 for face in faces if face.codim == 1))
  print("Number of ridges:", sum(1 for face in faces if face.codim == 2))
  print("Number of peaks:", sum(1 for face in faces if face.codim == 3))
  print("Number of edges:", sum(1 for face in faces if face.codim == 4))
  print("Number of vertices:", sum(1 for face in faces if face.codim == 5))
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
  pairing_dict = polytope_exactness(wbs, faces)
  if pairing_dict is None:
      print("The Dirichlet-Selberg domain is not exact.")
  else:
      ridge_cycles = compute_ridge_cycles(wbs, faces)
      for i in range(len(ridge_cycles)):
          ridge_cycle = ridge_cycles[i]
          print("The indices of ridges in the", i, "th cycle:", ridge_cycle.ridge)
          angle_sum_quotient = angle_sum(wbs, faces, ridge_cycle)
          if angle_sum_quotient is None:
              print("The", i, "th ridge cycle does not satisfy the angle sum condition.")
          else:
              print("The angle sum divisor for the", i, "th ridge cycle equals", angle_sum_quotient)
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
      if is_word_recovered(wbs, faces, matrix):
          print("the", ind, "th generator is recovered by the product of facet pairings with indices:", find_path_word(wbs, faces, matrix.T @ center @ matrix))
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
