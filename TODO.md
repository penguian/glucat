TODO for GluCat 0.98a2 with PyClical
====================================

## Roadmap

### 0.98a2 (Current)
* [x] Achieve 100% C++ function coverage across the core library via `test_doctest`.
* [x] Reinforce representation-independence by defining `nbr_rows()` and `nbr_cols()` as protected in `matrix_multi`.

### 0.98a3: Defense in Depth
* **C++ Quality**: Implement `clang-format` and `cppcheck` linting.
* **Cython/Python Quality**: Implement `ruff` and `pylint` for PyClical.
* **Automation**: Integrate `pre-commit` hooks for all quality checks.
* **CI/CD**: Establish GitHub Actions matrix builds (GCC/Clang, multiple Boost versions).
* **Environment Automation**: Automate Jupyter and Mayavi setup as much as possible; fully document manual steps for non-automatable constraints.
* Further harden the stochastic identity tests in `test_doctest`.
* Complete the transition of all legacy tests to the modern unit testing framework.
* Improve the packaging of the example and test programs.

### 0.98a4: Raytracer
* Integrate `rayconfglucat` as `raytracer` subdirectory under GPL v3.
* Implement `--enable-raytracer` configure option and dependency checks (FLTK, OpenGL, libpng, zlib).
* Update coordinate transformations in `scene.cpp` to use `operator|` and `versor_exp`.

### 0.98a5: Documentation
* **ReadTheDocs**: Set up RTD formatting and automated deployment.
* **User Guide**: Develop a comprehensive User Guide using Sphinx + MyST.
* **Maintenance Guide**: Develop a technical Maintenance Guide using Doxygen + Sphinx (Breathe).
* **Case Studies**: Write a programmer's guide with descriptions of usage via use cases.

## Future Stabilisation (0.99bx series)
* Completely remove legacy `icc`/`icpc` compile flags, quirks, and workarounds in the next major version release.

## Post-1.0 Possibilities

### Sage and Boost integration
* Expand the Cython-based Python extension module PyClical into a Sage interface.
* Try defining concepts and more numeric traits so that GluCat can eventually become a Boost library.

### Machine Learning and Clifford Neural Networks
* Develop a unified PyTorch-PyClical Autograd Extension to integrate high-dimensional Clifford algebra operations into deep learning training pipelines.
* Implement custom PyTorch `torch.autograd.Function` wrapper layers that offload geometric computation to pre-compiled GluCat templates in the forward pass, and compute analytical bilinear gradients using Clifford conjugation in the backward pass.
* Establish integrations with geometric deep learning frameworks (such as Microsoft's *cliffordlayers*) to solve high-dimensional optimization problems (e.g. searching for LPS generators on $S^7$ and $S^{15}$) without automatic differentiation graph-tracking bloat.

### Transcendental functions
* Devise better algorithms and better implementations of existing algorithms for
  the transcendental functions. In particular, pay more attention to radius of
  convergence, condition number of matrices and poking out of the subalgebra.
* Investigate the use of matrix decompositions in the evaluation of
  transcendental functions. See N. J. Higham, Functions of Matrices: Theory and
  Computation, 2008.

### Experimental investigations
* Expand the use of numeric type promotion and demotion from transcendental
  functions to operations such as division.
* Investigate the use of expression templates.
* Try refactoring the relationship between `matrix_multi`, `framed_multi` and
  `clifford_algebra` to allow more flexibility with template parameters.
  Possibly use enable_if and SFINAE to do this.
* Try removing the template parameters `LO` and `HI` from `framed_multi` and
  `matrix_multi`, and using `DEFAULT_LO` and `DEFAULT_HI` where these are needed in
  `framed_multi.h`, `matrix_multi.h`, etc.
* Add convenience constructors to `index_set<>`: `index_set<>(int, int)`,
  `index_set<>(int, int, int)`, ... etc.
* Try replacing multiplication by +/-1 within inner products by addition and
  subtraction.
* Try creating a class, `vector_multi<>`, which uses `std::vector` rather than
  `std::map`. This should be faster than `framed_multi<>`, if tuned properly, for
  full multivectors. For sparse multivectors, it may be slower.
