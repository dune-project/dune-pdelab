PDELab
======

This is the 2.4.0-rc1 version of PDELab, a PDE discretization toolkit built
on top of the [DUNE][] framework. License information can be found in the file
[LICENSE.md][].

PDELab 2.4 is a major release with many changes. For details and an overview of the bug fixes in
this release, see the changelog below.

If you need help, please ask on our [mailing list][]. Bugs can also be submitted
to the [PDELab bugtracker][] instead.

Changes
=======

PDELab 2.4
----------

-   PDELab has updated its minimum compiler requirements. You now need a compiler that is at
    least compatible with GCC 4.7 in C++11 mode.

-   The PDELab build system now uses the dune_enable_all_packages() feature and thus requires
    at least CMake 2.8.12 to build.

-   PDELab 2.4.0 requires at least version 2.4.0 of the core modules.

-   In stride with the changes to the core modules, a lot of backwards compatibility code for
    older compilers was removed.

-   There has been a **major** cleanup of the local operators included in PDELab. There were lots of
    duplicate implementations that had similar features, but often used very different interfaces. A lot
    of the older, badly maintained operators have been deprecated or removed, and there is now a much smaller
    set of operators that you should use. Moreover, different operators for the same type of problem (e.g.
    versions for Continuous Galerkin, DG and Finite Volumes) now typically use the same parameter interfaces,
    making it much easier to test different discretizations. In particular:

    -   For convection-diffusion-reaction problems, there are now three files [convectiondiffusionfem.hh][],
        [convectiondiffusiondg.hh][] and [convectiondiffusionccfv.hh][], with a unified parameter interface in
        [convectiondiffusionparameter.hh][]. The implementation in [diffusionmixed.hh][] also uses this parameter
        interface whereas [convectiondiffusion.hh][] which has been renamed to [nonlinearconvectiondiffusionfem.hh][]
        uses its own parameter interface. All other diffusion-type and convection-type operators have been
        deprecated and will be removed after PDELab 2.4.

    -   The parameter class `ConvectionDiffusion_Diffusion_Adapter` has been deprecated and moved from
        [convectiondiffusionparameter.hh][] to the deprecated parameter interface in [diffusionparam.hh][].
        **Note** that the usage of this old parameter interface is now strongly discouraged since it leads
        to a hard **compile error** in order to avoid a deprecation warning whenever [convectiondiffusionparameter.hh][]
        is being included.

    -   New Darcy velocity adapters in [darcy_CCFV.hh][] and [darcy_FEM.hh][] as well as a permeability adapter in
        [permeability_adapter.hh][].

    -   There has been a massive reorganization of the (Navier-)Stokes code. All implementations now use a common
        parameter class, and we now have three implementations: one based on standard Taylor-Hood elements (in
        [taylorhoodnavierstokes.hh][], renamed from [cg_stokes.hh][]), a similar implementation using a DG
        discretization (in [dgnavierstokes.hh][]) and a discontinuous Galerkin discretization that uses a vector-valued
        finite element for the velocity (in [dgnavierstokesvelvecfem.hh][]). All of these implementations now also
        have fully analytic jacobians.

    -   `vectorwave.hh` has been broken for a long time, and as we could not identify any users, it was removed
        without a deprecation period.

    -   The method `g()` in the parameter class for the convection diffusion operators now expects an unwrapped entity
        (that means you have to call it with `eg.entity()` instead of just `eg()`. The version of `g()` that can be
            called with an intersection has been deprecated, please always call the version taking an entity.

    -   All of the operators were updated to new standards in the 2.4 release of the core modules (copyable entities and
        intersections, renamed dimension constants, range-based for loops etc.).

    -   All remaining traces of mimetic finite differences support (which has been broken since at least PDELab 1.1)
        have been removed.

    -   The header [instationary/onestep.hh][] has been split into separate headers for the implicit and explicit one
        step methods and the parameter classes with the Butcher tableaus. The implementation of the class `FilenameHelper`
        has been moved to [common/instationaryfilenamehelper.hh][].

-   The linear algebra backends have also seen a large cleanup:

    -   The code for the different backends has been moved into separate subdirectories like [backend/istl][],
        [backend/eigen][] etc. As part of this change, many files have had their naming improved. As an example,
        [istlvectorbackend.hh][] is now simply [istl/vector.hh][].

    -   Similarly, all of the classes have been moved to corresponding subnamespaces of `Dune::PDELab`. The old classes
        in the namespace `Dune::PDELab` have all been deprecated. Note that when you switch from
        `Dune::PDELab::ISTLVectorBackend` to `Dune::PDELab::istl::VectorBackend`, the type of the `enum` used to describe
        the desired blocking also changes to `Dune::PDELab::istl::Blocking`. The values of the old and the new `enum` can
        be mapped using the following table:

        |`Dune::PDELab::ISTLParameters::Blocking`|`Dune::PDELab::istl::Blocking`|
        |----------------------------------------|------------------------------|
        | `no_blocking`                          | `none`                       |
        | `dynamic_blocking`                     | `bcrs`                       |
        | `static_blocking`                      | `fixed`                      |

        The new identifiers are hopefully better at conveying the actual result of choosing a particular type of
        blocking.

    -   There are now alias templates that provide a much more readable way of extracting vector and matrix types from
        function spaces that feels a lot more like simply using a class template:

        -   Vectors

            ```c++
            // old
            typedef typename Dune::PDELab::BackendVectorSelector<GFS,Field>::Type Vec;
            // new
            using Vec = Dune::PDELab::Backend::Vector<GFS,Field>;
            ```

        -   Matrices

            ```c++
            // old
            typedef typename Dune::PDELab::BackendMatrixSelector<Backend,ColumnVector,RowVector,Field>::Type Mat;
            // new
            using Mat = Dune::PDELab::Backend::Matrix<Backend,ColumnVector,RowVector,Field>;
            ```

    -   The backend infrastructure now has a common method for extracting the native vectors and matrices for all backends
        from the PDELab-specific wrappers. In PDELab 2.0, you could use the `raw()` functions and `raw_type` metafunctions
        for the ISTL backends, but the other backends didn't have any comparable feature. In 2.4, there is now a
        function `native()` and an alias template `Native<>` for this purpose:

        -   Types

            ```c++
            using NativeVector = Dune::PDELab::Backend::Native<Vec>; // a native ISTL BlockVector
            using IdemPotentVector = Dune::PDELab::Backend::Native<NativeVector>; // the functionality is idempotent
            ```

        -   Objects

            ```c++
            auto& native_vector = native(vec);
            ```

        The `native()` function can typically be found using ADL, so you don't have to specify a namespace (like the
        entity iteration functions in `dune-grid`).

    -   The older, ISTL-specific mechanism using `raw()` and `raw_type<>` has been deprecated and will be removed after
        PDELab 2.4.

    -   The old `ISTLMatrixBackend` has been deprecated and will be removed after PDELab 2.4. Please switch to the new
        `istl::BCRSMatrixBackend`, which is much faster during pattern construction. Note, however, that the new
        construction method requires you to provide an approximate guess for the average number of non-zero matrix entries
        per row. A wrong guess will slow down the program, but it will not cause any fatal errors. After matrix
        construction, the matrix provides you with some statistics about the quality of your guess through the member
        function `BCRSMatrix::patternStatistics()`.

        In the case of structured grids it is possible to derive a reasonable estimate for the number of non-zeros that
        is in the most cases even exact. Considering the simplest case of a continuous Galerkin discretization of a
        scalar valued quantity with polynomial degree *deg* in *dim* dimensions the number of non-zeros can be set to
        *(2*deg+1)^dim* which corresponds to the stencil of a discretization with cubic elements.

        The number of non-zeros also depends on the blocking of the unknowns. For a discontinuous Galerkin discretization
        with `Dune::PDELab::istl::Blocking` set to `fixed` this number is independent on the polynomial degree.
        It only depends on the number of faces the mesh elements have which leads to *2*dim+1* on cubic grids and
        *dim+2* on simplicial grids. When considering discretizations involving vector-valued quantities
        the number of non-zeros depend both on the blocking and on the ordering of the unknows. The following table
        summarizes a reasonable choice for the pattern construction in common cases:

        -   Scalar valued quantities

            |Discretization                            |Number of non-zeros|
            |------------------------------------------|-------------------|
            | Continuous Galerkin                      | (2 deg + 1)^dim   |
            | DG on cubic grids, blocking enabled      | 2 dim + 1         |
            | DG on simplicial grids, blocking enabled | dim + 2           |
            | DG, blocking enabled, mass matrix        | 1                 |

        -   Vector valued quantities, using the ordering `Dune::PDELab::EntityBlockedOrderingTag`

            |Discretization                            |Number of non-zeros|
            |------------------------------------------|-------------------|
            | DG on cubic grids, blocking enabled      | 2 dim + 1         |
            | DG on simplicial grids, blocking enabled | dim + 2           |
            | DG, blocking enabled, mass matrix        | 1                 |


            Note that the number of non-zeros is equal to the scalar case.

        -   Vector valued quantities, using the ordering `Dune::PDELab::LexicographicOrderingTag`

            |Discretization                            |Number of non-zeros     |
            |------------------------------------------|------------------------|
            | DG on cubic grids, blocking enabled      | # Childs * (2 dim + 1) |
            | DG on simplicial grids, blocking enabled | # Childs * (dim + 2)   |
            | DG, blocking enabled, mass matrix        | # Childs               |

-   Tests for PDELab are now created using a new CMake function `pdelab_add_test()`, which makes it possible to have
    tests that run on multiple MPI ranks as well as a number of other interesting features. If you are interested, you
    can also use this function in your own modules -- take a look at `cmake/modules/DunePdelabTestMacros.cmake` for
    the documentation. DUNE 3.0 will contain a similar feature in the core modules. Moreover, `make test` will **not**
    build the PDELab tests anymore before running them, you have to explicitly build them using `make build_tests`.
    The manual approach avoids lots of dark CMake magic and makes it possible to build multiple tests in parallel.

-   The support for nonoverlapping parallel computations has been completely rewritten due to changes in the upstream
    modules and is now much more robust. This rewrite does, however, fundamentally change a number of PDELab internals:

    -   `GridFunctionSpace`s and everything built on top of them (`GridOperator`s, constraints etc.) are now defined on an
    `EntitySet` instead of a `GridView`. This `EntitySet` can span a set of parallel partitions that is smaller than
    `Partitions::all`. Its interface is very similar to that of a `GridView`. The most important difference is that it
    provides an `IndexSet` that is restricted to the underlying parallel partition set, i.e. the indices of that set are
    consecutive on that partition set. Users can still create a `GridFunctionSpace` on top of a `GridView`, which will
    automatically be wrapped during construction of the space, but internally, all of PDELab now expects an `EntitySet`.

    -   For nonoverlapping computations, users now **must** manually construct a correct `NonOverlappingEntitySet` and
        pass it to the `GridFunctionSpace`. The `EntitySet` has value semantics, but it is important to only create it
        once and then copy it afterwards, as all copies will share a single `IndexSet`, and this index set can be
        expensive in terms of both setup time and memory usage.

    -   Nonoverlapping computations should now use `ConformingDirichletConstraints` instead of
        `NonoverlappingConformingDirichletConstraints`, as there are no ghost DOFs that need to be constrained anymore.

    -   The template parameter `nonoverlapping_mode` of the `GridOperator` is deprecated, the correct parallelization
        model is extracted from the function spaces and their `EntitySet`.

    -   Take a look at `test/testnonoverlappingsinglephaseflow.cc` for an example of how to port your existing programs.

    -   The `FiniteElementMap` API has been extended with a method `bool hasDOFs(int codim)`. All `FiniteElementMap`
        implementations must support this method, which has to return `true` if it is possible that DOFs might be attached
        to the given codimension.

-   The `StationaryMatrixLinearSolver` has been deprecated. Please use the `StationaryLinearProblemSolver` instead.

-   The `PermutationOrdering` has been deprecated. Please use `PermutedOrdering` instead.

-   There is a new ordering decorator that chunks the index spaces of its children according to a simple list of chunk
    sizes. Take a look at `test/testchunkedblockordering.cc` for an example of how to use this decorator.

-   The deprecated and broken support for multi step methods has been removed.

-   [gridfunctionspace/gridfunctionspaceutilities.hh][] now contains grid functions for the divergence and curl of a
    vector field, even for `VectorGridFunctionSpace`s as in the Taylor-Hood case.

-   The Newton solver implementation now defaults to **not** reallocating the matrix for each iteration, which will
    significantly speed up the solver in many cases. If this setting is problematic for your program, it can be overridden
    using either a method on the `Newton` class and the `ParameterTree` interface.

-   (Hopefully) all of the APIs deprecated in PDELab 2.0 have been removed.

-   The `PermutationOrderingTag` and its implementation have been deprecated. If you need to permute an ordering, apply
    the `Permuted<>` decorator instead.

-   There are probably some additional APIs that have been deprecated for removal after the release of PDELab 2.4.

-   We have added a few additional tests and fixed some of the existing ones by either removing clearly broken tests or
    updating them to work again.

-   A lot of existing code has been updated to take advantage of C++11 features like range-based for loops and `auto`.

-   Lots and lots of bug fixes.

### Release history

###### PDELab 2.4.0-rc1 ######

-   Initial release candidate

PDELab 2.0
----------

-   The TypeTree library is now an external dependency. See README for information on
    how to obtain it.

-   PDELab now supports building with CMake in addition to autotools.

-   The DOF handling and the linear algebra backends have been completely rewritten to allow
    for more flexibility (nested vectors / matrices, more elaborate blocking / DOF reordering).
    Most of these changes are transparent to the average user, but you will have to change the
    typedef of your vector and matrix backends. Changing the matrix backend might also speed up
    matrix pattern construction (see below for further details). If you have been working with
    some of the internals of the GridFunctionSpace (e.g. by writing your own discrete grid functions
    that need to load data from a solution vector), you will have to adapt your code to the new
    index mapping and DOF access structure. See the file [doc/README.Changes-2.0][] for further details.

-   There has been an important change to the way nested `GridFunctionSpace`s can be used: It is no longer
    allowed to use a `GridFunctionSpace` on its own and then integrate it into a larger system of spaces at
    a later time. Trying to do so will result in an exception at run time. This change was necessary as
    part of the move to the new DOF handling.

-   There is a new, vastly improved pattern building algorithm for ISTL matrices, which is a lot faster
    than the old version and also avoid the problem of requiring twice as much memory as the final matrix
    during pattern construction. You can select the new algorithm by using the new class
    `Dune::PDELab::istl::BCRSMatrixBackend`. The new algorithm requires you to supply an approximate guesstimate
    of the average number of matrix entries per row, which you pass to the constructor of the matrix backend.
    That matrix backend is then passed to the constructor of the `GridOperator`. After matrix construction, you can
    query the matrix container for statistics about the pattern creation process, in particular about the quality
    of your guess. If your guess was bad, the performance (both memory and run time wise) will slowly degrade, so
    you should try to pass in a good value, even though it does not have to be optimal. All of the examples in the
    Howto have been ported to the new backend and can be used as a template for porting user programs.

-   There are a number of new finite element maps, e.g. Brezzi-Douglas-Marini, a larger number of
    Raviart-Thomas elements, and new QkDG finite element maps.

-   General overhaul of many `FiniteElementMap`s: Instead of having different classes depending
    on the order and dimension, there is now usually a single class that can be parameterized
    using template parameters. The old headers and classes have been deprecated and will be
    removed after the 2.0 release.

-   As part of the new DOF handling, `FiniteElementMap`s now need to export some information about
    the DOF layout, in particular whether they are fixed size (i.e. always attach the same number
    of DOFs to an entity of given GeometryType) and if they do, the number of those DOFs. This
    extra information allows us to avoid a grid traversal for fixed size FEMs.

-   The `ExplicitOneStepMethod` now supports limiters. There is a rather minimal interface (the limiter
    gets called at the beginning of each stage with the old solution and at the end of each stage with
    the new one). It uses a virtual interface, allowing you to switch limiters at run time.

-   Nonoverlapping grids can now be used without creating DOFs on the ghost partition, improving surface-to-
    volume ratio and simplifying the AMG implementation for nonoverlapping grids. In order to use this
    mode, you have to pass a special `OrderingTag` to the leaf `GridFunctionSpace`s and define them on a `GridView`
    that is restricted to the `InteriorBorder` partition. See [dune/pdelab/test/testnonoverlapping.cc][] for
    an example.

-   The constraints files were cleaned up. Now all constraints files have been moved to the directory
    [dune/pdelab/constraints][] (instead of the previous split between [dune/pdelab/constraints][] and
    [dune/pdelab/finiteelementmap][]). The infrastructure headers are now in [dune/pdelab/constraints/common][].
    The old headers are still there in this release, but they have been deprecated and will be removed after
    2.0.

-   There is new support infrastructure for automatically disassembling the `GridFunctionSpace` of a system
    and outputting all of the components to a VTK file. See [doc/README.Changes-2.0][] for further information.

-   The new `VectorGridFunctionSpace` can be used to represent vector-valued fields. It acts like a combination
    of a scalar `GridFunctionSpace` and a `PowerGridFunctionSpace` on top of it. As an added bonus, it will
    automatically be detected by the new VTK output code and will be written as a vector-valued solution to
    VTK.

-   The adaptivity support has been fixed and now also works for systems of variables. As long as there are
    no hanging nodes, it should also work for arbitrary discretizations.

-   The new `PermutedOrdering` is a decorator that can be used to wrap an existing ordering tag and perform
    arbitrary permutations of the DOFs in that underlying ordering.

-   The ISTL vector backend now provides an iterator for flat iteration over the scalar entries in the vector,
    automatically traversing a possible block structure. This can be very useful when you want to simply dump
    or load the contents of a vector for debug or restart purposes.

-   The `GridOperatorSpace`, which had been broken for a long time, has now been completely removed.

-   The `GridFunctionSpace` data handles have been completely rewritten and now support communication per
    entity or per DOF. To avoid code duplication, they have been split into data handles and gather/scatter
    objects responsible for the actual data handling.

-   There have been numerous improvements and bugfixes in the local operators.

-   The method `Dune::PDELab::istl::raw(x)` provides idempotent access to the naked ISTL objects from the
    PDELab wrappers. They will also work correctly if passed an unpacked ISTL objects by returning the
    object itself.

-   The `StdVectorBackend` is gone, but you can use the new simple backend instead. This new backend even
    provides a basic CRS matrix.

-   The constructors of the `StationaryLinearProblemSolver` have been modified to use same order of parameters
    as the one of the `Newton` solver. Moreover, there are now also constructors that load the parameters from
    a `Dune::ParameterTree`. The old constructors have been deprecated and will be removed after PDELab 2.0.

-   The Eigen backend is mostly functional now and works correctly with the new ordering infrastructure.

-   A number of optimizations to the `GridOperator` and its engines make for important performance improvements
    when working with higher-order methods that have a large number of DOFs per element.

-   It is now possible to use a diagonal local matrix for the jacobian assembly to reduce the required memory
    from N^2 to N if N is the number of DOFs per element.

Releasy history

  -   PDELab 2.0.0
      -   Fix instructions for obtaining dune-typetree.
      -   Link `libdunepdelab` against `libdunecommon`.
      -   Documentation fixes.

  -   PDELab 2.0.0-rc2
      -   Buildsystem fixes to include missing headers etc.
      -   Improved handling of constraints engine storage in `GridFunctionSpace`.
      -   Fix PkFEM for k >= 3 in 2D, which was broken due to a bug in the variant
          selection. Also extended interpolation tests to cover this case.
      -   Improvements to coding style and correctness in some places.
      -   Check for correct cell `GeometryType` in PkFEM and VariableMonomFEM.
      -   Fix broken testordering.
      -   Add infrastructure support required for parallel computations with the
          dune-multidomain extension module.

  -   PDELab 2.0.0-rc1
      -   initial release candidate.


PDELab 1.1
----------

This is the first release of PDELab after the move from Subversion to Git for
version control. You CANNOT find this release on the Subversion server. If you prefer
to directly check out the sources instead of using the provided tarballs, please see
our website for information on how to access the Git repository.

-   In the directory boilerplate/ you can find a lot of useful infrastructure classes
    to simplify writing most PDELab programs. For examples on how to use the boilerplate
    infrastructure, take a look at the Howto.

-   There is now a Jacobi preconditioner for scalar nonoverlapping problems, along
    with a BiCGStab-based backend using it.

-   Improved support for calculations with complex field type.

-   The parameter class interface for Stokes problems has been redesigned to increase
    commonality between CG and DG versions.

-   Working adaptivity, including hanging nodes. This functionality is currently restricted
    to scalar problems.

-   Reimplemented support for matrix-free methods using GridOperator::jacobian_apply().

-   We fixed most of the deprecation warnings due to deprecated APIs in the core modules.

-   Fix for builds from repository with recent versions of autotools.

-   Numerous bug fixes.

    -   PDELab 1.1.0
        -   Improved documentation for nonoverlapping AMG solvers

    -   PDELab 1.1-rc2
        -   fix for compilation problem with boilerplate examples in Howto
        -   build tarballs using GNU tar and make sure they work without automake installed

    -   PDELab 1.1-rc1
        -   initial release candidate


Caveats
=======

The following list is a non-exhaustive overview of possible problems you might
encounter with this release.

Assembler
--------

Solvers
-------

*   Both the Newton solver and the linear solve currently allocate a new matrix on
    each call to apply(), which can incur a significant overhead if they are
    called frequently (e.g. for instationary problems). This will be fixed in a
    later release.

Linear Algebra Backends
-----------------------

*   Due to changes in the ISTL backend, users who construct their own solvers
    directly from ISTL primitives will have to make sure to use native ISTL types
    and variables for this. These can be accessed by the nested typedef ::BaseT
    and the method .base() for both vectors and matrices. For an example, see
    src/convection-diffusion/poisson.cc in dune-pdelab-howto. In general, we
    encourage usage of the predefined solver backends, though.

*   The [PETSc][] backend is currently broken.

*   The alternative backend for [Eigen][] is not as mature as the ISTL backend yet.

General
-------

*   Compile times can be really long for non-trivial problems. Some developers
    have had good success with using the clang compiler instead of GCC during
    development and bug-testing to reduce compile times.

*   After PDELab 2.0, the minimum compiler requirement of PDELab will be increased
    to GCC 4.5. Please be aware of this change in minimum requirements.

Links
=====

[DUNE]: https://www.dune-project.org
[mailing list]: http://lists.dune-project.org/mailman/listinfo/dune-pdelab
[PDELab bugtracker]: https://gitlab.dune-project.org/pdelab/dune-pdelab/issues
[PETSc]: http://www.mcs.anl.gov/petsc/
[Eigen]: http://eigen.tuxfamily.org
[LICENSE.md]: LICENSE.md
[convectiondiffusionfem.hh]: dune/pdelab/localoperator/convectiondiffusionfem.hh
[convectiondiffusiondg.hh]: dune/pdelab/localoperator/convectiondiffusiondg.hh
[convectiondiffusionccfv.hh]: dune/pdelab/localoperator/convectiondiffusionccfv.hh
[convectiondiffusionparameter.hh]: dune/pdelab/localoperator/convectiondiffusionparameter.hh
[diffusionparam.hh]: dune/pdelab/localoperator/diffusionparam.hh
[convectiondiffusion.hh]: dune/pdelab/localoperator/convectiondiffusion.hh
[nonlinearconvectiondiffusionfem.hh]: dune/pdelab/localoperator/nonlinearconvectiondiffusionfem.hh
[diffusionmixed.hh]: dune/pdelab/localoperator/diffusionmixed.hh
[darcy_CCFV.hh]: dune/pdelab/localoperator/darcy_CCFV.hh
[darcy_FEM.hh]: dune/pdelab/localoperator/darcy_FEM.hh
[permeability_adapter.hh]: dune/pdelab/localoperator/permeability_adapter.hh
[taylorhoodnavierstokes.hh]: dune/pdelab/localoperator/taylorhoodnavierstokes.hh
[cg_stokes.hh]: dune/pdelab/localoperator/cg_stokes.hh
[dgnavierstokes.hh]: dune/pdelab/localoperator/dgnavierstokes.hh
[dgnavierstokesvelvecfem.hh]: dune/pdelab/localoperator/dgnavierstokesvelvecfem.hh
[instationary/onestep.hh]: dune/pdelab/instationary/onestep.hh
[common/instationaryfilenamehelper.hh]: dune/pdelab/common/instationaryfilenamehelper.hh
[backend/istl]: dune/pdelab/backend/istl
[backend/eigen]: dune/pdelab/backend/eigen
[istlvectorbackend.hh]: dune/pdelab/backend/istlvectorbackend.hh
[istl/vector.hh]: dune/pdelab/backend/istl/vector.hh
[gridfunctionspace/gridfunctionspaceutilities.hh]: dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh
[doc/README.Changes-2.0]: doc/README.Changes-2.0
[dune/pdelab/test/testnonoverlapping.cc]: dune/pdelab/test/testnonoverlapping.cc
[dune/pdelab/constraints]: dune/pdelab/constraints
[dune/pdelab/finiteelementmap]: dune/pdelab/finiteelementmap
[dune/pdelab/constraints/common]: dune/pdelab/constraints/common
