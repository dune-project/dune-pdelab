// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_TREEPATH_HH
#define DUNE_PDELAB_COMMON_TREEPATH_HH

#include <cstddef>

#include <dune/common/documentation.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/typetraits.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup TreePath
    //! \{

    //! number used a dummy child number, similar to Nil
    /**
     * \note This should be used directly, it is an implementation detail of
     *       class TreePath.
     */
    static const std::size_t noChildIndex = ~std::size_t(0);

    //! class used as a key to denote nodes in a tree
    /**
     * Class TreePath is a collection of indices.  It is like a tuple for
     * integral types instead of types.  Like a tuple, the number of entries
     * is not prescribed.
     *
     * Given a tree, class TreePath statically denotes one particular node
     * (subtree) in that tree.  That makes it possible for node visitors to
     * know exactly where in the tree they are, as opposed to only knowing the
     * type of the node they operate on.
     *
     * \note Due to limitations in the C++ language, the length of a tree path
     *       is arbitrarily limited to an implementation defined value
     *       (currently 10).  If this value is too restrictive, patches that
     *       extend it are welcome.
     */
    template<std::size_t i0 = noChildIndex, std::size_t i1 = noChildIndex,
             std::size_t i2 = noChildIndex, std::size_t i3 = noChildIndex,
             std::size_t i4 = noChildIndex, std::size_t i5 = noChildIndex,
             std::size_t i6 = noChildIndex, std::size_t i7 = noChildIndex,
             std::size_t i8 = noChildIndex, std::size_t i9 = noChildIndex>
    class TreePath {
      dune_static_assert(i0 == noChildIndex ? i1 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i1 == noChildIndex ? i2 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i2 == noChildIndex ? i3 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i3 == noChildIndex ? i4 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i4 == noChildIndex ? i5 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i5 == noChildIndex ? i6 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i6 == noChildIndex ? i7 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i7 == noChildIndex ? i8 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
      dune_static_assert(i8 == noChildIndex ? i9 == noChildIndex : true,
                         "Only trailing indices my be noChildIndex");
    };

    //
    // The following classes operate on the front of a TreePath: They extract,
    // remove or prepend the first element.
    //

    //! remove first element of a tree path
    template<class TP>
    class TreePathPopFront {
      dune_static_assert(AlwaysTrue<TP>::value,
                         "TreePathPopFront works on TreePaths only");
    public:
      //! The tree path with the first element removed
      typedef ImplementationDefined type;
    };

#ifndef DOXYGEN
    template<std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3,
             std::size_t i4, std::size_t i5, std::size_t i6, std::size_t i7,
             std::size_t i8, std::size_t i9>
    struct TreePathPopFront<TreePath<i0, i1, i2, i3, i4, i5, i6, i7, i8, i9> >
    {
      dune_static_assert(i0 != noChildIndex,
                         "Can't pop first element from an empty TreePath");
      typedef TreePath<i1, i2, i3, i4, i5, i6, i7, i8, i9> type;
    };
#endif // DOXYGEN

    //! get first element of a tree path
    template<class TP>
    class TreePathFront {
      dune_static_assert(AlwaysTrue<TP>::value,
                         "TreePathPopFront works on TreePaths only");
    public:
      //! value of the first element
      static const std::size_t value = implementationDefined;
    };

#ifndef DOXYGEN
    template<std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3,
             std::size_t i4, std::size_t i5, std::size_t i6, std::size_t i7,
             std::size_t i8, std::size_t i9>
    class TreePathFront<TreePath<i0, i1, i2, i3, i4, i5, i6, i7, i8, i9> > :
      public integral_constant<std::size_t, i0>
    {
      dune_static_assert(i0 != noChildIndex,
                         "Can't take first element of an empty TreePath");
    };
#endif // DOXYGEN

    //! prepend one element to a tree path
    template<class TP, std::size_t i>
    class TreePathPushFront {
      dune_static_assert(AlwaysTrue<TP>::value,
                         "TreePathPushFront works on TreePaths only");
    public:
      //! The tree path with i prepended
      typedef ImplementationDefined type;
    };

#ifndef DOXYGEN
    template<std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3,
             std::size_t i4, std::size_t i5, std::size_t i6, std::size_t i7,
             std::size_t i8, std::size_t i9, std::size_t i>
    struct TreePathPushFront<TreePath<i0, i1, i2, i3, i4, i5, i6, i7, i8, i9>,
                             i>
    {
      dune_static_assert(i0 ? false : false, "TreePathPushFront: exceeded "
                         "implementation limit on TreePath size");
    };
    template<std::size_t i0, std::size_t i1, std::size_t i2, std::size_t i3,
             std::size_t i4, std::size_t i5, std::size_t i6, std::size_t i7,
             std::size_t i8, std::size_t i>
    struct TreePathPushFront<TreePath<i0, i1, i2, i3, i4, i5, i6, i7, i8>, i> {
      typedef TreePath<i, i0, i1, i2, i3, i4, i5, i6, i7, i8> type;
    };
#endif // DOXYGEN

    //
    // Using the front operations from above, is is easy to reverse an
    // arbitrary TreePath
    //

#ifdef DOXYGEN
    //! reverse a tree path
    template<class TP>
    class TreePathReverse {
    public:
      //! type of the reversed tree path
      typedef ImplementationDefined type;
    };
#else // !DOXYGEN
    template<class TP, class TP2 = TreePath<> >
    struct TreePathReverse :
      public TreePathReverse<
        typename TreePathPopFront<TP>::type,
        typename TreePathPushFront<TP2, TreePathFront<TP>::value>::type>
    { };
    template<class TP2>
    struct TreePathReverse<TreePath<>, TP2>
    { typedef TP2 type; };
#endif // DOXYGEN

    //
    // The following classes operate on the back of a TreePath: They extract,
    // remove or append the last element.  They are easy to implement using
    // the reverse operation.
    //

#ifdef DOXYGEN
    //! remove last element of a tree path
    template<class TP>
    struct TreePathPopBack {
      //! type of the tree path with the last element removed
      typedef ImplementationDefined type;
    };
#else // !DOXYGEN
    template<class TP>
    struct TreePathPopBack :
      public TreePathReverse<
        typename TreePathPopFront<typename TreePathReverse<TP>::type>::type
      >
    { };
#endif // DOXYGEN

#ifdef DOXYGEN
    //! get last element of a tree path
    template<class TP>
    struct TreePathBack {
      //! value of the last element
      static const std::size_t value = implementationDefined;
    };
#else // !DOXYGEN
    template<class TP>
    struct TreePathBack :
      public TreePathFront<typename TreePathReverse<TP>::type>
    { };
#endif // DOXYGEN

#ifdef DOXYGEN
    //! append one element to a tree path
    template<class TP, std::size_t i>
    struct TreePathPushBack {
      //! type of the tree path with i appended
      typedef ImplementationDefined type;
    };
#else // !DOXYGEN
    template<class TP, std::size_t i>
    struct TreePathPushBack :
      public TreePathReverse<
        typename TreePathPushFront<typename TreePathReverse<TP>::type, i>::type
      >
    { };
#endif // DOXYGEN

    //
    // Finally some tuple-like classes to determine the length of the TreePath
    // and to extract particular elements.
    //

#ifdef DOXYGEN
    //! determine size of a tree path
    template<class TP>
    struct TreePathSize {
      //! size of the tree path
      static const std::size_t value = implementationDefined;
    };
#else // !DOXYGEN
    template<class TP>
    struct TreePathSize :
      public integral_constant<
        std::size_t,
        TreePathSize<typename TreePathPopFront<TP>::type>::value+1
      >
    { };
    template<>
    struct TreePathSize<TreePath<> > :
      public integral_constant<std::size_t, 0>
    { };
#endif // DOXYGEN

#ifdef DOXYGEN
    //! extract a certain element of a path
    template<class TP>
    struct TreePathElement {
      //! value of the element
      static const std::size_t value = implementationDefined;
    };
#else // !DOXYGEN
    template<std::size_t pos, class TP>
    struct TreePathElement :
      public TreePathElement<pos-1, typename TreePathPopFront<TP>::type>
    { };
    template<class TP>
    struct TreePathElement<0, TP> :
      public TreePathFront<TP>
    { };
#endif // DOXYGEN

    //! \} group TreePath

  } // namespace PDELab
} //namespace Dune

#endif // DUNE_PDELAB_COMMON_TREEPATH_HH

