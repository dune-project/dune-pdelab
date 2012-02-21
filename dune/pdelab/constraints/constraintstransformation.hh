// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH

#include <unordered_map>
#include <dune/common/tuples.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! \brief a class holding transformation for constrained spaces
    template<typename DI, typename F>
    class ConstraintsTransformation
      : public std::unordered_map<DI,std::unordered_map<DI,F> >
    {
    public:
      //! export ElementType
      typedef F ElementType;
      typedef F Field;

      class LocalTransformation
        : public std::unordered_map<typename DI::size_type,std::unordered_map<typename DI::size_type,F> >
      {

      public:

        typedef F ElementType;
        typedef F Field;

        typedef std::unordered_map<typename DI::size_type,F> RowType;

      };

      //! export RowType
      typedef typename ConstraintsTransformation::mapped_type RowType;
    };

    class EmptyTransformation : public ConstraintsTransformation<char,char>
    {
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH
