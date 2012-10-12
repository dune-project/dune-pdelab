// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH

#include <dune/common/tuples.hh>
#include <dune/pdelab/common/unordered_map.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! \brief a class holding transformation for constrained spaces
    template<typename DI, typename CI, typename F>
    class ConstraintsTransformation
      : public unordered_map<CI,unordered_map<CI,F> >
    {
    public:
      //! export ElementType
      typedef F ElementType;
      typedef F Field;

      class LocalTransformation
        : public unordered_map<DI,unordered_map<DI,F> >
      {

      public:

        typedef F ElementType;
        typedef F Field;

        typedef unordered_map<DI,F> RowType;

      };

      //! export RowType
      typedef typename ConstraintsTransformation::mapped_type RowType;

      template<typename IndexCache>
      void import_local_transformation(const LocalTransformation& local_transformation, const IndexCache& index_cache)
      {
        typedef typename LocalTransformation::const_iterator LocalConstraintIterator;
        typedef typename LocalTransformation::mapped_type::const_iterator LocalEntryIterator;
        typedef std::unordered_map<CI,F> GlobalConstraint;

        for (LocalConstraintIterator lc_it = local_transformation.begin(),
               lc_end = local_transformation.end();
             lc_it != lc_end;
             ++lc_it)
          {
            GlobalConstraint& global_constraint = (*this)[index_cache.container_index(lc_it->first)];
            for (LocalEntryIterator le_it = lc_it->second.begin(),
                   le_end = lc_it->second.end();
                 le_it != le_end;
                 ++le_it)
              {
                global_constraint[index_cache.container_index(le_it->first)] = le_it->second;
              }
          }
      }

    };

    class EmptyTransformation : public ConstraintsTransformation<char,char,char>
    {
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_CONSTRAINTSTRANSFORMATION_HH
