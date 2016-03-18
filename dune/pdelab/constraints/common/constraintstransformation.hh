// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:
#ifndef DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSTRANSFORMATION_HH
#define DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSTRANSFORMATION_HH

#include <algorithm>
#include <unordered_map>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //! \brief a class holding transformation for constrained spaces
    template<typename DI, typename CI, typename F>
    class ConstraintsTransformation
      : public std::unordered_map<CI,std::unordered_map<CI,F> >
    {

      typedef std::unordered_map<CI,std::unordered_map<CI,F> > BaseT;

    public:
      //! export ElementType
      typedef F ElementType;
      typedef F Field;

      class LocalTransformation
        : public std::unordered_map<DI,std::unordered_map<DI,F> >
      {

      public:

        typedef F ElementType;
        typedef F Field;

        typedef std::unordered_map<DI,F> RowType;

        bool containsNonDirichletConstraints() const
        {
          return std::any_of(
            this->begin(), this->end(),
            [] (const std::pair<DI,RowType> & c)
                -> bool { return (!c.second.empty()); });
        }

      };

      //! export RowType
      typedef typename ConstraintsTransformation::mapped_type RowType;

      ConstraintsTransformation()
        : _contains_non_dirichlet_constraints(false)
      {}

      void clear()
      {
        BaseT::clear();
        _contains_non_dirichlet_constraints = false;
      }

      template<typename IndexCache>
      void import_local_transformation(const LocalTransformation& local_transformation, const IndexCache& index_cache)
      {
        typedef typename IndexCache::ContainerIndex ContainerIndex;
        typedef typename ConstraintsTransformation::iterator GlobalConstraintIterator;
        typedef typename ConstraintsTransformation::mapped_type GlobalConstraint;

        for (const auto& local_constraint : local_transformation)
          {
            const ContainerIndex& ci = index_cache.containerIndex(local_constraint.first);

            std::pair<GlobalConstraintIterator,bool> r =
              this->insert(make_pair(ci,GlobalConstraint()));

            GlobalConstraint& global_constraint = r.first->second;

            // Don't modify an existing Dirichlet constraint
            if (!r.second && global_constraint.empty())
              continue;

            // The new constraint is a Dirichlet constraint
            // Clear out any existing entries in the global constraint and stop
            if (local_constraint.second.empty())
              {
                global_constraint.clear();
                continue;
              }

            // We have a non-Dirichlet constraint
            _contains_non_dirichlet_constraints = true;

            // Accumulate new entries into global constraint
            for (const auto& local_entry : local_constraint.second)
              global_constraint[index_cache.containerIndex(local_entry.first)] = local_entry.second;
          }
      }

      bool containsNonDirichletConstraints() const
      {
        return _contains_non_dirichlet_constraints;
      }

    private:

      bool _contains_non_dirichlet_constraints;

    };

    class EmptyTransformation : public ConstraintsTransformation<char,char,char>
    {

    public:

      bool containsNonDirichletConstraints() const
      {
        return false;
      }

    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTSTRANSFORMATION_HH
