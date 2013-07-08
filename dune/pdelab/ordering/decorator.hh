// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_DECORATOR_HH
#define DUNE_PDELAB_ORDERING_DECORATOR_HH

#include <dune/pdelab/common/typetree.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/orderingbase.hh>

namespace Dune {
  namespace PDELab {

    namespace ordering {

#ifndef DOXYGEN

      // Forward declaration for is_decorated specialization
      template<typename D, typename U>
      struct decorated_ordering_tag;

      namespace impl {

        struct is_decorated
        {};

        template<typename DecoratedOrderingTag, std::size_t level>
        struct basetag_access_provider
        {};

        // provide access to the base tag for disambiguation of member
        // accesses in case those are overloaded by a decorator.
        template<typename D, typename U>
        struct basetag_access_provider<decorated_ordering_tag<D,U>,0>
          : public is_decorated
        {
          typedef U BaseTag;

          const BaseTag& baseTag() const
          {
            return static_cast<decorated_ordering_tag<D,U>& >(*this);
          }

          BaseTag& baseTag()
          {
            return static_cast<decorated_ordering_tag<D,U>& >(*this);
          }

        };

        template<typename T, T v>
        struct lazy_constant
        {
          typedef integral_constant<T,v> type;
        };

        template<typename D>
        struct lazy_level
        {
          typedef integral_constant<std::size_t,D::level> type;
        };

        template<typename D>
        struct decoration_level
          : public conditional<
                     is_base_of<is_decorated, D>::value,
                     lazy_level<D>,
                     lazy_constant<std::size_t,0>
                   >::type::type
        {};

      } // namespace impl

#endif // DOXYGEN

      template<typename D, typename U>
      struct decorated_ordering_tag
        : U
        , impl::basetag_access_provider<decorated_ordering_tag<D,U>,impl::decoration_level<U>::value>
      {

        typedef D Decorator;
        typedef U Undecorated;

        static const std::size_t level = impl::decoration_level<U>::value + 1;

        decorated_ordering_tag()
        {}

        decorated_ordering_tag(const Undecorated& u)
          : Undecorated(u)
        {}

#if HAVE_RVALUE_REFERENCES

        decorated_ordering_tag(Undecorated&& u)
          : Undecorated(std::move(u))
        {}

#endif // HAVE_RVALUE_REFERENCES

      };


      template<typename GFS,typename Transformation,typename Undecorated,typename GlueTag, typename Tag>
      struct gfs_to_decorator_descriptor
        : public meta_function
      {
        typedef DUNE_DECLTYPE(
          register_gfs_to_decorator_descriptor(
            declptr<GFS>(),             // the source GridFunctionSpace
            declptr<Transformation>(),  // the full transformation descriptor
            declptr<Undecorated>(),     // the type of the undecorated Ordering to be wrapped in the decorator
            declptr<GlueTag>(),         // the decorated_ordering_tag for the current decoration nesting level
            declptr<Tag>()              // the decorator tag
            )
          ) type;
      };


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct leaf_gfs_to_decorated
      {

        static const bool recursive = false;

        typedef typename leaf_gfs_to_ordering_descriptor<
          GFS,
          Transformation,
          typename OrderingTag::Undecorated
          >::type undecorated_descriptor;

        typedef typename undecorated_descriptor::transformed_type undecorated_type;

        typedef typename gfs_to_decorator_descriptor<
          GFS,
          Transformation,
          undecorated_type,
          OrderingTag,
          typename OrderingTag::Decorator
          >::type decorator_descriptor;

        typedef typename decorator_descriptor::transformed_type transformed_type;
        typedef typename decorator_descriptor::transformed_storage_type transformed_storage_type;

        static transformed_type transform(const GFS& gfs, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs,t,make_shared<undecorated_type>(undecorated_descriptor::transform(gfs,t)));
        }

        static transformed_storage_type transform(shared_ptr<const GFS>& gfs_pointer, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs_pointer,t,undecorated_descriptor::transform(gfs_pointer,t));
        }

      };

      template<typename GFS, typename Transformation, typename D, typename U>
      leaf_gfs_to_decorated<GFS,Transformation,decorated_ordering_tag<D,U> >
      register_leaf_gfs_to_ordering_descriptor(GFS*,Transformation*,decorated_ordering_tag<D,U>*);


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct recursive_power_gfs_to_decorated
      {

        static const bool recursive = true;

        template<typename TC>
        struct result
        {

          typedef typename power_gfs_to_ordering_descriptor<
            GFS,
            Transformation,
            typename OrderingTag::Undecorated
            >::type undecorated_descriptor;

          typedef typename undecorated_descriptor::template result<TC>::type undecorated_type;
          typedef typename gfs_to_decorator_descriptor<
            GFS,
            Transformation,
            undecorated_type,
            OrderingTag,
            typename OrderingTag::Decorator
            >::type decorator_descriptor;

          typedef typename decorator_descriptor::transformed_type type;
          typedef typename decorator_descriptor::transformed_storage_type storage_type;

        };

        template<typename TC>
        static typename result<TC>::type transform(const GFS& gfs, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
        {
          return result<TC>::decorator_descriptor::transform(gfs,t,make_shared<typename result<TC>::undecorated_type>(result<TC>::undecorated_descriptor::transform(gfs,t,children)));
        }

        template<typename TC>
        static typename result<TC>::storage_type transform_storage(shared_ptr<const GFS> gfs_pointer, const Transformation& t, const array<shared_ptr<TC>,GFS::CHILDREN>& children)
        {
          return result<TC>::decorator_descriptor::transform(gfs_pointer,t,result<TC>::undecorated_descriptor::transform(gfs_pointer,t,children));
        }

      };


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct nonrecursive_power_gfs_to_decorated
      {

        static const bool recursive = false;

        typedef typename power_gfs_to_ordering_descriptor<
          GFS,
          Transformation,
          typename OrderingTag::Undecorated
          >::type undecorated_descriptor;

        typedef typename undecorated_descriptor::transformed_type undecorated_type;

        typedef typename gfs_to_decorator_descriptor<
          GFS,
          Transformation,
          undecorated_type,
          OrderingTag,
          typename OrderingTag::Decorator
          >::type decorator_descriptor;

        typedef typename decorator_descriptor::transformed_type transformed_type;
        typedef typename decorator_descriptor::transformed_storage_type transformed_storage_type;

        static transformed_type transform(const GFS& gfs, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs,t,make_shared<undecorated_type>(undecorated_descriptor::transform(gfs,t)));
        }

        static transformed_storage_type transform(shared_ptr<const GFS>& gfs_pointer, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs_pointer,t,undecorated_descriptor::transform(gfs_pointer,t));
        }

      };


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct power_gfs_to_decorated
        : public std::conditional<
            power_gfs_to_ordering_descriptor<
              GFS,
              Transformation,
              typename OrderingTag::Undecorated
              >::type::recursive,
            recursive_power_gfs_to_decorated<
              GFS,
              Transformation,
              OrderingTag
              >,
            nonrecursive_power_gfs_to_decorated<
              GFS,
              Transformation,
              OrderingTag>
          >::type
      {};

      template<typename GFS, typename Transformation, typename D, typename U>
      power_gfs_to_decorated<
        GFS,
        Transformation,
        decorated_ordering_tag<D,U>
        >
      register_power_gfs_to_ordering_descriptor(GFS*,Transformation*,decorated_ordering_tag<D,U>*);






      template<typename GFS, typename Transformation, typename OrderingTag>
      struct recursive_composite_gfs_to_decorated
      {

        static const bool recursive = true;

        template<typename... TC>
        struct result
        {

          typedef typename composite_gfs_to_ordering_descriptor<
            GFS,
            Transformation,
            typename OrderingTag::Undecorated
            >::type undecorated_descriptor;

          typedef typename undecorated_descriptor::template result<TC...>::type undecorated_type;
          typedef typename gfs_to_decorator_descriptor<
            GFS,
            Transformation,
            undecorated_type,
            OrderingTag,
            typename OrderingTag::Decorator
            >::type decorator_descriptor;

          typedef typename decorator_descriptor::transformed_type type;
          typedef typename decorator_descriptor::transformed_storage_type storage_type;

        };

        template<typename... TC>
        static typename result<TC...>::type transform(const GFS& gfs, const Transformation& t, shared_ptr<TC>... children)
        {
          return result<TC...>::decorator_descriptor::transform(gfs,t,make_shared<typename result<TC...>::undecorated_type>(result<TC...>::undecorated_descriptor::transform(gfs,t,children...)));
        }

        template<typename... TC>
        static typename result<TC...>::storage_type transform_storage(shared_ptr<const GFS> gfs_pointer, const Transformation& t, shared_ptr<TC>... children)
        {
          return result<TC...>::decorator_descriptor::transform(gfs_pointer,t,result<TC...>::undecorated_descriptor::transform(gfs_pointer,t,children...));
        }

      };


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct nonrecursive_composite_gfs_to_decorated
      {

        static const bool recursive = false;

        typedef typename composite_gfs_to_ordering_descriptor<
          GFS,
          Transformation,
          typename OrderingTag::Undecorated
          >::type undecorated_descriptor;

        typedef typename undecorated_descriptor::transformed_type undecorated_type;

        typedef typename gfs_to_decorator_descriptor<
          GFS,
          Transformation,
          undecorated_type,
          OrderingTag,
          typename OrderingTag::Decorator
          >::type decorator_descriptor;

        typedef typename decorator_descriptor::transformed_type transformed_type;
        typedef typename decorator_descriptor::transformed_storage_type transformed_storage_type;

        static transformed_type transform(const GFS& gfs, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs,t,make_shared<undecorated_type>(undecorated_descriptor::transform(gfs,t)));
        }

        static transformed_storage_type transform(shared_ptr<const GFS>& gfs_pointer, const Transformation& t)
        {
          return decorator_descriptor::transform(gfs_pointer,t,undecorated_descriptor::transform(gfs_pointer,t));
        }

      };


      template<typename GFS, typename Transformation, typename OrderingTag>
      struct composite_gfs_to_decorated
        : public std::conditional<
            composite_gfs_to_ordering_descriptor<
              GFS,
              Transformation,
              typename OrderingTag::Undecorated
              >::type::recursive,
            recursive_composite_gfs_to_decorated<
              GFS,
              Transformation,
              OrderingTag
              >,
            nonrecursive_composite_gfs_to_decorated<
              GFS,
              Transformation,
              OrderingTag>
          >::type
      {};


      template<typename GFS, typename Transformation, typename D, typename U>
      composite_gfs_to_decorated<GFS,Transformation,decorated_ordering_tag<D,U> >
      register_composite_gfs_to_ordering_descriptor(GFS*,Transformation*,decorated_ordering_tag<D,U>*);

    } // namespace ordering

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_DECORATOR_HH
