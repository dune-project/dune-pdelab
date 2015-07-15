// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTS_HH
#define DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTS_HH

#include<dune/common/exceptions.hh>
#include<dune/common/float_cmp.hh>

#include<dune/pdelab/common/function.hh>
#include<dune/pdelab/common/geometrywrapper.hh>
#include<dune/pdelab/common/typetraits.hh>
#include<dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include"constraintstransformation.hh"
#include"constraintsparameters.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    namespace { // hide internals
      // do method invocation only if class has the method

      template<typename C, bool doIt>
      struct ConstraintsCallBoundary
      {
        template<typename P, typename IG, typename LFS, typename T>
        static void boundary (const C& c, const P& p, const IG& ig, const LFS& lfs, T& trafo)
        {
        }
      };
      template<typename C, bool doIt>
      struct ConstraintsCallProcessor
      {
        template<typename IG, typename LFS, typename T>
        static void processor (const C& c, const IG& ig, const LFS& lfs, T& trafo)
        {
        }
      };
      template<typename C, bool doIt>
      struct ConstraintsCallSkeleton
      {
        template<typename IG, typename LFS, typename T>
        static void skeleton (const C& c,  const IG& ig,
                              const LFS& lfs_e, const LFS& lfs_f,
                              T& trafo_e, T& trafo_f)
        {
        }
      };
      template<typename C, bool doIt>
      struct ConstraintsCallVolume
      {
        template<typename P, typename EG, typename LFS, typename T>
        static void volume (const C& c, const P&, const EG& eg, const LFS& lfs, T& trafo)
        {
        }
      };


      template<typename C>
      struct ConstraintsCallBoundary<C,true>
      {
        template<typename P, typename IG, typename LFS, typename T>
        static void boundary (const C& c, const P& p, const IG& ig, const LFS& lfs, T& trafo)
        {
          if (lfs.size())
            c.boundary(p,ig,lfs,trafo);
        }
      };
      template<typename C>
      struct ConstraintsCallProcessor<C,true>
      {
        template<typename IG, typename LFS, typename T>
        static void processor (const C& c, const IG& ig, const LFS& lfs, T& trafo)
        {
          if (lfs.size())
            c.processor(ig,lfs,trafo);
        }
      };
      template<typename C>
      struct ConstraintsCallSkeleton<C,true>
      {
        template<typename IG, typename LFS, typename T>
        static void skeleton (const C& c, const IG& ig,
                              const LFS& lfs_e, const LFS& lfs_f,
                              T& trafo_e, T& trafo_f)
        {
          if (lfs_e.size() || lfs_f.size())
            c.skeleton(ig, lfs_e, lfs_f, trafo_e, trafo_f);
        }
      };
      template<typename C>
      struct ConstraintsCallVolume<C,true>
      {
        template<typename P, typename EG, typename LFS, typename T>
        static void volume (const C& c, const P& p, const EG& eg, const LFS& lfs, T& trafo)
        {
          if (lfs.size())
            c.volume(p,eg,lfs,trafo);
        }
      };


      struct ParameterizedConstraintsBase
        : public TypeTree::TreePairVisitor
      {
        // This acts as a catch-all for unsupported leaf- / non-leaf combinations in the two
        // trees. It is necessary because otherwise, the visitor would fall back to the default
        // implementation in TreeVisitor, which simply does nothing. The resulting bugs would
        // probably be hell to find...
        template<typename P, typename LFS, typename TreePath>
        void leaf(const P& p, const LFS& lfs, TreePath treePath) const
        {
          static_assert((AlwaysFalse<P>::Value),
                        "unsupported combination of function and LocalFunctionSpace");
        }
      };


      template<typename P, typename IG, typename CL>
      struct BoundaryConstraintsForParametersLeaf
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {

          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;

          // iterate over boundary, need intersection iterator
          ConstraintsCallBoundary<C,C::doBoundary>::boundary(lfs.constraints(),p,ig,lfs,cl);
        }

        BoundaryConstraintsForParametersLeaf(const P& p_, const IG& ig_, CL& cl_)
          : p(p_)
          , ig(ig_)
          , cl(cl_)
        {}

        const P& p;
        const IG& ig;
        CL& cl;

      };


      template<typename IG, typename CL>
      struct BoundaryConstraints
        : public ParameterizedConstraintsBase
        , public TypeTree::DynamicTraversal
      {

        // standard case - leaf in both trees
        template<typename P, typename LFS, typename TreePath>
        typename enable_if<P::isLeaf && LFS::isLeaf>::type
        leaf(const P& p, const LFS& lfs, TreePath treePath) const
        {
          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;

          // iterate over boundary, need intersection iterator
          ConstraintsCallBoundary<C,C::doBoundary>::boundary(lfs.constraints(),p,ig,lfs,cl);
        }

        // reuse constraints parameter information from p for all LFS children
        template<typename P, typename LFS, typename TreePath>
        typename enable_if<P::isLeaf && (!LFS::isLeaf)>::type
        leaf(const P& p, const LFS& lfs, TreePath treePath) const
        {
          // traverse LFS tree and reuse parameter information
          TypeTree::applyToTree(lfs,BoundaryConstraintsForParametersLeaf<P,IG,CL>(p,ig,cl));
        }

        BoundaryConstraints(const IG& ig_, CL& cl_)
          : ig(ig_)
          , cl(cl_)
        {}

      private:
        const IG& ig;
        CL& cl;

      };


      template<typename IG, typename CL>
      struct ProcessorConstraints
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;

          // iterate over boundary, need intersection iterator
          ConstraintsCallProcessor<C,C::doProcessor>::processor(lfs.constraints(),ig,lfs,cl);
        }

        ProcessorConstraints(const IG& ig_, CL& cl_)
          : ig(ig_)
          , cl(cl_)
        {}

      private:
        const IG& ig;
        CL& cl;

      };


      template<typename IG, typename CL>
      struct SkeletonConstraints
        : public TypeTree::TreePairVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs_e, const LFS& lfs_f, TreePath treePath) const
        {
          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;

          // as LFS::constraints() just returns the constraints of the
          // GridFunctionSpace, lfs_e.constraints() is equivalent to
          // lfs_f.constraints()
          const C & c = lfs_e.constraints();

          // iterate over boundary, need intersection iterator
          ConstraintsCallSkeleton<C,C::doSkeleton>::skeleton(c,ig,lfs_e,lfs_f,cl_e,cl_f);
        }

        SkeletonConstraints(const IG& ig_, CL& cl_e_, CL& cl_f_)
          : ig(ig_)
          , cl_e(cl_e_)
          , cl_f(cl_f_)
        {}

      private:
        const IG& ig;
        CL& cl_e;
        CL& cl_f;

      };


      template<typename P, typename EG, typename CL>
      struct VolumeConstraintsForParametersLeaf
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        template<typename LFS, typename TreePath>
        void leaf(const LFS& lfs, TreePath treePath) const
        {
          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;
          const C & c = lfs.constraints();

          // iterate over boundary, need intersection iterator
          ConstraintsCallVolume<C,C::doVolume>::volume(c,p,eg,lfs,cl);
        }

        VolumeConstraintsForParametersLeaf(const P& p_, const EG& eg_, CL& cl_)
          : p(p_)
          , eg(eg_)
          , cl(cl_)
        {}

      private:
        const P& p;
        const EG& eg;
        CL& cl;

      };


      template<typename EG, typename CL>
      struct VolumeConstraints
        : public ParameterizedConstraintsBase
        , public TypeTree::DynamicTraversal
      {

        // standard case - leaf in both trees
        template<typename P, typename LFS, typename TreePath>
        typename enable_if<P::isLeaf && LFS::isLeaf>::type
        leaf(const P& p, const LFS& lfs, TreePath treePath) const
        {
          // allocate local constraints map
          CL cl;

          // extract constraints type
          typedef typename LFS::Traits::ConstraintsType C;
          const C & c = lfs.constraints();

          // iterate over boundary, need intersection iterator
          ConstraintsCallVolume<C,C::doVolume>::volume(c,p,eg,lfs,cl);
        }

        // reuse constraints parameter information from p for all LFS children
        template<typename P, typename LFS, typename TreePath>
        typename enable_if<P::isLeaf && (!LFS::isLeaf)>::type
        leaf(const P& p, const LFS& lfs, TreePath treePath) const
        {
          // traverse LFS tree and reuse parameter information
          TypeTree::applyToTree(lfs,VolumeConstraintsForParametersLeaf<P,EG,CL>(p,eg,cl));
        }

        VolumeConstraints(const EG& eg_, CL& cl_)
          : eg(eg_)
          , cl(cl_)
        {}

      private:
        const EG& eg;
        CL& cl;

      };


    } // anonymous namespace

    template<typename... Children>
    struct CompositeConstraintsOperator
      : public TypeTree::CompositeNode<Children...>
    {
      typedef TypeTree::CompositeNode<Children...> BaseT;

      CompositeConstraintsOperator(Children&... children)
        : BaseT(children...)
      {}

      CompositeConstraintsOperator(std::shared_ptr<Children>... children)
        : BaseT(children...)
      {}

      // aggregate flags

      // forward methods to childs

    };

    template<typename... Children>
    struct CompositeConstraintsParameters
      : public TypeTree::CompositeNode<Children...>
    {
      typedef TypeTree::CompositeNode<Children...> BaseT;

      CompositeConstraintsParameters(Children&... children)
        : BaseT(children...)
      {}

      CompositeConstraintsParameters(std::shared_ptr<Children>... children)
        : BaseT(children...)
      {}
    };

    template<typename T, std::size_t k>
    struct PowerConstraintsParameters
      : public TypeTree::PowerNode<T,k>
    {
      typedef TypeTree::PowerNode<T,k> BaseT;

      PowerConstraintsParameters()
        : BaseT()
      {}

      PowerConstraintsParameters(T& c)
        : BaseT(c)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1)
        : BaseT(c0,c1)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2)
        : BaseT(c0,c1,c2)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3)
        : BaseT(c0,c1,c2,c3)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {}

      PowerConstraintsParameters (T& c0,
                              T& c1,
                              T& c2,
                              T& c3,
                              T& c4,
                              T& c5,
                              T& c6,
                              T& c7,
                              T& c8,
                              T& c9)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8,c9)
      {}

      PowerConstraintsParameters (const array<std::shared_ptr<T>,k>& children)
        : BaseT(children)
      {}
    };

#ifndef DOXYGEN

    //! wraps a BoundaryGridFunction in a OldStyleConstraintsParameter class
    template<typename F>
    class OldStyleConstraintsWrapper
      : public TypeTree::LeafNode
    {
      std::shared_ptr<const F> _f;
      unsigned int _i;
    public:

      template<typename Transformation>
      OldStyleConstraintsWrapper(std::shared_ptr<const F> f, const Transformation& t, unsigned int i=0)
        : _f(f)
        , _i(i)
      {}

      template<typename Transformation>
      OldStyleConstraintsWrapper(const F & f, const Transformation& t, unsigned int i=0)
        : _f(stackobject_to_shared_ptr(f))
        , _i(i)
      {}

      template<typename I>
      bool isDirichlet(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        typename F::Traits::RangeType bctype;
        _f->evaluate(intersection,coord,bctype);
        return bctype[_i] > 0;
      }

      template<typename I>
      bool isNeumann(const I & intersection, const FieldVector<typename I::ctype, I::dimension-1> & coord) const
      {
        typename F::Traits::RangeType bctype;
        _f->evaluate(intersection,coord,bctype);
        return bctype[_i] == 0;
      }
    };

    // Tag to name trafo GridFunction -> OldStyleConstraintsWrapper
    struct gf_to_constraints {};

    // register trafos GridFunction -> OldStyleConstraintsWrapper
    template<typename F, typename Transformation>
    struct MultiComponentOldStyleConstraintsWrapperDescription
    {

      static const bool recursive = false;

      enum { dim = F::Traits::dimRange };
      typedef OldStyleConstraintsWrapper<F> node_type;
      typedef PowerConstraintsParameters<node_type, dim> transformed_type;
      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const F& s, const Transformation& t)
      {
        std::shared_ptr<const F> sp = stackobject_to_shared_ptr(s);
        array<std::shared_ptr<node_type>, dim> childs;
        for (int i=0; i<dim; i++)
          childs[i] = std::make_shared<node_type>(sp,t,i);
        return transformed_type(childs);
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const F> s, const Transformation& t)
      {
        array<std::shared_ptr<node_type>, dim> childs;
        for (int i=0; i<dim; i++)
          childs[i] = std::make_shared<node_type>(s,t,i);
        return std::make_shared<transformed_type>(childs);
      }

    };
    // trafos for leaf nodes
    template<typename GridFunction>
    typename conditional<
      (GridFunction::Traits::dimRange == 1),
      // trafo for scalar leaf nodes
      TypeTree::GenericLeafNodeTransformation<GridFunction,gf_to_constraints,OldStyleConstraintsWrapper<GridFunction> >,
      // trafo for multi component leaf nodes
      MultiComponentOldStyleConstraintsWrapperDescription<GridFunction,gf_to_constraints>
      >::type
    registerNodeTransformation(GridFunction*, gf_to_constraints*, GridFunctionTag*);

    // trafo for power nodes
    template<typename PowerGridFunction>
    TypeTree::SimplePowerNodeTransformation<PowerGridFunction,gf_to_constraints,PowerConstraintsParameters>
    registerNodeTransformation(PowerGridFunction*, gf_to_constraints*, PowerGridFunctionTag*);

    // trafos for composite nodes
    template<typename CompositeGridFunction>
    TypeTree::SimpleCompositeNodeTransformation<CompositeGridFunction,gf_to_constraints,CompositeConstraintsParameters>
    registerNodeTransformation(CompositeGridFunction*, gf_to_constraints*, CompositeGridFunctionTag*);

    //! construct constraints
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \tparam P   Type implementing a constraints parameter tree
     * \tparam GFS Type implementing the model GridFunctionSpace
     * \tparam CG  Type implementing the model
     *             GridFunctionSpace::ConstraintsContainer::Type
     * \tparam isFunction bool to identify old-style parameters, which were implemented the Dune::PDELab::FunctionInterface
     */
    template<typename P, typename GFS, typename GV, typename CG, bool isFunction>
    struct ConstraintsAssemblerHelper
    {
      //! construct constraints from given boundary condition function
      /**
       * \code
       * #include <dune/pdelab/constraints/common/constraints.hh>
       * \endcode
       * \tparam P   Type implementing a constraints parameter tree
       * \tparam GFS Type implementing the model GridFunctionSpace
       * \tparam CG  Type implementing the model
       *             GridFunctionSpace::ConstraintsContainer::Type
       *
       * \param p       The parameter object
       * \param gfs     The gridfunctionspace
       * \param cg      The constraints container
       * \param verbose Print information about the constaints at the end
       */
      static void
      assemble(const P& p, const GFS& gfs, const GV& gv, CG& cg, const bool verbose)
      {
        // get some types
        typedef typename GV::Traits::template Codim<0>::Entity Element;
        typedef typename GV::Intersection Intersection;

        // make local function space
        typedef LocalFunctionSpace<GFS> LFS;
        LFS lfs_e(gfs);
        LFSIndexCache<LFS> lfs_cache_e(lfs_e);
        LFS lfs_f(gfs);
        LFSIndexCache<LFS> lfs_cache_f(lfs_f);

        // get index set
        const typename GV::IndexSet& is=gv.indexSet();

        // helper to compute offset dependent on geometry type
        const int chunk=1<<28;
        int offset = 0;
        std::map<Dune::GeometryType,int> gtoffset;

        // loop once over the grid
        for (const auto& cell : elements(gv))
        {
          // assign offset for geometry type;
          if (gtoffset.find(cell.type())==gtoffset.end())
          {
            gtoffset[cell.type()] = offset;
            offset += chunk;
          }

          const typename GV::IndexSet::IndexType id = is.index(cell)+gtoffset[cell.type()];

          // bind local function space to element
          lfs_e.bind(cell);

          typedef typename CG::LocalTransformation CL;

          CL cl_self;

          typedef ElementGeometry<Element> ElementWrapper;
          TypeTree::applyToTreePair(p,lfs_e,VolumeConstraints<ElementWrapper,CL>(ElementWrapper(cell),cl_self));

          // iterate over intersections and call metaprogram
          unsigned int intersection_index = 0;
          for (const auto& intersection : intersections(gv,cell))
          {
            if (intersection.boundary())
            {
              typedef IntersectionGeometry<Intersection> IntersectionWrapper;
              TypeTree::applyToTreePair(p,lfs_e,BoundaryConstraints<IntersectionWrapper,CL>(IntersectionWrapper(intersection,intersection_index),cl_self));
            }

            // ParallelStuff: BEGIN support for processor boundaries.
            if ((!intersection.boundary()) && (!intersection.neighbor()))
            {
              typedef IntersectionGeometry<Intersection> IntersectionWrapper;
              TypeTree::applyToTree(lfs_e,ProcessorConstraints<IntersectionWrapper,CL>(IntersectionWrapper(intersection,intersection_index),cl_self));
            }
            // END support for processor boundaries.

            if (intersection.neighbor()){

              auto outside_cell = intersection.outside();
              Dune::GeometryType gtn = outside_cell.type();
              const typename GV::IndexSet::IndexType idn = is.index(outside_cell)+gtoffset[gtn];

              if(id>idn){
                // bind local function space to element in neighbor
                lfs_f.bind(outside_cell);

                CL cl_neighbor;

                typedef IntersectionGeometry<Intersection> IntersectionWrapper;
                TypeTree::applyToTreePair(lfs_e,lfs_f,SkeletonConstraints<IntersectionWrapper,CL>(IntersectionWrapper(intersection,intersection_index),cl_self,cl_neighbor));

                if (!cl_neighbor.empty())
                  {
                    lfs_cache_f.update();
                    cg.import_local_transformation(cl_neighbor,lfs_cache_f);
                  }

              }
            }
          }

          if (!cl_self.empty())
            {
              lfs_cache_e.update();
              cg.import_local_transformation(cl_self,lfs_cache_e);
            }

        }

        // print result
        if(verbose){
          std::cout << "constraints:" << std::endl;

          std::cout << cg.size() << " constrained degrees of freedom" << std::endl;

          for (const auto& col : cg)
          {
            std::cout << col.first << ": "; // col index
            for (const auto& row : col.second)
              std::cout << "(" << row.first << "," << row.second << ") "; // row index , value
            std::cout << std::endl;
          }
        }
      }
    }; // end ConstraintsAssemblerHelper



    // Disable constraints assembly for empty transformation
    template<typename F, typename GFS, typename GV>
    struct ConstraintsAssemblerHelper<F, GFS, GV, EmptyTransformation, true>
    {
      static void assemble(const F& f, const GFS& gfs, const GV& gv, EmptyTransformation& cg, const bool verbose)
      {}
    };

    // Disable constraints assembly for empty transformation
    template<typename F, typename GFS, typename GV>
    struct ConstraintsAssemblerHelper<F, GFS, GV, EmptyTransformation, false>
    {
      static void assemble(const F& f, const GFS& gfs, const GV& gv, EmptyTransformation& cg, const bool verbose)
      {}
    };



    // Backwards compatibility shim
    template<typename F, typename GFS, typename GV, typename CG>
    struct ConstraintsAssemblerHelper<F, GFS, GV, CG, true>
    {
      static void
      assemble(const F& f, const GFS& gfs, const GV& gv, CG& cg, const bool verbose)
      {
        // type of transformed tree
        typedef typename TypeTree::TransformTree<F,gf_to_constraints> Transformation;
        typedef typename Transformation::Type P;
        // transform tree
        P p = Transformation::transform(f);
        // call parameter based implementation
        ConstraintsAssemblerHelper<P, GFS, GV, CG, IsGridFunction<P>::value>::assemble(p,gfs,gv,cg,verbose);
      }
    };
#endif

    //! construct constraints
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \tparam GFS Type implementing the model GridFunctionSpace
     * \tparam CG  Type implementing the model
     *             GridFunctionSpace::ConstraintsContainer::Type
     *
     * \param gfs     The gridfunctionspace
     * \param cg      The constraints container
     * \param verbose Print information about the constaints at the end
     */
    template<typename GFS, typename CG>
    void constraints(const GFS& gfs, CG& cg,
                     const bool verbose = false)
    {
      typedef typename GFS::Traits::GridViewType GV;
      NoConstraintsParameters p;
      ConstraintsAssemblerHelper<NoConstraintsParameters, GFS, GV, CG, false>::assemble(p,gfs,gfs.gridView(),cg,verbose);
    }

    //! construct constraints from given constraints parameter tree
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \tparam P   Type implementing a constraits parameter tree
     * \tparam GFS Type implementing the model GridFunctionSpace
     * \tparam CG  Type implementing the model
     *             GridFunctionSpace::ConstraintsContainer::Type
     *
     * \param p       The condition parameters
     * \param gfs     The gridfunctionspace
     * \param cg      The constraints container
     * \param verbose Print information about the constaints at the end
     *
     * \note For backwards compatibility you can implement the parameter tree as an Implemention Dune::PDELab::FunctionInterface
     *
     */
    template<typename P, typename GFS, typename CG>
    void constraints(const P& p, const GFS& gfs, CG& cg,
                     const bool verbose = false)
    {
      typedef typename GFS::Traits::GridViewType GV;
      // clear global constraints
      cg.clear();
      ConstraintsAssemblerHelper<P, GFS, GV, CG, IsGridFunction<P>::value>::assemble(p,gfs,gfs.gridView(),cg,verbose);
    }

    //! construct constraints from given boundary condition function
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \tparam CG Type of ConstraintsContainer
     * \tparam XG Type of coefficients container
     *
     * \param cg The ConstraintsContainer
     * \param x  The value to assign
     * \param xg The container with the coefficients
     */
    template<typename CG, typename XG>
    void set_constrained_dofs(const CG& cg,
                              typename XG::ElementType x,
                              XG& xg)
    {
      for (const auto& col : cg)
        xg[col.first] = x;
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void set_constrained_dofs(const EmptyTransformation& cg,
                              typename XG::ElementType x,
                              XG& xg)
    {}

#endif // DOXYGEN


    //! check that constrained dofs match a certain value
    /**
     * as if they were set by set_constrained_dofs()
     *
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     *
     * \tparam CG  Type of ConstraintsContainer
     * \tparam XG  Type of coefficients container
     * \tparam Cmp Type of Comparison object to use (something with the
     *             interface of Dune::FloatCmpOps)
     *
     * \param cg  The ConstraintsContainer
     * \param x   The value to compare with
     * \param xg  The container with the coefficients
     * \param cmp The comparison object to use.
     *
     * \returns true if all constrained dofs match the given value, false
     *          otherwise.
     */
    template<typename CG, typename XG, typename Cmp>
    bool check_constrained_dofs(const CG& cg, typename XG::ElementType x,
                                XG& xg, const Cmp& cmp = Cmp())
    {
      for (const auto& col : cg)
        if(cmp.ne(xg[col.first], x))
          return false;
      return true;
    }

    //! check that constrained dofs match a certain value
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     *
     * This just calls check_constrained_dofs(cg, x, xg, cmp) with a default
     * comparison object FloatCmpOps<typename XG::ElementType>().
     *
     * \tparam CG  Type of ConstraintsContainer
     * \tparam XG  Type of coefficients container
     *
     * \param cg  The ConstraintsContainer
     * \param x   The value to compare with
     * \param xg  The container with the coefficients
     *
     * \returns true if all constrained dofs match the given value, false
     *          otherwise.
     */
    template<typename CG, typename XG>
    bool check_constrained_dofs(const CG& cg, typename XG::ElementType x,
                                XG& xg)
    {
      return check_constrained_dofs(cg, x, xg,
                                    FloatCmpOps<typename XG::ElementType>());
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG, typename Cmp>
    bool check_constrained_dofs(const EmptyTransformation& cg, typename XG::ElementType x,
                                XG& xg, const Cmp& cmp = Cmp())
    {
      return true;
    }

    // Specialized version for unconstrained spaces
    template<typename XG>
    bool check_constrained_dofs(const EmptyTransformation& cg, typename XG::ElementType x,
                                XG& xg)
    {
      return true;
    }

#endif // DOXYGEN


    //! transform residual into transformed basis: r -> r~
    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     */
    template<typename CG, typename XG>
    void constrain_residual (const CG& cg, XG& xg)
    {
      for (const auto& col : cg)
        for (const auto& row : col.second)
          xg[row.first] += row.second * xg[col.first];

      // extra loop because constrained dofs might have contributions
      // to constrained dofs
      for (const auto& col : cg)
        xg[col.first] = typename XG::ElementType(0);
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void constrain_residual (const EmptyTransformation& cg, XG& xg)
    {}

#endif // DOXYGEN

    //! Modify coefficient vector based on constrained dofs as given
    //! in the constraints container
    //! @{

    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     */
    template<typename CG, typename XG>
    void copy_constrained_dofs (const CG& cg, const XG& xgin, XG& xgout)
    {
      for (const auto& col : cg)
        xgout[col.first] = xgin[col.first];
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void copy_constrained_dofs (const EmptyTransformation& cg, const XG& xgin, XG& xgout)
    {}

#endif // DOXYGEN


    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \todo reimplement; this is horribly inefficient
     */
    template<typename CG, typename XG>
    void set_nonconstrained_dofs (const CG& cg, typename XG::ElementType x, XG& xg)
    {
      XG tmp(xg);
      xg = x;
      copy_constrained_dofs(cg,tmp,xg);
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void set_nonconstrained_dofs (const EmptyTransformation& cg, typename XG::ElementType x, XG& xg)
    {
      xg = x;
    }

#endif // DOXYGEN


    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \todo reimplement; this is horribly inefficient
     */
    template<typename CG, typename XG>
    void copy_nonconstrained_dofs (const CG& cg, const XG& xgin, XG& xgout)
    {
      XG tmp(xgin);
      copy_constrained_dofs(cg,xgout,tmp);
      xgout = tmp;
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void copy_nonconstrained_dofs (const EmptyTransformation& cg, const XG& xgin, XG& xgout)
    {
      xgout = xgin;
    }

#endif // DOXYGEN


    /**
     * \code
     * #include <dune/pdelab/constraints/common/constraints.hh>
     * \endcode
     * \todo reimplement; this is horribly inefficient
     */
    template<typename CG, typename XG>
    void set_shifted_dofs (const CG& cg, typename XG::ElementType x, XG& xg)
    {
      XG tmp(xg);
      tmp = x;

      for (const auto& col : cg)
      {
        // this is our marker for Dirichlet constraints
        if (col.second.size() == 0)
        {
          tmp[col.first] = xg[col.first];
        }
      }

      xg = tmp;
    }


#ifndef DOXYGEN

    // Specialized version for unconstrained spaces
    template<typename XG>
    void set_shifted_dofs (const EmptyTransformation& cg, typename XG::ElementType x, XG& xg)
    {}

#endif // DOXYGEN

    //! @}

    //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_CONSTRAINTS_COMMON_CONSTRAINTS_HH
