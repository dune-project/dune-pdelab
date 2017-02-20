//-*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FUNCTION_HH
#define DUNE_PDELAB_FUNCTION_HH

#include <iostream>
#include <sstream>

#include <dune/common/deprecated.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/typetraits.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dune/grid/utility/hierarchicsearch.hh>

#include <dune/typetree/typetree.hh>

#include "vtkexport.hh"
#include "geometrywrapper.hh"
#include "typetraits.hh"

namespace Dune {
  namespace PDELab {

    //! \addtogroup PDELab_Function Function
    //! \ingroup PDELab
    //! \{

    //! traits class holding function signature, same as in local function
    //! \tparam DF The numeric type of the field representing the domain.
    //! \tparam dimension of the domain.
    //! \tparam D The type of the domain.
    //! \tparam m The dimension of the range.
    //! \tparam RF The numeric type of the field representing the range.
    //! \tparam R The type of the range.
    template<class DF, int n, class D, class RF, int m, class R>
    struct FunctionTraits
    {
      //! \brief Export type for domain field
      typedef DF DomainFieldType;

      //! \brief Enum for domain dimension
      enum {
        //! \brief dimension of the domain
        dimDomain = n
      };

      //! \brief domain type in dim-size coordinates
      typedef D DomainType;

      //! \brief Export type for range field
      typedef RF RangeFieldType;

      //! \brief Enum for range dimension
      enum {
        //! \brief dimension of the range
        dimRange = m
      };

      //! \brief range type
      typedef R RangeType;
    };

    //! \brief a Function that maps x in DomainType to y in RangeType
    //! \tparam T The type of the function traits
    //! \tparam Imp The type implementing the interface.
    template<class T, class Imp>
    class FunctionInterface
    {
    public:
      //! \brief Export type traits
      typedef T Traits;

      /** \brief Evaluate all basis function at given position

          Evaluates all shape functions at the given position and returns
          these values in a vector.
      */
      inline void evaluate (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        asImp().evaluate(x,y);
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //! \brief Default class for additional methods in instationary functions
    class InstationaryFunctionDefaults
    {
    public:
      //! set time for subsequent evaluation
      /**
       * This method sets the time for subsequent calls to any of the
       * evaluation methods.
       *
       * \note This default method does nothing, it just ensures setTime() can
       *       be called without ill effects.
       * \note Function implementation are free to restrict the types of
       *       acceptable parameters.  This should be noted in the function
       *       classes documentation.
       */
      template<typename Time>
      inline void setTime(Time t)
      { }
    };

    //! \brief GV The type of the grid view the function lives on.
    template<typename GV>
    struct PowerCompositeGridFunctionTraits
    {
      //! \brief The type of the grid view the function lives on.
      typedef GV GridViewType;

      //! \brief codim 0 entity
      typedef typename GV::Traits::template Codim<0>::Entity ElementType;

    };


    //! Mixin base class for specifying output hints to I/O routines like VTK.
    class GridFunctionOutputParameters
    {

    public:

      //! Namespace for output-related data types and enums.
      struct Output
      {
        //! The type of the data set.
        /**
         * This information can be used by a VTKWriter to pick the correct
         * VTK data set type.
         */
        enum DataSetType
          {
            vertexData, //!< A data set with vertex values.
            cellData    //!< A data set with cell values.
          };
      };

      //! Standard constructor.
      /**
       * \param dataSetType The type of the data set represented by this function.
       */
      GridFunctionOutputParameters(Output::DataSetType dataSetType = Output::vertexData)
        : _dataSetType(dataSetType)
      {}

      //! Return the data set type of this function.
      Output::DataSetType dataSetType() const
      {
        return _dataSetType;
      }

      //! Set the data set type of this function.
      void setDataSetType(Output::DataSetType dataSetType)
      {
        _dataSetType = dataSetType;
      }

    private:

      Output::DataSetType _dataSetType;

    };

    //! \brief traits class holding the function signature, same as in local function
    //! \brief GV The type of the grid view the function lives on.
    //! \brief RF The numeric type used in the range of the function.
    //! \brief m The dimension of the range.
    //! \tparam R The numeric type of the field representing the range.
    template<class GV, class RF, int m, class R>
    struct GridFunctionTraits
      : public FunctionTraits<typename GV::Grid::ctype, GV::dimension,
                              Dune::FieldVector<typename GV::Grid::ctype,
                                                GV::dimension>,
                              RF, m, R>
      , public PowerCompositeGridFunctionTraits<GV>
    {
    };

    //! \brief a GridFunction maps x in DomainType to y in RangeType
    template<class T, class Imp>
    class GridFunctionInterface
      : public GridFunctionOutputParameters
    {
    public:
      //! \brief Export type traits
      typedef T Traits;

      GridFunctionInterface(Output::DataSetType dataSetType = Output::vertexData)
        : GridFunctionOutputParameters(dataSetType)
      {}

      /** \brief Evaluate the GridFunction at given position

          Evaluates components of the grid function at the given position and
          returns these values in a vector.

          \param[in]  e The entity to evaluate on
          \param[in]  x The position in entity-local coordinates
          \param[out] y The result of the evaluation
      */
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        asImp().evaluate(e,x,y);
      }

      //! \brief get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return asImp().getGridView();
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //! \brief traits class holding function signature, same as in local function
    //! \tparam GV The type of the grid view the function lives on.
    //! \tparam RF The numeric type of the field representing the range.
    //! \tparam m The dimension of the range.
    //! \tparam R The type of the range.
    template<class GV, class RF, int m, class R>
    struct BoundaryGridFunctionTraits
      : public FunctionTraits<typename GV::Grid::ctype, GV::dimension-1,
                              Dune::FieldVector<typename GV::Grid::ctype,
                                                GV::dimension-1>,
                              RF, m, R>
    {
      //! \brief Export grid view type in addition
      typedef GV GridViewType;
    };


    //! \brief A BoundaryGridFunction allows evaluation on boundary intersections
    // \tparam T The type of the BoundaryGridFunctionTraits.
    // \tparam Imp The type of the implementing class.
    template<class T, class Imp>
    class BoundaryGridFunctionInterface
    {
    public:
      //! \brief Export type traits of the boundary grid function.
      typedef T Traits;

      /** \brief Evaluate the GridFunction at given position

          Evaluates components of the grid function at the given position and
          returns these values in a vector.

          \param[in]  ig geometry of intersection with boundary
          \param[in]  x The position in entity-local coordinates
          \param[out] y The result of the evaluation
      */
      template<typename I>
      inline void evaluate (const IntersectionGeometry<I>& ig,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        asImp().evaluate(ig,x,y);
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return asImp().getGridView();
      }

    private:
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };

    //============================
    // Function tree
    //============================

    //! \addtogroup GridFunctionTree
    //! \{

    struct GridFunctionTag {};

    /** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
    template<class T, class Imp>
    class GridFunctionBase
      : public GridFunctionInterface<T,Imp>
      , public TypeTree::LeafNode
    {
      using Base = GridFunctionInterface<T,Imp>;
    public:
      typedef GridFunctionTag ImplementationTag;
      //! Type of the GridView
      typedef typename T::GridViewType GridViewType;

      using Output = typename Base::Output;

      GridFunctionBase(typename Output::DataSetType dataSetType = Output::vertexData)
        : Base(dataSetType)
      {}
    };


    /** \brief leaf of a function tree
     *
     *  Classes derived from this class implement a \ref GridFunctionTree.
     *
     *  \tparam T   Traits class holding the functions signature
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              GridFunctionBase in some way (Barton-Nackman-Trick).
     */
    template<class T, class Imp>
    class BoundaryGridFunctionBase
      : public BoundaryGridFunctionInterface<T,Imp>
      , public TypeTree::LeafNode
    {
    public:
      typedef GridFunctionTag ImplementationTag;
      //! Type of the GridView
      typedef typename T::GridViewType GridViewType;
    };


    /** \brief Visitor for Power- and CompositeGridFunctions calling
        the setTime() method on the leafs of the corresponding
        function trees.

        \tparam Scalar type representing time.
    */
    template<typename TT>
    struct PowerCompositeSetTimeVisitor
      : public TypeTree::TreeVisitor, public TypeTree::DynamicTraversal
    {
      TT time;
      PowerCompositeSetTimeVisitor(const TT time_) : time(time_) {}

      template<typename LeafNode, typename TreePath>
      void leaf(LeafNode& node, TreePath treePath) const
      {
        node.setTime(time);
      }
    };

    struct PowerGridFunctionTag {};

    /** \brief product of identical functions
     *
     *  This collects k instances of T in a \ref GridFunctionTree.
     *
     *  \tparam T The type of the children of this node in the tree.
     *  \tparam k The number of children this node has.
     */
    template<class T, std::size_t k>
    class PowerGridFunction
      : public TypeTree::PowerNode<T,k>
    {

      typedef TypeTree::PowerNode<T,k> BaseT;

    public:

      typedef PowerCompositeGridFunctionTraits<typename T::GridViewType> Traits;

      typedef PowerGridFunctionTag ImplementationTag;

      //! record the GridView
      typedef typename T::GridViewType GridViewType;

      //! Set the time in all leaf nodes of this function tree
      template <typename TT>
      void setTime(TT time){
        PowerCompositeSetTimeVisitor<TT> visitor(time);
        TypeTree::applyToTree(*this,visitor);
      }

      PowerGridFunction()
      {}

      //! Construct a PowerGridFunction with k clones of the function t
      PowerGridFunction (T& t)
        : BaseT(t) {}

      /** \brief Initialize all children with different function objects
       *
       *  This constructor is only available in the non-specialized version
       *
       *  \param t Points to an array of pointers to function objects of type
       *           T.  The function pointed to by the first pointer will be
       *           used to initialize the first child, the second pointer for
       *           the second child and so on.
       */
      // TODO: PowerGridFunction (T** t) : ...

#ifdef DOXYGEN
      /** \brief Initialize all children with different function objects
       *
       *  Currently there exist specializations for 2 <= k <= 9.  Each
       *  specialization has a constructor which takes the initializers for
       *  its children as arguments.
       *
       *  @param t0 The initializer for the first child.
       *  @param t1 The initializer for the second child.
       *  @param ... more initializers
       */
      PowerGridFunction (T& t0, T& t1, ...)
      {
      }

#else

      PowerGridFunction (T& c0,
                         T& c1)
        : BaseT(c0,c1)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2)
        : BaseT(c0,c1,c2)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3)
        : BaseT(c0,c1,c2,c3)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4)
        : BaseT(c0,c1,c2,c3,c4)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5)
        : BaseT(c0,c1,c2,c3,c4,c5)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6)
        : BaseT(c0,c1,c2,c3,c4,c5,c6)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6,
                         T& c7)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7)
      {
      }

      PowerGridFunction (T& c0,
                         T& c1,
                         T& c2,
                         T& c3,
                         T& c4,
                         T& c5,
                         T& c6,
                         T& c7,
                         T& c8)
        : BaseT(c0,c1,c2,c3,c4,c5,c6,c7,c8)
      {
      }

      PowerGridFunction (T& c0,
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
      {
      }

#endif // DOXYGEN
    };

    struct CompositeGridFunctionTag {};

    /** \brief composite functions
     *
     *  Collect instances of possibly different function types Tn within a
     *  \ref GridFunctionTree.  This impolements a \ref GridFunctionTree
     *
     *  \tparam Tn The base types.  Tn==EmptyChild means that slot n is
     *             unused.  Currently, up to 9 slots are supported, making 8
     *             the maximum n.
     */
    template<typename... Children>
    class CompositeGridFunction
      : public TypeTree::CompositeNode<Children...>
    {

      typedef TypeTree::CompositeNode<Children...> BaseT;

    public:

      typedef CompositeGridFunctionTag ImplementationTag;

      typedef PowerCompositeGridFunctionTraits<typename BaseT::template Child<0>::Type::GridViewType> Traits;

      //! record the GridView
      typedef typename BaseT::template Child<0>::Type::GridViewType GridViewType;

      CompositeGridFunction()
      {}

      CompositeGridFunction (Children&... children)
        : BaseT(TypeTree::assertGridViewType<typename BaseT::template Child<0>::Type>(children)...)
      {
      }

      //! Set the time in all leaf nodes of this function tree
      template <typename TT>
      void setTime(TT time){
        PowerCompositeSetTimeVisitor<TT> visitor(time);
        TypeTree::applyToTree(*this,visitor);
      }

#ifdef DOXYGEN
      /** \brief Initialize all children
       *
       *  @param t0 The initializer for the first child.
       *  @param t1 The initializer for the second child.
       *  @param ... more initializers
       *
       *  The actual number of arguments for this constructor corresponds to
       *  the number of slots used in the template parameter list of the class.
       */
      CompositeGridFunction (T0& t0, T1& t1, ...) {}
#endif //DOXYGEN
    };

    //========================================================
    // helper template to turn an ordinary GridFunction into a
    // GridFunctionTree leaf
    //========================================================
    //! Turn an ordinary GridFunction into a GridFunctionTree leaf
    /**
     *  \tparam Imp Class implementing the function.
     */
    template<class Imp>
    class GridFunctionBaseAdapter
      : public GridFunctionBase<typename Imp::Traits,
                                GridFunctionBaseAdapter<Imp> >
    {
      const Imp &imp;

    public:
      //! construct a GridFunctionBaseAdapter
      /**
       * \param imp_ The underlying ordinary GridFunction.  A reference to
       *             this Object is stored, so the object must be valid for as
       *             long as this GridFunctionBaseAdapter is used.
       */
      GridFunctionBaseAdapter(const Imp& imp_)
        : imp(imp_)
      { }

      //! Evaluate the GridFunction at given position
      /**
       * Evaluates components of the grid function at the given position and
       * returns these values in a vector.
       *
       * \param[in]  e The entity to evaluate on
       * \param[in]  x The position in entity-local coordinates
       * \param[out] y The result of the evaluation
       */
      inline void evaluate (const typename Imp::Traits::ElementType& e,
                            const typename Imp::Traits::DomainType& x,
                            typename Imp::Traits::RangeType& y) const
      {
        imp.evaluate(e,x,y);
      }

      //! get a reference to the GridView
      inline const typename Imp::Traits::GridViewType& getGridView () const
      {
        return imp.getGridView();
      }
    };

    //=======================================
    // helper template for analytic functions
    //=======================================

    //! function signature for analytic functions on a grid
    template<typename GV, typename RF, int m>
    struct AnalyticGridFunctionTraits
      : public GridFunctionTraits<GV, RF, m, Dune::FieldVector<RF,m> >
    {
    };

    /** \brief an analytic grid function
     *
     *  This is a convenience class which eases the creation of analytic
     *  GridFunctions.  Classes derived from it need only implement a method
     *  evaluateGlobal(const Dune::FieldVector<typename Traits::DomainFieldType,GV::dimensionworld> &x_global, RangeType &y) to have a
     *  full-fledged GridFunction.
     *
     *  \tparam T   The Traits class
     *  \tparam Imp Class implementing the function.  Imp must be derived from
     *              AnalyticGridFunctionBase in some way
     *              (Barton-Nackman-Trick).
     */
    template<typename T, typename Imp>
    class AnalyticGridFunctionBase
      : public GridFunctionBase<T,AnalyticGridFunctionBase<T,Imp> >
    {
    public:
      typedef T Traits;

      //! Construct an Analytic GridFunctionBase given a GridView g_
      AnalyticGridFunctionBase (const typename Traits::GridViewType& g_) : g(g_) {}

      //! \copydoc GridFunctionBase::evaluate()
      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        asImp().evaluateGlobal(e.geometry().global(x),y);
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return g;
      }

    private:
      typename Traits::GridViewType g;
      Imp& asImp () {return static_cast<Imp &> (*this);}
      const Imp& asImp () const {return static_cast<const Imp &>(*this);}
    };


    // Adapter takes a vector-valued grid function and provides evaluation
    // of normal flux on the interior of faces.
    template<typename T>
    class NormalFluxGridFunctionAdapter
      : public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,
                                                                                    typename T::Traits::RangeFieldType,
                                                                                    1,
                                                                                    Dune::FieldVector<typename T::Traits::RangeFieldType,1>
                                                                                    >,
                                                   NormalFluxGridFunctionAdapter<T> >
    {
    public:
      typedef Dune::PDELab::GridFunctionTraits<typename T::Traits::GridViewType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > Traits;
      typedef Dune::PDELab::GridFunctionInterface<Traits,NormalFluxGridFunctionAdapter<T> > BaseT;

      NormalFluxGridFunctionAdapter (const T& t_) : t(stackobject_to_shared_ptr(t_)) {}


      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // ensure correct size
        static_assert((static_cast<int>(T::Traits::GridViewType::dimension)==static_cast<int>(T::Traits::dimRange)),"number of components must equal dimension");

        // evaluate velocity
        typename T::Traits::RangeType v;
        t->evaluate(e,x,v);

        // implementation only handles triangles so far
        if (!e.geometry().type().isTriangle())
          DUNE_THROW(Dune::NotImplemented, "only implemented for triangles");

        // start and end corner in local numbering
        int n0, n1;

        typename Traits::DomainType nu;

        // determine outer unit normal
        if (std::abs(x[0])<1E-10)
          {
            // edge 1
            n0 = 2;
            n1 = 0;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        if (std::abs(x[1])<1E-10)
          {
            // edge 2
            n0 = 0;
            n1 = 1;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        if (std::abs(x[0]+x[1]-1.0)<1E-10)
          {
            // edge 0
            n0 = 1;
            n1 = 2;

            nu = e.geometry().corner(n1);
            nu -= e.geometry().corner(n0);
            typename Traits::DomainFieldType temp = nu[0];
            nu[0] = nu[1];
            nu[1] = -temp;
            nu /= nu.two_norm();
            y = v[0]*nu[0]+v[1]*nu[1];
            return;
          }

        DUNE_THROW(Dune::Exception, "x needs to be on an edge");
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return t->getGridView();
      }

    private:
      shared_ptr<T const> t;
    };

    // Adapter takes a vector-valued grid function and applies
    // backward Piola transformation on each element
    template<typename T>
    class PiolaBackwardAdapter
      : public Dune::PDELab::GridFunctionBase<typename T::Traits,PiolaBackwardAdapter<T> >
    {
    public:
      typedef typename T::Traits::GridViewType GridViewType;
      typedef typename T::Traits Traits;
      typedef Dune::PDELab::GridFunctionBase<Traits,PiolaBackwardAdapter<T> > BaseT;
      // typedef GridFunctionTag ImplementationTag;

      PiolaBackwardAdapter (const T& t_) : t(stackobject_to_shared_ptr(t_)) {}


      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        // evaluate velocity
        typename T::Traits::RangeType v;
        t->evaluate(e,x,v);

        // apply Piola transformation
        typename Traits::ElementType::Geometry::JacobianInverseTransposed
          J = e.geometry().jacobianInverseTransposed(x);
        y = 0;
        J.umtv(v,y);
        y *= e.geometry().integrationElement(x);
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return t->getGridView();
      }

    private:
      shared_ptr<T const> t;
    };


    //==========================
    // template metaprograms
    //==========================

    namespace {

      //! implement VisitingFunctor for vtkwriter_tree_addvertexdata
      template<typename VTKWriter>
      struct AddGridFunctionsToVTKWriter
        : public TypeTree::TreeVisitor
        , public TypeTree::DynamicTraversal
      {

        VTKWriter& w;
        const std::string s;

        AddGridFunctionsToVTKWriter(VTKWriter& w_, const std::string & s_) :
          w(w_), s(s_) {}

        template<typename T, typename TreePath>
        void leaf(const T& t, TreePath treePath) {
          std::stringstream name;
          name << s;
          for (std::size_t i=0; i < treePath.size(); ++i)
            name << "_" << treePath.element(i);
          w.addVertexData(make_shared< VTKGridFunctionAdapter<T> >(t,name.str()));
        }
      };

    } // anonymous namespace

    /** \brief add vertex data from a \ref GridFunctionTree to a VTKWriter
     *
     *  \tparam GV The GridView for the VTKWriter
     *  \tparam T  The \ref GridFunctionTree
     */
    template<typename GV, typename T>
    void vtkwriter_tree_addvertexdata (Dune::VTKWriter<GV>& w, const T& t, std::string s = "data")
    {
      AddGridFunctionsToVTKWriter<Dune::VTKWriter<GV> > visitor(w,s);
      TypeTree::applyToTree(t,visitor);
    }

    //! \} GridFunctionTree

    //! \addtogroup PDELab_FunctionAdapters Function Adapters
    //! \{

    /** \brief make a GridFunction from a Function
     *
     *  \tparam G The GridView type
     *  \tparam T The function type
     */
    template<typename G, typename T>
    class FunctionToGridFunctionAdapter :
      public GridFunctionBase<GridFunctionTraits<
                                     G,
                                     typename T::Traits::RangeFieldType,
                                     T::Traits::dimRange,
                                     typename T::Traits::RangeType>,
                                   FunctionToGridFunctionAdapter<G,T> >
    {
    public:
      typedef GridFunctionTraits<G,
                                 typename T::Traits::RangeFieldType,
                                 T::Traits::dimRange,
                                 typename T::Traits::RangeType> Traits;
      static_assert(
                    (std::is_same<typename T::Traits::DomainFieldType,
                     typename Traits::DomainFieldType>::value),
                    "GridView's and wrapped Functions DomainFieldType don't match");
      static_assert(
                    T::Traits::dimDomain==Traits::dimDomain,
                    "GridView's and wrapped Functions dimDomain don't match");
      static_assert(
                    (std::is_same<typename T::Traits::DomainType,
                     typename Traits::DomainType>::value),
                    "GridView's and wrapped Functions DomainType don't match");

      /** \brief Create a FunctionToGridFunctionAdapter
       *
       *  \param g_ The GridView
       *  \param t_ The function
       */
      FunctionToGridFunctionAdapter (const G& g_, const T& t_) : g(g_), t(t_) {}

      inline void evaluate (const typename Traits::ElementType& e,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        t.evaluate(e.geometry().global(x),y);
      }

      inline const typename Traits::GridViewType& getGridView () const
      {
        return g;
      }

    private:
      G g;
      const T& t;
    };

    /** \brief make a Function from a GridFunction
     *
     *  \tparam GF The GridFunction type
     */
    template<typename GF>
    class GridFunctionToFunctionAdapter
      : public FunctionInterface<FunctionTraits<typename GF::Traits::GridViewType::ctype,
                                                GF::Traits::GridViewType::dimensionworld,
                                                Dune::FieldVector<typename GF::Traits::GridViewType::ctype,
                                                                  GF::Traits::GridViewType::dimensionworld
                                                                  >,
                                                typename GF::Traits::RangeFieldType,
                                                GF::Traits::dimRange,
                                                Dune::FieldVector<typename GF::Traits::RangeFieldType,
                                                                  GF::Traits::dimRange>
                                                >,
                                 GridFunctionToFunctionAdapter<GF> >
    {
    public:
      //! \brief Export type traits
      typedef FunctionTraits<typename GF::Traits::GridViewType::ctype,
                             GF::Traits::GridViewType::dimensionworld,
                             Dune::FieldVector<typename GF::Traits::GridViewType::ctype,
                                               GF::Traits::GridViewType::dimensionworld
                                               >,
                             typename GF::Traits::RangeFieldType,
                             GF::Traits::dimRange,
                             Dune::FieldVector<typename GF::Traits::RangeFieldType,
                                               GF::Traits::dimRange>
                             > Traits;

      //! make a GridFunctionToFunctionAdapter
      GridFunctionToFunctionAdapter(const GF &gf_)
        : gf(gf_)
        , hsearch(gf.getGridView().grid(), gf.getGridView().indexSet())
      { }

      /** \brief Evaluate all basis function at given position

          Evaluates all shape functions at the given position and returns
          these values in a vector.
      */
      inline void evaluate (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        typename GF::Traits::GridViewType::Grid::Traits::template Codim<0>::EntityPointer
          ep = hsearch.findEntity(x);
        gf.evaluate(*ep, ep->geometry().local(x), y);
      }

    private:
      const GF &gf;
      const Dune::HierarchicSearch<typename GF::Traits::GridViewType::Grid,
                                   typename GF::Traits::GridViewType::IndexSet> hsearch;
    };


    /** \brief make a Function in local coordinates from a Function in global coordinates
     *
     *  \tparam T Type of the global function
     *  \tparam E Type of the grid's element
     */
    template<typename T, typename E>
    class GlobalFunctionToLocalFunctionAdapter :
      public FunctionInterface<typename T::Traits,
                               GlobalFunctionToLocalFunctionAdapter<T,E> >
    {
    public:
      typedef typename T::Traits Traits;

      /** \brief Create a GlobalFunctionToLocalFunctionAdapter
       *
       *  \param t_ Global function
       *  \param e_ Grid's element where the local function is defined
       */
      GlobalFunctionToLocalFunctionAdapter (const T& t_, const E& e_) : t(t_), e(e_) {}

      /** \brief Evaluate the local function at the given position

          \param[in]  x The position in local coordinates
          \param[out] y The result of the evaluation
      */
      inline void evaluate (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        t.evaluate(e.geometry().global(x),y);
      }

    private:
      const T& t;
      const E& e;
    };


    /** \brief make a LocalFunction from a GridFunction using local coordinates
     *
     *  \tparam T type of the GridFunction
     */
    template<typename T> // T: GridFunction, E: Entity
    class GridFunctionToLocalFunctionAdapter :
      public FunctionInterface<typename T::Traits,
                               GridFunctionToLocalFunctionAdapter<T> >
    {
    public:
      typedef typename T::Traits Traits;

      /** \brief Create a GridFunctionToLocalFunctionAdapter
       *
       *  \param t_ GridFunction
       *  \param e_ Grid's element where the local function is defined
       */
      GridFunctionToLocalFunctionAdapter (const T& t_,
                                          const typename Traits::ElementType& e_)
        : t(t_), e(e_) {}

      /** \brief Evaluate the local function at the given position

          \param[in]  x The position in local coordinates
          \param[out] y The result of the evaluation
      */
      inline void evaluate (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        t.evaluate(e,x,y);
      }

    private:
      const T& t;
      const typename Traits::ElementType& e;
    };


    //! a Function maps x in DomainType to y in RangeType
    template<class T>
    class SelectComponentAdapter : public FunctionInterface<FunctionTraits<typename T::Traits::DomainFieldType,T::Traits::dimDomain,typename T::Traits::DomainType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > , SelectComponentAdapter<T> >
    {
      typedef FunctionInterface<FunctionTraits<typename T::Traits::DomainFieldType,T::Traits::dimDomain,typename T::Traits::DomainType,typename T::Traits::RangeFieldType,1,Dune::FieldVector<typename T::Traits::RangeFieldType,1> > , SelectComponentAdapter<T> > BaseT;
    public:
      //! \brief Export type traits
      typedef typename BaseT::Traits Traits;

      SelectComponentAdapter (const T& t_, int k_) : t(t_), k(k_) {}

      /** \brief Evaluate all basis function at given position

          Evaluates all shape functions at the given position and returns
          these values in a vector.
      */
      inline void evaluate (const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        typename T::Traits::RangeType Y;
        t.evaluate(x,Y);
        y = Y[k];
      }

      //! set component to be selected
      void select (int k_)
      {
        k = k_;
      }

    private:
      const T& t;
      int k;
    };

    //! Takes a BoundaryGridFunction and acts as a single component
    template<class T>
    class BoundaryGridFunctionSelectComponentAdapter
      : public BoundaryGridFunctionInterface<BoundaryGridFunctionTraits<typename T::Traits::GridViewType,
                                                                        typename T::Traits::RangeFieldType,1,
                                                                        Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                             BoundaryGridFunctionSelectComponentAdapter<T> >
    {
      typedef BoundaryGridFunctionInterface<BoundaryGridFunctionTraits<typename T::Traits::GridViewType,
                                                                       typename T::Traits::RangeFieldType,1,
                                                                       Dune::FieldVector<typename T::Traits::RangeFieldType,1> > ,
                                            BoundaryGridFunctionSelectComponentAdapter<T> > BaseT;
    public:
      //! \brief Export type traits
      typedef typename BaseT::Traits Traits;

      BoundaryGridFunctionSelectComponentAdapter (const T& t_, int k_) : t(t_), k(k_) {}

      /** \brief Evaluate all basis function at given position

          Evaluates all shape functions at the given position and returns
          these values in a vector.
      */
      template<typename I>
      inline void evaluate (const IntersectionGeometry<I>& ig,
                            const typename Traits::DomainType& x,
                            typename Traits::RangeType& y) const
      {
        typename T::Traits::RangeType Y;
        t.evaluate(ig,x,Y);
        y = Y[k];
      }

      //! get a reference to the GridView
      inline const typename Traits::GridViewType& getGridView () const
      {
        return t.getGridView();
      }


      //! set component to be selected
      void select (int k_)
      {
        k = k_;
      }

    private:
      const T& t;
      int k;
    };

    //! \}

    //! \} Function

  } // namespace PDELab
} // namespace Dune

#endif
