// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_COMMON_VTKEXPORT_HH
#define DUNE_PDELAB_COMMON_VTKEXPORT_HH

#include<cstddef>
#include<string>
#include<vector>

#include <dune/common/shared_ptr.hh>

#include<dune/grid/io/file/vtk/vtkwriter.hh>

#include<dune/pdelab/common/range.hh>

namespace Dune {
  namespace PDELab {

    //! wrap a GridFunction so it can be used with the VTKWriter from dune-grid.
    template<typename T> // T is a grid function
    class VTKGridFunctionAdapter
      : public Dune::VTKFunction<typename T::Traits::GridViewType>
    {
      typedef typename T::Traits::GridViewType::Grid::ctype DF;
      enum {n=T::Traits::GridViewType::dimension};
      typedef typename T::Traits::GridViewType::Grid::template Codim<0>::Entity Entity;

    public:
      //! construct a VTKGridFunctionAdapter
      /**
       * \param t_     GridFunction object to wrap.  A reference to the grid
       *               function object is stored internally and the
       *               constructed object becomes invalid as soon as that
       *               reference becomes invalid.
       * \param s_     Name of the field as returned by name().
       * \param remap_ How components are remapped between PDELab and VTK.
       *               The default value yields the identity map with entries
       *               (0,1,...,T::Traits::dimRange-1).
       *
       * The resulting VTKFunction will have \c remap_.size() components.
       * None of the elements of \c remap_ should be greater or equal to \c
       * T::Traits::dimRange.  When component \c c is requested by the
       * VTKWriter, it will be mapped to the component \c remap_[c] of the
       * GridFunction.
       */
      VTKGridFunctionAdapter(const T& t_, std::string s_,
                             const std::vector<std::size_t> &remap_ =
                               rangeVector(std::size_t(T::Traits::dimRange)))
        : t(stackobject_to_shared_ptr(t_)), s(s_), remap(remap_)
      {}

      //! construct a VTKGridFunctionAdapter
      /**
       * \param t_     Shared pointer to a GridFunction object to wrap.
       * \param s_     Name of the field as returned by name().
       * \param remap_ How components are remapped between PDELab and VTK.
       *               The default value yields the identity map with entries
       *               (0,1,...,T::Traits::dimRange-1).
       *
       * The resulting VTKFunction will have \c remap_.size() components.
       * None of the elements of \c remap_ should be greater or equal to \c
       * T::Traits::dimRange.  When component \c c is requested by the
       * VTKWriter, it will be mapped to the component \c remap_[c] of the
       * GridFunction.
       */
      VTKGridFunctionAdapter(const std::shared_ptr<const T>& t_, std::string s_,
                             const std::vector<std::size_t> &remap_ =
                               rangeVector(std::size_t(T::Traits::dimRange)))
        : t(t_), s(s_), remap(remap_)
      { }

      virtual int ncomps () const
      {
        return remap.size();;
      }

      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,n>& xi) const
      {
        typename T::Traits::DomainType x;
        typename T::Traits::RangeType y;

        for (int i=0; i<n; i++)
          x[i] = xi[i];
        t->evaluate(e,x,y);
        return y[remap[comp]];
      }

      virtual std::string name () const
      {
        return s;
      }

    private:
      std::shared_ptr<const T> t;
      std::string s;
      std::vector<std::size_t> remap;
    };

    //! construct a VTKGridFunctionAdapter
    /**
     * \param gf    GridFunction object to wrap.  A reference to the grid
     *              function object is stored internally and the constructed
     *              object becomes invalid as soon as that reference becomes
     *              invalid.
     * \param name  Name of the field as returned by name().
     * \param remap How components are remapped between PDELab and VTK.  The
     *              default value yields the identity map with entries
     *              (0,1,...,T::Traits::dimRange-1).
     *
     * The resulting VTKFunction will have \c remap_.size() components.  None
     * of the elements of \c remap should be greater or equal to \c
     * T::Traits::dimRange.  When component \c c is requested by the
     * VTKWriter, it will be mapped to the component \c remap[c] of the
     * GridFunction.
     */
    template<class GF>
    std::shared_ptr<VTKGridFunctionAdapter<GF> > makeVTKGridFunctionAdapter
    ( const std::shared_ptr<GF> &gf, const std::string &name,
      const std::vector<std::size_t> &remap =
        rangeVector(std::size_t(GF::Traits::dimRange)))
    { return std::make_shared<VTKGridFunctionAdapter<GF> >(gf, name, remap); }

    //! construct a VTKGridFunctionAdapter
    /**
     * \param gf    GridFunction object to wrap.  A reference to the grid
     *              function object is stored internally and the constructed
     *              object becomes invalid as soon as that reference becomes
     *              invalid.
     * \param name  Name of the field as returned by name().
     * \param remap How components are remapped between PDELab and VTK.  The
     *              default value yields the identity map with entries
     *              (0,1,...,T::Traits::dimRange-1).
     *
     * The resulting VTKFunction will have \c remap_.size() components.  None
     * of the elements of \c remap should be greater or equal to \c
     * T::Traits::dimRange.  When component \c c is requested by the
     * VTKWriter, it will be mapped to the component \c remap[c] of the
     * GridFunction.
     */
    template<class GF>
    std::shared_ptr<VTKGridFunctionAdapter<GF> > makeVTKGridFunctionAdapter
    ( const GF &gf, const std::string &name,
      const std::vector<std::size_t> &remap =
        rangeVector(std::size_t(GF::Traits::dimRange)))
    { return std::make_shared<VTKGridFunctionAdapter<GF> >(stackobject_to_shared_ptr(gf), name, remap); }

    //! construct a VTKGridFunctionAdapter
    /**
     * \param gf    Shared pointer to a GridFunction object to wrap.
     * \param name  Name of the field as returned by name().

     * \param remap How components are remapped between PDELab and VTK. The
     *              default value yields the identity map with entries
     *              (0,1,...,T::Traits::dimRange-1).
     *
     * The resulting VTKFunction will have \c remap.size() components.  None
     * of the elements of \c remap should be greater or equal to \c
     * T::Traits::dimRange.  When component \c c is requested by the
     * VTKWriter, it will be mapped to the component \c remap[c] of the
     * GridFunction.
     */
    template<class GF>
    std::shared_ptr<VTKGridFunctionAdapter<GF> > makeVTKGridFunctionAdapter
    ( const std::shared_ptr<const GF> &gf, const std::string &name,
      const std::vector<std::size_t> &remap =
        rangeVector(std::size_t(GF::Traits::dimRange)))
    { return std::make_shared<VTKGridFunctionAdapter<GF> >(gf, name, remap); }

    /** vizualize order of a hp-FiniteElementMap so it can be used with the VTKWriter from dune-grid.
        @tparam GV GridView to vizualize on
        @tparam FEM FiniteElementMap to vizualize
     */
    template<typename G, typename FEM>
    class VTKFiniteElementMapAdapter
      : public Dune::VTKFunction<G>
    {
      typedef typename G::ctype DF;
      enum {n=G::dimension};
      typedef typename G::template Codim<0>::Entity Entity;

    public:
      VTKFiniteElementMapAdapter (const FEM& fem_, std::string s_)
        : fem(fem_), s(s_)
      {}

      virtual int ncomps () const
      {
        return 1;
      }

      virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,n>& xi) const
      {
        return fem.getOrder(e);
      }

      virtual std::string name () const
      {
        return s;
      }

    private:
      const FEM& fem;
      std::string s;
    };

  }
}

#endif // DUNE_PDELAB_COMMON_VTKEXPORT_HH
