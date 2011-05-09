// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_VTKEXPORT_HH
#define DUNE_PDELAB_VTKEXPORT_HH

#include<cstddef>
#include<vector>

#include<dune/grid/io/file/vtk/vtkwriter.hh>

namespace Dune {
  namespace PDELab {

    //! wrap a GridFunction so it can be used with the VTKWriter from dune-grid.
	template<typename T> // T is a grid function
	class VTKGridFunctionAdapter
	  : public Dune::VTKFunction<typename T::Traits::GridViewType::Grid>
    {
	  typedef typename T::Traits::GridViewType::Grid::ctype DF;
	  enum {n=T::Traits::GridViewType::dimension};
	  typedef typename T::Traits::GridViewType::Grid::template Codim<0>::Entity Entity;

	public:
	  VTKGridFunctionAdapter (const T& t_, std::string s_)
        : t(t_), s(s_), remap(T::Traits::dimRange)
      {
        for(std::size_t c = 0; c < T::Traits::dimRange; ++c)
          remap[c] = c;
      }

      //! construct a VTKGridFunctionAdapter
      /**
       * \param t_     GridFunction object to wrap.
       * \param s_     Name of the field as returned by name().
       * \param remap_ How components are remapped between PDELab and VTK.
       *
       * The resulting VTKFunction will have \c remap_.size() components.
       * None of the elements of \c remap_ should be greater or equal to \c
       * T::Traits::dimRange.  When component \c c is requested by the
       * VTKWriter, it will be mapped to the component \c remap_[c] of the
       * GridFunction.
       */
      VTKGridFunctionAdapter (const T& t_, std::string s_,
                              const std::vector<std::size_t> &remap_)
        : t(t_), s(s_), remap(remap_)
	  {}

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
		t.evaluate(e,x,y);
        return y[remap[comp]];
	  }
	  
	  virtual std::string name () const
	  {
		return s;
	  }

	private:
	  const T& t;
	  std::string s;
      std::vector<std::size_t> remap;
	};

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

#endif
