// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_VTKEXPORT_HH
#define DUNE_PDELAB_VTKEXPORT_HH

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
		: t(t_), s(s_)
	  {}

	  virtual int ncomps () const
	  {
		return T::Traits::dimRange;
	  }

	  virtual double evaluate (int comp, const Entity& e, const Dune::FieldVector<DF,n>& xi) const
	  {
		typename T::Traits::DomainType x;
		typename T::Traits::RangeType y;

		for (int i=0; i<n; i++) 
		  x[i] = xi[i];
		t.evaluate(e,x,y);
		return y[comp];
	  }
	  
	  virtual std::string name () const
	  {
		return s;
	  }

	private:
	  const T& t;
	  std::string s;
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
