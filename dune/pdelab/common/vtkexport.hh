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
	  : public Dune::VTKWriter<typename T::Traits::GridViewType>::VTKFunction
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

  }
}

#endif
