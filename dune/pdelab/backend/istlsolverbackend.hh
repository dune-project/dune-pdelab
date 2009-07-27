#ifndef DUNE_ISTLSOLVERBACKEND_HH
#define DUNE_ISTLSOLVERBACKEND_HH

#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/io.hh>

#include "../gridfunctionspace/genericdatahandle.hh"

namespace Dune {
  namespace PDELab {

	template<typename X, typename Y, typename GOS>
	class OnTheFlyOperator : public Dune::LinearOperator<X,Y>
	{
	public:
	  typedef X domain_type;
	  typedef Y range_type;
	  typedef typename X::field_type field_type;

	  enum {category=Dune::SolverCategory::sequential};

	  OnTheFlyOperator (GOS& gos_)
		: gos(gos_)
	  {}

	  virtual void apply (const X& x, Y& y) const
	  {
		y = 0.0;
		gos.jacobian_apply(x,y);
	  }

	  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
	  {
		Y temp(y);
		temp = 0.0;
		gos.jacobian_apply(x,temp);
		y.axpy(alpha,temp);
	  }

	private:
	  GOS& gos;
	};


	// operator that resets result to zero at constrained DOFS
	template<typename GFS>
	class NonoverlappingHelper
	{
	  class InteriorBorderGatherScatter
	  {
	  public:
		InteriorBorderGatherScatter (int rank_)
		  : rank(rank_)
		{
		}

		template<class MessageBuffer, class EntityType, class DataType>
		void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  buff.write(data);
		}
  
		template<class MessageBuffer, class EntityType, class DataType>
		void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  DataType x; buff.read(x);
		  if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
			data = -1.0;
		}
	  private:
		int rank;
	  };

	  typedef typename GFS::template VectorContainer<double>::Type V;
	  
	public:

	  NonoverlappingHelper (const GFS& gfs_)
		: gfs(gfs_), v(gfs_,(double)gfs.gridview().comm().rank())
	  {
		// fill interior/border DOFS with rank
		InteriorBorderGatherScatter gs(gfs.gridview().comm().rank());
		Dune::PDELab::GenericDataHandle2<GFS,V,InteriorBorderGatherScatter> dh(gfs,v,gs);
		gfs.gridview().communicate(dh,Dune::All_All_Interface,Dune::ForwardCommunication);

		// compute minimum rank for each interior/border DOF
		Dune::PDELab::MinDataHandle<GFS,V> mindh(gfs,v);
		gfs.gridview().communicate(mindh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

		// convert vector into mask vector
		for (typename V::size_type i=0; i<v.N(); ++i)
		  for (typename V::size_type j=0; j<v[i].N(); ++j)
			if (v[i][j]==gfs.gridview().comm().rank())
			  v[i][j] = 1.0;
			else
			  v[i][j] = 0.0;

// 		for (typename V::size_type i=0; i<v.N(); ++i)
// 		  for (typename V::size_type j=0; j<v[i].N(); ++j)
// 			std::cout << "/" << gfs.gridview().comm().rank() 
// 					  << "/: i=" << i
// 					  << " j=" << j
// 					  << " v=" << v[i][j]
// 					  << std::endl;
	  }

	  // keep only DOFs assigned to this processor
	  template<typename W>
	  void mask (W& w) const
	  {
		for (typename V::size_type i=0; i<v.N(); ++i)
		  for (typename V::size_type j=0; j<v[i].N(); ++j)
			  w[i][j] *= v[i][j];
	  }

	  // access to mask vector
	  double mask (typename V::size_type i, typename V::size_type j) const
	  {
		return v[i][j];
	  }

	private:
	  const GFS& gfs;
	  V v;
	};

	// operator that resets result to zero at nonowned DOFS
	template<class GFS, class M, class X, class Y>
	class NonoverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
	{
	public:
	  //! export types
	  typedef M matrix_type;
	  typedef X domain_type;
	  typedef Y range_type;
	  typedef typename X::field_type field_type;

	  //redefine the category, that is the only difference
	  enum {category=Dune::SolverCategory::nonoverlapping};

	  NonoverlappingOperator (const GFS& gfs_, const M& A, const NonoverlappingHelper<GFS>& helper_) 
		: gfs(gfs_), _A_(A), helper(helper_)
	  {
	  }

	  //! apply operator to x:  \f$ y = A(x) \f$
	  virtual void apply (const X& x, Y& y) const
	  {
		// apply local operator; now we have sum y_p = sequential y
		_A_.mv(x,y);

		// accumulate y on border
		Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
		gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }
  
	  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
	  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
	  {
		// apply local operator; now we have sum y_p = sequential y
		_A_.usmv(alpha,x,y);

		// accumulate y on border
		Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
		gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }
  
	  //! get matrix via *
	  virtual const M& getmat () const
	  {
		return _A_;
	  }
  
	private:
	  const GFS& gfs;
	  const M& _A_;
	  const NonoverlappingHelper<GFS>& helper;
	};


	// parallel scalar product assuming no overlap
	template<class GFS, class X>
	class NonoverlappingScalarProduct : public Dune::ScalarProduct<X>
	{
	public:
	  //! export types
	  typedef X domain_type;
	  typedef typename X::ElementType field_type;

	  //! define the category
	  enum {category=Dune::SolverCategory::nonoverlapping};

	  /*! \brief Constructor needs to know the grid function space
	   */
	  NonoverlappingScalarProduct (const GFS& gfs_, const NonoverlappingHelper<GFS>& helper_)
		: gfs(gfs_), helper(helper_)
	  {}

	  /*! \brief Dot product of two vectors. 
		It is assumed that the vectors are consistent on the interior+border
		partition.
	  */
	  virtual field_type dot (const X& x, const X& y)
	  {
		// do local scalar product on unique partition
		field_type sum = 0;
		for (typename X::size_type i=0; i<x.N(); ++i)
		  for (typename X::size_type j=0; j<x[i].N(); ++j)
			sum += (x[i][j]*y[i][j])*helper.mask(i,j);

		// do global communication
		return gfs.gridview().comm().sum(sum);
	  }

	  /*! \brief Norm of a right-hand side vector. 
		The vector must be consistent on the interior+border partition
	  */
	  virtual double norm (const X& x)
	  {
		return sqrt(static_cast<double>(this->dot(x,x)));
	  }

	private:
	  const GFS& gfs;
	  const NonoverlappingHelper<GFS>& helper;
	};

	// parallel Richardson preconditioner
	template<class GFS, class X, class Y>
	class NonoverlappingRichardson : public Dune::Preconditioner<X,Y> 
	{
	public:
	  //! \brief The domain type of the preconditioner.
	  typedef X domain_type;
	  //! \brief The range type of the preconditioner.
	  typedef Y range_type;
	  //! \brief The field type of the preconditioner.
	  typedef typename X::ElementType field_type;

	  // define the category
	  enum {
		//! \brief The category the preconditioner is part of.
		category=Dune::SolverCategory::nonoverlapping
	  };

	  /*! \brief Constructor.
      
		Constructor gets all parameters to operate the prec.
		\param A The matrix to operate on.
		\param n The number of iterations to perform.
		\param w The relaxation factor.
	  */
	  NonoverlappingRichardson (const GFS& gfs_, const NonoverlappingHelper<GFS>& helper_)
		: gfs(gfs_), helper(helper_)
	  {
	  }

	  /*!
		\brief Prepare the preconditioner.
      
		\copydoc Preconditioner::pre(X&,Y&)
	  */
	  virtual void pre (X& x, Y& b) {}

	  /*!
		\brief Apply the precondioner.
      
		\copydoc Preconditioner::apply(X&,const Y&)
	  */
	  virtual void apply (X& v, const Y& d)
	  {
		v = d;
		helper.mask(v);
		Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
		gfs.gridview().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	  }

	  /*!
		\brief Clean up.
      
		\copydoc Preconditioner::post(X&)
	  */
	  virtual void post (X& x) {}

	private:
	  const GFS& gfs;
	  const NonoverlappingHelper<GFS>& helper;
	};



  } // namespace PDELab
} // namespace Dune

#endif
