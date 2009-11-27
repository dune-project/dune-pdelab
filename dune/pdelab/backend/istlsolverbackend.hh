#ifndef DUNE_ISTLSOLVERBACKEND_HH
#define DUNE_ISTLSOLVERBACKEND_HH

#include <dune/istl/solvercategory.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/scalarproducts.hh>
#include <dune/istl/io.hh>
#include<dune/istl/superlu.hh>

#include "../gridfunctionspace/constraints.hh"
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


	//========================================================
	// A parallel helper class providing a nonoverlapping 
	// decomposition of all degrees of freedom
	//========================================================

	// operator that resets result to zero at constrained DOFS
	template<typename GFS>
	class ParallelISTLHelper
	{
	  class GhostGatherScatter
	  {
	  public:
		template<class MessageBuffer, class EntityType, class DataType>
		void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
			data = (1<<24);
		  buff.write(data);
		}
  
		template<class MessageBuffer, class EntityType, class DataType>
		void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  DataType x; 
		  buff.read(x);
		  if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
			data = (1<<24);
		}
	  };

	  class InteriorBorderGatherScatter
	  {
	  public:
		template<class MessageBuffer, class EntityType, class DataType>
		void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
			data = (1<<24);
		  buff.write(data);
		}
  
		template<class MessageBuffer, class EntityType, class DataType>
		void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
		{
		  DataType x; 
		  buff.read(x);
		  if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
			data = x;
		  else
			data = std::min(data,x);
		}
	  };

	  typedef typename GFS::template VectorContainer<double>::Type V;

	public:

	  ParallelISTLHelper (const GFS& gfs_)
		: gfs(gfs_), v(gfs,(double)gfs.gridview().comm().rank())
	  {
		// find out about ghosts
		Dune::PDELab::GenericDataHandle2<GFS,V,GhostGatherScatter> gdh(gfs,v,GhostGatherScatter());
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(gdh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

		// partition interior/border
		Dune::PDELab::GenericDataHandle2<GFS,V,InteriorBorderGatherScatter> dh(gfs,v,InteriorBorderGatherScatter());
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(dh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);

		// convert vector into mask vector
		for (typename V::size_type i=0; i<v.N(); ++i)
		  for (typename V::size_type j=0; j<v[i].N(); ++j)
			if (v[i][j]==gfs.gridview().comm().rank())
			  v[i][j] = 1.0;
			else
			  v[i][j] = 0.0;
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


	//========================================================
	// Generic support for nonoverlapping grids
	//========================================================

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

	  NonoverlappingOperator (const GFS& gfs_, const M& A, const ParallelISTLHelper<GFS>& helper_) 
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
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }
  
	  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
	  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
	  {
		// apply local operator; now we have sum y_p = sequential y
		_A_.usmv(alpha,x,y);

		// accumulate y on border
		Dune::PDELab::AddDataHandle<GFS,Y> adddh(gfs,y);
		if (gfs.gridview().comm().size()>1)
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
	  const ParallelISTLHelper<GFS>& helper;
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
	  NonoverlappingScalarProduct (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
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

	  /*! \brief make additive vector consistent
	  */
	  void make_consistent (X& x) const
	  {
 		Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,x);
 		if (gfs.gridview().comm().size()>1)
 		  gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }

	private:
	  const GFS& gfs;
	  const ParallelISTLHelper<GFS>& helper;
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
	  NonoverlappingRichardson (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
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
		// no communication is necessary here because defect is already consistent!
// 		helper.mask(v);
// 		Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
// 		if (gfs.gridview().comm().size()>1)
// 		  gfs.gridview().communicate(adddh,Dune::InteriorBorder_InteriorBorder_Interface,Dune::ForwardCommunication);
	  }

	  /*!
		\brief Clean up.
      
		\copydoc Preconditioner::post(X&)
	  */
	  virtual void post (X& x) {}

	private:
	  const GFS& gfs;
	  const ParallelISTLHelper<GFS>& helper;
	};


	//========================================================
	// Generic support for overlapping grids
	// (need to be used with appropriate constraints)
	//========================================================

	// operator that resets result to zero at constrained DOFS
	template<class CC, class M, class X, class Y>
	class OverlappingOperator : public Dune::AssembledLinearOperator<M,X,Y>
	{
	public:
	  //! export types
	  typedef M matrix_type;
	  typedef X domain_type;
	  typedef Y range_type;
	  typedef typename X::field_type field_type;

	  //redefine the category, that is the only difference
	  enum {category=Dune::SolverCategory::overlapping};

	  OverlappingOperator (const CC& cc_, const M& A) 
		: cc(cc_), _A_(A) 
	  {}

	  //! apply operator to x:  \f$ y = A(x) \f$
	  virtual void apply (const X& x, Y& y) const
	  {
		_A_.mv(x,y);
		Dune::PDELab::set_constrained_dofs(cc,0.0,y);
	  }
  
	  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
	  virtual void applyscaleadd (field_type alpha, const X& x, Y& y) const
	  {
		_A_.usmv(alpha,x,y);
		Dune::PDELab::set_constrained_dofs(cc,0.0,y);
	  }
  
	  //! get matrix via *
	  virtual const M& getmat () const
	  {
		return _A_;
	  }
  
	private:
	  const CC& cc;
	  const M& _A_;
	};

	// new scalar product assuming at least overlap 1
	// uses unique partitioning of nodes for parallelization
	template<class GFS, class X>
	class OverlappingScalarProduct : public Dune::ScalarProduct<X>
	{
	public:
	  //! export types
	  typedef X domain_type;
	  typedef typename X::ElementType field_type;

	  //! define the category
	  enum {category=Dune::SolverCategory::overlapping};

	  /*! \brief Constructor needs to know the grid function space
	   */
	  OverlappingScalarProduct (const GFS& gfs_, const ParallelISTLHelper<GFS>& helper_)
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
	  const ParallelISTLHelper<GFS>& helper;
	};

	// wrapped sequential preconditioner
	template<class CC, class GFS, class P>
	class OverlappingWrappedPreconditioner 
	  : public Dune::Preconditioner<typename P::domain_type,typename P::range_type> 
	{
	public:
	  //! \brief The domain type of the preconditioner.
	  typedef typename P::domain_type domain_type;
	  //! \brief The range type of the preconditioner.
	  typedef typename P::range_type range_type;

	  // define the category
	  enum {
		//! \brief The category the preconditioner is part of.
		category=Dune::SolverCategory::overlapping
	  };

	  /*! \brief Constructor.
      
		Constructor gets all parameters to operate the prec.
		\param A The matrix to operate on.
		\param n The number of iterations to perform.
		\param w The relaxation factor.
	  */
	  OverlappingWrappedPreconditioner (const GFS& gfs_, P& prec_, const CC& cc_, 
										const ParallelISTLHelper<GFS>& helper_)
		: gfs(gfs_), prec(prec_), cc(cc_), helper(helper_)
	  {}

	  /*!
		\brief Prepare the preconditioner.
      
		\copydoc Preconditioner::pre(domain_type&,range_type&)
	  */
	  virtual void pre (domain_type& x, range_type& b) 
	  {
		prec.pre(x,b);
	  }

	  /*!
		\brief Apply the precondioner.
      
		\copydoc Preconditioner::apply(domain_type&,const range_type&)
	  */
	  virtual void apply (domain_type& v, const range_type& d)
	  {
		range_type dd(d);
		set_constrained_dofs(cc,0.0,dd);
		prec.apply(v,dd);
		Dune::PDELab::AddDataHandle<GFS,domain_type> adddh(gfs,v);
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
	  }

	  /*!
		\brief Clean up.
      
		\copydoc Preconditioner::post(domain_type&)
	  */
	  virtual void post (domain_type& x) 
	  {
		prec.post(x);
	  }

	private:
	  const GFS& gfs;
	  P& prec;
	  const CC& cc;
	  const ParallelISTLHelper<GFS>& helper;
	};


#ifdef HAVE_SUPERLU
	// exact subdomain solves with SuperLU as preconditioner
	template<class GFS, class M, class X, class Y>
	class SuperLUSubdomainSolver : public Dune::Preconditioner<X,Y> 
	{
	  typedef typename M::BaseT ISTLM;

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
		category=Dune::SolverCategory::overlapping
	  };

	  /*! \brief Constructor.
      
		Constructor gets all parameters to operate the prec.
		\param A The matrix to operate on.
		\param n The number of iterations to perform.
		\param w The relaxation factor.
	  */
	  SuperLUSubdomainSolver (const GFS& gfs_, const M& A_)
		: gfs(gfs_), A(A_), solver(A_,false) // this does the decomposition
	  {}

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
		Dune::InverseOperatorResult stat;
		Y b(d); // need copy, since solver overwrites right hand side
		solver.apply(v,b,stat);
		Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
	  }

	  /*!
		\brief Clean up.
      
		\copydoc Preconditioner::post(X&)
	  */
	  virtual void post (X& x) {}

	private:
	  const GFS& gfs;
	  const M& A;
	  Dune::SuperLU<ISTLM> solver;
	};

	// exact subdomain solves with SuperLU as preconditioner
	template<class GFS, class M, class X, class Y>
	class RestrictedSuperLUSubdomainSolver : public Dune::Preconditioner<X,Y> 
	{
	  typedef typename M::BaseT ISTLM;

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
		category=Dune::SolverCategory::overlapping
	  };

	  /*! \brief Constructor.
      
		Constructor gets all parameters to operate the prec.
		\param A The matrix to operate on.
		\param n The number of iterations to perform.
		\param w The relaxation factor.
	  */
	  RestrictedSuperLUSubdomainSolver (const GFS& gfs_, const M& A_,  
										const ParallelISTLHelper<GFS>& helper_)
		: gfs(gfs_), A(A_), solver(A_,false), helper(helper_) // this does the decomposition
	  {}

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
		Dune::InverseOperatorResult stat;
		Y b(d); // need copy, since solver overwrites right hand side
		solver.apply(v,b,stat);
		helper.mask(v);
		Dune::PDELab::AddDataHandle<GFS,X> adddh(gfs,v);
		if (gfs.gridview().comm().size()>1)
		  gfs.gridview().communicate(adddh,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);
	  }

	  /*!
		\brief Clean up.
      
		\copydoc Preconditioner::post(X&)
	  */
	  virtual void post (X& x) {}

	private:
	  const GFS& gfs;
	  const M& A;
	  Dune::SuperLU<ISTLM> solver;
	  const ParallelISTLHelper<GFS>& helper;
	};

#endif



  } // namespace PDELab
} // namespace Dune

#endif
