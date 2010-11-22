// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH
#define DUNE_PDELAB_GRIDOPERATORSPACEUTILITIES_HH

#include <vector>

#include <dune/common/tuples.hh>

#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridoperatorspace/localmatrix.hh>

namespace Dune {
  namespace PDELab {

	//! collect types exported by a grid operator space
	template<typename GFSU, typename GFSV, typename B, 
			 typename CU, typename CV>
	struct GridOperatorSpaceTraits
	{
	  typedef GFSU TrialGridFunctionSpace;

	  typedef CU TrialConstraintsType;

	  typedef GFSV TestGridFunctionSpace;

	  typedef CV TestConstraintsType;

	  //! \brief the grid view where grid function is defined upon
	  typedef typename GFSU::Traits::GridViewType GridViewType;

	  //! \brief vector backend
	  typedef B BackendType;

	  //! \brief short cut for size type exported by Backend
	  typedef typename B::size_type SizeType;
	};


    class EmptyTransformation : public ConstraintsTransformation<int,float>
    {
    };


    class NoSubTriangulationImp
    {
    public:
      template<class E,class EG>
      void create_geometries(const E &, const EG &, const bool d = false)
      {DUNE_THROW(Dune::NotImplemented,"This should never be called.");}
      template<class E,class EG>
      void create_edges(const E &, const EG &)
      {DUNE_THROW(Dune::NotImplemented,"This should never be called.");}
      template<class E,class EG>
      void create_boundaries(const E &, const EG &)
      {DUNE_THROW(Dune::NotImplemented,"This should never be called.");}
    };
    
    template<typename GV>
    class NoSubTriangulation
    {
    public:
      static const bool hasSubTriangulation = false;
	  typedef typename GV::Traits::template Codim<0>::Entity Entity;
	  typedef typename GV::IntersectionIterator IntersectionIterator;
	  typedef typename IntersectionIterator::Intersection Intersection;

	  typedef ElementGeometry<Entity> SubEntity;
      typedef std::list<SubEntity> SubEntityList;
      typedef typename SubEntityList::iterator SubEntityIterator;
	  
      typedef IntersectionGeometry<Intersection>  SubIntersection;
      typedef std::list<SubIntersection> SubIntersectionList;

      struct BindSubEntity{
        template<typename LFS, typename SE>
        inline static void rebind(const LFS &, const SE &){}
      };

      struct BindInsideSubIntersection{
        template<typename LFS, typename SE>
        inline static void rebind(const LFS &, const SE &){}
      };

      struct BindOutsideSubIntersection{
        template<typename LFS, typename SE>
        inline static void rebind(const LFS &, const SE &){}
      };

      struct BindSubIntersection{
        template<typename LFS, typename SE>
        inline static void rebind(const LFS &, const SE &){}
      };

      // The iterator has to wrap the true intersection
      // iterator. Otherwise a list of entity pointers would have to
      // be stored.
      class SubIntersectionIterator
      {
      public:
        SubIntersectionIterator(const IntersectionIterator & it_, const IntersectionIterator & eit_)
          : intersection_index(0), it(it_), eit(eit_),
            sub_intersection(new SubIntersection(*it,intersection_index))
        {}

        SubIntersectionIterator(const SubIntersectionIterator & sit_)
          : intersection_index(0),
            it(sit_.it), eit(sit_.eit), sub_intersection(sit_.sub_intersection)
        {}

        SubIntersectionIterator & operator++()
        { 
          ++it;  ++intersection_index; 
          if(it != eit)
            sub_intersection.reset(new SubIntersection(*it,intersection_index));
          return *this; 
        }

        SubIntersection & operator*()
        {
          return *sub_intersection;
        }

        SubIntersection* operator->()
        {
          return sub_intersection.get();
        }
        
        bool operator==(const SubIntersectionIterator & at) const
        {
          return it == at.it;
        }

        bool operator!=(const SubIntersectionIterator & at) const
        {
          return it != at.it;
        }

      private: 
        int intersection_index;
        IntersectionIterator it;
        IntersectionIterator eit;
        mutable std::auto_ptr<SubIntersection> sub_intersection;
      };
      
      NoSubTriangulation(const GV & gv_, const NoSubTriangulationImp &) 
        : gv(gv_)
      {}

      NoSubTriangulation(const NoSubTriangulation & c_) 
        : gv(c_.gv), sub_entities(c_.sub_entities)
      {}

      void create(const Entity & e) const
      {
        sub_entities.clear();
        sub_entities.push_back(SubEntity(e));
      }

      SubEntityIterator begin() const
      {
        return sub_entities.begin();
      }

      SubEntityIterator end() const
      {
        return sub_entities.end();
      }

      SubIntersectionIterator ibegin() const
      {
        const Entity & e = sub_entities.front().entity();
        return SubIntersectionIterator(gv.ibegin(e),gv.iend(e));
      }

      SubIntersectionIterator iend() const
      {
        const Entity & e = sub_entities.front().entity();
        return SubIntersectionIterator(gv.iend(e),gv.iend(e));
      }

    private:
      const GV & gv;
      mutable SubEntityList sub_entities;
    };

	/**@ingroup FlatOperatorSpaceGroup
	   \brief Entry in sparsity pattern

	   The sparsity pattern of a linear operator is described by by connecting
	   degrees of freedom in one element with degrees of freedom in the
	   same element (intra) or an intersecting element (inter). 

	   This numbering is with respect to the depth-first canonical order of the 
	   degrees of freedom of an entity.

	   \nosubgrouping
	*/
	class SparsityLink : public Dune::tuple<int,int>
	{
	public:
	  //! \brief Standard constructor for uninitialized local index
	  SparsityLink ()
	  {}

	  //! \brief Initialize all components
	  SparsityLink (int i, int j)
		: Dune::tuple<int,int>(i,j)
	  {}

	  //! \brief Return first component
	  inline int i () const
	  {
		return Dune::get<0>(*this);
	  } 

	  //! \brief Return second component
	  inline int j () const
	  {
		return Dune::get<1>(*this);
	  } 

	  //! \brief Set both components
	  void set (int i, int j)
	  {
		Dune::get<0>(*this) = i;
		Dune::get<1>(*this) = j;
	  } 
	};

	/**@ingroup FlatOperatorSpaceGroup
	   \brief Layout description for a sparse linear operator
       \see SparsityLink
       
	   \nosubgrouping
	*/
	class LocalSparsityPattern : public std::vector<SparsityLink>
	{};

	//================================================
	// Default matrix backend
	//================================================

	// Simple Backend for std::vector
	class StdVectorFlatMatrixBackend
	{
	public:
	  // Matrix construction
	  template<typename T, typename E>
	  class Matrix : public std::vector<E>
	  {
		typedef std::vector<E> BaseT;
	  public:
		friend class StdVectorFlatMatrixBackend; // for access to line length

		typedef E ElementType;
        typedef StdVectorFlatMatrixBackend Backend;

		Matrix (const T& t) 
		  : n(t.globalSizeU()), 
			BaseT(t.globalSizeU()*t.globalSizeV()) 
		{}

	  private:
		std::vector<int>::size_type line () const
		{
		  return n;
		}

		std::vector<int>::size_type n;
	  };

	  // extract type of matrix element from a Matrix
	  template<class C>
	  struct Value
	  {
		typedef typename C::value_type Type;
	  };

	  //! The size type
	  typedef std::vector<int>::size_type size_type;

	  // get const_reference to container element
	  template<typename C>
	  static const typename C::value_type& access (const C& c, size_type i, size_type j)
	  {
		return c.operator[](i*c.line()+j);
	  }

	  // type to store sparsity pattern
	  class Pattern
	  {
	  public:
		void add_link (size_type i, size_type j)
		{
		}
	  };

	  // get non const_reference to container element 
	  // note: this method does not depend on T!
	  template<typename C>
	  static typename C::value_type& access (C& c, size_type i, size_type j)
	  {
		return c.operator[](i*c.line()+j);
	  }

// 	  // read a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void read (const C& c, 
// 						const RI& row_index, const CI& col_index, 
// 						LocalMatrix<T>& submatrix)
// 	  {
// 		submatrix.resize(row_index.size(),col_index.size());
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			submatrix(i,j) = c.operator[](row_index[i]*c.line()+col_index[j]);
// 	  }

// 	  // write a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void write (const RI& row_index, const CI& col_index, 
// 						 const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			c.operator[](row_index[i]*c.line()+col_index[j]) = submatrix(i,j);
// 	  }

// 	  // write a submatrix given by global indices
// 	  template<typename C, typename RI, typename CI, typename T>
// 	  static void add (const RI& row_index, const CI& col_index, 
// 					   const LocalMatrix<T>& submatrix, C& c)
// 	  {
// 		for (int j=0; j<col_index.size(); j++)
// 		  for (int i=0; i<row_index.size(); i++)
// 			c.operator[](row_index[i]*c.line()+col_index[j]) += submatrix(i,j);
// 	  }

	  // clear one row of the matrix
	  template<typename C, typename RI>
	  static void clear_row (RI row, C& c)
	  {
		for (int j=0; j<c.line(); j++)
		  c.operator[](row*c.line()+j) = 0;
	  }
	};

    template<typename GV>
    class MultiGeomUniqueIDMapper
    {
    private:
      static const int chunk = 1<<28;
      int offset;
      const typename GV::IndexSet & is;
      std::map<Dune::GeometryType,int> gtoffset;
      typedef typename GV::template Codim<0>::Entity Entity;
      typedef typename GV::IndexSet::IndexType SizeType;
      
    public:
      MultiGeomUniqueIDMapper(const GV & gv) 
        : offset(0), is(gv.indexSet())
      {}

      SizeType map(const Entity & e)
      {
        const Dune::GeometryType & gtype = e.type();

        // assign offset for geometry type;
        if (gtoffset.find(gtype)==gtoffset.end())
          {
            gtoffset[gtype] = offset;
            offset += chunk;
          }

        // compute unique id
        const SizeType id = is.index(e) + gtoffset[gtype];
        return id;
      }
    };

	// compile time switching of function call
    template<typename LA, bool doIt>
    struct LocalAssemblerCallSwitch
    {
      template<typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s, 
                                  const LFSU& lfsu_n, const LFSV& lfsv_n, 
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns)
      {
      }
      template<typename LFSU, typename LFSV>
      static void pattern_boundary(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   LocalSparsityPattern& pattern_ss)
      {
      }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_boundary (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s)
      {
      }

	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
                                  const LFSV& lfsv_s, const LFSV& lfsv_n,
                                  R& r_s, R& r_n)
      {
      }
 	  template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           Y& y_s, Y& y_n)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_boundary (const LA& la, const IG& ig, 
                                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                           Y& y_s)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_skeleton (const LA& la, const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn, 
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
      {
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_boundary (const LA& la, const IG& ig, 
                                     const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                     LocalMatrix<R>& mat_ss)
      {
      }
    };
    template<typename LA>
    struct LocalAssemblerCallSwitch<LA,true>
    {
      template<typename LFSU, typename LFSV>
      static void pattern_volume (const LA& la, const LFSU& lfsu, const LFSV& lfsv, LocalSparsityPattern& pattern)
      {
        la.pattern_volume(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_volume_post_skeleton
      ( const LA& la,
        const LFSU& lfsu, const LFSV& lfsv,
        LocalSparsityPattern& pattern)
      {
        la.pattern_volume_post_skeleton(lfsu,lfsv,pattern);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_skeleton (const LA& la, const LFSU& lfsu_s, const LFSV& lfsv_s, 
                                  const LFSU& lfsu_n, const LFSV& lfsv_n, 
                                   LocalSparsityPattern& pattern_sn,
                                   LocalSparsityPattern& pattern_ns)
      {
        la.pattern_skeleton(lfsu_s,lfsv_s,lfsu_n,lfsv_n,
                            pattern_sn, pattern_ns);
      }
      template<typename LFSU, typename LFSV>
      static void pattern_boundary(const LA& la,
                                   const LFSU& lfsu_s, const LFSV& lfsv_s,
                                   LocalSparsityPattern& pattern_ss)
      {
        la.pattern_boundary(lfsu_s,lfsv_s,pattern_ss);
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume(eg,lfsu,x,lfsv,r);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r)
	  {
		la.alpha_volume_post_skeleton(eg,lfsu,x,lfsv,r);
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           R& r_s, R& r_n)
      {
        la.alpha_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,r_s,r_n);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void alpha_boundary (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           R& r_s)
      {
        la.alpha_boundary(ig,lfsu_s,x_s,lfsv_s,r_s);
      }

	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume(eg,lfsv,r);
      }
	  template<typename EG, typename LFSV, typename R>
      static void lambda_volume_post_skeleton (const LA& la, const EG& eg, const LFSV& lfsv, R& r)
      {
        la.lambda_volume_post_skeleton(eg,lfsv,r);
      }
      template<typename IG, typename LFSV, typename R>
      static void lambda_skeleton(const LA& la, const IG& ig,
                                  const LFSV& lfsv_s, const LFSV& lfsv_n,
                                  R& r_s, R& r_n)
      {
        la.lambda_skeleton(ig, lfsv_s, lfsv_n, r_s, r_n);
      }
 	  template<typename IG, typename LFSV, typename R>
      static void lambda_boundary (const LA& la, const IG& ig, const LFSV& lfsv, R& r)
      {
        la.lambda_boundary(ig,lfsv,r);
      }

	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_apply_volume(eg,lfsu,x,lfsv,y);
	  }
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, Y& y)
	  {
		la.jacobian_apply_volume_post_skeleton(eg,lfsu,x,lfsv,y);
	  }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_skeleton (const LA& la, const IG& ig, 
                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                           const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                           Y& y_s, Y& y_n)
      {
        la.jacobian_apply_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,y_s,y_n);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
	  static void jacobian_apply_boundary (const LA& la, const IG& ig, 
                                           const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                           Y& y_s)
      {
        la.jacobian_apply_boundary(ig,lfsu_s,x_s,lfsv_s,y_s);
      }

      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
        la.jacobian_volume(eg,lfsu,x,lfsv,mat);
      }
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_volume_post_skeleton (const LA& la, const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, LocalMatrix<R>& mat)
      {
        la.jacobian_volume_post_skeleton(eg,lfsu,x,lfsv,mat);
      }
 	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_skeleton (const LA& la, const IG& ig, 
                              const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                              const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n, 
                              LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn, 
                              LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn)
      {
        la.jacobian_skeleton(ig,lfsu_s,x_s,lfsv_s,lfsu_n,x_n,lfsv_n,
                             mat_ss, mat_sn, mat_ns, mat_nn);
      }
	  template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	  static void jacobian_boundary (const LA& la, const IG& ig, 
                                     const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                     LocalMatrix<R>& mat_ss)
      {
        la.jacobian_boundary(ig,lfsu_s,x_s,lfsv_s,mat_ss);
      }
   };

    /**
       \brief Base class for grid operators.

       This class provides some generic behavior required for most
       grid operators. This includes the access of the global vectors
       and matrices via local indices and local function spaces with
       regard to the constraint mappings.
       
       \tparam GFSU GridFunctionSpace for ansatz functions
       \tparam GFSV GridFunctionSpace for test functions
       \tparam CU   Constraints maps for the individual dofs (trial space)
       \tparam CV   Constraints maps for the individual dofs (test space)
       \tparam B The vector backend used for the coefficient vector
       
     */
	template<typename GFSU, typename GFSV,
			 typename CU=EmptyTransformation,
			 typename CV=EmptyTransformation,
			 typename B=StdVectorFlatMatrixBackend>
    class GridOperatorBase{
    public:

	  typedef GridOperatorSpaceTraits<GFSU,GFSV,B,CU,CV> Traits;
      
      //! construct GridOperatorSpace
	  GridOperatorBase (const GFSU& gfsu_, const GFSV& gfsv_) 
		: gfsu(gfsu_), gfsv(gfsv_),
          pconstraintsu(&emptyconstraintsu), pconstraintsv(&emptyconstraintsv),
          lfsu(gfsu), lfsv(gfsv), lfsun(gfsu), lfsvn(gfsv)
	  {}

      //! construct GridOperatorSpace, with constraints
	  GridOperatorBase (const GFSU& gfsu_, const CU& cu,
						 const GFSV& gfsv_, const CV& cv) 
		: gfsu(gfsu_), gfsv(gfsv_),
          pconstraintsu(&cu), pconstraintsv(&cv),
          lfsu(gfsu), lfsv(gfsv), lfsun(gfsu), lfsvn(gfsv)
	  {}

      //! get dimension of space u
	  typename GFSU::Traits::SizeType globalSizeU () const
	  {
		return gfsu.globalSize();
	  }

      //! get dimension of space v
	  typename GFSV::Traits::SizeType globalSizeV () const
	  {
		return gfsv.globalSize();
	  }

      //! get the trial grid function space
      const GFSU& trialGridFunctionSpace() const
      {
        return gfsu;
      }

      //! get the test grid function space
      const GFSV& testGridFunctionSpace() const
      {
        return gfsv;
      }

      //! get the constraints on the trial grid function space
      const CU& trialConstraints() const
      {
        return *pconstraintsu;
      }

      //! get the constraints on the test grid function space
      const CV& testConstraints() const
      {
        return *pconstraintsv;
      }


      /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
          V\f$ to \f$ V'\f$. If postrestrict == true then
          \f$\boldsymbol{R}^T_{\boldsymbol{\tilde U}', \boldsymbol{U}'}
          \boldsymbol{S}_{\boldsymbol{\tilde V}}\f$ is applied
           instead of the full transformation.  */
      template<typename X>
      void forwardtransform(X & x, const bool postrestrict = false)
      {
        typedef typename CV::const_iterator global_col_iterator;	  
        for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
          typedef typename global_col_iterator::value_type::first_type GlobalIndex;
          const GlobalIndex & contributor = cit->first;

          typedef typename global_col_iterator::value_type::second_type ContributedMap;
          typedef typename ContributedMap::const_iterator global_row_iterator;
          const ContributedMap & contributed = cit->second;
          global_row_iterator it  = contributed.begin();
          global_row_iterator eit = contributed.end();
          
          for(;it!=eit;++it)
            x[it->first] += it->second * x[contributor];
        }

        if(postrestrict)
          for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit)
            x[cit->first]=0.;
      }

      /** \brief Transforms a vector \f$ \boldsymbol{x} \f$ from \f$
          V'\f$ to \f$ V\f$. If prerestrict == true then
          \f$\boldsymbol{S}^T_{\boldsymbol{\tilde U}}\f$ is applied
           instead of the full transformation.  */
      template<typename X>
      void backtransform(X & x, const bool prerestrict = false)
      {
        typedef typename CV::const_iterator global_col_iterator;  
        for (global_col_iterator cit=pconstraintsv->begin(); cit!=pconstraintsv->end(); ++cit){
          typedef typename global_col_iterator::value_type::first_type GlobalIndex;
          const GlobalIndex & contributor = cit->first;

          typedef typename global_col_iterator::value_type::second_type ContributedMap;
          typedef typename ContributedMap::const_iterator global_row_iterator;
          const ContributedMap & contributed = cit->second;
          global_row_iterator it  = contributed.begin();
          global_row_iterator eit = contributed.end();
          
          if(prerestrict)
            x[contributor] = 0.;

          for(;it!=eit;++it)
            x[contributor] += it->second * x[it->first];
        }
      }

    protected:

      /** \brief read local stiffness matrix for entity */  
      template<typename LFSV, typename LFSU, typename GC, typename T>
      void eread (const LFSV& lfsv, const LFSU& lfsu, const GC& globalcontainer, 
                  LocalMatrix<T>& localcontainer) const
      {
        for (int i=0; i<lfsv.size(); i++)
          for (int j=0; j<lfsu.size(); j++)
            localcontainer(i,j) = B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j));
      }

      /** \brief write local stiffness matrix for entity */  
      template<typename LFSV, typename LFSU, typename T, typename GC>
      void ewrite (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
      {
        for (int i=0; i<lfsv.size(); i++)
          for (int j=0; j<lfsu.size(); j++)
            B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j)) = localcontainer(i,j);
      }

      /** \brief write local stiffness matrix for entity */  
      template<typename LFSV, typename LFSU, typename T, typename GC>
      void eadd (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
      {
        for (size_t i=0; i<lfsv.size(); i++)
          for (size_t j=0; j<lfsu.size(); j++)
            B::access(globalcontainer,lfsv.globalIndex(i),lfsu.globalIndex(j)) += localcontainer(i,j);
      }

      /** \brief Add local matrix \f$m\f$ to global Jacobian \f$J\f$
          and apply constraints transformation. Hence we perform: \f$
          \boldsymbol{J} := \boldsymbol{J} + \boldsymbol{S}_{
          \boldsymbol{\tilde V}} m \boldsymbol{S}^T_{
          \boldsymbol{\tilde U}} \f$*/  
      template<typename LFSV, typename LFSU, typename T, typename GC>
      void etadd (const LFSV& lfsv, const LFSU& lfsu, const LocalMatrix<T>& localcontainer, GC& globalcontainer) const
      {

        for (size_t i=0; i<lfsv.size(); i++)
          for (size_t j=0; j<lfsu.size(); j++){
            typename Traits::SizeType gi = lfsv.globalIndex(i);
            typename Traits::SizeType gj = lfsu.globalIndex(j);
            
            // Get global constraints containers for test and ansatz space
            const CV & cv = *pconstraintsv;
            const CU & cu = *pconstraintsu;

            typedef typename CV::const_iterator global_vcol_iterator;
            typedef typename global_vcol_iterator::value_type::second_type global_vrow_type;
            typedef typename global_vrow_type::const_iterator global_vrow_iterator;

            typedef typename CU::const_iterator global_ucol_iterator;
            typedef typename global_ucol_iterator::value_type::second_type global_urow_type;
            typedef typename global_urow_type::const_iterator global_urow_iterator;

            // Check whether the global indices are constrained indices
            global_vcol_iterator gvcit = cv.find(gi);
            global_ucol_iterator gucit = cu.find(gj);

            // Set constrained_v true if gi is constrained dof
            bool constrained_v(false);
            global_vrow_iterator gvrit;
            if(gvcit!=cv.end()){
              gvrit = gvcit->second.begin();              
              constrained_v = true;
            }

            T vf = 1;
            do{
              // if gi is index of constrained dof
              if(constrained_v){

                if(gvrit == gvcit->second.end())
                  break;

                // otherwise set gi to an index to a contributed dof
                // and set vf to the contribution weight
                gi = gvrit->first;
                vf = gvrit->second;
              }

            // Set constrained_u true if gj is constrained dof
              bool constrained_u(false);
              global_urow_iterator gurit;
              if(gucit!=cu.end()){
                gurit = gucit->second.begin();
                constrained_u = true;
                if(gurit == gucit->second.end()){
                  T t = localcontainer(i,j) * vf;
                  if(t != 0.0)                 // entry might not be present in the matrix
                    B::access(globalcontainer,gi,gj) += t;
                }
              }

              T uf = 1;
              do{
                // if gj is index of constrained dof
                if(constrained_u){

                  if(gurit == gucit->second.end())
                    break;

                  // otherwise set gj to an index to a contributed dof
                  // and set uf to the contribution weight
                  gj = gurit->first;
                  uf = gurit->second;
                }

                // add weighted local entry to global matrix
                T t = localcontainer(i,j) * uf * vf;
                if (t != 0.0)                 // entry might not be present in the matrix
                  B::access(globalcontainer,gi,gj) += t;

                if(constrained_u && gurit != gucit->second.end())
                  ++gurit;
                else 
                  break;

              }while(true);

              if(constrained_v && gvrit != gvcit->second.end())
                ++gvrit;
              else
                break;

            }while(true);

          }
      }

      /** \brief Adding matrix entry to pattern with respect to the
       constraints contributions. This assembles the entries addressed
       by etadd(..). See the documentation there for more information
       about the matrix pattern. */
      template<typename GI, typename P>
      void add_entry(P & globalpattern, GI gi, GI gj) const
      {
        const CV & cv = *pconstraintsv;
        const CU & cu = *pconstraintsu;

        typedef typename CV::const_iterator global_vcol_iterator;
        typedef typename global_vcol_iterator::value_type::second_type global_vrow_type;
        typedef typename global_vrow_type::const_iterator global_vrow_iterator;

        typedef typename CU::const_iterator global_ucol_iterator;
        typedef typename global_ucol_iterator::value_type::second_type global_urow_type;
        typedef typename global_urow_type::const_iterator global_urow_iterator;
            
        global_vcol_iterator gvcit = cv.find(gi);
        global_ucol_iterator gucit = cu.find(gj);

        if(gi==gj)
          globalpattern.add_link(gi,gj);

        bool constrained_v(false);
        global_vrow_iterator gvrit;
        if(gvcit!=cv.end()){
          gvrit = gvcit->second.begin();              
          constrained_v = true;
          if(gvrit == gvcit->second.end())
            globalpattern.add_link(gi,gj);
        }

        do{
          if(constrained_v){
            if(gvrit == gvcit->second.end())
              break;
            gi = gvrit->first;
          }

          bool constrained_u(false);
          global_urow_iterator gurit;
          if(gucit!=cu.end()){
            gurit = gucit->second.begin();
            constrained_u = true;
            if(gurit == gucit->second.end())
              globalpattern.add_link(gi,gj);
          }

          do{
            if(constrained_u){
              if(gurit == gucit->second.end())
                break;

              gj = gurit->first;
            }
                
            globalpattern.add_link(gi,gj);

            if(constrained_u && gurit != gucit->second.end())
              ++gurit;
            else 
              break;

          }while(true);

          if(constrained_v && gvrit != gvcit->second.end())
            ++gvrit;
          else
            break;

        }while(true);

      }

      /** \brief insert dirichlet constraints for row and assemble
          T^T_U in constrained rows
      */  
      template<typename GI, typename GC, typename CG>
      void set_trivial_row (GI i, const CG & cv_i, GC& globalcontainer) const
      {
        //std::cout << "clearing row " << i << std::endl;
        // set all entries in row i to zero
        B::clear_row(i,globalcontainer);

        // set diagonal element to 1
        B::access(globalcontainer,i,i) = 1;
      }

      /* global function spaces */
      const GFSU& gfsu;
      const GFSV& gfsv;
      /* constraints */
      const CU* pconstraintsu;
      const CV* pconstraintsv;
      static CU emptyconstraintsu;
      static CV emptyconstraintsv;
      /* local function spaces */
      typedef LocalFunctionSpace<GFSU, TrialSpaceTag> LFSU;
      typedef LocalFunctionSpace<GFSV, TestSpaceTag> LFSV;
      // local function spaces in local cell
      mutable LFSU lfsu;
      mutable LFSV lfsv;
      // local function spaces in neighbor
      mutable LFSU lfsun;
      mutable LFSV lfsvn;
    };

	template<typename GFSU, typename GFSV, typename CU, typename CV, typename B>
    CU GridOperatorBase<GFSU,GFSV,CU,CV,B>::emptyconstraintsu;
    template<typename GFSU, typename GFSV, typename CU, typename CV, typename B>
    CV GridOperatorBase<GFSU,GFSV,CU,CV,B>::emptyconstraintsv;
    

  } // namespace PDELab
} // namespace Dune

#endif
