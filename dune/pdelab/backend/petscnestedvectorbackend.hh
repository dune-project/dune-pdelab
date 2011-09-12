// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PETSCNESTEDVECTORBACKEND_HH
#define DUNE_PETSCNESTEDVECTORBACKEND_HH

#include<vector>

#include<dune/common/fvector.hh>

#include <dune/pdelab/backend/petscutility.hh>
#include <dune/pdelab/backend/petscvectorbackend.hh>

#include <petscvec.h>

#include "backendselector.hh"

namespace Dune {
  namespace PDELab {

    class PetscNestedVectorBackend;

    class PetscNestedVectorContainer
    {

      typedef std::vector<shared_ptr<PetscVectorContainer> > ChildContainer;

    public:
      typedef PetscScalar ElementType;
      typedef ElementType E;
      typedef Vec ContainerType;

      typedef E block_type;
      typedef block_type field_type;

      typedef std::size_t size_type;
      typedef PetscVectorBackend Backend;

      struct create_child_vectors
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::DynamicTraversal
      {
        create_child_vectors(ChildContainer& children)
          : _children(children)
        {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(T&& t, Child&& child, TreePath treePath, ChildIndex childIndex)
        {
          _children[childIndex] = make_shared<PetscVectorContainer>(child);
        }

        ChildContainer& _children;
      };

      template<typename GFS>
      PetscNestedVectorContainer (const GFS& gfs)
        : _data(NULL)
      {
        const size_t _blockCount = GFS::CHILDREN;
        ChildContainer children(_blockCount);
        create_child_vectors ccv(children);
        TypeTree::applyToTree(gfs,ccv);
        Vec child_array[_blockCount];
        for (int i = 0; i < _blockCount; ++i)
          child_array[i] = children[i]->_v;
        PETSC_CALL(VecCreateNest(PETSC_COMM_SELF,_blockCount,PETSC_NULL,child_array,&_v));
      }

      template<typename GFS>
      PetscNestedVectorContainer (const GFS& gfs, const E& e)
        : _data(NULL)
      {
        const size_t _blockCount = GFS::CHILDREN;
        ChildContainer children(_blockCount);
        create_child_vectors ccv(children);
        TypeTree::applyToTree(gfs,ccv);
        Vec child_array[_blockCount];
        for (int i = 0; i < _blockCount; ++i)
          child_array[i] = *(children[i]);
        PETSC_CALL(VecCreateNest(PETSC_COMM_SELF,_blockCount,PETSC_NULL,child_array,&_v));
        this->operator=(e);
      }

      PetscNestedVectorContainer (const PetscNestedVectorContainer& rhs)
        : _data(NULL)
      {
        PETSC_CALL(VecDuplicate(rhs._v,&_v));
        PETSC_CALL(VecCopy(rhs._v,_v));
      }

      PetscNestedVectorContainer& operator= (const PetscNestedVectorContainer& rhs)
      {
        checkin();
        PETSC_CALL(VecCopy(rhs._v,_v));
        _data = NULL;
      }

      ~PetscNestedVectorContainer()
      {
        checkin();
        PETSC_CALL(VecDestroy(&_v));
      }

      size_type N() const
      {
        int size = 0;
        PETSC_CALL(VecGetLocalSize(_v,&size));
        return size;
      }


      PetscNestedVectorContainer& operator= (const E& e)
      {
        checkin();
        PETSC_CALL(VecSet(_v,e));
        return *this;
      }

      PetscNestedVectorContainer& operator*= (const E& e)
      {
        checkin();
        PETSC_CALL(VecScale(_v,e));
        return *this;
      }


      PetscNestedVectorContainer& operator+= (const E& e)
      {
        checkin();
        PETSC_CALL(VecShift(_v,e));
        return *this;
      }

      PetscNestedVectorContainer& operator+= (const PetscNestedVectorContainer& rhs)
      {
        return axpy(1.0,rhs);
      }

      PetscNestedVectorContainer& operator-= (const PetscNestedVectorContainer& rhs)
      {
        return axpy(-1.0,rhs);
      }

      field_type& operator[](std::size_t i)
      {
        checkout();
        return _data[i];
      }

      const field_type& operator[](std::size_t i) const
      {
        checkout();
        return _data[i];
      }

      E two_norm() const
      {
        checkin();
        double norm = 0;
        PETSC_CALL(VecNorm(_v,NORM_2,&norm));
        return norm;
      }

      E one_norm() const
      {
        checkin();
        double norm = 0;
        PETSC_CALL(VecNorm(_v,NORM_1,&norm));
        return norm;
      }

      E infinity_norm() const
      {
        checkin();
        double norm = 0;
        PETSC_CALL(VecNorm(_v,NORM_INFINITY,&norm));
        return norm;
      }

      E operator*(const PetscNestedVectorContainer& y) const
      {
        checkin();
        double dotproduct = 0;
        PETSC_CALL(VecDot(_v,y._v,&dotproduct));
        return dotproduct;
      }

      PetscNestedVectorContainer& axpy(const E& a, const PetscNestedVectorContainer& y)
      {
        checkin();
        PETSC_CALL(VecAXPY(_v,a,y._v));
        return *this;
      }

      // for debugging and AMG access
      ContainerType& base ()
      {
        checkin();
        return _v;
      }

      const ContainerType& base () const
      {
        checkin();
        return _v;
      }

      operator ContainerType&()
      {
        checkin();
        return _v;
      }

      operator const ContainerType&() const
      {
        checkin();
        return _v;
      }

      /*
      iterator begin()
      {
        return container.begin();
      }


      const_iterator begin() const
      {
        return container.begin();
      }

      iterator end()
      {
        return container.end();
      }


      const_iterator end() const
      {
        return container.end();
      }
      */

      size_t flatsize() const
      {
        return N();
      }

      template<typename X>
      void std_copy_to (std::vector<X>& x) const
      {
        size_t n = flatsize();
        x.resize(n);
        checkout();
        for (size_t i=0; i<n; i++)
          x[i] = _data[i];
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        //test if x has the same size as the container
        assert (x.size() == flatsize());
        checkout();
        for (size_t i=0; i<flatsize(); i++)
          _data[i] = x[i];
      }

    private:
      Vec _v;
      mutable E* _data;

      void checkin() const
      {
        if (_data)
          {
            std::cout << "restoring " << _data << std::endl;
            PETSC_CALL(VecRestoreArray(_v,&_data));
            _data = NULL;
          }
      }

      void checkout() const
      {
        if (!_data)
          {
            PETSC_CALL(VecGetArray(_v,&_data));
            std::cout << "got " << _data << std::endl;
          }
      }


    };


    //! ISTL backend for FunctionSpace
    class PetscNestedVectorBackend
    {
    public:
      /*      enum{
      //! \brief export the block size
      BlockSize = BLOCKSIZE
      };*/

      //export Matrix Backend Type
      //typedef ISTLBCRSMatrixBackend<BLOCKSIZE,BLOCKSIZE> MatrixBackend;

      //! container construction

      // extract type of container element
      template<class C>
      struct Value
      {
        typedef typename C::ElementType Type;
      };

      //! The size type
      typedef typename std::size_t size_type;

      // get const_reference to container element
      // note: this method does not depend on T!
      template<typename C>
      static const typename C::field_type& access (const C& c, size_type i)
      {
        return c[i];
      }

      // get non const_reference to container element
      // note: this method does not depend on T!
      template<typename C>
      static typename C::field_type& access (C& c, size_type i)
      {
        return c[i];
      }
    };

    template<typename T, typename E>
    struct BackendVectorSelectorHelper<PetscNestedVectorBackend,T,E>
    {
      dune_static_assert((is_same<E,double>::value),"Petsc only supports double for now");
      typedef PetscNestedVectorContainer Type;
    };



  } // namespace PDELab
} // namespace Dune

#endif
