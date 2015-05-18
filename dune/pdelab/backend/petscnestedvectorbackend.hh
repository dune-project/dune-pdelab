// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PETSCNESTEDVECTORBACKEND_HH
#define DUNE_PETSCNESTEDVECTORBACKEND_HH

#if HAVE_PETSC

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

      typedef std::vector<std::shared_ptr<PetscVectorContainer> > ChildContainer;

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
          _children[childIndex] = std::make_shared<PetscVectorContainer>(child);
        }

        ChildContainer& _children;
      };

      template<typename GFS>
      PetscNestedVectorContainer (const GFS& gfs)
        : _data(gfs.globalSize())
        , _sub_data(GFS::CHILDREN)
        , _checkedIn(true)
      {
        const size_t _blockCount = GFS::CHILDREN;
        ChildContainer children(_blockCount);
        create_child_vectors ccv(children);
        TypeTree::applyToTree(gfs,ccv);
        Vec child_array[_blockCount];
        for (size_t i = 0; i < _blockCount; ++i)
          child_array[i] = children[i]->_v;
        PETSC_CALL(VecCreateNest(PETSC_COMM_SELF,_blockCount,PETSC_NULL,child_array,&_v));
      }

      template<typename GFS>
      PetscNestedVectorContainer (const GFS& gfs, const E& e)
        : _data(gfs.globalSize())
        , _sub_data(GFS::CHILDREN)
        , _checkedIn(true)
      {
        const size_t _blockCount = GFS::CHILDREN;
        ChildContainer children(_blockCount);
        create_child_vectors ccv(children);
        TypeTree::applyToTree(gfs,ccv);
        Vec child_array[_blockCount];
        for (size_t i = 0; i < _blockCount; ++i)
          child_array[i] = *(children[i]);
        PETSC_CALL(VecCreateNest(PETSC_COMM_SELF,_blockCount,PETSC_NULL,child_array,&_v));
        this->operator=(e);
      }

      PetscNestedVectorContainer (const PetscNestedVectorContainer& rhs)
        : _data(rhs._data.size())
        , _sub_data(rhs._sub_data.size())
        , _checkedIn(true)
      {
        rhs.checkin();
        PETSC_CALL(VecDuplicate(rhs._v,&_v));
        PETSC_CALL(VecCopy(rhs._v,_v));
      }

      PetscNestedVectorContainer& operator= (const PetscNestedVectorContainer& rhs)
      {
        checkin();
        rhs.checkin();
        PETSC_CALL(VecCopy(rhs._v,_v));
        _checkedIn = true;
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
        PetscScalar dotproduct = 0;
        PETSC_CALL(VecTDot(_v,y._v,&dotproduct));
        return dotproduct;
      }

      E dot(const PetscNestedVectorContainer& y) const
      {
         checkin();
         PetscScalar dotproduct = 0;
         PETSC_CALL(VecDot(_v,y._v,&dotproduct));
         return dotproduct;
      }


      PetscNestedVectorContainer& axpy(const E& a, const PetscNestedVectorContainer& y)
      {
        checkin();
        y.checkin();
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
        checkout();
        std::copy(_data.begin(),_data.end(),x.begin());
      }

      template<typename X>
      void std_copy_from (const std::vector<X>& x)
      {
        checkout();
        std::copy(x.begin(),x.end(),_data.begin());
      }

      PetscVectorContainer subVector(size_type i)
      {
        checkin();
        Vec sub_vec;
        PETSC_CALL(VecNestGetSubVec(_v,i,&sub_vec));
        return PetscVectorContainer(sub_vec,false);
      }

    private:
      Vec _v;
      mutable std::vector<PetscScalar> _data;
      mutable std::vector<PetscScalar*> _sub_data;
      mutable bool _checkedIn;

      void checkin() const
      {
        if (!_checkedIn)
          {
            int N;
            Vec* sv;
            PETSC_CALL(VecNestGetSubVecs(_v,&N,&sv));
            auto it = _data.begin();
            for (int n = 0; n < N; ++n)
              {
                int CN;
                PETSC_CALL(VecGetLocalSize(sv[n],&CN));
                auto endit = it + CN;
                for(auto childit = _sub_data[n]; it != endit; ++it, ++childit)
                  *childit = *it;
                PETSC_CALL(VecRestoreArray(sv[n],&(_sub_data[n])));
              }
            _checkedIn = true;
          }
      }

      void checkout() const
      {
        if (_checkedIn)
          {
            int N;
            Vec* sv;
            PETSC_CALL(VecNestGetSubVecs(_v,&N,&sv));
            auto it = _data.begin();
            for (int n = 0; n < N; ++n)
              {
                int CN;
                PETSC_CALL(VecGetArray(sv[n],&(_sub_data[n])));
                PETSC_CALL(VecGetLocalSize(sv[n],&CN));
                auto endit = it + CN;
                for(auto childit = _sub_data[n]; it != endit; ++it, ++childit)
                  *it = *childit;
              }
            _checkedIn = false;
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
      typedef PetscNestedVectorContainer Type;
    };



  } // namespace PDELab
} // namespace Dune

#endif // HAVE_PETSC

#endif // DUNE_PETSCNESTEDVECTORBACKEND_HH
