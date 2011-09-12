// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PETSCMATRIXBACKEND_HH
#define DUNE_PETSCMATRIXBACKEND_HH

#include<vector>

#include <dune/pdelab/backend/petscutility.hh>
#include <dune/pdelab/backend/petscvectorbackend.hh>
#include <dune/pdelab/backend/petscnestedvectorbackend.hh>

#include <petscmat.h>

namespace Dune {
  namespace PDELab {

    class PetscMatrixBackend;
    class PetscMatrixContainer;

    class PetscNestedMatrixBackend;
    class PetscNestedMatrixContainer;

    class PetscMatrixAccessorBase;

    template<typename LFSV, typename LFSU>
    class PetscMatrixAccessor;

    template<typename LFSV, typename LFSU>
    class PetscNestedMatrixAccessor;

    template<typename T>
    class Pattern : public std::vector< std::set<T> >
    {

      typedef std::vector< std::set<T> > BaseT;

    public:

      Pattern (T m_, T n_)
        : BaseT(m_)
      {}

      void add_link (T i, T j)
      {
        (*this)[i].insert(j);
      }
    };


    struct petsc_types
    {
      typedef std::size_t size_type;
      typedef Dune::PDELab::Pattern<size_type> Pattern;
      typedef shared_ptr<PetscMatrixContainer> MatrixPtr;
      typedef std::vector<MatrixPtr> MatrixArray;
    };

    template<typename GFSV, typename GFSVB, typename GFSU, typename GFSUB>
    struct petsc_matrix_builder
    {
      dune_static_assert((AlwaysFalse<GFSV>::value),"Unsupported combination of trial and test space blocking");
    };

    template<typename GFSV, typename GFSU>
    petsc_types::MatrixPtr build_matrix(const GFSV& gfsv, const GFSU& gfsu, const petsc_types::Pattern& pattern)
    {
      return petsc_matrix_builder<GFSV,typename GFSV::Traits::BackendType,GFSU,typename GFSU::Traits::BackendType>::build(gfsv,gfsu,pattern);
    }

    template<typename GFSV,typename GFSU>
    struct petsc_matrix_builder<GFSV,PetscVectorBackend,GFSU,PetscVectorBackend>
      : public petsc_types
    {
      static MatrixPtr build(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern)
      {
        const size_type M = gfsv.size();
        const size_type N = gfsu.size();

        std::vector<int> nnz(M);
        std::vector<int>::iterator oit = nnz.begin();
        for (Pattern::const_iterator it = pattern.begin(); it != pattern.end(); ++it, ++oit)
          (*oit) = it->size();

        return make_shared<PetscMatrixContainer>(M,N,nnz);
      }
    };

    template<typename GFSV, typename GFSU, petsc_types::size_type row = 0, petsc_types::size_type col = 0>
    struct recursive_matrix_builder;

    template<typename GFSV,typename GFSU>
    struct petsc_matrix_builder<GFSV,PetscNestedVectorBackend,GFSU,PetscNestedVectorBackend>
      : public petsc_types
    {
      static MatrixPtr build(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern)
      {
        const size_type M = GFSV::CHILDREN;
        const size_type N = GFSU::CHILDREN;

        MatrixArray matrices(N*M);
        recursive_matrix_builder<GFSV,GFSU>::build_children(gfsv,gfsu,pattern,matrices);
        return make_shared<PetscNestedMatrixContainer>(M,N,matrices);
      }
    };

    struct recursive_matrix_builder_base
      : public petsc_types
    {

      static Pattern extract_pattern(const Pattern& p, const size_type r_l, const size_type r_h, const size_type c_l, const size_type c_h)
      {
        Pattern out(r_h-r_l,c_h-c_l);
        Pattern::const_iterator endit = p.begin() + r_h;
        Pattern::iterator outit = out.begin();
        typedef Pattern::value_type::const_iterator ColIterator;
        for (Pattern::const_iterator rit = p.begin() + r_l; rit != endit; ++rit, ++outit)
          {
            for (ColIterator cit = rit->begin(); cit != rit->end(); ++cit)
              {
                if (*cit >= c_l && *cit < c_h)
                  outit->insert(*cit-c_l);
              }
          }
        return std::move(out);
      }

    };

    template<typename GFSV, typename GFSU, petsc_types::size_type row, petsc_types::size_type col>
    struct recursive_matrix_builder
      : public recursive_matrix_builder_base
    {

      static const size_type M = GFSV::CHILDREN;
      static const size_type N = GFSU::CHILDREN;

      static MatrixPtr build_child(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern)
      {
        auto cgfsv = gfsv.template child<row>();
        auto cgfsu = gfsu.template child<col>();
        const size_type row_l = gfsv.ordering().subMap(row,0);
        const size_type row_h = gfsv.ordering().subMap(row,cgfsv.ordering().size()-1) + 1;
        const size_type col_l = gfsu.ordering().subMap(col,0);
        const size_type col_h = gfsu.ordering().subMap(col,cgfsu.ordering().size()-1) + 1;
        Pattern child_pattern = extract_pattern(pattern,row_l,row_h,col_l,col_h);
        return build_matrix(cgfsv,cgfsu,child_pattern);
      }

      template<std::size_t r = row>
      static typename enable_if<(r < M) && (col < N-1)>::type
      build_children(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern, MatrixArray& matrices)
      {
        matrices[row*N + col] = build_child(gfsv,gfsu,pattern);
        recursive_matrix_builder<GFSV,GFSU,row,col+1>::build_children(gfsv,gfsu,pattern,matrices);
      }

      template<std::size_t r = row>
      static typename enable_if<(r < M) && (col == N-1)>::type
      build_children(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern, MatrixArray& matrices)
      {
        matrices[row*N + col] = build_child(gfsv,gfsu,pattern);
        recursive_matrix_builder<GFSV,GFSU,row+1,0>::build_children(gfsv,gfsu,pattern,matrices);
      }

      template<std::size_t r = row>
      static typename enable_if<(r == M)>::type
      build_children(const GFSV& gfsv, const GFSU& gfsu, const Pattern& pattern, MatrixArray& matrices)
      {
        // end recursion
      }

    };

    class PetscMatrixContainer
      : public petsc_types
    {

      friend class PetscMatrixBackend;

      friend class PetscMatrixAccessorBase;

      template<typename LFSV, typename LFSU>
      friend class PetscMatrixAccessor;

      friend class PetscNestedMatrixContainer;

      enum AccessorState {
        clean,
        readValues,
        setValues,
        addValues,
      };

    public:
      typedef PetscScalar ElementType;
      typedef ElementType E;
      typedef Mat ContainerType;

      typedef E block_type;
      typedef block_type field_type;

      typedef PetscMatrixBackend Backend;



      //! construct container
      template<typename GO>
      PetscMatrixContainer (const GO& go)
        : _m(NULL)
        , _accessorState(clean)
        , _rowsToClear()
        , _managed(true)
      {
        dune_static_assert((is_same<typename GO::Traits::MatrixBackend,PetscMatrixBackend>::value),"Wrong matrix backend type");

        Pattern pattern(go.globalSizeV(),go.globalSizeU());
        go.fill_pattern(pattern);
        MatrixPtr matrix(build_matrix(go.testGridFunctionSpace(),go.trialGridFunctionSpace(),pattern));
        matrix->_managed = false;
        _m = matrix->_m;
      }

      PetscMatrixContainer (const size_type M, const size_type N, const std::vector<int>& nnz)
        : _accessorState(clean)
        , _rowsToClear()
        , _managed(true)
      {
        PETSC_CALL(MatCreateSeqAIJ(M,N,&(nnz[0]),&_m));
        PETSC_CALL(MatSetOption(_m,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE));
        PETSC_CALL(MatSetOption(_m,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE)); // we only ever zero our own rows
      }

      PetscMatrixContainer (const PetscMatrixContainer& rhs)
        : _accessorState(clean)
        , _rowsToClear()
        , _managed(true)
      {
        PETSC_CALL(MatDuplicate(rhs._m,MAT_COPY_VALUES,&_m));
        PETSC_CALL(MatSetOption(_m,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE));
        PETSC_CALL(MatSetOption(_m,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE)); // we only ever zero our own rows
      }

      PetscMatrixContainer (Mat m, bool managed = true)
        : _m(m)
        , _accessorState(clean)
        , _rowsToClear()
        , _managed(managed)
      {}

    protected:

      PetscMatrixContainer (bool managed = true)
        : _m(NULL)
        , _accessorState(clean)
        , _rowsToClear()
        , _managed(managed)
      {}

    public:

      virtual ~PetscMatrixContainer()
      {
        if (_managed)
          PETSC_CALL(MatDestroy(&_m));
      }

      PetscMatrixContainer& operator=(const PetscMatrixContainer& rhs)
      {
        assert(_accessorState == clean);
        assert(_rowsToClear.empty());
        PETSC_CALL(MatCopy(rhs._m,_m, DIFFERENT_NONZERO_PATTERN));
        _managed = true;
      }

      Mat& base ()
      {
        return _m;
      }

      const Mat& base () const
      {
        return _m;
      }

      operator Mat&()
      {
        return _m;
      }

      operator const Mat&() const
      {
        return _m;
      }

      size_type M() const
      {
        int size = 0;
        PETSC_CALL(MatGetLocalSize(_m,&size,PETSC_NULL));
        return size;
      }

      size_type N() const
      {
        int size = 0;
        PETSC_CALL(MatGetLocalSize(_m,PETSC_NULL,&size));
        return size;
      }

      PetscMatrixContainer& operator*= (const E& e)
      {
        PETSC_CALL(MatScale(_m,e));
        return *this;
      }

    protected:
      Mat _m;
      AccessorState _accessorState;

      std::map<double,std::vector<int> > _rowsToClear;

      bool _managed;

    private:

      void enqueue_row_clear(size_type i, double diagonal_value)
      {
        _rowsToClear[diagonal_value].push_back(i);
      }

      void access(AccessorState state)
      {
        if (_accessorState == clean)
          _accessorState = state;
        else if (_accessorState != state)
          {
            DUNE_THROW(InvalidStateException,"PETSc matrices do not allow mixing different access types between calls to flush()");
          }
      }

      virtual void flush(MatAssemblyType assemblyType)
      {
        if (_accessorState == addValues || _accessorState == setValues || assemblyType == MAT_FINAL_ASSEMBLY)
          {
            PETSC_CALL(MatAssemblyBegin(_m,assemblyType));
            PETSC_CALL(MatAssemblyEnd(_m,assemblyType));
          }
        _accessorState = clean;
        if (_rowsToClear.size() > 0) // TODO: communicate in MPI case
          {
            for (auto it = _rowsToClear.begin(); it != _rowsToClear.end(); ++it)
              {
                PETSC_CALL(MatZeroRows(_m,it->second.size(),&(it->second[0]),it->first,PETSC_NULL,PETSC_NULL));
              }
            _rowsToClear.clear();
          }
      }


      void checkin() const
      {
        /*
        if (_data)
          {
            std::cout << "restoring " << _data << std::endl;
            PETSC_CALL(VecRestoreArray(_v,&_data));
          }
        */
      }

      void checkout() const
      {
      /*
        if (!_data)
          {
            PETSC_CALL(VecGetArray(_v,&_data));
            std::cout << "got " << _data << std::endl;
          }
        */
      }


    };



    class PetscMatrixAccessorBase
      : petsc_types
    {

      template<typename,typename>
      friend class PetscNestedMatrixAccessor;

      void setup(PetscMatrixContainer::AccessorState state)
      {
        if (_state == PetscMatrixContainer::clean)
          {
            _m.access(state);
            if (state == PetscMatrixContainer::readValues)
              {
                PETSC_CALL(MatGetValues(_m.base(),_M,&(_rows[0]),_N,&(_cols[0]),&(_vals[0])));
              }
            _state = state;
          } else if (_state != state)
          {
            DUNE_THROW(InvalidStateException,"Attempt to mix different access patterns within accessor");
          }
      }

    public:

      PetscMatrixAccessorBase(PetscMatrixContainer& m, size_type M, size_type N, size_type row_offset = 0, size_type col_offset = 0)
        : _m(m)
        , _M(M)
        , _N(N)
        , _rows(_M)
        , _cols(_N)
        , _vals(_M * _N)
        , _state(PetscMatrixContainer::clean)
        , _row_offset(row_offset)
        , _col_offset(col_offset)
      {}

      ~PetscMatrixAccessorBase()
      {
        switch (_state)
          {
          case PetscMatrixContainer::setValues:
            PETSC_CALL(MatSetValues(_m.base(),_M,&(_rows[0]),_N,&(_cols[0]),&(_vals[0]),INSERT_VALUES));
          case PetscMatrixContainer::addValues:
            PETSC_CALL(MatSetValues(_m.base(),_M,&(_rows[0]),_N,&(_cols[0]),&(_vals[0]),ADD_VALUES));
          default:
            break;
          }
      }

      double get(size_type i, size_type j)
      {
        setup(PetscMatrixContainer::readValues);
        return _vals[i * _cols.size() + j];
      }

      void set(size_type i, size_type j, double v)
      {
        setup(PetscMatrixContainer::setValues);
        _vals[i * _cols.size() + j] = v;
      }

      void add(size_type i, size_type j, double v)
      {
        setup(PetscMatrixContainer::addValues);
        _vals[i * _cols.size() + j] += v;
      }

    private:

      PetscMatrixContainer& _m;

    protected:
      const size_type _M;
      const size_type _N;
      std::vector<int> _rows;
      std::vector<int> _cols;

    private:
      std::vector<double> _vals;
      PetscMatrixContainer::AccessorState _state;
      const size_type _row_offset;
      const size_type _col_offset;

    };



    template<typename LFSV,typename LFSU>
    class PetscMatrixAccessor
      : public PetscMatrixAccessorBase
    {

      friend class PetscMatrixBackend;

    protected:

      PetscMatrixAccessor(PetscMatrixContainer& m, const LFSV& lfsv, const LFSU& lfsu, size_type row_offset = 0, size_type col_offset = 0)
        : PetscMatrixAccessorBase(m,lfsv.localVectorSize(),lfsu.localVectorSize(),row_offset,col_offset)
      {
        int m_offset = 0;
        int n_offset = 0;
        PETSC_CALL(MatGetOwnershipRange(m.base(),&m_offset,PETSC_NULL));
        PETSC_CALL(MatGetOwnershipRangeColumn(m.base(),&n_offset,PETSC_NULL));
        m_offset -= row_offset;
        n_offset -= col_offset;
        for (size_type i = 0; i < _M; ++i)
          _rows[i] = m_offset + lfsv.globalIndex(i);
        for (size_type i = 0; i < _N; ++i)
          _cols[i] = n_offset + lfsu.globalIndex(i);
      }

    };


    class PetscNestedMatrixContainer
      : public PetscMatrixContainer
    {

      friend class PetscNestedMatrixBackend;

      template<typename LFSV, typename LFSU>
      friend class PetscMatrixAccessor;

    public:
      typedef PetscNestedMatrixBackend Backend;



      PetscNestedMatrixContainer(const size_type M, const size_type N, const MatrixArray& children)
      {
        Mat matrices[N*M];
        for (size_type i = 0; i < N*M; ++i)
          matrices[i] = children[i]->base();

        PETSC_CALL(MatCreateNest(PETSC_COMM_SELF,M,PETSC_NULL,N,PETSC_NULL,matrices,&_m));
        PETSC_CALL(MatSetOption(_m,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE));
        PETSC_CALL(MatSetOption(_m,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE)); // we only ever zero our own rows
      }

      //! construct container
      template<typename GO>
      PetscNestedMatrixContainer (const GO& go)
      {
        dune_static_assert((is_same<typename GO::Traits::MatrixBackend,PetscNestedMatrixBackend>::value),"Wrong matrix backend type");

        Pattern pattern(go.globalSizeV(),go.globalSizeU());
        go.fill_pattern(pattern);
        MatrixPtr matrix(build_matrix(go.testGridFunctionSpace(),go.trialGridFunctionSpace(),pattern));
        matrix->_managed = false;
        _m = matrix->_m;
        PETSC_CALL(MatSetOption(_m,MAT_IGNORE_ZERO_ENTRIES,PETSC_TRUE));
        PETSC_CALL(MatSetOption(_m,MAT_NO_OFF_PROC_ZERO_ROWS,PETSC_TRUE)); // we only ever zero our own rows
      }

      virtual ~PetscNestedMatrixContainer()
      {
      }

      PetscMatrixContainer subMatrix(size_type i, size_type j)
      {
        Mat sub_mat;
        PETSC_CALL(MatNestGetSubMat(_m,i,j,&sub_mat));
        return PetscMatrixContainer(sub_mat,false);
      }

      virtual void flush(MatAssemblyType assemblyType)
      {
        if (_accessorState == addValues || _accessorState == setValues || assemblyType == MAT_FINAL_ASSEMBLY)
          {
            PETSC_CALL(MatAssemblyBegin(_m,assemblyType));
            PETSC_CALL(MatAssemblyEnd(_m,assemblyType));
          }
        _accessorState = clean;
        if (_rowsToClear.size() > 0) // TODO: communicate in MPI case
          {
            int M;
            int N;
            Mat** sm;
            PETSC_CALL(MatNestGetSubMats(_m,&M,&N,&sm));
            int row_offsets[M+1];
            int col_offsets[N+1];
            row_offsets[0] = col_offsets[0] = 0;

            for (int i = 0; i < M; ++i)
              {
                int MM;
                PETSC_CALL(MatGetLocalSize(sm[i][0],&MM,PETSC_NULL));
                row_offsets[i+1] = row_offsets[i] + MM;
              }

            for (int j = 0; j < N; ++j)
              {
                int NN;
                PETSC_CALL(MatGetLocalSize(sm[0][j],PETSC_NULL,&NN));
                col_offsets[j+1] = col_offsets[j] + NN;
              }

            for (auto it = _rowsToClear.begin(); it != _rowsToClear.end(); ++it)
              {
                auto rows = it->second;
                std::sort(rows.begin(),rows.end());
                auto rowit = rows.begin();
                auto rowend = rowit;
                int row_block = 0;
                do {
                  rowend = std::find_if(rowit,rows.end(),std::bind1st(std::less_equal<int>(),row_offsets[row_block+1]));
                  for (; rowit != rowend; ++rowit)
                    {
                      int local_row = (*rowit) - row_offsets[row_block];
                      for (int j = 0; j < N; ++j)
                        {
                          PETSC_CALL(MatZeroRows(sm[row_block][j],1,&local_row,0.0,PETSC_NULL,PETSC_NULL));
                          if (col_offsets[j] <= *rowit && *rowit < col_offsets[j+1])
                            {
                              PETSC_CALL(MatSetValue(sm[row_block][j],local_row,*rowit-col_offsets[j],it->first,INSERT_VALUES));
                              PETSC_CALL(MatAssemblyBegin(sm[row_block][j],MAT_FINAL_ASSEMBLY));
                              PETSC_CALL(MatAssemblyEnd(sm[row_block][j],MAT_FINAL_ASSEMBLY));
                            }
                        }
                    }
                  ++row_block;
                } while (rowit != rows.end());
              }
            PETSC_CALL(MatAssemblyBegin(_m,MAT_FINAL_ASSEMBLY));
            PETSC_CALL(MatAssemblyEnd(_m,MAT_FINAL_ASSEMBLY));
            _rowsToClear.clear();
          }
      }

    };



    class PetscMatrixBackend
    {
    public:

      template<typename LFSV, typename LFSU>
      class Accessor
        : public PetscMatrixAccessor<LFSV,LFSU>
      {

      public:

        Accessor(PetscMatrixContainer& m, const LFSV& lfsv, const LFSU& lfsu)
          : PetscMatrixAccessor<LFSV,LFSU>(m,lfsv,lfsu)
        {}

      };

      template<typename E>
      class Matrix
        : public PetscMatrixContainer
      {

        dune_static_assert((is_same<E,double>::value),"The PETSc backend currently only supports double as a field type");

      public:

        template<typename GridOperator>
        Matrix(const GridOperator& go)
        : PetscMatrixContainer(go)
        {}
      };

      typedef petsc_types::Pattern Pattern;

      template<class C>
      struct Value
      {
        typedef typename C::ElementType Type;
      };

      //! The size type
      typedef typename std::size_t size_type;

      static void clear_row (size_type i, PetscMatrixContainer& c, double diagonal_value = 0.0)
      {
        c.enqueue_row_clear(i,diagonal_value);
      }

      static void flush(PetscMatrixContainer& c)
      {
        c.flush(MAT_FLUSH_ASSEMBLY);
      }

      static void finalize(PetscMatrixContainer& c)
      {
        c.flush(MAT_FINAL_ASSEMBLY);
      }

    };



    template<typename LFSV,typename LFSU>
    class PetscNestedMatrixAccessor
      : public petsc_types
    {

      struct get_child_offsets
        : public TypeTree::DirectChildrenVisitor
        , public TypeTree::DynamicTraversal
      {
        get_child_offsets(std::vector<std::size_t>& local_offsets, std::vector<std::size_t>& global_offsets, std::vector<std::vector<int> >& indices)
          : _local_offsets(local_offsets)
          , _global_offsets(global_offsets)
          , _indices(indices)
        {}

        template<typename T, typename Child, typename TreePath, typename ChildIndex>
        void beforeChild(const T& t, const Child& child, TreePath treePath, ChildIndex childIndex)
        {
          _local_offsets[childIndex+1] = _local_offsets[childIndex] + child.size();
          _global_offsets[childIndex+1] = _global_offsets[childIndex] + child.gridFunctionSpace().size();
          _indices[childIndex].resize(child.size());
          for (std::size_t i = 0; i < child.size(); ++i)
            _indices[childIndex][i] = child.globalIndex(i) - _global_offsets[childIndex];
        }

        std::vector<std::size_t>& _local_offsets;
        std::vector<std::size_t>& _global_offsets;
        std::vector<std::vector<int> >& _indices;
      };

      friend class PetscNestedMatrixBackend;

      static const std::size_t M = LFSV::CHILDREN;
      static const std::size_t N = LFSU::CHILDREN;

      PetscNestedMatrixAccessor(PetscNestedMatrixContainer& m, const LFSV& lfsv, const LFSU& lfsu)
        : _m(m)
        , _lfsv(lfsv)
        , _lfsu(lfsu)
        , _accessors(N*M)
        , _matrices(N*M)
        , _local_row_offsets(M+1)
        , _local_col_offsets(N+1)
        , _global_row_offsets(M+1)
        , _global_col_offsets(N+1)
        , _row_indices(M)
        , _col_indices(N)
      {
        get_child_offsets rows(_local_row_offsets,_global_row_offsets,_row_indices);
        TypeTree::applyToTree(lfsv,rows);
        get_child_offsets cols(_local_col_offsets,_global_col_offsets,_col_indices);
        TypeTree::applyToTree(lfsu,cols);
      }

      PetscMatrixAccessorBase& accessor(size_type& i, size_type& j)
      {
        auto row = std::find_if(_local_row_offsets.begin(),_local_row_offsets.end(),std::bind1st(std::less<size_type>(),i));
        --row;
        auto col = std::find_if(_local_col_offsets.begin(),_local_col_offsets.end(),std::bind1st(std::less<size_type>(),j));
        --col;
        const size_type mi = row - _local_row_offsets.begin();
        const size_type mj = col - _local_col_offsets.begin();

        shared_ptr<PetscMatrixAccessorBase> a;
        if (!(a = _accessors[mi*N + mj]))
          {
            Mat submat;
            PETSC_CALL(MatNestGetSubMat(_m.base(),mi,mj,&submat));
            shared_ptr<PetscMatrixContainer> c(make_shared<PetscMatrixContainer>(submat,false));
            a = make_shared<PetscMatrixAccessorBase>(*c,(*(row+1))-(*row),(*(col+1))-(*col),_global_row_offsets[mi],_global_row_offsets[mj]);
            swap(a->_rows,_row_indices[mi]);
            swap(a->_cols,_col_indices[mj]);
            _accessors[mi*N + mj] = a;
            _matrices[mi*N + mj] = c;
          }
        i -= *row;
        j -= *col;
        return *a;
      }

    public:

      typedef PetscNestedMatrixContainer::size_type size_type;

      double get(size_type i, size_type j)
      {
        return accessor(i,j).get(i,j);
      }

      void set(size_type i, size_type j, double v)
      {
        accessor(i,j).set(i,j,v);
      }

      void add(size_type i, size_type j, double v)
      {
        accessor(i,j).add(i,j,v);
      }

    private:

      PetscNestedMatrixContainer& _m;
      const LFSV& _lfsv;
      const LFSU& _lfsu;
      std::vector<shared_ptr<PetscMatrixAccessorBase> > _accessors;
      std::vector<shared_ptr<PetscMatrixContainer> > _matrices;
      std::vector<std::size_t> _local_row_offsets;
      std::vector<std::size_t> _local_col_offsets;
      std::vector<std::size_t> _global_row_offsets;
      std::vector<std::size_t> _global_col_offsets;
      std::vector<std::vector<int> > _row_indices;
      std::vector<std::vector<int> > _col_indices;

    };




    class PetscNestedMatrixBackend
      : public PetscMatrixBackend
    {
    public:

      template<typename E>
      class Matrix
        : public PetscNestedMatrixContainer
      {

        dune_static_assert((is_same<E,double>::value),"The PETSc backend currently only supports double as a field type");

      public:

        template<typename GridOperator>
        Matrix(const GridOperator& go)
        : PetscNestedMatrixContainer(go)
        {}
      };


      template<typename LFSV, typename LFSU>
      class Accessor
        : public PetscNestedMatrixAccessor<LFSV,LFSU>
      {

      public:

        Accessor(PetscNestedMatrixContainer& m, const LFSV& lfsv, const LFSU& lfsu)
          : PetscNestedMatrixAccessor<LFSV,LFSU>(m,lfsv,lfsu)
        {}

      };


    };

  } // namespace PDELab
} // namespace Dune

#endif
