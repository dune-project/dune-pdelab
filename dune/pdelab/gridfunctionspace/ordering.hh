// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERING_HH
#define DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERING_HH

#include <vector>

namespace Dune {
  namespace PDELab {

    //! \addtogroup GridFunctionSpace
    //! \ingroup PDELab
    //! \{

    //===============================================================
    // Utilities for the power and composite gfs
    // ===============================================================

    //! \brief Indicates lexicographics ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    struct GridFunctionSpaceLexicographicMapper {};

    //! \brief Indicates using block-wise ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! The exact blocking structure can be passed as template parameters
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    template<int s0 = 1, int s1 = 1, int s2 = 1, int s3 = 1, int s4 = 1, int s5 = 1, int s6 = 1, int s7 = 1, int s8 = 1, int s9 = 1>
    struct GridFunctionSpaceComponentBlockwiseMapper
    {
      static const int size[];
      static const int offset[];
    };
    template<int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    const int GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>::
    size[] = { s0, s1, s2, s3, s4, s5, s6, s7, s8, s9 };
    template<int s0, int s1, int s2, int s3, int s4, int s5, int s6, int s7, int s8, int s9>
    const int GridFunctionSpaceComponentBlockwiseMapper<s0,s1,s2,s3,s4,s5,s6,s7,s8,s9>::
    offset[] = { 0, s0, s0+s1, s0+s1+s2, s0+s1+s2+s3, s0+s1+s2+s3+s4,
                 s0+s1+s2+s3+s4+s5, s0+s1+s2+s3+s4+s5+s6, s0+s1+s2+s3+s4+s5+s6+s7,
                 s0+s1+s2+s3+s4+s5+s6+s7+s8, s0+s1+s2+s3+s4+s5+s6+s7+s8+s9 };

    //! \brief Indicates using block-wise ordering of the unknowns for composite
    //! grid function spaces.
    //!
    //! this class may be used to pass compile-time
    //! parameters to the implementation of
    //! \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or
    //! \link CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    struct GridFunctionSpaceBlockwiseMapper : GridFunctionSpaceComponentBlockwiseMapper<> {};


    namespace DynamicBlockwiseMapperImp{
      //! Meta program to compute the individual children's offsets
      //! for the dynamic blockwise mapper. Interior node of template
      //! tree.
      template<typename T, int n, int i>
      struct GetChildOffsetsMetaProgram
      {
        template<class EntityType, class CO>
        static void getChildOffsets(const T& t, const EntityType& e, CO & childOffsets)
        {
          // Fill vector with global coefficients
          std::vector<typename T::Traits::SizeType> global;
          t.template child<i>().dataHandleGlobalIndices(e,global);

          // Update new offset for child
          typename T::Traits::SizeType entity_index = t.gridview().indexSet().index(e);
          childOffsets[i][entity_index+1] = childOffsets[i][entity_index] + global.size();

          // Continue metaprogram loop
          GetChildOffsetsMetaProgram<T,n,i+1>::getChildOffsets(t,e,childOffsets);
        }
      };

      //! Meta program to compute the individual children's offsets
      //! for the dynamic blockwise mapper. Leaf node of template
      //! tree.
      template<typename T, int n>
      struct GetChildOffsetsMetaProgram<T,n,n>
      {
        template<class EntityType, class CO>
        static void getChildOffsets (const T& t, const EntityType& e, CO & childOffsets)
        {
        }
      };
    } // namespace DynamicBlockwiseMapperImp


    /**
       \brief Indicates using block-wise ordering of the unknowns for
       composite grid function spaces. However, in this "dynamic"
       version, the number of DOFs may vary for each node of the local
       function space tree and each grid cell. Notice, that this
       behavior can only be achieved by storing individual offset values
       for each grid cell and each node treated with this mapping.

       This class may be used to pass compile-time parameters to the
       implementation of

       \link PowerGridFunctionSpace PowerGridFunctionSpace \endlink or \link
       CompositeGridFunctionSpace CompositeGridFunctionSpace \endlink
    */
    struct GridFunctionSpaceDynamicBlockwiseMapper
    {
    };

   //! \} group GridFunctionSpace
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_GRIDFUNCTIONSPACE_ORDERING_HH
