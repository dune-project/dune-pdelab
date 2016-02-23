// -*- tab-width: 8; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=8 sw=2 sts=2:

#ifndef DUNE_PDELAB_ORDERING_LEAFLOCALORDERING_HH
#define DUNE_PDELAB_ORDERING_LEAFLOCALORDERING_HH

#include <dune/typetree/leafnode.hh>

#include <dune/geometry/referenceelements.hh>
#include <dune/localfunctions/common/interfaceswitch.hh>
#include <dune/localfunctions/common/localkey.hh>
#include <dune/pdelab/ordering/localorderingbase.hh>

namespace Dune {
  namespace PDELab {

    //! \addtogroup Ordering
    //! \{

    template<typename OrderingTag, typename FEM, typename ES, typename DI, typename CI>
    class LeafLocalOrdering
      : public TypeTree::LeafNode
      , public LocalOrderingBase<ES,DI,CI>
    {

      template<typename>
      friend struct collect_used_geometry_types_from_cell_visitor;

      template<typename>
      friend struct extract_per_entity_sizes_from_cell_visitor;

      typedef LocalOrderingBase<ES,DI,CI> BaseT;

    public:

      typedef typename BaseT::Traits Traits;

      LeafLocalOrdering(const std::shared_ptr<const FEM>& fem, const ES& es, bool backend_blocked, typename BaseT::GFSData* gfs_data)
        : BaseT(*this,backend_blocked,gfs_data)
        , _fem(fem)
        , _es(es)
      {}

      const typename Traits::EntitySet& entitySet() const
      {
        return _es;
      }

      const typename Traits::GridView& gridView() const
      {
        return _es.gridView();
      }

      const FEM& finiteElementMap() const
      {
        return *_fem;
      }

      template<typename CodimMask>
      void collect_used_codims(CodimMask& codims) const
      {
        for (typename ES::dim_type codim = 0; codim <= ES::dimension; ++codim)
          if (_fem->hasDOFs(codim))
            codims.set(codim);
      }

      void update_a_priori_fixed_size()
      {
        this->_fixed_size = _fem->fixedSize();
      }

      void setup_fixed_size_possible()
      {
        this->_fixed_size_possible = true;
      }

    private:

      typedef FiniteElementInterfaceSwitch<
      typename FEM::Traits::FiniteElement
      > FESwitch;

      void collect_used_geometry_types_from_cell(const typename Traits::EntitySet::Element& cell)
      {
        FESwitch::setStore(_pfe,_fem->find(cell));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_pfe);

        this->_max_local_size = std::max(this->_max_local_size,coeffs.size());

        const auto& ref_el =
          ReferenceElements<typename Traits::EntitySet::Traits::CoordinateField,Traits::EntitySet::dimension>::general(cell.type());

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            Dune::GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            this->_gt_used[GlobalGeometryTypeIndex::index(gt)] = true;
            this->_codim_used.set(key.codim());
          }
      }


      void extract_per_entity_sizes_from_cell(const typename Traits::EntitySet::Element& cell,
                                              std::vector<typename Traits::SizeType>& gt_sizes)
      {
        if (this->_fixed_size_possible)
          std::fill(gt_sizes.begin(),gt_sizes.end(),0);

        FESwitch::setStore(_pfe,_fem->find(cell));

        const typename FESwitch::Coefficients& coeffs =
          FESwitch::coefficients(*_pfe);

        this->_max_local_size = std::max(this->_max_local_size,coeffs.size());

        typedef typename Traits::SizeType size_type;

        const auto& ref_el =
          ReferenceElements<typename Traits::EntitySet::Traits::CoordinateField,Traits::EntitySet::dimension>::general(cell.type());

        for (std::size_t i = 0; i < coeffs.size(); ++i)
          {
            const LocalKey& key = coeffs.localKey(i);
            GeometryType gt = ref_el.type(key.subEntity(),key.codim());
            const size_type geometry_type_index = GlobalGeometryTypeIndex::index(gt);

            const size_type entity_index = _es.indexSet().subIndex(cell,key.subEntity(),key.codim());
            const size_type index = this->_gt_entity_offsets[geometry_type_index] + entity_index;
            gt_sizes[geometry_type_index] = this->_entity_dof_offsets[index] = std::max(this->_entity_dof_offsets[index],static_cast<size_type>(key.index() + 1));
          }

        if (this->_fixed_size_possible)
          {
            for (size_type i = 0; i < gt_sizes.size(); ++i)
              {
                if (gt_sizes[i] > 0)
                  {
                    if (this->_gt_dof_offsets[i] == 0)
                      this->_gt_dof_offsets[i] = gt_sizes[i];
                    else if (this->_gt_dof_offsets[i] != gt_sizes[i])
                      {
                        this->_fixed_size_possible = false;
                        break;
                      }
                  }
                else
                  {
                    // catch special case where some grid entities sometimes don't have any DOFs associated with
                    // them, but are const size if there are any DOFs, i.e. a pattern of #(DOFs) == 0 || #(DOFs) == c > 0
                    this->_fixed_size_possible = this->_fixed_size_possible && (this->_gt_dof_offsets[i] == 0);
                  }
              }
          }
      }

      std::shared_ptr<const FEM> _fem;
      ES _es;
      typename FESwitch::Store _pfe;

    };

    template<typename GFS, typename Transformation, typename Params>
    struct leaf_gfs_to_local_ordering_descriptor<GFS,Transformation,LeafOrderingTag<Params> >
    {

      static const bool recursive = false;

      typedef LeafLocalOrdering<
        typename GFS::Traits::OrderingTag,
        typename GFS::Traits::FiniteElementMap,
        typename GFS::Traits::EntitySet,
        typename Transformation::DOFIndex,
        typename Transformation::ContainerIndex
        > transformed_type;

      typedef std::shared_ptr<transformed_type> transformed_storage_type;

      static transformed_type transform(const GFS& gfs, const Transformation& t)
      {
        return transformed_type(gfs.finiteElementMapStorage(),gfs.entitySet(),false,&const_cast<GFS*>(gfs));
      }

      static transformed_storage_type transform_storage(std::shared_ptr<const GFS> gfs, const Transformation& t)
      {
        return std::make_shared<transformed_type>(gfs->finiteElementMapStorage(),gfs->entitySet(),false,const_cast<GFS*>(gfs.get()));
      }

    };

    //! \}

  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_ORDERING_LEAFLOCALORDERING_HH
