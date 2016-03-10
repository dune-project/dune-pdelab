// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_PDELAB_LOCALOPERATOR_LAPLACE_HH
#define DUNE_PDELAB_LOCALOPERATOR_LAPLACE_HH
#warning This file is deprecated and will be removed after the Dune-PDELab 2.4 release! Use the ConvectionDiffusionFEM local operator from dune/pdelab/localoperator/convectiondiffusionfem.hh instead!

#include <vector>

#include <dune/common/fvector.hh>

#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>

#include <dune/localfunctions/common/interfaceswitch.hh>

#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/flags.hh>

namespace Dune {
  namespace PDELab {
    //! \addtogroup LocalOperator
    //! \ingroup PDELab
    //! \{

    /** a local operator for solving the Laplace equation
     *
     * \f{align*}{
     *           - \Delta u &=& 0 \mbox{ in } \Omega,          \\
     *  -\nabla u \cdot \nu &=& 0 \mbox{ on } \partial\Omega_N \\
     * \f}
     * with conforming finite elements on all types of grids in any dimension.
     *
     * In other words, it only assembles the Laplace matrix.
     *
     */
    class Laplace
    : public FullVolumePattern,
      public LocalOperatorDefaultFlags
    {
    public:
      // pattern assembly flags
      enum { doPatternVolume = true };

      // residual assembly flags
      enum { doAlphaVolume = true };

      /** \brief Constructor
       *
       * \param quadOrder Order of the quadrature rule used for integrating over the element
       */
      DUNE_DEPRECATED_MSG("Deprecated in Dune-PDELab 2.4, use the local operator ConvectionDiffusionFEM instead!")
      Laplace (unsigned int quadOrder)
      : quadOrder_(quadOrder)
      {}

      /** \brief Compute Laplace matrix times a given vector for one element
       *
       * This is used for matrix-free algorithms for the Laplace equation
       *
       * \param [in]  eg The grid element we are assembling on
       * \param [in]  lfsu Local ansatz function space basis
       * \param [in]  lfsv Local test function space basis
       * \param [in]  x Input vector
       * \param [out] r The product of the Laplace matrix times x
       */
      template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
      void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
      {
        // domain and range field type
        typedef FiniteElementInterfaceSwitch<
        typename LFSU::Traits::FiniteElementType
        > FESwitch;
        typedef BasisInterfaceSwitch<
        typename FESwitch::Basis
        > BasisSwitch;
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;

        // dimensions
        static const int dimLocal = EG::Geometry::mydimension;
        static const int dimGlobal = EG::Geometry::coorddimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dimLocal>& rule =
        Dune::QuadratureRules<DF,dimLocal>::rule(gt,quadOrder_);

        // loop over quadrature points
        for(typename Dune::QuadratureRule<DF,dimLocal>::const_iterator it =
          rule.begin(); it!=rule.end(); ++it)
        {
          // evaluate gradient of shape functions
          // (we assume Galerkin method lfsu=lfsv)
          std::vector<Dune::FieldMatrix<RF,1,dimGlobal> >
          gradphiu(lfsu.size());
          BasisSwitch::gradient(FESwitch::basis(lfsu.finiteElement()),
                                eg.geometry(), it->position(), gradphiu);
          std::vector<Dune::FieldMatrix<RF,1,dimGlobal> >
          gradphiv(lfsv.size());
          BasisSwitch::gradient(FESwitch::basis(lfsv.finiteElement()),
                                eg.geometry(), it->position(), gradphiv);

          // compute gradient of u
          Dune::FieldVector<RF,dimGlobal> gradu(0.0);
          for (size_t i=0; i<lfsu.size(); i++)
            gradu.axpy(x(lfsu,i),gradphiu[i][0]);

          // integrate grad u * grad phi_i
          RF factor = r.weight() * it->weight() * eg.geometry().integrationElement(it->position());
          for (size_t i=0; i<lfsv.size(); i++)
            r.rawAccumulate(lfsv,i,(gradu*gradphiv[i][0])*factor);
        }
      }

      /** \brief Compute the Laplace stiffness matrix for the element given in 'eg'
       *
       * \tparam M Type of the element stiffness matrix
       *
       * \param [in]  eg The grid element we are assembling on
       * \param [in]  lfsu Local ansatz function space basis
       * \param [in]  lfsv Local test function space basis
       * \param [in]  x Current configuration; gets ignored for linear problems like this one
       * \param [out] matrix Element stiffness matrix
       */
      template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
      void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, M & matrix) const
      {
        // Switches between local and global interface
        typedef FiniteElementInterfaceSwitch<
          typename LFSU::Traits::FiniteElementType
          > FESwitch;
        typedef BasisInterfaceSwitch<
          typename FESwitch::Basis
          > BasisSwitch;

        // domain and range field type
        typedef typename BasisSwitch::DomainField DF;
        typedef typename BasisSwitch::RangeField RF;
        typedef typename LFSU::Traits::SizeType size_type;

        // dimensions
        const int dim = EG::Entity::dimension;

        // select quadrature rule
        Dune::GeometryType gt = eg.geometry().type();
        const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,quadOrder_);

        // loop over quadrature points
        for (typename QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
        {
          std::vector<Dune::FieldMatrix<RF,1,dim> > gradphi(lfsu.size());
          BasisSwitch::gradient(FESwitch::basis(lfsu.finiteElement()),
                                eg.geometry(), it->position(), gradphi);

          // geometric weight
          RF factor = it->weight() * eg.geometry().integrationElement(it->position());

          for (size_type i=0; i<lfsu.size(); i++)
          {
            for (size_type j=0; j<lfsv.size(); j++)
            {
              // integrate grad u * grad phi
              matrix.accumulate(lfsv,j,lfsu,i, gradphi[i][0] * gradphi[j][0] * factor);
            }
          }
        }
      }

    protected:
      // Quadrature rule order
      unsigned int quadOrder_;
    };

     //! \} group LocalOperator
  } // namespace PDELab
} // namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_LAPLACE_HH
