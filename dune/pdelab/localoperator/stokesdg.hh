// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_STOKESDG_HH
#define DUNE_PDELAB_STOKESDG_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/common/configparser.hh>
#include <dune/grid/common/quadraturerules.hh>

#include "../common/geometrywrapper.hh"
#include "../gridoperatorspace/gridoperatorspace.hh"
#include "defaultimp.hh"
#include "pattern.hh"
#include "flags.hh"

#ifndef VBLOCK
#define VBLOCK 0
#endif
#define PBLOCK (- VBLOCK + 1)

namespace Dune {
    namespace PDELab {

        /**
           These are the boundary condition types as to be returned by
           the employed boundary type function. 

           Possible types:

           <ul>

           <li>\a DoNothing : Do not evaluate boundary integrals.

           <li>\a VelocityDirichlet : Dirichlet conditions for velocity.

           <li>\a PressureDirichlet : Natural Neumann conditions for the
           impulse flux. These are equivalent to a fixed pressure
           condition \b if \f$ \forall i : n \cdot \nabla v_i = 0 \f$.

           </ul>
         */
        struct StokesBoundaryCondition {
            enum Type {
                DoNothing = 0,
                VelocityDirichlet = 1,
                PressureDirichlet = 2
            };
        };

        /** 
            \brief This is the default implementation for the interior
            penalty factor.

            It computes the factor according to \f$
            \frac{\sigma}{|e|^\beta} \f$ for each face \f$ e \f$. It
            assumes that the intersection geometries passed to the
            local assembler allow to compute \f$|e|\f$ via \a
            ig.geometry().volume().
        */
        template <typename RF>
        class DefaultInteriorPenalty
        {
        private:
            RF beta;
            RF sigma;
            RF mu;
        public:

            DefaultInteriorPenalty(const std::string method, const RF mu_)
                : mu(mu_)
            {
                std::string s = method;
                std::transform(s.begin(), s.end(), s.begin(), tolower);

                // nipg (epsilon=1) 2d p1 -> Klaus sagt sollte auch sigma 1 klappen
                if (s.find("nipg") != std::string::npos)
                {
                    beta = 1;
                    if (sscanf(s.c_str(), "nipg %lg", &sigma) != 1)
                        sigma = 3.9;
                    return;
                }
                // sipg (epsilon=-1) 2d p1 -> Klaus sagt sigma=3.9irgendwas
                if (s.find("sipg") != std::string::npos)
                {
                    beta = 1;
                    if (sscanf(s.c_str(), "sipg %lg", &sigma) != 1)
                        sigma = 3.9;
                    return;
                }
                // obb sigma = 0, epsilon = 
                if (s == "obb")
                {
                    beta = 1;
                    sigma = 0;
                    return;
                }
                // extract parameters
                {
                    int epsilon;
                    if (3 == sscanf(s.c_str(), "%d %lg %lg", &epsilon, &sigma, &beta))
                        return;
                }
                DUNE_THROW(Dune::Exception, "Unknown DG type " << method);
            }

            DefaultInteriorPenalty(const Dune::ParameterTree & config, const RF mu_)
                : mu(mu_)
            {
                beta = config.get<double>("beta");
                sigma = config.get<double>("ip_sigma");
            }

            template<typename I>
            RF getFaceIP(const I & ig) const
            {
                return mu * sigma / std::pow(ig.geometry().volume(),beta);
            }
        };

        
        /** \brief a local operator for solving the stokes equation using a DG discretization
            
            \tparam F velocity source term function
            \tparam B boundary condition function
            \tparam V dirichlet velocity boundary condition function
            \tparam P dirichlet pressure boundary condition function
            \tparam IP a class providing the interior penalty factor for each face
         */
        template<typename F, typename B, typename V, typename P, 
                 typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
        class StokesDG :
            public LocalOperatorDefaultFlags,
            public FullSkeletonPattern, public FullVolumePattern
            //
            ,public JacobianBasedAlphaVolume< StokesDG<F,B,V,P,IP> >
            ,public JacobianBasedAlphaSkeleton< StokesDG<F,B,V,P,IP> >
            ,public JacobianBasedAlphaBoundary< StokesDG<F,B,V,P,IP> >
        {
            typedef StokesBoundaryCondition BC;
            typedef typename V::Traits::RangeFieldType RF;
        public:
            typedef IP InteriorPenaltyFactor;

            // pattern assembly flags
            enum { doPatternVolume = true };
            enum { doPatternSkeleton = true };

            // call the assembler for each face only once
            enum { doSkeletonTwoSided = false };

            // residual assembly flags
            enum { doAlphaVolume    = true };
            enum { doAlphaSkeleton  = true };
            enum { doAlphaBoundary  = true };
            enum { doLambdaVolume   = true };
            enum { doLambdaBoundary = true };

            StokesDG (const std::string & method, const IP & ip_factor_, const RF mu_,
                      const F & _f, const B & _b, const V & _v, const P & _p, int _qorder=4) :
                f(_f), b(_b), v(_v), p(_p), qorder(_qorder), mu(mu_), ip_factor(ip_factor_)
            {
                std::string s = method;
                std::transform(s.begin(), s.end(), s.begin(), tolower);

                // nipg (epsilon=1) 2d p1 -> Klaus sagt sollte auch sigma 1 klappen
                if (s.find("nipg") != std::string::npos)
                {
                    epsilon = 1;
                    return;
                }
                // sipg (epsilon=-1) 2d p1 -> Klaus sagt sigma=3.9irgendwas
                if (s.find("sipg") != std::string::npos)
                {
                    epsilon = -1;
                    return;
                }
                // obb sigma = 0, epsilon = 
                if (s == "obb")
                {
                    epsilon = 1;
                    return;
                }
                // extract parameters
                {
                    double sigma, beta;
                    if (3 == sscanf(s.c_str(), "%d %lg %lg", &epsilon, &sigma, &beta))
                        return;
                }
                DUNE_THROW(Dune::Exception, "Unknown DG type " << method);
            }

            StokesDG (const Dune::ParameterTree & configuration,const IP & ip_factor_, const RF mu_,
                      const F & _f, const B & _b, const V & _v, const P & _p, int _qorder=4) :
                f(_f), b(_b), v(_v), p(_p), qorder(_qorder), mu(mu_), ip_factor(ip_factor_)
            {
                epsilon = configuration.get<int>("epsilon");
            }

            // volume integral depending only on test functions,
            // contains f on the right hand side
            template<typename EG, typename LFSV, typename R>
            void lambda_volume (const EG& eg, const LFSV& lfsv, R& r) const
            {
                // dimensions
                static const unsigned int dim = EG::Geometry::dimension;
                static const unsigned int dimw = EG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert
                  ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

                typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
                const LFSV_PFS_V& lfsv_pfs_v = lfsv.template getChild<VBLOCK>();

                dune_static_assert
                  ((LFSV_PFS_V::CHILDREN == dim),"You seem to use the wrong function space for StokesDG");

                // we assume all velocity components are the same type
                typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
                const LFSV_V& lfsv_v = lfsv_pfs_v.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();

                // domain and range field type
                typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
                typedef typename BasisSwitch_V::DomainField DF;
                typedef typename BasisSwitch_V::Range RT;
                typedef typename BasisSwitch_V::RangeField RF;
                typedef typename BasisSwitch_V::Range Range_V;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gt = eg.geometry().type();
                const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
                
                // loop over quadrature points
                for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    const Dune::FieldVector<DF,dim> local = it->position();
                    const Dune::FieldVector<DF,dimw> global = eg.geometry().global(local);
                    
                    // values of velocity shape functions
                    std::vector<RT> phi_v(vsize);
                    FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);

                    const RF weight = it->weight() * eg.geometry().integrationElement(it->position());

                    // evaluate source term
                    typename F::Traits::RangeType fval;
                    f.evaluateGlobal(global,fval);
                    
                    //================================================//
                    // \int (f*v)
                    //================================================//
                    const RF factor = mu * weight;
                    for (unsigned int d=0; d<dim; d++)
                    {
                        const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);

                        // and store for each velocity component
                        for (size_type i=0; i<vsize; i++)
                        {
                            RF val = phi_v[i]*factor;
                            r[lfsv_v.localIndex(i)] += fval[d] * val;
                        }
                    }
                }
            }

            // boundary integral independent of ansatz functions,
            // Neumann and Dirichlet boundary conditions, DG penalty term's right hand side
            template<typename IG, typename LFSV, typename R>
            void lambda_boundary (const IG& ig, const LFSV& lfsv, R& r) const
            {
                // dimensions
                static const unsigned int dim = IG::Geometry::dimension;
                static const unsigned int dimw = IG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert
                    ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

                typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
                const LFSV_PFS_V& lfsv_pfs_v = lfsv.template getChild<VBLOCK>();

                dune_static_assert
                    ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

                // ... we assume all velocity components are the same
                typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
                const LFSV_V& lfsv_v = lfsv_pfs_v.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
                const LFSV_P& lfsv_p = lfsv.template getChild<PBLOCK>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
                typedef typename BasisSwitch_V::DomainField DF;
                typedef typename BasisSwitch_V::Range RT;
                typedef typename BasisSwitch_V::RangeField RF;
                typedef typename BasisSwitch_V::Range Range_V;
                typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
                typedef typename BasisSwitch_P::Range Range_P;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                const RF penalty_factor = ip_factor.getFaceIP(ig);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim-1> flocal = it->position();
                    Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(flocal);
                    Dune::FieldVector<DF,dimw> global = ig.geometry().global(flocal);
                    
                    // value of velocity shape functions
                    std::vector<RT> phi_v(vsize);
                    FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p(psize);
                    FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

                    std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
                    BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                          ig.inside()->geometry(), local, grad_phi_v);

                    const Dune::FieldVector<DF,dim> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());

                    // evaluate boundary condition type
                    typename B::Traits::RangeType bctype;
                    b.evaluate(ig,flocal,bctype);
                    
                    if (bctype == BC::VelocityDirichlet)
                    {
                        typename V::Traits::RangeType u0;
                        v.evaluateGlobal(global,u0);
                        
                        //================================================//
                        // \mu \int \nabla u_0 \cdot v \cdot n
                        //================================================//
                        RF factor = mu * weight;
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            const RF val = (grad_phi_v[i][0]*normal) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                r[lfsv_v.localIndex(i)] -= val * u0[d];
                            }
                        }
                        //================================================//
                        // \mu \int \sigma / |\gamma|^\beta v u_0
                        //================================================//
                        factor = penalty_factor * weight;
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            const RF val = phi_v[i] * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                r[lfsv_v.localIndex(i)] -= val * u0[d];
                            }
                        }
                        //================================================//
                        // \int q u_0 n
                        //================================================//
                        for (unsigned int i=0;i<psize;++i) // test
                        {
                            RF val = phi_p[i]*(u0 * normal) * weight;
                            r[lfsv_p.localIndex(i)] -= val;
                        }
                    }
                    if (bctype == BC::PressureDirichlet)
                    {
                        typename P::Traits::RangeType p0;
                        p.evaluateGlobal(global,p0);
                    
                        //std::cout << "Pdirichlet\n";
                        //================================================//
                        // \int p u n
                        //================================================//            
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                RF val = p0*normal[d]*phi_v[i] * weight;
                                r[lfsv_v.localIndex(i)] += val;
                            }
                        }
                    }
                }
            }

            // jacobian of volume term
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void jacobian_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv,
                LocalMatrix<R>& mat) const
            {
                // dimensions
                static const unsigned int dim = EG::Geometry::dimension;

                // subspaces
                dune_static_assert
                    ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
                const LFSV_PFS_V& lfsv_pfs_v = lfsv.template getChild<VBLOCK>();
                dune_static_assert
                    ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

                // ... we assume all velocity components are the same
                typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
                const LFSV_V& lfsv_v = lfsv_pfs_v.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
                const LFSV_P& lfsv_p = lfsv.template getChild<PBLOCK>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
                typedef typename BasisSwitch_V::DomainField DF;
                typedef typename BasisSwitch_V::Range RT;
                typedef typename BasisSwitch_V::RangeField RF;
                typedef typename BasisSwitch_V::Range Range_V;
                typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
                typedef typename BasisSwitch_P::Range Range_P;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gt = eg.geometry().type();
                const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
                
                // loop over quadrature points
                for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    const Dune::FieldVector<DF,dim> local = it->position();
                    
                    // and value of pressure shape functions
                    std::vector<RT> phi_p(psize);
                    FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

                    // compute gradients
                    std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
                    BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                            eg.geometry(), local, grad_phi_v);

                    const RF weight = it->weight() * eg.geometry().integrationElement(it->position());
                    
                    //================================================//
                    // \int (mu*grad_u*grad_v)
                    //================================================//
                    const RF factor = mu * weight;
                    for (size_type j=0; j<vsize; j++)
                    {
                        for (size_type i=0; i<vsize; i++)
                        {
                            // grad_phi_j*grad_phi_i
                            RF val = (grad_phi_v[j][0]*grad_phi_v[i][0])*factor;
                            // and store for each velocity component
                            for (unsigned int d=0; d<dim; d++)
                            {
                                const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                mat(lfsv_v.localIndex(i),lfsv_v.localIndex(j)) += val;
                            }
                        }
                    }

                    //================================================//
                    // - q * div u
                    // - p * div v
                    //================================================//            
                    for (size_type j=0; j<psize; j++) // test (q)
                    {
                        RF val = -1.0 * phi_p[j]*weight;
                        for (size_type i=0; i<vsize; i++) // ansatz (u)
                        {
                            for (unsigned int d=0; d<dim; d++)
                            {
                                const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                mat(lfsv_p.localIndex(j),lfsv_v.localIndex(i)) += val*grad_phi_v[i][0][d];
                                mat(lfsv_v.localIndex(i),lfsv_p.localIndex(j)) += val*grad_phi_v[i][0][d];
                            }
                        }
                    }
                }
            }

            // jacobian of skeleton term
            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void jacobian_skeleton (const IG& ig,
                const LFSU& lfsu_s, const X&, const LFSV& lfsv_s,
                const LFSU& lfsu_n, const X&, const LFSV& lfsv_n,
                LocalMatrix<R>& mat_ss, LocalMatrix<R>& mat_sn,
                LocalMatrix<R>& mat_ns, LocalMatrix<R>& mat_nn) const
            {
                // dimensions
                static const unsigned int dim = IG::Geometry::dimension;
                static const unsigned int dimw = IG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert
                    ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

                typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
                const LFSV_PFS_V& lfsv_s_pfs_v = lfsv_s.template getChild<VBLOCK>();
                const LFSV_PFS_V& lfsv_n_pfs_v = lfsv_n.template getChild<VBLOCK>();
                dune_static_assert
                    ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

                // ... we assume all velocity components are the same
                typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.template getChild<0>();
                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.template getChild<0>();
                const unsigned int vsize_s = lfsv_s_v.size();
                const unsigned int vsize_n = lfsv_n_v.size();
                typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
                const LFSV_P& lfsv_s_p = lfsv_s.template getChild<PBLOCK>();
                const LFSV_P& lfsv_n_p = lfsv_n.template getChild<PBLOCK>();
                const unsigned int psize_s = lfsv_s_p.size();
                const unsigned int psize_n = lfsv_n_p.size();

                // domain and range field type
                typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
                typedef typename BasisSwitch_V::DomainField DF;
                typedef typename BasisSwitch_V::Range RT;
                typedef typename BasisSwitch_V::RangeField RF;
                typedef typename BasisSwitch_V::Range Range_V;
                typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
                typedef typename BasisSwitch_P::Range Range_P;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                const RF penalty_factor = ip_factor.getFaceIP(ig);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
                    Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(it->position());

                    // value of velocity shape functions
                    std::vector<RT> phi_v_s(vsize_s);
                    std::vector<RT> phi_v_n(vsize_n);
                    FESwitch_V::basis(lfsv_s_v.finiteElement()).evaluateFunction(local_s,phi_v_s);
                    FESwitch_V::basis(lfsv_n_v.finiteElement()).evaluateFunction(local_n,phi_v_n);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p_s(psize_s);
                    std::vector<RT> phi_p_n(psize_n);
                    FESwitch_P::basis(lfsv_s_p.finiteElement()).evaluateFunction(local_s,phi_p_s);
                    FESwitch_P::basis(lfsv_n_p.finiteElement()).evaluateFunction(local_n,phi_p_n);

                    // compute gradients
                    std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_s(vsize_s);
                    BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_s_v.finiteElement()),
                                            ig.inside()->geometry(), local_s, grad_phi_v_s);

                    std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v_n(vsize_n);
                    BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_n_v.finiteElement()),
                                            ig.outside()->geometry(), local_n, grad_phi_v_n);

                    const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());
                    
                    //================================================//
                    // - (\mu \int < \nabla u > . normal . [v])  
                    //================================================//
                    assert(vsize_s == vsize_n);
                    RF factor = mu * weight;
                    for (unsigned int i=0;i<vsize_s;++i)
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            RF val = (0.5*(grad_phi_v_s[i][0]*normal)*phi_v_s[j]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                mat_ss(lfsv_s_v.localIndex(j),lfsv_s_v.localIndex(i)) -= val;
                                mat_ss(lfsv_s_v.localIndex(i),lfsv_s_v.localIndex(j)) += epsilon*val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            // the normal vector flipped, thus the sign flips
                            RF val = (-0.5*(grad_phi_v_s[i][0]*normal)*phi_v_n[j]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_ns(lfsv_n_v.localIndex(j),lfsv_s_v.localIndex(i)) -= val;
                                mat_sn(lfsv_s_v.localIndex(i),lfsv_n_v.localIndex(j)) += epsilon*val;
                            }
                        }
                    }
                    for (unsigned int i=0;i<vsize_n;++i) 
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            RF val = (0.5*(grad_phi_v_n[i][0]*normal)*phi_v_s[j]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_sn(lfsv_s_v.localIndex(j),lfsv_n_v.localIndex(i)) -= val;
                                mat_ns(lfsv_n_v.localIndex(i),lfsv_s_v.localIndex(j)) += epsilon*val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            // the normal vector flipped, thus the sign flips
                            RF val = (-0.5*(grad_phi_v_n[i][0]*normal)*phi_v_n[j]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_nn(lfsv_n_v.localIndex(j),lfsv_n_v.localIndex(i)) -= val;
                                mat_nn(lfsv_n_v.localIndex(i),lfsv_n_v.localIndex(j)) += epsilon*val;
                            }
                        }
                    }
                    //================================================//
                    // \mu \int \sigma / |\gamma|^\beta v u
                    //================================================//
                    factor = penalty_factor * weight;
                    for (unsigned int i=0;i<vsize_s;++i)
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            RF val = phi_v_s[i]*phi_v_s[j] * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                mat_ss(lfsv_s_v.localIndex(j),lfsv_s_v.localIndex(i)) += val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            RF val = phi_v_s[i]*phi_v_n[j] * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_ns(lfsv_n_v.localIndex(j),lfsv_s_v.localIndex(i)) -= val;
                            }
                        }
                    }
                    for (unsigned int i=0;i<vsize_n;++i) 
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            RF val = phi_v_n[i]*phi_v_s[j] * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_sn(lfsv_s_v.localIndex(j),lfsv_n_v.localIndex(i)) -= val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            RF val = phi_v_n[i]*phi_v_n[j] * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);
                                mat_nn(lfsv_n_v.localIndex(j),lfsv_n_v.localIndex(i)) += val;
                            }
                        }
                    }
                    //================================================//
                    // \int <q> [u] n
                    // \int <p> [v] n
                    //================================================//            
                    for (unsigned int i=0;i<vsize_s;++i) 
                    {
                        for (unsigned int j=0;j<psize_s;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                RF val = 0.5*(phi_p_s[j]*normal[d]*phi_v_s[i]) * weight;
                                mat_ss(lfsv_s_v.localIndex(i),lfsv_s_p.localIndex(j)) += val;
                                mat_ss(lfsv_s_p.localIndex(j),lfsv_s_v.localIndex(i)) += val;
                            }
                        }
                        for (unsigned int j=0;j<psize_n;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_s_v = lfsv_s_pfs_v.getChild(d);
                                RF val = 0.5*(phi_p_n[j]*normal[d]*phi_v_s[i]) * weight;
                                mat_sn(lfsv_s_v.localIndex(i),lfsv_n_p.localIndex(j)) += val;
                                mat_ns(lfsv_n_p.localIndex(j),lfsv_s_v.localIndex(i)) += val;
                            }
                        }
                    }
                    for (unsigned int i=0;i<vsize_n;++i) 
                    {
                        for (unsigned int j=0;j<psize_s;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);

                                // the normal vector flipped, thus the sign flips
                                RF val = -0.5*(phi_p_s[j]*normal[d]*phi_v_n[i]) * weight;
                                mat_ns(lfsv_n_v.localIndex(i),lfsv_s_p.localIndex(j)) += val;
                                mat_sn(lfsv_s_p.localIndex(j),lfsv_n_v.localIndex(i)) += val;
                            }
                        }
                        for (unsigned int j=0;j<psize_n;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                const LFSV_V& lfsv_n_v = lfsv_n_pfs_v.getChild(d);

                                // the normal vector flipped, thus the sign flips
                                RF val = -0.5*(phi_p_n[j]*normal[d]*phi_v_n[i]) * weight;
                                mat_nn(lfsv_n_v.localIndex(i),lfsv_n_p.localIndex(j)) += val;
                                mat_nn(lfsv_n_p.localIndex(j),lfsv_n_v.localIndex(i)) += val;
                            }
                        }
                    }
                }
            }

            // jacobian of boundary term
            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void jacobian_boundary (const IG& ig,
                const LFSU& lfsu, const X& x, const LFSV& lfsv,
                LocalMatrix<R>& mat) const
            {
                // dimensions
                static const unsigned int dim = IG::Geometry::dimension;
                static const unsigned int dimw = IG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert
                    ((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");

                typedef typename LFSV::template Child<VBLOCK>::Type LFSV_PFS_V;
                const LFSV_PFS_V& lfsv_pfs_v = lfsv.template getChild<VBLOCK>();
                dune_static_assert
                    ((LFSV_PFS_V::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");

                // ... we assume all velocity components are the same
                typedef typename LFSV_PFS_V::template Child<0>::Type LFSV_V;
                const LFSV_V& lfsv_v = lfsv_pfs_v.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<PBLOCK>::Type LFSV_P;
                const LFSV_P& lfsv_p = lfsv.template getChild<PBLOCK>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef FiniteElementInterfaceSwitch<typename LFSV_V::Traits::FiniteElementType > FESwitch_V;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_V::Basis > BasisSwitch_V;
                typedef typename BasisSwitch_V::DomainField DF;
                typedef typename BasisSwitch_V::Range RT;
                typedef typename BasisSwitch_V::RangeField RF;
                typedef typename BasisSwitch_V::Range Range_V;
                typedef FiniteElementInterfaceSwitch<typename LFSV_P::Traits::FiniteElementType > FESwitch_P;
                typedef PDELab::BasisInterfaceSwitch<typename FESwitch_P::Basis > BasisSwitch_P;
                typedef typename BasisSwitch_P::Range Range_P;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                // evaluate boundary condition type
                typename B::Traits::RangeType bctype;
                b.evaluate(ig,rule.begin()->position(),bctype);

                const RF penalty_factor = ip_factor.getFaceIP(ig);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
                    
                    // value of velocity shape functions
                    std::vector<RT> phi_v(vsize);
                    FESwitch_V::basis(lfsv_v.finiteElement()).evaluateFunction(local,phi_v);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p(psize);
                    FESwitch_P::basis(lfsv_p.finiteElement()).evaluateFunction(local,phi_p);

                    std::vector<Dune::FieldMatrix<RF,1,dim> > grad_phi_v(vsize);
                    BasisSwitch_V::gradient(FESwitch_V::basis(lfsv_v.finiteElement()),
                                          ig.inside()->geometry(), local, grad_phi_v);

                    const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());
                    
                    // velocity boundary condition
                    if (bctype == BC::VelocityDirichlet)
                    {
                        //================================================//
                        // - (\mu \int \nabla u. normal . v)  
                        //================================================//
                        const RF factor = - mu * weight;
                        for (unsigned int i=0;i<vsize;++i) // ansatz
                        {
                            for (unsigned int j=0;j<vsize;++j) // test
                            {
                                RF val = ((grad_phi_v[j][0]*normal)*phi_v[i]) * factor;
                                for (unsigned int d=0;d<dim;++d)
                                {
                                    const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                    mat(lfsv_v.localIndex(i),lfsv_v.localIndex(j)) += epsilon * val;
                                    mat(lfsv_v.localIndex(j),lfsv_v.localIndex(i)) += - val;
                                }
                            }
                        }
                        //================================================//
                        // \int q u n
                        // \int p v n
                        //================================================//
                        for (unsigned int i=0;i<vsize;++i) // ansatz
                        {
                            for (unsigned int j=0;j<psize;++j) // test
                            {
                                for (unsigned int d=0;d<dim;++d)
                                {
                                    const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                    RF val = (phi_p[j]*normal[d]*phi_v[i]) * weight;
                                    mat(lfsv_p.localIndex(j),lfsv_v.localIndex(i)) += val; // q u n
                                    mat(lfsv_v.localIndex(i),lfsv_p.localIndex(j)) += val; // p v n
                                }
                            }
                        }
                        //================================================//
                        // \mu \int \sigma / |\gamma|^\beta v u
                        //================================================//
                        const RF p_factor = penalty_factor * weight;
                        for (unsigned int i=0;i<vsize;++i)
                        {
                            for (unsigned int j=0;j<vsize;++j) 
                            {
                                RF val = phi_v[i]*phi_v[j] * p_factor;
                                for (unsigned int d=0;d<dim;++d)
                                {
                                    const LFSV_V& lfsv_v = lfsv_pfs_v.getChild(d);
                                    mat(lfsv_v.localIndex(j),lfsv_v.localIndex(i)) += val;
                                }
                            }
                        }

                    }
                }
            }

        private:
            const F& f;
            const B& b;
            const V& v;
            const P& p;
            // values for NIPG / NIPG
            int    epsilon;
            int    qorder;
            // physical parameters
            double mu;
            const IP & ip_factor;
        };

        //! \} group GridFunctionSpace
    } // namespace PDELab

  template<>
  struct FieldTraits<PDELab::StokesBoundaryCondition::Type>
  {
    typedef PDELab::StokesBoundaryCondition::Type field_type;
    typedef PDELab::StokesBoundaryCondition::Type real_type;
};

} // namespace Dune

#endif
