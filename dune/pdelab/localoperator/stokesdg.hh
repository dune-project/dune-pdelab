// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_STOKESDG_HH
#define DUNE_PDELAB_STOKESDG_HH

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/static_assert.hh>
#include <dune/grid/common/quadraturerules.hh>
#include <dune/grid/common/referenceelements.hh>

#include "../common/geometrywrapper.hh"
#include "../gridoperatorspace/gridoperatorspace.hh"
#include "../gridoperatorspace/gridoperatorspaceutilities.hh"
#include "pattern.hh"
#include "flags.hh"

namespace Dune {
    namespace PDELab {

        struct StokesBoundaryCondition {
            enum Type {
                DoNothing = 0,
                VelocityDirichlet = 1,
                PressureDirichlet = 2
            };
        };
        
        // a local operator for solving the stokes equation
        template<typename F, typename B, typename V, typename P>
        class StokesDG :
            public LocalOperatorDefaultFlags,
            public FullSkeletonPattern, public FullVolumePattern
            //
            ,public JacobianBasedAlphaVolume< StokesDG<F,B,V,P> >
            ,public JacobianBasedAlphaSkeleton< StokesDG<F,B,V,P> >
            ,public JacobianBasedAlphaBoundary< StokesDG<F,B,V,P> >
        {
            typedef StokesBoundaryCondition BC;
        public:
            // pattern assembly flags
            enum { doPatternVolume = true };
            enum { doPatternSkeleton = true };

            //
            enum { doSkeletonTwoSided = false };

            // residual assembly flags
            enum { doAlphaVolume    = true };
            enum { doAlphaSkeleton  = true };
            enum { doAlphaBoundary  = true };
            enum { doLambdaVolume   = true };
            enum { doAlphaVolumePostSkeleton = true };
            enum { doLambdaBoundary = true };

            StokesDG (const std::string & method,
                const F & _f, const B & _b, const V & _v, const P & _p) :
                f(_f), b(_b), v(_v), p(_p), qorder(4), mu(1)
            {
                std::string s = method;
                std::transform(s.begin(), s.end(), s.begin(), tolower);

                // nipg (epsilon=1) 2d p1 -> Klaus sagt sollte auch sigma 1 klappen
                if (s.find("nipg") != std::string::npos)
                {
                    epsilon = 1;
                    beta = 1;
                    if (sscanf(s.c_str(), "nipg %lg", &sigma) != 1)
                        sigma = 3.9;
                    return;
                }
                // sipg (epsilon=-1) 2d p1 -> Klaus sagt sigma=3.9irgendwas
                if (s.find("sipg") != std::string::npos)
                {
                    epsilon = -1;
                    beta = 1;
                    if (sscanf(s.c_str(), "sipg %lg", &sigma) != 1)
                        sigma = 3.9;
                    return;
                }
                // obb sigma = 0, epsilon = 
                if (s == "obb")
                {
                    epsilon = 1;
                    beta = 1;
                    sigma = 0;
                    return;
                }
                // extract parameters
                {
                    if (3 == sscanf(s.c_str(), "%d %lg %lg", &epsilon, &sigma, &beta))
                        return;
                }
                DUNE_THROW(Dune::Exception, "Unknown DG type " << method);
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
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_vel = lfsv.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_v = lfsv_vel.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();

                // domain and range field type
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RT;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::JacobianType JacobianType_v;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_v;
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
                    lfsv_v.localFiniteElement().localBasis().evaluateFunction(local,phi_v);

                    const RF weight = it->weight() * eg.geometry().integrationElement(it->position());

                    // evaluate source term
                    typename F::Traits::RangeType fval;
                    f.evaluateGlobal(global,fval);
                    
                    //================================================//
                    // \int (f*v)
                    //================================================//
                    const RF factor = mu * weight;
                    for (size_type i=0; i<vsize; i++)
                    {
                        // f*phi_i
                        RF val = phi_v[i]*factor;
                        // and store for each velocity component
                        for (unsigned int d=0; d<dim; d++)
                        {
                            r[i+d*vsize] += fval[d] * val;
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
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_vel = lfsv.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_v = lfsv_vel.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<1>::Type LFSV_p;
                const LFSV_p& lfsv_p = lfsv.template getChild<1>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RT;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::JacobianType JacobianType_v;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_v;
                typedef typename LFSV_p::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_p;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim-1> flocal = it->position();
                    Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(flocal);
                    Dune::FieldVector<DF,dimw> global = ig.geometry().global(flocal);
                    
                    // evaluate gradient of velocity shape functions (we assume Galerkin method lfsu=lfsv)
                    std::vector<JacobianType_v> jac_v_s(vsize);
                    lfsv_v.localFiniteElement().localBasis().evaluateJacobian(local,jac_v_s);
                    // value of velocity shape functions
                    std::vector<RT> phi_v(psize);
                    lfsv_v.localFiniteElement().localBasis().evaluateFunction(local,phi_v);

                    // transform gradient to real element
                    const Dune::FieldMatrix<DF,dimw,dim> jInvT =
                        ig.inside()->geometry().jacobianInverseTransposed(local);
                    std::vector<Dune::FieldVector<RF,dim> > grad_phi_v(vsize);
                    for (typename LFSV::Traits::SizeType i=0; i<vsize; i++)
                    {
                        grad_phi_v[i] = 0.0;
                        jInvT.umv(jac_v_s[i][0],grad_phi_v[i]);
                    }

                    const Dune::FieldVector<DF,dim> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());

                    // evaluate boundary condition type
                    typename B::Traits::RangeType bctype;
                    b.evaluate(ig,flocal,bctype);
                    
                    // get bc value
                    typename V::Traits::RangeType u0;
                    v.evaluateGlobal(global,u0);
                    typename P::Traits::RangeType p0;
                    p.evaluateGlobal(global,p0);
                    
                    if (bctype == BC::VelocityDirichlet)
                    {
                        //================================================//
                        // TERM: 4
                        // \mu \int \nabla u_0 \cdot v \cdot n
                        //================================================//
                        const RF factor = mu * weight;
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            const RF val = (grad_phi_v[i]*normal) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                r[i+d*vsize] += val * u0[d];
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
                static const unsigned int dimw = EG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_vel = lfsv.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_v = lfsv_vel.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<1>::Type LFSV_p;
                const LFSV_p& lfsv_p = lfsv.template getChild<1>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RT;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::JacobianType JacobianType_v;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_v;
                typedef typename LFSV_p::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_p;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gt = eg.geometry().type();
                const Dune::QuadratureRule<DF,dim>& rule = Dune::QuadratureRules<DF,dim>::rule(gt,qorder);
                
                // loop over quadrature points
                for (typename Dune::QuadratureRule<DF,dim>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    const Dune::FieldVector<DF,dim> local = it->position();
                    
                    // evaluate gradient of velocity shape functions (we assume Galerkin method lfsu=lfsv)
                    std::vector<JacobianType_v> jac_v_s(vsize);
                    lfsv_v.localFiniteElement().localBasis().evaluateJacobian(local,jac_v_s);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p(psize);
                    lfsv_p.localFiniteElement().localBasis().evaluateFunction(local,phi_p);

                    // transform gradient to real element
                    const Dune::FieldMatrix<DF,dimw,dim> jInvT =
                        eg.geometry().jacobianInverseTransposed(it->position());
                    std::vector<Dune::FieldVector<RF,dim> > grad_phi_v(vsize);
                    for (typename LFSV::Traits::SizeType i=0; i<vsize; i++)
                    {
                        grad_phi_v[i] = 0.0;
                        jInvT.umv(jac_v_s[i][0],grad_phi_v[i]);
                    }

                    const RF weight = it->weight() * eg.geometry().integrationElement(it->position());
                    
                    //================================================//
                    // TERM: 1
                    // \int (mu*grad_u*grad_v)
                    //================================================//
                    const RF factor = mu * weight;
                    for (size_type j=0; j<vsize; j++)
                    {
                        for (size_type i=0; i<vsize; i++)
                        {
                            // grad_phi_j*grad_phi_i
                            RF val = (grad_phi_v[j]*grad_phi_v[i])*factor;
                            // and store for each velocity component
                            for (unsigned int d=0; d<dim; d++)
                            {
                                mat(i+d*vsize,j+d*vsize) += val;
                            }
                        }
                    }

                    //================================================//
                    // TERM: 8, 11
                    // - p * div v
                    // - q * div u
                    //================================================//            
                    for (size_type j=0; j<psize; j++)
                    {
                        RF val = -1.0 * phi_p[j]*weight;
                        for (size_type i=0; i<vsize; i++)
                        {
                            for (unsigned int d=0; d<dim; d++)
                            {
                                mat(j+dim*vsize,i+d*vsize) += val*grad_phi_v[i][d];
                                mat(i+d*vsize,j+dim*vsize) += val*grad_phi_v[i][d];
                            }
                        }
                    }
                }
            }

            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                            const LFSV& lfsv, R& r) const
            {}
            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void jacobian_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                               const LFSV& lfsv, LocalMatrix<R>& mat) const
            {
                // dimensions
                static const unsigned int dim = EG::Geometry::dimension;

                // subspaces
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_vel = lfsv.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_v = lfsv_vel.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<1>::Type LFSV_p;
                const LFSV_p& lfsv_p = lfsv.template getChild<1>();
                const unsigned int psize = lfsv_p.size();

                // fix one pressure DOF
                static int cnt = 0;
                if (cnt == 0)
                {
                    for (unsigned int i=0; i<dim*vsize+psize; i++)
                        mat(dim*vsize, i) = 0;
                    mat(dim*vsize, dim*vsize) = 1;
                }
                cnt++;
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
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_s_vel = lfsv_s.template getChild<0>();
                const LFSV_vel& lfsv_n_vel = lfsv_n.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_s_v = lfsv_s_vel.template getChild<0>();
                const LFSV_v& lfsv_n_v = lfsv_n_vel.template getChild<0>();
                const unsigned int vsize_s = lfsv_s_v.size();
                const unsigned int vsize_n = lfsv_n_v.size();
                typedef typename LFSV::template Child<1>::Type LFSV_p;
                const LFSV_p& lfsv_s_p = lfsv_s.template getChild<1>();
                const LFSV_p& lfsv_n_p = lfsv_n.template getChild<1>();
                const unsigned int psize_s = lfsv_s_p.size();
                const unsigned int psize_n = lfsv_n_p.size();

                // domain and range field type
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RT;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::JacobianType JacobianType_v;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_v;
                typedef typename LFSV_p::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_p;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim> local_s = ig.geometryInInside().global(it->position());
                    Dune::FieldVector<DF,dim> local_n = ig.geometryInOutside().global(it->position());

                    std::cout << ig.geometry().global(it->position()) << " / " << local_s << " / " << local_n << "\n";
                    
                    // evaluate gradient of velocity shape functions (we assume Galerkin method lfsu=lfsv)
                    std::vector<JacobianType_v> jac_v_s_s(vsize_s);
                    std::vector<JacobianType_v> jac_v_n_s(vsize_n);
                    lfsv_s_v.localFiniteElement().localBasis().evaluateJacobian(local_s,jac_v_s_s);
                    lfsv_n_v.localFiniteElement().localBasis().evaluateJacobian(local_n,jac_v_n_s);
                    // value of velocity shape functions
                    std::vector<RT> phi_v_s(psize_s);
                    std::vector<RT> phi_v_n(psize_n);
                    lfsv_s_v.localFiniteElement().localBasis().evaluateFunction(local_s,phi_v_s);
                    lfsv_n_v.localFiniteElement().localBasis().evaluateFunction(local_n,phi_v_n);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p_s(psize_s);
                    std::vector<RT> phi_p_n(psize_n);
                    lfsv_s_p.localFiniteElement().localBasis().evaluateFunction(local_s,phi_p_s);
                    lfsv_n_p.localFiniteElement().localBasis().evaluateFunction(local_n,phi_p_n);

                    // transform gradient to real element
                    const Dune::FieldMatrix<DF,dimw,dim> jInvT_s =
                        ig.inside()->geometry().jacobianInverseTransposed(local_s);
                    const Dune::FieldMatrix<DF,dimw,dim> jInvT_n =
                        ig.outside()->geometry().jacobianInverseTransposed(local_n);
                    std::vector<Dune::FieldVector<RF,dim> > grad_phi_v_s(vsize_s);
                    std::vector<Dune::FieldVector<RF,dim> > grad_phi_v_n(vsize_n);
                    for (typename LFSV::Traits::SizeType i=0; i<vsize_s; i++)
                    {
                        grad_phi_v_s[i] = 0.0;
                        jInvT_s.umv(jac_v_s_s[i][0],grad_phi_v_s[i]);
                    }
                    for (typename LFSV::Traits::SizeType i=0; i<vsize_n; i++)
                    {
                        grad_phi_v_n[i] = 0.0;
                        jInvT_n.umv(jac_v_n_s[i][0],grad_phi_v_n[i]);
                    }

                    const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());
                    
                    //================================================//
                    // TERM: 4
                    // - (\mu \int < \nabla u > . normal . [v])  
                    //================================================//
                    assert(vsize_s == vsize_n);
                    const RF factor = mu * weight;
                    for (unsigned int i=0;i<vsize_s;++i) 
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            RF val = (0.5*(grad_phi_v_s[j]*normal)*phi_v_s[i]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                mat_ss(i+d*vsize_s,j+d*vsize_s) -= val;
                                mat_ss(j+d*vsize_s,i+d*vsize_s) += epsilon*val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            RF val = (0.5*(grad_phi_v_n[j]*normal)*phi_v_s[i]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                mat_sn(i+d*vsize_s,j+d*vsize_n) -= val;
                                mat_ns(j+d*vsize_n,i+d*vsize_s) += epsilon*val;
                            }
                        }
                    }
                    for (unsigned int i=0;i<vsize_n;++i) 
                    {
                        for (unsigned int j=0;j<vsize_s;++j) 
                        {
                            // the normal vector flipped, thus the sign flips
                            RF val = (-0.5*(grad_phi_v_s[j]*normal)*phi_v_n[i]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                mat_ns(i+d*vsize_n,j+d*vsize_s) -= val;
                                mat_sn(j+d*vsize_s,i+d*vsize_n) += epsilon*val;
                            }
                        }
                        for (unsigned int j=0;j<vsize_n;++j) 
                        {
                            // the normal vector flipped, thus the sign flips
                            RF val = (-0.5*(grad_phi_v_n[j]*normal)*phi_v_n[i]) * factor;
                            for (unsigned int d=0;d<dim;++d)
                            {
                                mat_nn(i+d*vsize_n,j+d*vsize_n) -= val;
                                mat_nn(j+d*vsize_n,i+d*vsize_n) += epsilon*val;
                            }
                        }
                    }
                    //================================================//
                    // TERM: 10
                    // \int <q> [u] n
                    //================================================//            
                    for (unsigned int i=0;i<vsize_s;++i) 
                    {
                        for (unsigned int j=0;j<psize_s;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                RF val = 0.5*(phi_p_s[j]*normal[d]*phi_v_s[i]) * weight;
                                mat_ss(i+d*vsize_s,j+dim*vsize_s) += val;
                                mat_ss(j+dim*vsize_s,i+d*vsize_s) += val;
                            }
                        }
                        for (unsigned int j=0;j<psize_n;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                RF val = 0.5*(phi_p_n[j]*normal[d]*phi_v_s[i]) * weight;
                                mat_sn(i+d*vsize_s,j+dim*vsize_n) += val;
                                mat_ns(j+dim*vsize_n,i+d*vsize_s) += val;
                            }
                        }
                    }
                    for (unsigned int i=0;i<vsize_n;++i) 
                    {
                        for (unsigned int j=0;j<psize_s;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                // the normal vector flipped, thus the sign flips
                                RF val = -0.5*(phi_p_s[j]*normal[d]*phi_v_n[i]) * weight;
                                mat_ns(i+d*vsize_n,j+dim*vsize_s) += val;
                                mat_sn(j+dim*vsize_s,i+d*vsize_n) += val;
                            }
                        }
                        for (unsigned int j=0;j<psize_n;++j) 
                        {
                            for (unsigned int d=0;d<dim;++d)
                            {
                                // the normal vector flipped, thus the sign flips
                                RF val = -0.5*(phi_p_n[j]*normal[d]*phi_v_n[i]) * weight;
                                mat_nn(i+d*vsize_n,j+dim*vsize_n) += val;
                                mat_nn(j+dim*vsize_n,i+d*vsize_n) += val;
                            }
                        }
                    }
                }
            }

            // jacobian of volume term
            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void jacobian_boundary (const IG& ig,
                const LFSU& lfsu, const X& x, const LFSV& lfsv,
                LocalMatrix<R>& mat) const
            {
                // dimensions
                static const unsigned int dim = IG::Geometry::dimension;
                static const unsigned int dimw = IG::Geometry::dimensionworld;

                // subspaces
                dune_static_assert((LFSV::CHILDREN == 2), "You seem to use the wrong function space for StokesDG");
                typedef typename LFSV::template Child<0>::Type LFSV_vel;
                const LFSV_vel& lfsv_vel = lfsv.template getChild<0>();
                dune_static_assert((LFSV_vel::CHILDREN == dim), "You seem to use the wrong function space for StokesDG");
                // ... we assume all velocity components are the same
                typedef typename LFSV_vel::template Child<0>::Type LFSV_v;
                const LFSV_v& lfsv_v = lfsv_vel.template getChild<0>();
                const unsigned int vsize = lfsv_v.size();
                typedef typename LFSV::template Child<1>::Type LFSV_p;
                const LFSV_p& lfsv_p = lfsv.template getChild<1>();
                const unsigned int psize = lfsv_p.size();

                // domain and range field type
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::DomainFieldType DF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RT;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeFieldType RF;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::JacobianType JacobianType_v;
                typedef typename LFSV_v::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_v;
                typedef typename LFSV_p::Traits::LocalFiniteElementType::
                    Traits::LocalBasisType::Traits::RangeType RangeType_p;
                typedef typename LFSV::Traits::SizeType size_type;

                // select quadrature rule
                Dune::GeometryType gtface = ig.geometryInInside().type();
                const Dune::QuadratureRule<DF,dim-1>& rule = Dune::QuadratureRules<DF,dim-1>::rule(gtface,qorder);

                // evaluate boundary condition type
                typename B::Traits::RangeType bctype;
                b.evaluate(ig,rule.begin()->position(),bctype);

                // loop over quadrature points and integrate normal flux
                for (typename Dune::QuadratureRule<DF,dim-1>::const_iterator it=rule.begin(); it!=rule.end(); ++it)
                {
                    // position of quadrature point in local coordinates of element
                    Dune::FieldVector<DF,dim> local = ig.geometryInInside().global(it->position());
                    
                    // evaluate gradient of velocity shape functions (we assume Galerkin method lfsu=lfsv)
                    std::vector<JacobianType_v> jac_v_s(vsize);
                    lfsv_v.localFiniteElement().localBasis().evaluateJacobian(local,jac_v_s);
                    // value of velocity shape functions
                    std::vector<RT> phi_v(psize);
                    lfsv_v.localFiniteElement().localBasis().evaluateFunction(local,phi_v);
                    // and value of pressure shape functions
                    std::vector<RT> phi_p(psize);
                    lfsv_p.localFiniteElement().localBasis().evaluateFunction(local,phi_p);

                    // transform gradient to real element
                    const Dune::FieldMatrix<DF,dimw,dim> jInvT =
                        ig.inside()->geometry().jacobianInverseTransposed(local);
                    std::vector<Dune::FieldVector<RF,dim> > grad_phi_v(vsize);
                    for (typename LFSV::Traits::SizeType i=0; i<vsize; i++)
                    {
                        grad_phi_v[i] = 0.0;
                        jInvT.umv(jac_v_s[i][0],grad_phi_v[i]);
                    }

                    const Dune::FieldVector<DF,dimw> normal = ig.unitOuterNormal(it->position());
                    const RF weight = it->weight()*ig.geometry().integrationElement(it->position());
                    
                    // velocity boundary condition
                    if (bctype == BC::VelocityDirichlet) {
                        //================================================//
                        // TERM: 4
                        // - (\mu \int \nabla u. normal . v)  
                        //================================================//
                        const RF factor = - mu * weight;
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            for (unsigned int j=0;j<vsize;++j) 
                            {
                                RF val = ((grad_phi_v[j]*normal)*phi_v[i]) * factor;
                                for (unsigned int d=0;d<dim;++d)
                                {
                                    mat(i+d*vsize,j+d*vsize) += epsilon * val;
                                    mat(j+d*vsize,i+d*vsize) += - val;
                                }
                            }
                        }
                    }
                    if (bctype == BC::VelocityDirichlet || bctype == BC::PressureDirichlet) {
                        //================================================//
                        // TERM: 10
                        // \int p u n
                        //================================================//            
                        for (unsigned int i=0;i<vsize;++i) 
                        {
                            for (unsigned int j=0;j<psize;++j) 
                            {
                                for (unsigned int d=0;d<dim;++d)
                                {
                                    RF val = (phi_p[j]*normal[d]*phi_v[i]) * weight;
                                    mat(i+d*vsize,j+dim*vsize) += val;
//                                    mat(j+dim*vsize,i+d*vsize) += val;
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
            double sigma;
            double beta;
            int    qorder;
            // physical parameters
            double mu;
        };

        //! \} group GridFunctionSpace
    } // namespace PDELab
} // namespace Dune

#endif
