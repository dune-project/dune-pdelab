#ifndef DUNE_PDELAB_LOCALOPERATOR_STOKESDGPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_STOKESDGPARAMETER_HH

#include <dune/common/parametertreeparser.hh>

#include <dune/pdelab/common/geometrywrapper.hh>

#include "stokesparameter.hh"

namespace Dune {
    namespace PDELab {

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

        /** \brief Traits class for parameter classes based on the concept
         * StokesDGParameters
         */
        template<typename GV, typename RF>
        struct StokesDGParameterTraits
        {
            //! \brief the grid view
            typedef GV GridViewType;

            //! \brief Enum for domain dimension
            enum { 
                //! \brief dimension of the domain
                dimDomain = GV::dimension
            }; 

            //! \brief Export type for domain field
            typedef typename GV::Grid::ctype DomainFieldType;

            //! \brief domain type
            typedef Dune::FieldVector<DomainFieldType,dimDomain> DomainType;

            //! \brief domain type
            typedef Dune::FieldVector<DomainFieldType,dimDomain-1> IntersectionDomainType;

            //! \brief Export type for range field
            typedef RF RangeFieldType;

            //! \brief Export dimension of range
            enum { dimRange = GV::dimensionworld }; 

            //! \brief range type
            typedef Dune::FieldVector<RF,GV::dimensionworld> RangeType;

            //! Boundary condition indicator type
            typedef StokesBoundaryCondition::Type BoundaryConditionType;

            //! grid types
            //! @{
            typedef typename GV::Traits::template Codim<0>::Entity ElementType;
            typedef typename GV::Intersection IntersectionType;
            //! @}
        };


        /** \brief Parameter class for local operator StokesDG
            
            \tparam F Momentum source term function
            \tparam B Boundary condition function
            \tparam V Dirichlet velocity boundary condition function
            \tparam P Dirichlet pressure boundary condition function
            \tparam IP A class providing the interior penalty factor for each face
         */
        template<typename GV, typename F, typename B, typename V, typename P, 
                 typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
        class StokesDGParameters
        {
        private:

            void initFromString(const std::string & method)
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


        public:
            //! Common range field type
            typedef typename V::Traits::RangeFieldType RF;
            //! Traits class
            typedef StokesDGParameterTraits<GV,RF> Traits;


            StokesDGParameters(const std::string & method, const RF mu_,
                               F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
                : ff(ff_), bf(bf_), vf(vf_), pf(pf_), ip(ip_), viscosity(mu_)
            {
                initFromString(method);
            }

            StokesDGParameters(const Dune::ParameterTree & configuration, const RF mu_,
                               F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
                : ff(ff_), bf(bf_), vf(vf_), pf(pf_), ip(ip_), viscosity(mu_)
            {
                epsilon = configuration.get<int>("epsilon");
            }

            //! source term
            template<typename ElementGeometry>
            typename Traits::RangeType 
            f (const ElementGeometry& e, const typename Traits::DomainType& x) const
            {
                typename Traits::RangeType fvalue;
                ff.evaluate(e.entity(),x,fvalue);
                return fvalue;
            }

            //! boundary condition type from local intersection coordinate
            template<typename IntersectionGeometry>
            typename Traits::BoundaryConditionType
            bctype (const IntersectionGeometry& is, 
                    const typename Traits::IntersectionDomainType& x) const
            {
                typename B::Traits::RangeType ry;
                bf.evaluate(is,x,ry);
                typename Traits::BoundaryConditionType y(ry);
                return y;
            }

            //! Dynamic viscosity value from local cell coordinate
            template<typename ElementGeometry>
            typename Traits::RangeFieldType 
            mu (const ElementGeometry& e, const typename Traits::DomainType& x) const
            {
                return viscosity;
            }

            //! Dynamic viscosity value from local intersection coordinate
            template<typename IntersectionGeometry>
            typename Traits::RangeFieldType 
            mu (const IntersectionGeometry& ig, const typename Traits::IntersectionDomainType& x) const
            {
                return viscosity;
            }

            //! Dirichlet boundary condition value from local cell coordinate
            template<typename ElementGeometry>
            typename Traits::RangeType 
            g (const ElementGeometry& e, const typename Traits::DomainType& x) const
            {
                typename Traits::RangeType y;
                vf.evaluate(e.entity(),x,y);
                return y;
            }

            //! Dirichlet boundary condition value from local intersection coordinate
            template<typename IntersectionGeometry>
            typename Traits::RangeType 
            g (const IntersectionGeometry& ig, const typename Traits::IntersectionDomainType& x) const
            {
                typename IntersectionGeometry::EntityPointer ep = ig.inside();
                typename Traits::DomainType ex = ep->geometry().local(ig.geometry().global(x));
                const Dune::PDELab::ElementGeometry<typename IntersectionGeometry::EntityPointer::Entity>
                    eg(*ep);
                typename Traits::RangeType y(g(eg,ex));
                return y;
            }

            //! Neumann boundary condition (Pressure condition)
            template<typename IntersectionGeometry>
            typename Traits::RangeFieldType 
            j (const IntersectionGeometry& ig, const typename Traits::IntersectionDomainType& x) const
            {
                typename IntersectionGeometry::EntityPointer ep = ig.inside();
                typename Traits::DomainType ex = ep->geometry().local(ig.geometry().global(x));
                typename P::Traits::RangeType y;
                pf.evaluate(*(ep),ex,y);
                return y;
            }

            //! Rescaling factor for the incompressibility equation
            typename Traits::RangeFieldType 
            incompressibilityScaling ( typename Traits::RangeFieldType  dt ) const
            {
              typename P::Traits::RangeType y(1.0 / dt);
              return y;
            }

            /** \brief Interior penalty parameter.

                \return The coefficient of the interior penalty term.
            */
            template<typename I>
            typename Traits::RangeFieldType       
            getFaceIP(const I & ig) const
            {
                return ip.getFaceIP(ig);
            }

            //! Return the symmetry factor epsilon for this IP
            //! discretization
            int
            epsilonIPSymmetryFactor()
            {return epsilon;}


            //! Set current simulation time
            template<typename TT>
            void 
            setTime(TT time){
                ff.setTime(time);
                bf.setTime(time);
                vf.setTime(time);
                pf.setTime(time);
            }

        private:

            F & ff;               // Momentum source term function
            B & bf;               // Boundary condition type function
            V & vf;               // Velocity dirichlet function
            P & pf;               // Pressure dirichlet
            IP & ip;              // Interior penalty
            const RF viscosity;   // Dynamic viscosity
            int epsilon;          // IP symmetry factor
        };


        /** \brief Parameter class for local operator NavierStokesDG
            
            \tparam F Momentum source term function
            \tparam B Boundary condition function
            \tparam V Dirichlet velocity boundary condition function
            \tparam P Dirichlet pressure boundary condition function
            \tparam IP A class providing the interior penalty factor for each face
        */
        template<typename GV, typename F, typename B, typename V, typename P, 
                 typename IP = DefaultInteriorPenalty<typename V::Traits::RangeFieldType> >
        class NavierStokesDGParameters : public StokesDGParameters<GV,F,B,V,P,IP>
        {

        public:
            //! Type of base class
            typedef StokesDGParameters<GV,F,B,V,P,IP> Base;
            //! Common range field type
            typedef typename V::Traits::RangeFieldType RF;
            //! Traits class
            typedef StokesDGParameterTraits<GV,RF> Traits;

            NavierStokesDGParameters(const std::string & method, const RF rho_, const RF mu_,
                                     F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
                : Base(method,mu_,ff_,bf_,vf_,pf_,ip_), density(rho_)
            {}

            NavierStokesDGParameters(const Dune::ParameterTree configuration, const RF rho_, const RF mu_,
                               F & ff_, B & bf_, V & vf_, P & pf_, IP & ip_)
                : Base(configuration,mu_,ff_,bf_,vf_,pf_,ip_), density(rho_)
            {}

            //! Dynamic viscosity value from local cell coordinate
            template<typename ElementGeometry>
            typename Traits::RangeFieldType 
            rho (const ElementGeometry& eg, const typename Traits::DomainType& x) const
            {
                return density;
            }

            //! Dynamic viscosity value from local intersection coordinate
            template<typename IntersectionGeometry>
            typename Traits::RangeFieldType 
            rho (const IntersectionGeometry& ig, const typename Traits::IntersectionDomainType& x) const
            {
                return density;
            }

        private:
            RF density;         // Fluid density
        };
    }
}

#endif
