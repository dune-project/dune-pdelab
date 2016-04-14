#ifndef DUNE_PDELAB_LOCALOPERATOR_DGPARAMETER_HH
#define DUNE_PDELAB_LOCALOPERATOR_DGPARAMETER_HH

#warning This file is deprecated and will be removed. Use dune/pdelab/localoperator/dginteriorpenaltyparameter.hh instead!

#include <dune/common/parametertreeparser.hh>

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
        class DUNE_DEPRECATED_MSG("DefaultInteriorPenalty is deprecated. Please use the implementation from dune/pdelab/localoperator/dginteriorpenaltyparameter.hh instead!") DefaultInteriorPenalty
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

    } // end namespace PDELab
} // end namespace Dune

#endif // DUNE_PDELAB_LOCALOPERATOR_DGPARAMETER_HH
