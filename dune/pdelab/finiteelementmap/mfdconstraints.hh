#ifndef DUNE_PDELAB_MFDCONSTRAINTS_HH
#define DUNE_PDELAB_MFDCONSTRAINTS_HH

#include<dune/grid/common/genericreferenceelements.hh>
#include<dune/common/geometrytype.hh>

namespace Dune
{
    namespace PDELab
    {
        class MimeticConstraints
        {
        public:
            enum{doVolume=false};
            enum{doSkeleton=false};
            enum{doBoundary=true};
            enum{doProcessor=false};

            template<typename B, typename I, typename LFS, typename T>
            void boundary(const B& b, const I& ig, const LFS& lfs, T& trafo) const
            {
                static const unsigned int dimIntersection = B::Traits::dimDomain;
                typedef typename B::Traits::DomainFieldType ctype;

                GeometryType gt = ig.intersection().type();
                typename B::Traits::DomainType center
                    = GenericReferenceElements<ctype,dimIntersection>::general(gt).position(0,0);
                typename B::Traits::RangeType bctype;
                b.evaluate(ig, center, bctype);
                if (bctype == 1)
                {
                    typename T::RowType empty;
                    trafo[ig.intersectionIndex()] = empty;
                }
            }
        };
    }
}

#endif
