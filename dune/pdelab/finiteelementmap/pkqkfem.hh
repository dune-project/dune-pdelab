// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_PKQKFEM_HH
#define DUNE_PDELAB_PKQKFEM_HH

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/common/array.hh>
#include <dune/common/shared_ptr.hh>
#include "finiteelementmap.hh"
#include "qkfem.hh"
#include "pkfem.hh"

namespace Dune {
    namespace PDELab {

        namespace {
            template<class D, class R, int d, int k>
            struct InitPkQkLocalFiniteElementMap
            {
                template<typename C>
                static void init(C & c, unsigned int order)
                {
                    if (order == k)
                    {
                        typedef Dune::QkLocalFiniteElement<D,R,d,k> QkLFE;
                        typedef Dune::PkLocalFiniteElement<D,R,d,k> PkLFE;
                        typedef typename C::value_type ptr;
                        c[0] = ptr(new LocalFiniteElementVirtualImp<QkLFE>(QkLFE()));
                        c[1] = ptr(new LocalFiniteElementVirtualImp<PkLFE>(PkLFE()));
                    }
                    else
                        InitPkQkLocalFiniteElementMap<D,R,d,k-1>::init(c,order);
                }
            };
            template<class D, class R, int d>
            struct InitPkQkLocalFiniteElementMap<D,R,d,-1>
            {
                template<typename C>
                    static void init(C & c, unsigned int order)
                {
                    DUNE_THROW(Exception, "Sorry, but we failed to initialize a QkPk FiniteElementMap of order " << order);
                }
            };
        }

        //! FiniteElementMap which provides PkQkLocalFiniteElement instances, depending on the geometry type
        //! \ingroup FiniteElementMap
        template<class D, class R, int d, int maxP=6>
        class PkQkLocalFiniteElementMap
        {
            //! Type of finite element from local functions
            typedef LocalFiniteElementVirtualInterface<Dune::LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d>,0> > FiniteElementType;
        public:
            typedef FiniteElementMapTraits<FiniteElementType> Traits;

            PkQkLocalFiniteElementMap ()
            {
                InitPkQkLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_,maxP);
            }

            PkQkLocalFiniteElementMap (unsigned int order)
            {
                InitPkQkLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_,order);
            }

            //! \brief get local basis functions for entity
            template<class EntityType>
            const typename Traits::FiniteElementType& find (const EntityType& e) const
            {
                typename Dune::GeometryType geoType=e.type();
                return getFEM(geoType);
            }

            //! \brief get local basis functions for a given geometrytype
            const typename Traits::FiniteElementType& getFEM (Dune::GeometryType gt) const
            {
                if (gt.isCube())
                {
                    return *(finiteElements_[0]);
                }
                if (gt.isSimplex())
                {
                    return *(finiteElements_[1]);
                }
                DUNE_THROW(Exception, "We can only handle cubes and simplices");
            }

            bool fixedSize() const
            {
                return false;
            }

            std::size_t size(GeometryType gt) const
            {
                assert(false && "this method should never be called");
            }

            std::size_t maxLocalSize() const
            {
                return (1<<d);
            }

        private:
            Dune::array< Dune::shared_ptr<FiniteElementType>, 2 > finiteElements_;
        };
    }
}

#endif //DUNE_PDELAB_PKQKFEM_HH
