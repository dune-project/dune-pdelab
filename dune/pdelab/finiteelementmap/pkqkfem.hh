// -*- tab-width: 4; indent-tabs-mode: nil -*-
#ifndef DUNE_PDELAB_FINITEELEMENTMAP_PKQKFEM_HH
#define DUNE_PDELAB_FINITEELEMENTMAP_PKQKFEM_HH

#include <memory>

#include <dune/geometry/type.hh>

#include <dune/localfunctions/common/virtualwrappers.hh>
#include <dune/localfunctions/common/virtualinterface.hh>
#include <dune/common/array.hh>
#include "finiteelementmap.hh"
#include "qkfem.hh"
#include "pkfem.hh"

namespace Dune {
    namespace PDELab {

        namespace {

            /** \brief Construct LocalFiniteElement with a given run-time order
             *
             * \tparam k Loop variable in a static loop
             */
            template<class D, class R, int d, int k>
            struct InitPkQkLocalFiniteElementMap
            {
                /** \brief Set up LocalFiniteElement with order 'order' */
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

        /** \brief FiniteElementMap which provides PkQkLocalFiniteElement instances, depending on the geometry type
         * \ingroup FiniteElementMap
         *
         * \tparam D Type used for coordinates
         * \tparam R Type used for shape function values
         * \tparam d Grid dimension
         * \tparam maxP Approximation order: if you construct an object of this class with its default constructor,
         *    then this number is the approximation order that you get.  If you construct an object giving an order
         *    at run-time, then maxP is the maximal order that you can request.
         */
        template<class D, class R, int d, int maxP=6>
        class PkQkLocalFiniteElementMap
        {
            //! Type of finite element from local functions
            typedef LocalFiniteElementVirtualInterface<Dune::LocalBasisTraits<D,d,Dune::FieldVector<D,d>,R,1,Dune::FieldVector<R,1>,Dune::FieldMatrix<R,1,d> > > FiniteElementType;
        public:
            typedef FiniteElementMapTraits<FiniteElementType> Traits;

            //! The dimension of the finite elements returned by this map.
            static constexpr int dimension = d;

            /** \brief Default constructor.  Constructs a space of order maxP */
            PkQkLocalFiniteElementMap ()
              : order_(maxP)
            {
                InitPkQkLocalFiniteElementMap<D,R,d,maxP>::init(finiteElements_,maxP);
            }

            /** \brief Construct a space with a given order
             * \throw Dune::Exception if the requested order is larger than maxP
             */
            PkQkLocalFiniteElementMap (unsigned int order)
              : order_(order)
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

            static constexpr bool fixedSize()
            {
                return false;
            }

            bool hasDOFs(int codim) const
            {
              switch (codim)
                {
                case 0:
                  return order_ == 0 || order_ > 1;
                case d:
                  return order_ > 0;
                default:
                  return order_ > 1;
                }
            }

            std::size_t size(GeometryType gt) const
            {
                DUNE_THROW(NotImplemented, "PkQkLocalFiniteElement is not fixed-size!");
            }

            static constexpr std::size_t maxLocalSize()
            {
                return (1<<d);
            }

        private:
            std::array< std::shared_ptr<FiniteElementType>, 2 > finiteElements_;
            const std::size_t order_;
        };
    }
}

#endif // DUNE_PDELAB_FINITEELEMENTMAP_PKQKFEM_HH
