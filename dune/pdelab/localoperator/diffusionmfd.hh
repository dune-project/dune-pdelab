#ifndef DUNE_PDELAB_DIFFUSIONMFD_HH
#define DUNE_PDELAB_DIFFUSIONMFD_HH

#include "pattern.hh"
#include "flags.hh"
#include "mfdcommon.hh"
#include "diffusionparam.hh"

namespace Dune
{
    namespace PDELab
    {

        //! \addtogroup LocalOperator
        //! \ingroup PDELab
        //! \{

        /** a local operator for solving the diffusion equation
         *
         * \f{align*}{
         * - \nabla\cdot\{K(x) \nabla u\} + a_0 u &=& f \mbox{ in } \Omega,          \\
         *                                      u &=& g \mbox{ on } \partial\Omega_D \\
         *              -(K(x)\nabla u) \cdot \nu &=& j \mbox{ on } \partial\Omega_N \\
         * \f}
         * with a mimetic finite difference method on all types of grids in any dimension
         * \tparam Data contains methods K, a_0, f, g and j giving the equation data
         * and bcType to define the boundary condition type
         * \tparam WBuilder builds the local scalar product matrix
         */
        template<class Data, class WBuilder = MimeticBrezziW<typename Data::ctype,Data::dimension> >
        class DiffusionMFD
            : public FullVolumePattern
            , public LocalOperatorDefaultFlags
        {
            static const unsigned int dim = Data::dimension;
            typedef typename Data::ctype ctype;
            typedef typename Data::rtype rtype;

        public:
            // pattern assembly flags
            enum { doPatternVolume = true };

            // residual assembly flags
            enum { doAlphaVolume = true };
            enum { doSkeletonTwoSided = true };
            enum { doAlphaSkeleton = true };
            enum { doAlphaBoundary = true };
            enum { doAlphaVolumePostSkeleton = true };

            enum { doLambdaVolume = true };
            enum { doLambdaBoundary = true };

            DiffusionMFD(const Data& data_, const WBuilder& wbuilder_ = WBuilder())
                : data(data_), wbuilder(wbuilder_)
            {}

            // alpha ***********************************************

            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume(const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const
            {
                cell.init(eg.entity());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_skeleton(const IG& ig,
                                const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                R& r_s, R& r_n) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s,
                                const LFSV& lfsv_s, R& r_s) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
            void alpha_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                            const LFSV& lfsv, R& r) const
            {
                // extract subspaces
                typedef typename LFSU::template Child<0>::Type CellUnknowns;
                const CellUnknowns& cell_space = lfsu.template child<0>();
                typedef typename LFSU::template Child<1>::Type FaceUnknowns;
                const FaceUnknowns& face_space = lfsu.template child<1>();

                // get permeability for current cell
                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = ReferenceElements<ctype,dim>::general(gt).position(0,0);
                FieldMatrix<rtype,dim,dim> K = data.K(eg.entity(), localcenter);

                // build matrix W
                wbuilder.build_W(cell, K, W);

                // Compute residual
                for(unsigned int e = 0, m = 0; e < cell.num_faces; ++e)
                    for(unsigned int f = 0; f < cell.num_faces; ++f, ++m)
                    {
		      r.accumulate(cell_space,0,W[m]
				   * (x(cell_space,0) - x(face_space,f)));
		      r.accumulate(face_space,e,-W[m]
				   * (x(cell_space,0) - x(face_space,f)));
                    }

                // Helmholtz term
                r.accumulate(cell_space,0,data.a_0(eg.entity(), localcenter)
			     * x(cell_space,0));
            }

            // jacobian ********************************************

            template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
            void jacobian_volume(const EG& eg, const LFSU& lfsu, const X& x,
                                 const LFSV& lfsv, M& mat) const
            {
                cell.init(eg.entity());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
            void jacobian_skeleton(const IG& ig,
                                   const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                   const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                   M& mat_ss, M& mat_sn,
                                   M& mat_ns, M& mat_nn) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename M>
            void jacobian_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                   M& mat_ss) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename EG, typename LFSU, typename X, typename LFSV, typename M>
            void jacobian_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                               const LFSV& lfsv, M& mat) const
            {
                // extract subspaces
                typedef typename LFSU::template Child<0>::Type CellUnknowns;
                const CellUnknowns& cell_space = lfsu.template child<0>();
                typedef typename LFSU::template Child<1>::Type FaceUnknowns;
                const FaceUnknowns& face_space = lfsu.template child<1>();

                // get permeability for current cell
                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = ReferenceElements<ctype,dim>::general(gt).position(0,0);
                FieldMatrix<rtype,dim,dim> K = data.K(eg.entity(), localcenter);

                // build matrix W
                wbuilder.build_W(cell, K, W);

                // Compute residual
                for(unsigned int e = 0, m = 0; e < cell.num_faces; ++e)
                    for(unsigned int f = 0; f < cell.num_faces; ++f, ++m)
                    {
		      mat.accumulate(cell_space,0,cell_space,0,W[m]);
		      mat.accumulate(cell_space,0,face_space,f,-W[m]);
		      mat.accumulate(face_space,f,cell_space,0,-W[m]);
		      mat.accumulate(face_space,f,face_space,f,W[m]);
                    }

                // Helmholtz term
                mat.accumulate(cell_space,0,cell_space,0,
			       data.a_0(eg.entity(), localcenter));
            }

            // jacobian_apply **************************************

            template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
            void jacobian_apply_volume(const EG& eg, const LFSU& lfsu, const X& x,
                                       const LFSV& lfsv, Y& y) const
            {
                cell.init(eg.entity());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
            void jacobian_apply_skeleton(const IG& ig,
                                         const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                         const LFSU& lfsu_n, const X& x_n, const LFSV& lfsv_n,
                                         Y& y_s, Y& y_n) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename IG, typename LFSU, typename X, typename LFSV, typename Y>
            void jacobian_apply_boundary(const IG& ig, const LFSU& lfsu_s, const X& x_s, const LFSV& lfsv_s,
                                         Y& y_s) const
            {
                cell.add_face(ig.intersection());
            }

            template<typename EG, typename LFSU, typename X, typename LFSV, typename Y>
            void jacobian_apply_volume_post_skeleton(const EG& eg, const LFSU& lfsu, const X& x,
                                                     const LFSV& lfsv, Y& y) const
            {
                // extract subspaces
                typedef typename LFSU::template Child<0>::Type CellUnknowns;
                const CellUnknowns& cell_space = lfsu.template child<0>();
                typedef typename LFSU::template Child<1>::Type FaceUnknowns;
                const FaceUnknowns& face_space = lfsu.template child<1>();

                // get permeability for current cell
                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = ReferenceElements<ctype,dim>::general(gt).position(0,0);
                FieldMatrix<rtype,dim,dim> K = data.K(eg.entity(), localcenter);

                // build matrix W
                wbuilder.build_W(cell, K, W);

                // Compute residual
                for(int e = 0, m = 0; e < cell.num_faces; ++e)
                    for(int f = 0; f < cell.num_faces; ++f, ++m)
                    {
                        y[cell_space.localIndex(0)] += W[m] * x[cell_space.localIndex(0)];
                        y[cell_space.localIndex(0)] -= W[m] * x[face_space.localIndex(f)];
                        y[face_space.localIndex(f)] -= W[m] * x[cell_space.localIndex(0)];
                        y[face_space.localIndex(f)] += W[m] * x[face_space.localIndex(f)];
                    }

                // Helmholtz term
                y[cell_space.localIndex(0)] += data.a_0(eg.entity(), localcenter)
                    * x[cell_space.localIndex(0)];
            }

            // lambda **********************************************

            template<typename EG, typename LFSV, typename R>
            void lambda_volume(const EG& eg, const LFSV& lfsv, R& r) const
            {
                // extract subspaces
                typedef typename LFSV::template Child<0>::Type CellUnknowns;
                const CellUnknowns& cell_space = lfsv.template child<0>();

                GeometryType gt = eg.geometry().type();
                FieldVector<ctype,dim> localcenter = ReferenceElements<ctype,dim>::general(gt).position(0,0);
                r.accumulate(cell_space,0,-cell.volume * data.f(eg.entity(), localcenter));
            }

            template<typename IG, typename LFSV, typename R>
            void lambda_boundary(const IG& ig, const LFSV& lfsv, R& r) const
            {
                // local index of current face
                unsigned int e = cell.num_faces - 1;

                GeometryType gt = ig.intersection().type();
                FieldVector<ctype,dim-1> center = ReferenceElements<ctype,dim-1>::general(gt).position(0,0);
                FieldVector<ctype,dim> local_face_center = ig.geometryInInside().global(center);

                if (data.bcType(ig, center) == Data::bcNeumann)
                {
                    typedef typename LFSV::template Child<1>::Type FaceUnknowns;
                    const FaceUnknowns& face_space = lfsv.template child<1>();

		    r.accumulate(face_space,e,cell.face_areas[e]
				 * data.j(*(ig.inside()), local_face_center));
                }
            }

        private:
            const Data& data;
            const WBuilder wbuilder;
            mutable MimeticCellProperties<ctype,dim> cell;
            mutable std::vector<rtype> W;
        };

        //! \} group LocalOperator
    }
}

#endif
