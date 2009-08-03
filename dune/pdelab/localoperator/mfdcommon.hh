#ifndef DUNE_PDELAB_MFDCOMMON_HH
#define DUNE_PDELAB_MFDCOMMON_HH

#include <vector>
#include <dune/common/fvector.hh>
#include <dune/grid/common/genericreferenceelements.hh>

namespace Dune
{
    namespace PDELab
    {
        // Cache for cell information used by mimetic finite difference mthods
        template<class ctype,int dim>
        struct MimeticCellProperties
        {
            typedef FieldVector<ctype,dim> Vector;

            Vector center;
            ctype volume;
            unsigned int num_faces;

            std::vector<ctype> face_areas;
            std::vector<Vector> face_centers;
            std::vector<bool> boundary;
            std::vector<ctype> R;
            std::vector<ctype> N;

            // initialise MimeticCellProperties for given entity
            template<class Entity>
            void init(const Entity& ent)
            {
                const typename Entity::Geometry& cgeo = ent.geometry();
                const FieldVector<ctype,dim> center_local
                    = GenericReferenceElements<ctype,dim>::general(cgeo.type()).position(0,0);
                center = cgeo.global(center_local);
                volume = cgeo.volume();
                num_faces = 0;

                face_areas.clear();
                face_centers.clear();
                boundary.clear();
                R.clear();
                N.clear();
            }

            // Append information on face described by the intersection object
            template<class Intersection>
            void add_face(const Intersection& is)
            {
                const typename Intersection::Geometry& fgeo = is.geometry();
                const FieldVector<ctype,dim-1> face_center_local
                    = GenericReferenceElements<ctype,dim-1>::general(fgeo.type()).position(0,0);
                const Vector face_center = fgeo.global(face_center_local);
                const ctype face_area = fgeo.volume();

                face_areas.push_back(face_area);
                face_centers.push_back(face_center);
                boundary.push_back(is.boundary());

                ++num_faces;

                // construct matrix R
                for (int i = 0; i < dim; ++i)
                    R.push_back((face_center[i]-center[i]) * face_area);

                // construct matrix N
                Vector face_normal = is.unitOuterNormal(face_center_local);
                for (int i = 0; i < dim; ++i)
                    N.push_back(face_normal[i]);
            }

            void print() const
            {
                std::cout.precision(3);
                std::cout << "Cell center (" << center[0];
                for (int i = 1; i < dim; ++i)
                    std::cout << ", " << center[i];
                std::cout << "),   volume " << volume << std::endl;
                std::cout << "face areas " << std::setw(6) << face_areas[0];
                for (unsigned i = 1; i < face_areas.size(); ++i)
                    std::cout << ", " << std::setw(6) << face_areas[i];
                std::cout << "\n    is boundary  " << std::setw(4) << boundary[0];
                for (unsigned i = 1; i < boundary.size(); ++i)
                    std::cout << ", " << std::setw(4) << boundary[i];
                std::cout << std::endl;
                std::cout.precision(6);
            }
        };

        template<class ctype,int dim>
        class MimeticBrezziWBase
        {
        protected:
            mutable std::vector<ctype> W_indep;

            void build_W_indep(const MimeticCellProperties<ctype,dim>& cell) const
            {
                int R_max = cell.R.size();

                // orthonormalize R
                R_tilde = cell.R;
                sp.reserve(dim - 1);
                for(int i = 0; i < dim; ++i)
                {
                    // orthogonalize
                    sp.clear(); sp.resize(i, 0.0);
                    for(int j = 0; j < i; ++j)
                        for(int k = 0; k < R_max; k += dim)
                            sp[j] += R_tilde[k+i] * R_tilde[k+j];
                    for(int j = 0; j < i; ++j)
                        for(int k = 0; k < R_max; k += dim)
                            R_tilde[k+i] -= sp[j]*R_tilde[k+j];

                    // normalize
                    ctype norm = 0.0;
                    for(int k = i; k < R_max; k += dim)
                        norm += R_tilde[k] * R_tilde[k];
                    norm = 1.0 / std::sqrt(norm);
                    for(int k = i; k < R_max; k += dim)
                        R_tilde[k] *= norm;
                }

                // construct the part of matrix W that is independent of K
                W_indep.clear();
                W_indep.resize(cell.num_faces * cell.num_faces, 0.0);
                for (unsigned i = 0; i < W_indep.size(); i += cell.num_faces+1)
                    W_indep[i] = 1.0;
                for (int i = 0, m = 0; i < R_max; i += dim)
                    for (int j = 0; j < R_max; j += dim, ++m)
                        for (int k = 0; k < dim; ++k)
                            W_indep[m] -= R_tilde[i+k] * R_tilde[j+k];
            }

        private:
            // local variables of build_W_indep; defined here only for
            // performance reasons (no need to reallocate each time)
            mutable std::vector<ctype> R_tilde;
            mutable std::vector<ctype> sp;
        };

        template<class ctype,int dim>
        class MimeticBrezziW : public MimeticBrezziWBase<ctype,dim>
        {
            const ctype alpha;

        public:
            explicit MimeticBrezziW(ctype alpha_ = 1.0) : alpha(alpha_) {}

            void build_W(const MimeticCellProperties<ctype,dim>& cell,
                         const FieldMatrix<ctype,dim,dim>& K, std::vector<ctype>& W) const
            {
                this->build_W_indep(cell);

                W.resize(this->W_indep.size());
                ctype u = 0.0;
                for (int i = 0; i < dim; ++i)
                    u += K[i][i];
                u *= 2.0 * alpha / dim;
                ctype inv_volume = 1.0 / cell.volume;
                for (int i = 0, m = 0; i < cell.num_faces; ++i)
                    for (int j = 0; j < cell.num_faces; ++j, ++m)
                    {
                        W[m] = u * this->W_indep[m];
                        for (int k = 0; k < dim; ++k)
                            for (int l = 0; l < dim; ++l)
                                W[m] += cell.N[i*dim+k] * K[k][l] * cell.N[j*dim+l];
                        W[m] *= cell.face_areas[i] * cell.face_areas[j] * inv_volume;
                    }
            }
        };
    }
}

#endif
