/*
 * Modified version of the PDELab constraints assembler
 * Only assembles constraints on domain boundaries, not on processor boundaries
 */

namespace Dune {
  namespace PDELab {

    template<typename P, typename GFS, typename CG, bool isFunction>
    struct ConstraintsAssemblerHelperExterior
    {
      //! Modified version of the PDELab constraints assembler only assembling on domain boundaries
      static void
      assemble(const P& p, const GFS& gfs, CG& cg, const bool verbose)
      {
        // get some types
        using ES = typename GFS::Traits::EntitySet;
        using Element = typename ES::Traits::Element;
        using Intersection = typename ES::Traits::Intersection;

        ES es = gfs.entitySet();

        // make local function space
        using LFS = LocalFunctionSpace<GFS>;
        LFS lfs_e(gfs);
        LFSIndexCache<LFS> lfs_cache_e(lfs_e);
        LFS lfs_f(gfs);
        LFSIndexCache<LFS> lfs_cache_f(lfs_f);

        // get index set
        auto& is = es.indexSet();

        // loop once over the grid
        for (const auto& element : elements(es))
        {

          auto id = is.uniqueIndex(element);

          // bind local function space to element
          lfs_e.bind(element);

          using CL = typename CG::LocalTransformation;

          CL cl_self;

          using ElementWrapper = ElementGeometry<Element>;
          using IntersectionWrapper = IntersectionGeometry<Intersection>;

          // iterate over intersections and call metaprogram
          unsigned int intersection_index = 0;
          for (const auto& intersection : intersections(es,element))
          {

            auto intersection_data = classifyIntersection(es,intersection);
            auto intersection_type = std::get<0>(intersection_data);
            auto& outside_element = std::get<1>(intersection_data);

            switch (intersection_type) {

            case IntersectionType::skeleton:
            case IntersectionType::periodic:
              {
                break;
              }

            case IntersectionType::boundary:
              TypeTree::applyToTreePair(p,lfs_e,BoundaryConstraints<IntersectionWrapper,CL>(IntersectionWrapper(intersection,intersection_index),cl_self));
              break;

            case IntersectionType::processor:
              break;

            }
            ++intersection_index;
          }

          if (!cl_self.empty())
            {
              lfs_cache_e.update();
              cg.import_local_transformation(cl_self,lfs_cache_e);
            }

        }

        // print result
        if(verbose){
          std::cout << "constraints:" << std::endl;

          std::cout << cg.size() << " constrained degrees of freedom" << std::endl;

          for (const auto& col : cg)
          {
            std::cout << col.first << ": "; // col index
            for (const auto& row : col.second)
              std::cout << "(" << row.first << "," << row.second << ") "; // row index , value
            std::cout << std::endl;
          }
        }
      }
    }; // end ConstraintsAssemblerHelper

    /*! \brief Assemble constraints exclusively on domain boundary ignoring processor boundaries
     */
    template<typename P, typename GFS, typename CG>
    void constraints_exterior(const P& p, const GFS& gfs, CG& cg,
                     const bool verbose = false)
    {
      // clear global constraints
      cg.clear();
      ConstraintsAssemblerHelperExterior<P, GFS, CG, IsGridFunction<P>::value>::assemble(p,gfs,cg,verbose);
    }

  }
}
