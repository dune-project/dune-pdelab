#ifndef DUNE_PDELAB_GRID_WRITESUBDOMAINMESHTOFILE_HH
#define DUNE_PDELAB_GRID_WRITESUBDOMAINMESHTOFILE_HH

#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>
#include <iterator>

#include <dune/common/exceptions.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/referenceelements.hh>
#include <dune/grid/common/grid.hh>
#include <dune/grid/common/mcmgmapper.hh>

// #include <dune/grid/io/file/gmshwriter.hh>

namespace Dune {


  class Element
  {
    public:

    Element() {}
    int index;
    std::vector<int> nodes;
    std::size_t type;
  };

  template<class GridView>
  class SubDomainGmshWriter // : public GmshWriter<GridView>
  {
    private:
    const GridView gv;
    int precision;
    static const unsigned int dim = GridView::dimension;
    static const unsigned int dimWorld = GridView::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshWriter requires dimWorld <= 3." );


    /** \brief Returns index of i-th vertex of an element, plus 1 (for gmsh numbering) */
    template<typename Entity>
    std::size_t nodeIndexFromEntity(const Entity& entity, int i) const {
      return gv.indexSet().subIndex(entity, i, dim);
    }

    /** \brief Translate GeometryType to corresponding Gmsh element type number
      * \throws IOError if there is no equivalent type in Gmsh
      */
    static std::size_t translateDuneToGmshType(const GeometryType& type) {
      std::size_t element_type;

      if (type.isLine())
        element_type = 1;
      else if (type.isTriangle())
        element_type = 2;
      else if (type.isQuadrilateral())
        element_type = 3;
      else if (type.isTetrahedron())
        element_type = 4;
      else if (type.isHexahedron())
        element_type = 5;
      else if (type.isPrism())
        element_type = 6;
      else if (type.isPyramid())
        element_type = 7;
      else if (type.isVertex())
        element_type = 15;
      else
        DUNE_THROW(Dune::IOError, "GeometryType " << type << " is not supported by gmsh.");

      return element_type;
    }

    /** \brief Writes all the vertices of a grid line by line
     *
     * Each line has the format
     *  node-number x-coord y-coord z-coord
     * The node-numbers will most certainly not have the arrangement "1, 2, 3, ...".
     */
    void outputNodes(std::vector<std::ofstream>& files, std::vector< std::vector< int > >& NodesMapBack) const {

      int subdomains_number = files.size();
      std::vector<int> iterators(subdomains_number);
      for(int incrSD=0; incrSD< iterators.size(); incrSD++)
        iterators[incrSD]=0;

      for (const auto& vertex : vertices(gv)) {

        const auto globalCoord = vertex.geometry().center();
        const auto nodeIndex = gv.indexSet().index(vertex);

        for(int incrSD=0; incrSD< subdomains_number; incrSD++){
          // std::cout << "nodeIndex = " << nodeIndex << ", NodesToBeSavedForEachSubdomain[iterators[incrSD]] = " << NodesToBeSavedForEachSubdomain[incrSD][iterators[incrSD]] << std::endl;
          if(nodeIndex==NodesToBeSavedForEachSubdomain[incrSD][iterators[incrSD]]){
            if (1 == dimWorld)
              files[incrSD] << iterators[incrSD]+1 << " " << globalCoord[0] << " " << 0 << " " << 0 << std::endl;
            else if (2 == dimWorld)
              files[incrSD] << iterators[incrSD]+1 << " " << globalCoord[0] << " " << globalCoord[1] << " " << 0 << std::endl;
            else // (3 == dimWorld)
              files[incrSD] << iterators[incrSD]+1 << " " << globalCoord[0] << " " << globalCoord[1] << " " << globalCoord[2] << std::endl;

            NodesMapBack[incrSD][nodeIndex] = iterators[incrSD]+1;
            iterators[incrSD]+=1;
          }
        }
      }
    }

     /** \brief Writes all the elements of a grid line by line
     *
     * Each line has the format
     *    element-number element-type number-of-tags <tags> node-number-list
     * Counting of the element numbers starts by "1".
     *
     * If `physicalEntities` is not empty, each element has a tag representing its physical id.
     *
     * If `physicalBoundaries` is not empty, also the boundaries are written to the file with
     * the corresponding physical value.
     *
     * The physicalBoundaries vector need to be sorted according to the interesection
     * boundary segment index.
     */
    void outputElements(std::vector<std::ofstream>& files, std::vector< std::vector< int > >& NodesMapBack, const std::vector<int>& physicalEntities, const std::vector<int>& physicalBoundaries) const {

      for (int incrSD=0; incrSD< files.size(); incrSD++) {
        // Check whether the type is compatible. If not, close file and rethrow exception.

        int counter = 1;

        for (int incrE=0; incrE< ElementsToBeSavedForEachSubdomain[incrSD].size(); incrE++) {

          files[incrSD] << counter << " " << ElementsToBeSavedForEachSubdomain[incrSD][incrE].type;
          // If present, set the first tag to the physical entity
          if (!physicalEntities.empty())
            files[incrSD] << " " << 1 << " " << physicalEntities[ElementsToBeSavedForEachSubdomain[incrSD][incrE].index];
            // files[incrSD] << " " << 2 << " " << physicalEntities[ElementsToBeSavedForEachSubdomain[incrSD][incrE].index] << " " << ElementsToBeSavedForEachSubdomain[incrSD][incrE].index;
            // files[incrSD] << " " << 3 << " " << incrSD << " " << physicalEntities[ElementsToBeSavedForEachSubdomain[incrSD][incrE].index] << " " << ElementsToBeSavedForEachSubdomain[incrSD][incrE].index;
          else
            files[incrSD] << " " << 0;
            // files[incrSD] << " " << 1 << " " << ElementsToBeSavedForEachSubdomain[incrSD][incrE].index;
            // files[incrSD] << " " << 2 << " " << incrSD << " " << ElementsToBeSavedForEachSubdomain[incrSD][incrE].index;

          for (int k = 0; k < ElementsToBeSavedForEachSubdomain[incrSD][incrE].nodes.size(); ++k)
            files[incrSD] << " " << NodesMapBack[incrSD][ElementsToBeSavedForEachSubdomain[incrSD][incrE].nodes[k]];

          files[incrSD] << std::endl;
          counter++;

          // // Write boundaries
          // if (!physicalBoundaries.empty()) {
          //   auto refElement = referenceElement<typename GridView::ctype,dim>(entity.type());
          //   for(const auto& intersection : intersections(gv, entity)) {
          //     if(intersection.boundary()) {
          //       const auto faceLocalIndex(intersection.indexInInside());
          //       file << counter << " " << translateDuneToGmshType(intersection.type())
          //         << " " << 1 << " " << physicalBoundaries[intersection.boundarySegmentIndex()];
          //       for (int k = 0; k < intersection.geometry().corners(); ++k)
          //       {
          //         const auto vtxLocalIndex(refElement.subEntity(faceLocalIndex, 1, k, dim));
          //         file << " " << nodeIndexFromEntity(entity, vtxLocalIndex);
          //       }
          //       ++counter;
          //       file << std::endl;
          //     }
          //   }
          // }

        }
      }
    }

    public:

    /**
    * \brief Constructor expecting GridView of Grid to be written.
    * \param gridView GridView that will be written.
    * \param overlappingsubdomains std::vector<std::set<unsigned>> for each element corresp. a set of subomains containing this element.
    * \param partitions Number of subdomains.
    * \param numDigits Number of digits to use.
    */
    SubDomainGmshWriter(const GridView& gridView, const std::vector<std::set<unsigned>>& overlappingsubdomains, unsigned partitions, int numDigits=6)
    : gv(gridView),
      ElementToSubdomains_(overlappingsubdomains),
      partitions_(partitions),
      precision(numDigits)
      {
        ElementsToBeSavedForEachSubdomain.resize(partitions_);
        NodesToBeSavedForEachSubdomain.resize(partitions_);
        auto& indexset = gv.indexSet(); // to attach data to elements

        for (const auto& e : elements(gv)) {

          std::size_t element_type = translateDuneToGmshType(e.type());

          for(int i=0; i< ElementToSubdomains_[indexset.index(e)].size(); i++){
            auto it = ElementToSubdomains_[indexset.index(e)].begin();
            std::advance(it, i);

            Element tmpE;
            tmpE.index = indexset.index(e);
            tmpE.type = element_type;

              // Output list of nodes.
            // 3, 5 and 7 got different vertex numbering compared to Dune
            if (3 == element_type){
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 0));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 1));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 3));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 2));
              // NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 2));
              // NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 3));

              tmpE.nodes.push_back(nodeIndexFromEntity(e, 0));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 1));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 3));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 2));
              // tmpE.nodes.push_back(nodeIndexFromEntity(e, 2));
              // tmpE.nodes.push_back(nodeIndexFromEntity(e, 3));

            } else if (5 == element_type){
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 0));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 1));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 3));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 2));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 4));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 5));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 7));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 6));

              tmpE.nodes.push_back(nodeIndexFromEntity(e, 0));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 1));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 3));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 2));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 4));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 5));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 7));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 6));
            } else if (7 == element_type){
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 0));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 1));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 3));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 2));
              NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, 4));

              tmpE.nodes.push_back(nodeIndexFromEntity(e, 0));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 1));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 3));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 2));
              tmpE.nodes.push_back(nodeIndexFromEntity(e, 4));
            } else {
              for (int k = 0; k < e.geometry().corners(); ++k) {
                NodesToBeSavedForEachSubdomain[*it].push_back(nodeIndexFromEntity(e, k));

                tmpE.nodes.push_back(nodeIndexFromEntity(e, k));
              }
            }

            ElementsToBeSavedForEachSubdomain[*it].push_back(tmpE);
          }
        }

        for(int incrSD=0; incrSD< partitions_; incrSD++){
          std::sort(NodesToBeSavedForEachSubdomain[incrSD].begin(), NodesToBeSavedForEachSubdomain[incrSD].end());
          auto last = std::unique(NodesToBeSavedForEachSubdomain[incrSD].begin(), NodesToBeSavedForEachSubdomain[incrSD].end());
          NodesToBeSavedForEachSubdomain[incrSD].erase(last, NodesToBeSavedForEachSubdomain[incrSD].end());
        }
      }

    /**
     * \brief Set the number of digits to be used when writing the vertices. By default is 6.
     * \brief numDigits Number of digits to use.
     */
    void setPrecision(int numDigits) {
      precision = numDigits;
    }

    /**
    * \brief Write given grid in Gmsh 2.0 compatible ASCII file.
    * \param fileName Path of file. This method does not attach a ".msh"-extension by itself.
    * \param physicalEntities Physical entities for each element (optional).
    * \param physicalBoundaries Physical boundaries (optional).
    *
    * Opens the file with given name and path, stores the element data of the grid
    * and closes the file when done.
    *
    * If the optional parameter `physicalEntities` is provided, each element is written with
    * a tag representing its physical id.
    *
    * If the optional parameter `physicalBoundaries` is provided, also the boundaries
    * are written on file with the corresponding physical value.
    *
    * The physicalBoundaries vector need to be sorted according to the interesection boundary
    * segment index.
    *
    * Throws an IOError if file could not be opened or an unsupported element type is
    * encountered.
    */
    void write(const std::string& path_to_storage,
          const std::string& fileName,
          const std::vector<int>& physicalEntities=std::vector<int>(),
          const std::vector<int>& physicalBoundaries=std::vector<int>()) const {


      std::vector<std::vector<int>> NodesMapBack(partitions_);
      std::vector<std::ofstream> files(partitions_); // Open files

      for(int incrSD=0; incrSD<partitions_; incrSD++){

        NodesMapBack[incrSD].resize(gv.indexSet().size(dim), 0); // for each subdomain: need to save node of element in the new index set

        std::string localname = path_to_storage + std::to_string(incrSD)+"_"+fileName;
        files[incrSD].open(localname.c_str());
        if (!files[incrSD].is_open())
          DUNE_THROW(Dune::IOError, "Could not open " << fileName << " with write access.");

        // Set precision
        files[incrSD] << std::setprecision( precision );

        // Output Header
        files[incrSD] << "$MeshFormat" << std::endl
            << "2.2 0 " << sizeof(double) << std::endl // "2.2" for "version 2.2", "0" for ASCII
            << "$EndMeshFormat" << std::endl;

        // Output Nodes
        files[incrSD] << "$Nodes" << std::endl
            << NodesToBeSavedForEachSubdomain[incrSD].size() << std::endl;

      }

      outputNodes(files, NodesMapBack);

      for(int incrSD=0; incrSD<partitions_; incrSD++){
        files[incrSD] << "$EndNodes" << std::endl;

        // TODO : See how to deal with physical boundaries
        // Have a look at the initial gmshwriter.hh in dune grid
        // int boundariesSize(0);
        // if(!physicalBoundaries.empty())
        //   for(const auto& entity : elements(gv))
        //     for(const auto& intersection : intersections(gv, entity))
        //         if(intersection.boundary())
        //           ++boundariesSize;

        // Output Elements;
        files[incrSD] << "$Elements" << std::endl
            << ElementsToBeSavedForEachSubdomain[incrSD].size() << std::endl;
      }

      outputElements(files, NodesMapBack, physicalEntities, physicalBoundaries);

      for(int incrSD=0; incrSD<partitions_; incrSD++){
        files[incrSD] << "$EndElements" << std::endl;

        files[incrSD].close();
      }

      // Save local to global indices (for vertices and cells) out of the .msh file
      // because the gmshReader currently implemented (2.7) doesn't read multiple
      // tags for elements (only physical groups) and any tags for nodes
      for(int incrSD=0; incrSD<partitions_; incrSD++){
        Dune::BlockVector<Dune::FieldVector<int, 1>> LocalToGlobalELEMENT(ElementsToBeSavedForEachSubdomain[incrSD].size());
        for (int i=0; i<ElementsToBeSavedForEachSubdomain[incrSD].size();i++)
          LocalToGlobalELEMENT[i]=ElementsToBeSavedForEachSubdomain[incrSD][i].index;
        std::string filename_ltgE = path_to_storage  + std::to_string(incrSD) + "_LocalToGlobalElement.mm";
        Dune::storeMatrixMarket(LocalToGlobalELEMENT, filename_ltgE, 15);

        Dune::BlockVector<Dune::FieldVector<int, 1>> LocalToGlobalNODE(NodesToBeSavedForEachSubdomain[incrSD].size());
        for (int i=0; i<NodesToBeSavedForEachSubdomain[incrSD].size();i++)
          LocalToGlobalNODE[i]=NodesToBeSavedForEachSubdomain[incrSD][i];
        std::string filename_ltgN = path_to_storage  + std::to_string(incrSD) + "_LocalToGlobalNode.mm";
        Dune::storeMatrixMarket(LocalToGlobalNODE, filename_ltgN, 15);
      }
    }

    std::vector<std::set<unsigned>> ElementToSubdomains_;
    std::vector< std::vector< Element > > ElementsToBeSavedForEachSubdomain;
    std::vector< std::vector< int > > NodesToBeSavedForEachSubdomain;
    unsigned partitions_;

  };

}

#endif // DUNE_PDELAB_GRID_WRITESUBDOMAINMESHTOFILE_HH