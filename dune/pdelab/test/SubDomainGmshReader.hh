// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifndef DUNE_GMSHREADER_KEEPING_RENUMBER_HH
#define DUNE_GMSHREADER_KEEPING_RENUMBER_HH

#include <cstdarg>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/common/fvector.hh>
#include <dune/common/to_unique_ptr.hh>

#include <dune/geometry/type.hh>

#include <dune/grid/common/boundarysegment.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dune/grid/io/file/gmshreader.hh>

namespace Dune
{

  //! dimension independent parts for GmshReaderParser
  template<typename GridType>
  class SubDomainGmshReaderParser
  {
  protected:
    // private data
    Dune::GridFactory<GridType>& factory;
    bool verbose;
    bool insert_boundary_segments;
    unsigned int number_of_real_vertices;
    int boundary_element_count;
    int element_count;
    // read buffer
    char buf[512];
    std::string fileName;
    // exported data
    std::vector<int> boundary_id_to_physical_entity;
    std::vector<int> element_index_to_physical_entity;

    // static data
    static const int dim = GridType::dimension;
    static const int dimWorld = GridType::dimensionworld;
    static_assert( (dimWorld <= 3), "GmshReader requires dimWorld <= 3." );

    // typedefs
    typedef FieldVector< double, dimWorld > GlobalVector;

    // don't use something like
    //   readfile(file, 1, "%s\n", buf);
    // to skip the rest of of the line -- that will only skip the next
    // whitespace-separated word!  Use skipline() instead.
    void readfile(FILE * file, int cnt, const char * format,
                  void* t1, void* t2 = 0, void* t3 = 0, void* t4 = 0,
                  void* t5 = 0, void* t6 = 0, void* t7 = 0, void* t8 = 0,
                  void* t9 = 0, void* t10 = 0)
    {
      off_t pos = ftello(file);
      int c = fscanf(file, format, t1, t2, t3, t4, t5, t6, t7, t8, t9, t10);
      if (c != cnt)
        DUNE_THROW(Dune::IOError, "Error parsing " << fileName << " "
                   "file pos " << pos
                                                   << ": Expected '" << format << "', only read " << c << " entries instead of " << cnt << ".");
    }

    // skip over the rest of the line, including the terminating newline
    void skipline(FILE * file)
    {
      int c;
      do {
        c = std::fgetc(file);
      } while(c != '\n' && c != EOF);
    }

  public:

    SubDomainGmshReaderParser(Dune::GridFactory<GridType>& _factory, bool v, bool i) :
      factory(_factory), verbose(v), insert_boundary_segments(i) {}

    std::vector<int> & boundaryIdMap()
    {
      return boundary_id_to_physical_entity;
    }

    std::vector<int> & elementIndexMap()
    {
      return element_index_to_physical_entity;
    }

    void read_and_save_renumbering (const std::string& f, std::vector<int>& renumbering)
    {
      if (verbose) std::cout << "Reading " << dim << "d Gmsh grid..." << std::endl;

      // open file name, we use C I/O
      fileName = f;
      FILE* file = fopen(fileName.c_str(),"rb");
      if (file==0)
        DUNE_THROW(Dune::IOError, "Could not open " << fileName);

      //=========================================
      // Header: Read vertices into vector
      //         Check vertices that are needed
      //=========================================

      number_of_real_vertices = 0;
      boundary_element_count = 0;
      element_count = 0;

      // process header
      double version_number;
      int file_type, data_size;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$MeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $MeshFormat in first line");
      readfile(file,3,"%lg %d %d\n",&version_number,&file_type,&data_size);
      if( (version_number < 2.0) || (version_number > 2.2) )
        DUNE_THROW(Dune::IOError, "can only read Gmsh version 2 files");
      if (verbose) std::cout << "version " << version_number << " Gmsh file detected" << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndMeshFormat")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndMeshFormat");

      // node section
      int number_of_nodes;

      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Nodes")!=0)
        DUNE_THROW(Dune::IOError, "expected $Nodes");
      readfile(file,1,"%d\n",&number_of_nodes);
      if (verbose) std::cout << "file contains " << number_of_nodes << " nodes" << std::endl;

      // read nodes
      // The '+1' is due to the fact that gmsh numbers node starting from 1 rather than from 0
      std::vector< GlobalVector > nodes( number_of_nodes+1 );
      {
        int id;
        double x[ 3 ];
        for( int i = 1; i <= number_of_nodes; ++i )
        {
          readfile(file,4, "%d %lg %lg %lg\n", &id, &x[ 0 ], &x[ 1 ], &x[ 2 ] );

          if (id > number_of_nodes) {
            DUNE_THROW(Dune::IOError,
                       "Only dense sequences of node indices are currently supported (node index "
                       << id << " is invalid).");
          }

          // just store node position
          for( int j = 0; j < dimWorld; ++j )
            nodes[ id ][ j ] = x[ j ];
        }
        readfile(file,1,"%s\n",buf);
        if (strcmp(buf,"$EndNodes")!=0)
          DUNE_THROW(Dune::IOError, "expected $EndNodes");
      }

      // element section
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$Elements")!=0)
        DUNE_THROW(Dune::IOError, "expected $Elements");
      int number_of_elements;
      readfile(file,1,"%d\n",&number_of_elements);
      if (verbose) std::cout << "file contains " << number_of_elements << " elements" << std::endl;

      //=========================================
      // Pass 1: Select and insert those vertices in the file that
      //    actually occur as corners of an element.
      //=========================================

      off_t section_element_offset = ftello(file);
      std::map<int,unsigned int> renumber;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          // k == 1: physical entity (not used here)
          // k == 2: elementary entity (not used here either)
          // if version_number < 2.2:
          //   k == 3: mesh partition 0
          // else
          //   k == 3: number of mesh partitions
          //   k => 4: mesh partition k-4
        }
        pass1HandleElement(file, elm_type, renumber, nodes);
      }

      if (verbose) std::cout << "number of real vertices = " << number_of_real_vertices << std::endl;
      if (verbose) std::cout << "number of boundary elements = " << boundary_element_count << std::endl;
      if (verbose) std::cout << "number of elements = " << element_count << std::endl;
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");
      boundary_id_to_physical_entity.resize(boundary_element_count);
      element_index_to_physical_entity.resize(element_count);

      //==============================================
      // Pass 2: Insert boundary segments and elements
      //==============================================

      fseeko(file, section_element_offset, SEEK_SET);
      boundary_element_count = 0;
      element_count = 0;
      for (int i=1; i<=number_of_elements; i++)
      {
        int id, elm_type, number_of_tags;
        readfile(file,3,"%d %d %d ",&id,&elm_type,&number_of_tags);
        int physical_entity = -1;

        for (int k=1; k<=number_of_tags; k++)
        {
          int blub;
          readfile(file,1,"%d ",&blub);
          if (k==1) physical_entity = blub;
        }
        pass2HandleElement(file, elm_type, renumber, nodes, physical_entity);
      }
      readfile(file,1,"%s\n",buf);
      if (strcmp(buf,"$EndElements")!=0)
        DUNE_THROW(Dune::IOError, "expected $EndElements");

      fclose(file);

      renumbering.resize(number_of_nodes);
      for(int i=0; i<renumbering.size(); i++)
        renumbering[i] = renumber[i+1]; // +1 because renumber goes from 1 to n like in gmsh

    }

    /** \brief Process one element during the first pass through the list of all elements
     *
     * Mainly, the method inserts all vertices needed by the current element,
     * unless they have been inserted already for a previous element.
     */
    void pass1HandleElement(FILE* file, const int elm_type,
                            std::map<int,unsigned int> & renumber,
                            const std::vector< GlobalVector > & nodes)
    {
      // some data about gmsh elements
      const int nDofs[16]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10, -1, -1, -1, 1};
      const int nVertices[16]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4, -1, -1, -1, 1};
      const int elementDim[16] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3, -1, -1, -1, 0};

      // test whether we support the element type
      if ( not (elm_type > 0 && elm_type <= 15         // index in suitable range?
                && (elementDim[elm_type] == dim || elementDim[elm_type] == (dim-1) ) ) )         // real element or boundary element?
      {
        skipline(file);         // skip rest of line if element is unknown
        return;
      }

      // The format string for parsing is n times '%d' in a row
      std::string formatString = "%d";
      for (int i=1; i<nDofs[elm_type]; i++)
        formatString += " %d";
      formatString += "\n";

      // '10' is the largest number of dofs we may encounter in a .msh file
      std::vector<int> elementDofs(10);

      readfile(file,nDofs[elm_type], formatString.c_str(),
               &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]),
               &(elementDofs[3]),&(elementDofs[4]),&(elementDofs[5]),
               &(elementDofs[6]),&(elementDofs[7]),&(elementDofs[8]),
               &(elementDofs[9]));

      // insert each vertex if it hasn't been inserted already
      for (int i=0; i<nVertices[elm_type]; i++){
        if (renumber.find(elementDofs[i])==renumber.end())
        {
          renumber[elementDofs[i]] = number_of_real_vertices++;
          factory.insertVertex(nodes[elementDofs[i]]);
          // std::cout << "elementDofs[i]: " << elementDofs[i] << ", number_of_real_vertices++: " << number_of_real_vertices << std::endl;
          // std::cout << "renumber[elementDofs[i]]: " << renumber[elementDofs[i]] << std::endl;
        }
      }

      // count elements and boundary elements
      if (elementDim[elm_type] == dim)
        element_count++;
      else
        boundary_element_count++;

    }



    // generic-case: This is not supposed to be used at runtime.
    template <class E, class V, class V2>
    void boundarysegment_insert(
      const V& nodes,
      const E& elementDofs,
      const V2& vertices
      )
    {
      DUNE_THROW(Dune::IOError, "tried to create a 3D boundary segment in a non-3D Grid");
    }

    // 3d-case:
    template <class E, class V>
    void boundarysegment_insert(
      const std::vector<FieldVector<double, 3> >& nodes,
      const E& elementDofs,
      const V& vertices
      )
    {
      std::array<FieldVector<double,dimWorld>, 6> v;
      for (int i=0; i<6; i++)
        for (int j=0; j<dimWorld; j++)
          v[i][j] = nodes[elementDofs[i]][j];

      BoundarySegment<dim,dimWorld>* newBoundarySegment
        = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 3, 3 >( v[0], v[1], v[2],
                                                                                           v[3], v[4], v[5] );

      factory.insertBoundarySegment( vertices,
                                     std::shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment) );
    }



    /** \brief Process one element during the second pass through the list of all elements
     *
     * This method actually inserts the element into the grid factory.
     */
    virtual void pass2HandleElement(FILE* file, const int elm_type,
                                    std::map<int,unsigned int> & renumber,
                                    const std::vector< GlobalVector > & nodes,
                                    const int physical_entity)
    {
      // some data about gmsh elements
      const int nDofs[16]      = {-1, 2, 3, 4, 4, 8, 6, 5, 3, 6, -1, 10, -1, -1, -1, 1};
      const int nVertices[16]  = {-1, 2, 3, 4, 4, 8, 6, 5, 2, 3, -1, 4, -1, -1, -1, 1};
      const int elementDim[16] = {-1, 1, 2, 2, 3, 3, 3, 3, 1, 2, -1, 3, -1, -1, -1, 0};

      // test whether we support the element type
      if ( not (elm_type > 0 && elm_type <= 15         // index in suitable range?
                && (elementDim[elm_type] == dim || elementDim[elm_type] == (dim-1) ) ) )         // real element or boundary element?
      {
        skipline(file);         // skip rest of line if element is unknown
        return;
      }

      // The format string for parsing is n times '%d' in a row
      std::string formatString = "%d";
      for (int i=1; i<nDofs[elm_type]; i++)
        formatString += " %d";
      formatString += "\n";

      // '10' is the largest number of dofs we may encounter in a .msh file
      std::vector<int> elementDofs(10);

      readfile(file,nDofs[elm_type], formatString.c_str(),
               &(elementDofs[0]),&(elementDofs[1]),&(elementDofs[2]),
               &(elementDofs[3]),&(elementDofs[4]),&(elementDofs[5]),
               &(elementDofs[6]),&(elementDofs[7]),&(elementDofs[8]),
               &(elementDofs[9]));

      // correct differences between gmsh and Dune in the local vertex numbering
      switch (elm_type)
      {
      case 3 :          // 4-node quadrilateral
        std::swap(elementDofs[2],elementDofs[3]);
        break;
      case 5 :          // 8-node hexahedron
        std::swap(elementDofs[2],elementDofs[3]);
        std::swap(elementDofs[6],elementDofs[7]);
        break;
      case 7 :          // 5-node pyramid
        std::swap(elementDofs[2],elementDofs[3]);
        break;
      }

      // renumber corners to account for the explicitly given vertex
      // numbering in the file
      std::vector<unsigned int> vertices(nVertices[elm_type]);

      for (int i=0; i<nVertices[elm_type]; i++)
        vertices[i] = renumber[elementDofs[i]];

      // If it is an element, insert it as such
      if (elementDim[elm_type] == dim) {

        switch (elm_type)
        {
        case 1 :            // 2-node line
          factory.insertElement(Dune::GeometryTypes::line,vertices);
          break;
        case 2 :            // 3-node triangle
          factory.insertElement(Dune::GeometryTypes::triangle,vertices);
          break;
        case 3 :            // 4-node quadrilateral
          factory.insertElement(Dune::GeometryTypes::quadrilateral,vertices);
          break;
        case 4 :            // 4-node tetrahedron
          factory.insertElement(Dune::GeometryTypes::tetrahedron,vertices);
          break;
        case 5 :            // 8-node hexahedron
          factory.insertElement(Dune::GeometryTypes::hexahedron,vertices);
          break;
        case 6 :            // 6-node prism
          factory.insertElement(Dune::GeometryTypes::prism,vertices);
          break;
        case 7 :            // 5-node pyramid
          factory.insertElement(Dune::GeometryTypes::pyramid,vertices);
          break;
        case 9 :            // 6-node triangle
          factory.insertElement(Dune::GeometryTypes::triangle,vertices);
          break;
        case 11 :            // 10-node tetrahedron
          factory.insertElement(Dune::GeometryTypes::tetrahedron,vertices);
          break;
        }

      } else {
        // it must be a boundary segment then
        if (insert_boundary_segments) {

          switch (elm_type)
          {
          case 1 :              // 2-node line
            factory.insertBoundarySegment(vertices);
            break;

          case 2 :              // 3-node triangle
            factory.insertBoundarySegment(vertices);
            break;

          case 15 :             // 1-node point
            factory.insertBoundarySegment(vertices);
            break;

          case 8 : {              // 3-node line
            std::array<FieldVector<double,dimWorld>, 3> v;
            for (int i=0; i<dimWorld; i++) {
              v[0][i] = nodes[elementDofs[0]][i];
              v[1][i] = nodes[elementDofs[2]][i];                    // yes, the renumbering is intended!
              v[2][i] = nodes[elementDofs[1]][i];
            }
            BoundarySegment<dim,dimWorld>* newBoundarySegment
              = (BoundarySegment<dim,dimWorld>*) new GmshReaderQuadraticBoundarySegment< 2, dimWorld >(v[0], v[1], v[2]);
            factory.insertBoundarySegment(vertices,
                                          std::shared_ptr<BoundarySegment<dim,dimWorld> >(newBoundarySegment));
            break;
          }
          case 9 : {              // 6-node triangle
            boundarysegment_insert(nodes, elementDofs, vertices);
            break;
          }

          }

        }
      }

      // count elements and boundary elements
      if (elementDim[elm_type] == dim) {
        element_index_to_physical_entity[element_count] = physical_entity;
        element_count++;
      } else {
        boundary_id_to_physical_entity[boundary_element_count] = physical_entity;
        boundary_element_count++;
      }

    }

  };

  /**
     \ingroup Gmsh

     \brief Read Gmsh mesh file

     Read a .msh (version 2) file generated using Gmsh and construct a grid using the grid factory interface.

     The file format used by gmsh can hold grids that are more general than the simplex grids that
     the gmsh grid generator is able to construct.  We try to read as many grids as possible, as
     long as they are valid files.  You can test this by checking whether gmsh will load the file
     and display its content.

     All grids in a gmsh file live in three-dimensional Euclidean space.  If the world dimension
     of the grid type that you are reading the file into is less than three, the remaining coordinates
     are simply ignored.

     \note Recent versions of Gmsh introduced a new .msh file format (version 4) with a different syntax.
     This is currently not supported by GmshReader. One can export to an older .msh version as follows:
       - select File&rarr;Export (or CTRL+E)
       - select file format `.msh`
       - a dialog asks for options
       - select 'Version 2 ASCII' and mark 'Save all elements'

   */
  template<typename GridType>
  class SubDomainGmshReader
  {
    //! internal general reading method
    /**
     * This method does all the highlevel steering of the reader:
     * - it will register the SubGmshReader boundary segment implementation with
     *   the factory
     * - it will ensure the reader is called on all ranks (debug mode only)
     * - proceed to construct a parser (rank 0 only)
     * - use the parser to read the grid into the factory (rank 0 only)
     * - move entity and boundary data from the parser into the data vector
     *   arguments, or clear the data vector arguments, depending on rank
     *
     * \note That the parser always reads the data vectors from the files.
     *       However, if insertBoundarySegments is false, no boundary segments
     *       are inserted into the factory, and thus there will be no correct
     *       indexing of the boundarySegmentToPhysicalEntity vector possible.
     *       For this reason, this method is not exposed to the user, and the
     *       interface methods are responsible to ensure that
     *       boundarySegmentToPhysicalEntity is discarded if boundary segments
     *       are not inserted.
     */
    static std::vector<int> do_read_and_return(Dune::GridFactory<GridType> &factory,
                        const std::string &fileName,
                        std::vector<int>& boundarySegmentToPhysicalEntity,
                        std::vector<int>& elementToPhysicalEntity,
                        bool verbose, bool insertBoundarySegments)
    {
      // register boundary segment to boundary segment factory for possible load balancing
      // this needs to be done on all cores since the type might not be known otherwise
      GmshReaderQuadraticBoundarySegment< Grid::dimension, Grid::dimensionworld >::registerFactory();
      std::vector<int> renumbering;

#ifndef NDEBUG
      // check that this method is called on all cores
      factory.comm().barrier();
#endif

      // create parse object and read grid on process 0
      if (factory.comm().rank() == 0)
      {
        SubDomainGmshReaderParser<Grid> parser(factory,verbose,insertBoundarySegments);
        parser.read_and_save_renumbering(fileName, renumbering);

        boundarySegmentToPhysicalEntity = std::move(parser.boundaryIdMap());
        elementToPhysicalEntity = std::move(parser.elementIndexMap());
      }
      else
      {
        boundarySegmentToPhysicalEntity = {};
        elementToPhysicalEntity = {};
      }

      return renumbering;
    }

    //! pass a discarded lvalue argument to a function
    /**
     * This method is intended to be used in function calls that require
     * lvalue arguments, when the caller just wants to pass in temporary
     * variable that is immediately discarded after the return of the
     * function.  It expects an rvalue argument, that is turned into an
     * lvalue.  For instance:
     * ```c++
     * do_read(factory, fileName, discarded(std::vector<int>{}),
     *         discarded(std::vector<int>{}));
     * ```
     * Here, the vectors are constructed as rvalues, passed through
     * `discarded()` which turns them into lvalues, so they can be arguments
     * to `do_read()`.  `do_read()` will fill them with some data, and they
     * will be destroyed at the end of the full-expression containing the
     * function call.
     *
     * \note It is very likely an error to use this outside a function call
     *       argument.
     */
    template<class T>
    static T &discarded(T &&value) { return value; }

    struct DataArg {
      std::vector<int> *data_ = nullptr;
      DataArg(std::vector<int> &data) : data_(&data) {}
      DataArg(const decltype(std::ignore)&) {}
      DataArg() = default;
    };

    struct DataFlagArg : DataArg {
      bool flag_ = false;
      using DataArg::DataArg;
      DataFlagArg(bool flag) : flag_(flag) {}
    };

  public:
    typedef GridType Grid;

    // /** \todo doc me
    //  *
    //  * \return The return type is a special pointer type that casts into Grid*,
    //  *    std::unique_ptr<Grid>, and std::shared_ptr<Grid>.  It is scheduled
    //  *    to be replaced by std::unique_ptr<Grid> eventually.
    //  */
    // static ToUniquePtr<Grid> read (const std::string& fileName, bool verbose = true, bool insertBoundarySegments=true)
    // {
    //   // make a grid factory
    //   Dune::GridFactory<Grid> factory;

    //   read(factory, fileName, verbose, insertBoundarySegments);

    //   return factory.createGrid();
    // }

    // /** \todo doc me
    //  *
    //  * \return The return type is a special pointer type that casts into Grid*,
    //  *    std::unique_ptr<Grid>, and std::shared_ptr<Grid>.  It is scheduled
    //  *    to be replaced by std::unique_ptr<Grid> eventually.
    //  */
    // [[deprecated("Deprecated after Dune 2.7 (ca. 2020-02): This method is "
    //              "deprecated, see "
    //              "https://gitlab.dune-project.org/flyspray/FS/issues/1698, as "
    //              "it does not allow for relating the data to entities. Please "
    //              "use read(string, bool, bool) (if you don't use the data in "
    //              "the vectors, or use read(GridFactory&, string, "
    //              "vector<int>&, vector<int>&, bool, bool) if you do use "
    //              "them.")]]
    // static ToUniquePtr<Grid> read (const std::string& fileName,
    //                    std::vector<int>& boundarySegmentToPhysicalEntity,
    //                    std::vector<int>& elementToPhysicalEntity,
    //                    bool verbose = true, bool insertBoundarySegments=true)
    // {
    //   // make a grid factory
    //   Dune::GridFactory<Grid> factory;

    //   do_read(factory, fileName, boundarySegmentToPhysicalEntity,
    //           elementToPhysicalEntity, verbose, insertBoundarySegments);

    //   return factory.createGrid();
    // }

    // /** \todo doc me */
    // static void read (Dune::GridFactory<Grid>& factory, const std::string& fileName,
    //                   bool verbose = true, bool insertBoundarySegments=true)
    // {
    //   do_read(factory, fileName, discarded(std::vector<int>{}),
    //           discarded(std::vector<int>{}), verbose, insertBoundarySegments);
    // }

    /** \todo doc me */
    static std::vector<int> read_and_return (Dune::GridFactory<Grid>& factory, const std::string& fileName,
                      bool verbose = true, bool insertBoundarySegments=true)
    {
      return do_read_and_return(factory, fileName, discarded(std::vector<int>{}),
              discarded(std::vector<int>{}), verbose, insertBoundarySegments);
    }

  //   //! read Gmsh file, possibly with data
  //   /**
  //    * \param factory             The GridFactory to fill.
  //    * \param fileName            Name of the file to read from.
  //    * \param boundarySegmentData Container to fill with boundary segment
  //    *                            physical entity data, or `std::ignore`, or a
  //    *                            `bool` value.  Boundary segments are
  //    *                            inserted when a container or `true` is
  //    *                            given, otherwise they are not inserted.
  //    * \param elementData         Container to fill with element physical
  //    *                            entity data, or `std::ignore`.
  //    * \param verbose             Whether to be chatty.
  //    *
  //    * Containers to fill with data must be `std::vector<int>` lvalues.
  //    * Element data is indexed by the insertion index of the element,
  //    * boundarySegment data is indexed by the insertion index of the boundary
  //    * intersection.  These can be obtained from the `factory`, and are lost
  //    * once the grid gets modified (refined or load-balanced).
  //    *
  //    * \note At the moment the data containers are still filled internally,
  //    *       even if they are ignored.  So not having to pass them is more of
  //    *       a convenience feature and less of an optimization.  This may
  //    *       however change in the future.
  //    */
  //   static void read (Dune::GridFactory<Grid> &factory,
  //                     const std::string &fileName,
  //                     DataFlagArg boundarySegmentData,
  //                     DataArg elementData,
  //                     bool verbose=true)
  //   {
  //     do_read(factory, fileName,
  //             boundarySegmentData.data_
  //               ? *boundarySegmentData.data_ : discarded(std::vector<int>{}),
  //             elementData.data_
  //               ? *elementData.data_ : discarded(std::vector<int>{}),
  //             verbose,
  //             boundarySegmentData.flag_ || boundarySegmentData.data_);
  //   }

  //   /** \todo doc me */
  //   [[deprecated("Deprecated after Dune 2.7 (ca. 2020-02): The "
  //                "insertBoundarySegments argument for this overload of the "
  //                "read method is deprecated, please omit it.  See "
  //                "https://gitlab.dune-project.org/flyspray/FS/issues/1698.  "
  //                "When setting insertBoundarySegments=false there is no way "
  //                "to correctly use the data returned in "
  //                "boundarySegmentToPhysicalEntity, which makes using the old "
  //                "method error-prone.")]]
  //   static void read (Dune::GridFactory<Grid>& factory,
  //                     const std::string& fileName,
  //                     std::vector<int>& boundarySegmentToPhysicalEntity,
  //                     std::vector<int>& elementToPhysicalEntity,
  //                     bool verbose, bool insertBoundarySegments)
  //   {
  //     do_read(factory, fileName, boundarySegmentToPhysicalEntity,
  //             elementToPhysicalEntity, verbose, insertBoundarySegments);
  //   }
  };

  /** \} */

} // namespace Dune

#endif
