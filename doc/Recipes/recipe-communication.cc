// -*- tab-width: 4; indent-tabs-mode: nil -*-
// always include the config file
#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
// C++ includes
#include<math.h>
#include<iostream>
// dune-common includes
#include<dune/common/parallel/mpihelper.hh>
#include<dune/common/parametertreeparser.hh>
#include<dune/common/timer.hh>
// dune-geometry includes
#include<dune/geometry/referenceelements.hh>
#include<dune/geometry/quadraturerules.hh>
// dune-grid includes
#include<dune/grid/onedgrid.hh>
#include<dune/grid/yaspgrid.hh>
#include<dune/grid/utility/structuredgridfactory.hh>
#include<dune/grid/io/file/vtk/vtkwriter.hh>
#include<dune/grid/io/file/vtk/subsamplingvtkwriter.hh>
#include<dune/grid/io/file/gmshreader.hh>
#if HAVE_UG
#include<dune/grid/uggrid.hh>
#endif
#if HAVE_DUNE_ALUGRID
#include<dune/alugrid/grid.hh>
#include<dune/alugrid/dgf.hh>
#include<dune/grid/io/file/dgfparser/dgfparser.hh>
#endif

#include <dune/pdelab.hh>

/**
 * \page recipe-communication Communication in parallel programs
 *
 * This recipe explains two types of communication available in GridView:
 * neighbourwise communication designed for domain decomposition methods,
 * and CollectiveCommunication.
 *
 * Parallel solvers in DUNE already use the communication, so it is often
 * possible to run parallel models without communicating explicitly.
 *
 * For complete communication preview check section 4 in
 * <a href="https://www.dune-project.org/modules/dune-pdelab-tutorials/">tutorial06</a>.
 *
 * \section neighbourwise-communication Neighbourwise communication
 *
 * This type of communication is designed to communicate shared degrees
 * of freedom between domains. The communication happens only between
 * neighbouring domains, and only at specified parts they have in common.
 *
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~{.cpp}
 * template <class DataHandleImp, class DataType>
 * void Dune::GridView<ViewTraits>::communicate(CommDataHandleIF<DataHandleImp,DataType> &dh, InterfaceType iftype, CommunicationDirection dir) const
 * ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 *
 * The communicate function is accessible from <a href="https://dune-project.org/doxygen/master/group__GIGridView.html">GridView</a> object
 * \snippet recipe-communication.cc Define gv
 *
 * To tell the method what data to communicate, we provide data handle (dh)
 * encapsulating the data vector,
 * \snippet recipe-communication.cc Define DataHandle
 *
 * and InterfaceType describing which entities are sent and received.
 *
 * \snippet recipe-communication.cc Communication type
 * <table>
 * <caption>Table of interface types</caption>
 * <tr> <td> InteriorBorder_InteriorBorder_Interface <td> send/receive interior and border entities
 * <tr> <td> InteriorBorder_All_Interface            <td> send interior and border, receive all entities
 * <tr> <td> Overlap_OverlapFront_Interface          <td> send overlap, receive overlap and front entities
 * <tr> <td> Overlap_All_Interface                   <td> send overlap, receive all entities
 * <tr> <td> All_All_Interface                       <td> send all and receive all entities
 * </table>
 *
 * \section collective-communication Collective communication
 *
 * This type of communication shares data between all ranks. Offers
 * many MPI methods, for example
 *
 * <table>
 * <caption>Table of collective communication functions</caption>
 * <tr><th> Method name   <th>  Description
 * <tr><td>\lstinline rank       <td> obtain number (rank) of this process
 * <tr><td>\lstinline size       <td> obtain number of processes
 * <tr><td>\lstinline barrier    <td> wait until all process arrived at the barrier
 * <tr><td>\lstinline min        <td> global min of local values
 * <tr><td>\lstinline max        <td> global max of local values
 * <tr><td>\lstinline sum        <td> global sum of local values
 * <tr><td>\lstinline allreduce  <td> compute something over all processes for each component of \n
 *                                    an array and return result in every process
 * <tr><td>\lstinline broadcast  <td> broadcast from one process to all other processes
 * <tr><td>\lstinline scatter    <td> scatter individual data from root process to all other tasks
 * <tr><td>\lstinline gather, allgather  <td> gather data on root process (and distribute it to all other tasks)
 * </table>
 *
 * The communication object is a part of the GridView
 *
 * \snippet recipe-communication.cc Collective communication object
 *
 * Most methods take a constant reference and return that type. We need to use
 * a variable as an argument, and not forget to store the result.
 *
 * \snippet recipe-communication.cc Collective communication
 *
 * Full example code: @ref recipe-communication.cc
 * \example recipe-communication.cc
 * See explanation at @ref recipe-communication
*/


/**
 * Taken 6.6.2019 from tutorial06/exercise/solution/solution06-2.cc
 * Added commentaries, changed data[i] to store (100^rank) instead of (rank).
 * Removed ParameterTree and dependancy on .ini file.
 * Removed customized DataHandle, use default one.
*/
#ifndef COMMUNICATE_HH
#define COMMUNICATE_HH

template <typename GV>
void communicate(const GV& gv, int communicationType){

  using RF = int; // RangeField
  using CON = Dune::PDELab::NoConstraints;
  using VBE = Dune::PDELab::ISTL::VectorBackend<>;
  using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,RF,RF,1>;
  FEM fem(gv);
  using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM,CON,VBE>;
  GFS gfs(gv,fem);
  using Z = Dune::PDELab::Backend::Vector<GFS, RF>; // data type
  Z z(gfs);
  // [Define DataHandle]
  using DH = Dune::PDELab::AddDataHandle<GFS,Z>;
  DH dh(gfs,z);
  //! [Define DataHandle]

  // Create collective communication object
  // [Collective communication object]
  auto comm = gv.comm();
  //! [Collective communication object]
  // [Get rank]
  int myrank = comm.rank();
  //! [Get rank]
  // Store the 100^rank or rank of the current processor as data for each element.
  using Dune::PDELab::Backend::native;
  int size{0};
  for(auto& v : native(z)){
    v = 1;
    ++size;
    for(int j=0 ; j<myrank; ++j)
      // v *= 1000; // makes it easy to see which ranks communicated ths entity
      v += 1; // makes VTK outputs readable
  }

  // Different communication types for DataHandles:
  //
  // InteriorBorder_InteriorBorder_Interface: send/receive interior and border entities
  // InteriorBorder_All_Interface:            send interior and border, receive all entities
  // Overlap_OverlapFront_Interface:          send overlap, receive overlap and front entities
  // Overlap_All_Interface:                   send overlap, receive all entities
  // All_All_Interface:                       send all and receive all entities
  // [Communication type]
  switch (communicationType){
    case 1:  gv.communicate(dh, Dune::InteriorBorder_InteriorBorder_Interface ,Dune::ForwardCommunication); break;
    case 2:  gv.communicate(dh, Dune::InteriorBorder_All_Interface            ,Dune::ForwardCommunication); break;
    case 3:  gv.communicate(dh, Dune::Overlap_OverlapFront_Interface          ,Dune::ForwardCommunication); break;
    case 4:  gv.communicate(dh, Dune::Overlap_All_Interface                   ,Dune::ForwardCommunication); break;
    default: gv.communicate(dh, Dune::All_All_Interface                       ,Dune::ForwardCommunication);
  }
  //! [Communication type]

  // Calculate the sum of the data vector on this rank
  int sum = z.one_norm();

  // If we are on rank 0 print the results.
  if (myrank==0){
    std::cout << std::endl;
    std::cout << "== Output for rank " << myrank << std::endl;
    std::cout << std::endl;
    std::cout << "Each process stores values equal to 1000 powered to its rank, or only rank." << std::endl;
    std::cout << "The sum shows how many entities are communicated and from which ranks they are." << std::endl;
    std::cout << "The size of the data vector is equal to the number of all entities of this processor." << std::endl;
    std::cout << std::endl;
    std::cout << "Sum of the data vector: " << sum << std::endl;
    std::cout << "Size of the data vector: " << size << std::endl;
  }

  // Find the maximal and total sum on all ranks:
  int globmax{0};
  int globsum{0};
  // [Collective communication]
  globmax = comm.max(sum);
  globsum = comm.sum(sum);
  //! [Collective communication]
  if (myrank==0){
    std::cout << "Maximal sum on all ranks is " << globmax << std::endl;
    std::cout << "Total sum on all ranks is " << globsum << std::endl;
  }

  // Make a grid function out of z
  typedef Dune::PDELab::DiscreteGridFunction<GFS,Z> ZDGF;
  ZDGF zdgf(gfs,z);
  // prepare VTK writer and write the file,
  int subsampling{1};
  using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
  VTKWRITER vtkwriter(gv,Dune::refinementIntervals(subsampling));
  std::string filename="recipe-communication";
  std::string outputname="commType_"+std::to_string(communicationType);
  typedef Dune::PDELab::VTKGridFunctionAdapter<ZDGF> VTKF;
  vtkwriter.addVertexData(std::shared_ptr<VTKF>(new VTKF(zdgf,outputname)));
  vtkwriter.write(filename,Dune::VTK::appendedraw);
}

#endif


//===============================================================
// Main program with grid setup
//===============================================================
int main(int argc, char** argv)
{
  try{
    // Maybe initialize Mpi
    Dune::MPIHelper&
      helper = Dune::MPIHelper::instance(argc, argv);
    if(Dune::MPIHelper::isFake)
      std::cout<< "This is a sequential program." << std::endl;
    else
      std::cout << "Parallel code run on "
                << helper.size() << " process(es)" << std::endl;

    // read ini file
    const int overlap = 2;
    const int refinement = 0;

    // Create 2D YaspGrid
    constexpr int dim=2;
    typedef Dune::YaspGrid<dim> Grid;
    typedef Grid::ctype DF;
    Dune::FieldVector<DF,dim> L{1.,1.};
    std::array<int,dim> N{16,16};
    std::bitset<dim> B(false); // periodic boundary (left-right, up-bottom)
    std::shared_ptr<Grid> gridp = std::shared_ptr<Grid>(new Grid(L,N,B,overlap));//,Dune::MPIHelper::getCollectiveCommunication()));
    gridp->refineOptions(false); // keep overlap in cells
    gridp->globalRefine(refinement);
    //! [Define gv]
    typedef Grid::LeafGridView GV;
    GV gv=gridp->leafGridView();
    //! [Define gv]

    int communicationType = 1;
    communicate(gv,communicationType);
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
