// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#include "config.h"

#include <cassert>
#include <cstdlib>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab.hh>


/*
 * This test makes sure that the GenericDataHandle infrastructure correctly transmits
 * the sender rank if requested, both for per-DOF and per-entity communications, and
 * that the reported communication size is still correct.
 */

// This custom gather/scatter handler sends the rank and checks its value when receiving
struct GatherScatter
{

  template<typename Buffer, typename Entity, typename LocalView>
  bool gather(Buffer& buffer, const Entity& e, const LocalView& local_view) const
  {
    std::size_t size = _expected_size >= 0 ? _expected_size : local_view.size();
    for (std::size_t i = 0 ; i < size ; ++i)
      buffer.write(_rank);
    return false;
  }

  template<typename Buffer, typename Entity, typename LocalView>
  bool scatter(Buffer& buffer, std::size_t n, const Entity& e, const LocalView& local_view) const
  {
    std::size_t size = 0;
    if (_expected_size >= 0)
      size = _expected_size;
    else
      size = local_view.size();

    assert(size == n);

    for (std::size_t i = 0 ; i < local_view.size() ; ++i)
      {
        int rank;
        buffer.read(rank);
        assert(rank == buffer.senderRank());
      }
    return false;
  }

  template<typename Buffer, typename O1, typename O2, typename Entity, typename LocalView>
  bool scatter(Buffer& buffer, const O1& o1, const O2& o2, const Entity& e, const LocalView& local_view) const
  {
    assert(false && "Should never get here - this only exists to make the test compile");
    return false;
  }

  // If expected_size is < 0, take the size from the LocalView, otherwise use the prescribed value
  GatherScatter(int rank, int expected_size)
    : _rank(rank)
    , _expected_size(expected_size)
  {}

private:

  int _rank;
  int _expected_size;

};


template<typename GFS, typename V>
class DOFDataHandle
  : public Dune::PDELab::GFSDataHandle<
      GFS,V,GatherScatter,Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>>
{

  using Base = Dune::PDELab::GFSDataHandle<
    GFS,
    V,
    GatherScatter,
    Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>
    >;

public:

  DOFDataHandle(const GFS& gfs, V& v)
    : Base(gfs,v,GatherScatter(gfs.gridView().comm().rank(),-1))
  {}

};


template<typename GFS, typename V>
class EntityDataHandle
  : public Dune::PDELab::GFSDataHandle<
      GFS,V,GatherScatter,Dune::PDELab::EntityDataCommunicationDescriptor<typename V::ElementType,true>>
{

  using Base = Dune::PDELab::GFSDataHandle<
    GFS,
    V,
    GatherScatter,
    Dune::PDELab::DOFDataCommunicationDescriptor<typename V::ElementType,true>
    >;

public:

  EntityDataHandle(const GFS& gfs, V& v, std::size_t count)
    : Base(gfs,v,GatherScatter(gfs.gridView().comm().rank(),count))
  {}

};



int main(int argc, char** argv)
{

  try {
    auto& helper = Dune::MPIHelper::instance(argc,argv);

    if (helper.getCollectiveCommunication().size() == 1)
      {
        std::cerr << "Communication test requires at least 2 MPI ranks" << std::endl;
        std::exit(77);
      }

    Dune::YaspGrid<2> grid({{1.0,1.0}},{{16,16}},0,1);
    auto gv = grid.leafGridView();
    using GV = decltype(gv);
    using F = GV::ctype;

    using FEM = Dune::PDELab::QkLocalFiniteElementMap<GV,F,F,2>;
    FEM fem(gv);

    using GFS = Dune::PDELab::GridFunctionSpace<GV,FEM>;
    GFS gfs(gv,fem);

    using V = Dune::PDELab::Backend::Vector<GFS,int>;
    auto v = V(gfs);

    {
      // start with a standard DOF communication - no rank transmission
      auto handle = Dune::PDELab::CopyDataHandle<GFS,V>(gfs,v);
      gv.communicate(handle,Dune::All_All_Interface,Dune::ForwardCommunication);
    }

    {
      // DOF communication with rank transmission
      auto handle = DOFDataHandle<GFS,V>(gfs,v);
      gv.communicate(handle,Dune::All_All_Interface,Dune::ForwardCommunication);
    }

    {
      // standard per-entity communication - no rank transmission
      auto handle = Dune::PDELab::CopyDataHandle<GFS,V>(gfs,v);
      gv.communicate(handle,Dune::All_All_Interface,Dune::ForwardCommunication);
    }

    {
      // DOF communication with rank transmission
      auto handle = DOFDataHandle<GFS,V>(gfs,v);
      gv.communicate(handle,Dune::All_All_Interface,Dune::ForwardCommunication);
    }

    return 0;
  } catch (std::exception& e) {
    std::cerr << "error: caught exception: " << e.what() << std::endl;
    std::abort();
  } catch (...) {
    std::cerr << "error: caught unknown exception" << std::endl;
    std::abort();
  }
}
