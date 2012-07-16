// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GENERICDATAHANDLE_HH
#define DUNE_PDELAB_GENERICDATAHANDLE_HH

#include <vector>

#include<dune/common/exceptions.hh>
#include <dune/grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

#include <dune/pdelab/gridfunctionspace/entitycontainerindexcache.hh>

namespace Dune {
  namespace PDELab {


    /// \brief implement a data handle with a grid function space
    /// \tparam GFS a grid function space
    /// \tparam V   a vector container associated with the GFS
    /// \tparam T   gather/scatter methods with argumemts buffer, and data
    template<typename GFS, typename V, typename GatherScatter, typename E = typename V::ElementType>
    class GFSDataHandle
      : public Dune::CommDataHandleIF<GFSDataHandle<GFS,V,GatherScatter,E>,E>
    {

    public:

      typedef E DataType;
      typedef typename GFS::Traits::SizeType size_type;

      GFSDataHandle(const GFS& gfs, V& v, GatherScatter gather_scatter)
        : _gfs(gfs)
        , _index_cache(gfs)
        , _local_view(v)
        , _gather_scatter(gather_scatter)
      {}

      //! returns true if data for this codim should be communicated
      bool contains(int dim, int codim) const
      {
        return _gfs.dataHandleContains(dim,codim);
      }

      //!  \brief returns true if size per entity of given dim and codim is a constant
      bool fixedsize(int dim, int codim) const
      {
        return _gfs.dataHandleFixedSize(dim,codim);
      }

      /*!  \brief how many objects of type DataType have to be sent for a given entity

        Note: Only the sender side needs to know this size.
      */
      template<typename Entity>
      size_type size(Entity& e) const
      {
        return _gfs.dataHandleSize(e);
      }

      //! \brief pack data from user to message buffer
      template<typename MessageBuffer, typename Entity>
      void gather(MessageBuffer& buff, const Entity& e) const
      {
        _index_cache.update(e);
        _local_view.bind(_index_cache);
        _gather_scatter.gather(buff,e,_local_view);
        _local_view.unbind();
      }

      /*! \brief unpack data from message buffer to user

        n is the number of objects sent by the sender
      */
      template<typename MessageBuffer, typename Entity>
      void scatter(MessageBuffer& buff, const Entity& e, size_type n)
      {
        _index_cache.update(e);
        _local_view.bind(_index_cache);
        _gather_scatter.scatter(buff,n,e,_local_view);
        _local_view.commit();
        _local_view.unbind();
      }

    private:

      typedef EntityContainerIndexCache<GFS> IndexCache;
      typedef typename V::template LocalView<IndexCache> LocalView;

      const GFS& _gfs;
      mutable IndexCache _index_cache;
      mutable LocalView _local_view;
      mutable GatherScatter _gather_scatter;

    };


    template<typename GatherScatter>
    class DataGatherScatter
    {

    public:

      typedef std::size_t size_type;

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void gather(MessageBuffer& buff, const Entity& e, const LocalView& local_view) const
      {
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.gather(buff,local_view[i]);
      }

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void scatter(MessageBuffer& buff, size_type n, const Entity& e, LocalView& local_view) const
      {
        if (local_view.size() != n)
          DUNE_THROW(Exception,"size mismatch in generic data handle");
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.scatter(buff,local_view[i]);
      }

      DataGatherScatter(GatherScatter gather_scatter)
        : _gather_scatter(gather_scatter)
      {}

    private:

      GatherScatter _gather_scatter;

    };


    template<typename GatherScatter>
    class DataEntityGatherScatter
    {

    public:

      typedef std::size_t size_type;

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void gather(MessageBuffer& buff, const Entity& e, const LocalView& local_view) const
      {
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.gather(buff,e,local_view[i]);
      }

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void scatter(MessageBuffer& buff, size_type n, const Entity& e, LocalView& local_view) const
      {
        if (local_view.size() != n)
          DUNE_THROW(Exception,"size mismatch in generic data handle");
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.scatter(buff,e,local_view[i]);
      }

      DataEntityGatherScatter(GatherScatter gather_scatter)
        : _gather_scatter(gather_scatter)
      {}

    private:

      GatherScatter _gather_scatter;

    };


    template<typename GatherScatter>
    class DataContainerIndexGatherScatter
    {

    public:

      typedef std::size_t size_type;

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void gather(MessageBuffer& buff, const Entity& e, const LocalView& local_view) const
      {
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.gather(buff,local_view.cache().container_index(i),local_view[i]);
      }

      template<typename MessageBuffer, typename Entity, typename LocalView>
      void scatter(MessageBuffer& buff, size_type n, const Entity& e, LocalView& local_view) const
      {
        if (local_view.size() != n)
          DUNE_THROW(Exception,"size mismatch in generic data handle");
        for (std::size_t i = 0; i < local_view.size(); ++i)
          _gather_scatter.scatter(buff,local_view.cache().container_index(i),local_view[i]);
      }

      DataContainerIndexGatherScatter(GatherScatter gather_scatter)
        : _gather_scatter(gather_scatter)
      {}

    private:

      GatherScatter _gather_scatter;

    };


	class AddGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
      void gather (MessageBuffer& buff, DataType& data) const
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
      void scatter (MessageBuffer& buff, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
		data += x;
	  }
	};

    template<class GFS, class V>
	class AddDataHandle
      : public GFSDataHandle<GFS,V,DataGatherScatter<AddGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataGatherScatter<AddGatherScatter> > BaseT;

	public:

      AddDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataGatherScatter<AddGatherScatter>(AddGatherScatter()))
	  {}
	};

	class AddClearGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
      void gather (MessageBuffer& buff, DataType& data) const
	  {
		buff.write(data);
        data = (DataType) 0;
	  }
	  
	  template<class MessageBuffer, class DataType>
      void scatter (MessageBuffer& buff, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
		data += x;
	  }
	};
	
	template<class GFS, class V>
	class AddClearDataHandle
      : public GFSDataHandle<GFS,V,DataGatherScatter<AddClearGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataGatherScatter<AddClearGatherScatter> > BaseT;

	public:

      AddClearDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataGatherScatter<AddClearGatherScatter>(AddClearGatherScatter()))
	  {}
	};

	class CopyGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
      void gather (MessageBuffer& buff, DataType& data) const
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
      void scatter (MessageBuffer& buff, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
		data = x;
	  }
	};
	
	template<class GFS, class V>
	class CopyDataHandle
      : public GFSDataHandle<GFS,V,DataGatherScatter<CopyGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataGatherScatter<CopyGatherScatter >> BaseT;

	public:

      CopyDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataGatherScatter<CopyGatherScatter>(CopyGatherScatter()))
	  {}
	};

	class MinGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
      void gather (MessageBuffer& buff, DataType& data) const
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
      void scatter (MessageBuffer& buff, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
		data = std::min(data,x);
	  }
	};
	
	template<class GFS, class V>
	class MinDataHandle
      : public GFSDataHandle<GFS,V,DataGatherScatter<MinGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataGatherScatter<MinGatherScatter> > BaseT;

	public:

      MinDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataGatherScatter<MinGatherScatter>(MinGatherScatter()))
	  {}
	};

	class MaxGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
      void gather (MessageBuffer& buff, DataType& data) const
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
      void scatter (MessageBuffer& buff, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
		data = std::max(data,x);
	  }
	};
	
	template<class GFS, class V>
	class MaxDataHandle
      : public GFSDataHandle<GFS,V,DataGatherScatter<MaxGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataGatherScatter<MaxGatherScatter> > BaseT;

	public:

      MaxDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataGatherScatter<MaxGatherScatter>(MaxGatherScatter()))
	  {}
	};

    // assign degrees of freedoms to processors
    // owner is never a ghost
    class PartitionGatherScatter
    {
    public:
      template<class MessageBuffer, class EntityType, class DataType>
      void gather (MessageBuffer& buff, const EntityType& e, DataType& data) const
      {
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = (1<<24);
        buff.write(data);
      }
  
      template<class MessageBuffer, class EntityType, class DataType>
      void scatter (MessageBuffer& buff, const EntityType& e, DataType& data) const
      {
		DataType x; 
		buff.read(x);
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = x;
        else
          data = std::min(data,x);
      }
    };

	template<class GFS, class V>
	class PartitionDataHandle
      : public GFSDataHandle<GFS,V,DataEntityGatherScatter<PartitionGatherScatter> >
	{
      typedef GFSDataHandle<GFS,V,DataEntityGatherScatter<PartitionGatherScatter> > BaseT;

	public:

      PartitionDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataEntityGatherScatter<PartitionGatherScatter>(PartitionGatherScatter()))
	  {
        v_ = gfs_.gridView().comm().rank();
      }
	};

    // compute dofs assigned to ghost entities
	class GhostGatherScatter
	{
	public:
      template<class MessageBuffer, class EntityType, class DataType>
      void gather (MessageBuffer& buff, const EntityType& e, DataType& data) const
	  {
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = 1;
        buff.write(data);
	  }
	  
      template<class MessageBuffer, class EntityType, class DataType>
      void scatter (MessageBuffer& buff, const EntityType& e, DataType& data) const
	  {
		DataType x; 
		buff.read(x);
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = 1;
	  }
	};
	
	template<class GFS, class V>
	class GhostDataHandle
      : public Dune::PDELab::GFSDataHandle<GFS,V,DataEntityGatherScatter<GhostGatherScatter> >
	{
      typedef Dune::PDELab::GFSDataHandle<GFS,V,DataEntityGatherScatter<GhostGatherScatter> > BaseT;

	public:

      GhostDataHandle (const GFS& gfs_, V& v_)
        : BaseT(gfs_,v_,DataEntityGatherScatter<GhostGatherScatter>(GhostGatherScatter()))
	  {
        v_ = static_cast<typename V::ElementType>(0);
      }
	};


  }
}
#endif
