// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_PDELAB_GENERICDATAHANDLE_HH
#define DUNE_PDELAB_GENERICDATAHANDLE_HH

#include <vector>

#include<dune/common/exceptions.hh>
#include<grid/common/datahandleif.hh>
#include <dune/grid/common/gridenums.hh>

namespace Dune {
  namespace PDELab {

	// implement a data handle with a grid function space
	// GFS: a grid function space
	// V:   a vector container associated with the GFS
	// T:   gather/scatter methods
	template<class GFS, class V, class T>
	class GenericDataHandle
	  : public Dune::CommDataHandleIF<GenericDataHandle<GFS,V,T>,typename V::ElementType>
	{
	  typedef typename GFS::Traits::BackendType B;

	public:
	  GenericDataHandle (const GFS& gfs_, V& v_, T t_) 
		: gfs(gfs_), v(v_), t(t_), global(gfs.maxLocalSize())
	  {}

	  //! export type of data for message buffer
	  typedef typename V::ElementType DataType;

	  //! returns true if data for this codim should be communicated
	  bool contains (int dim, int codim) const
	  {
		return gfs.dataHandleContains(dim,codim);
	  }

	  //! returns true if size per entity of given dim and codim is a constant
	  bool fixedsize (int dim, int codim) const
	  {
		return gfs.dataHandleFixedSize(dim,codim);
	  }

	  /*! how many objects of type DataType have to be sent for a given entity

		Note: Only the sender side needs to know this size. 
	  */
	  template<class EntityType>
	  size_t size (EntityType& e) const
	  {
		return gfs.dataHandleSize(e);
	  }

	  //! pack data from user to message buffer
	  template<class MessageBuffer, class EntityType>
	  void gather (MessageBuffer& buff, const EntityType& e) const
	  {
		gfs.dataHandleGlobalIndices(e,global);
		for (size_t i=0; i<global.size(); ++i)
		  t.gather(buff,B::access(v,global[i]));
	  }

	  /*! unpack data from message buffer to user

		n is the number of objects sent by the sender
	  */
	  template<class MessageBuffer, class EntityType>
	  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	  {
		gfs.dataHandleGlobalIndices(e,global);
		if (global.size()!=n)
          DUNE_THROW(Exception,"size mismatch in generic data handle");
		for (size_t i=0; i<global.size(); ++i)
		  t.scatter(buff,B::access(v,global[i]));
	  }

	private:
	  const GFS& gfs;
	  V& v;
	  mutable T t;
	  mutable std::vector<typename GFS::Traits::SizeType> global;
	};


	// implement a data handle with a grid function space
	// GFS: a grid function space
	// V:   a vector container associated with the GFS
	// T:   gather/scatter methods
	template<class GFS, class V, class T>
	class GenericDataHandle2
	  : public Dune::CommDataHandleIF<GenericDataHandle2<GFS,V,T>,typename V::ElementType>
	{
	  typedef typename GFS::Traits::BackendType B;

	public:
	  GenericDataHandle2 (const GFS& gfs_, V& v_, T t_) 
		: gfs(gfs_), v(v_), t(t_), global(gfs.maxLocalSize())
	  {}

	  //! export type of data for message buffer
	  typedef typename V::ElementType DataType;

	  //! returns true if data for this codim should be communicated
	  bool contains (int dim, int codim) const
	  {
		return gfs.dataHandleContains(dim,codim);
	  }

	  //! returns true if size per entity of given dim and codim is a constant
	  bool fixedsize (int dim, int codim) const
	  {
		return gfs.dataHandleFixedSize(dim,codim);
	  }

	  /*! how many objects of type DataType have to be sent for a given entity

		Note: Only the sender side needs to know this size. 
	  */
	  template<class EntityType>
	  size_t size (EntityType& e) const
	  {
		return gfs.dataHandleSize(e);
	  }

	  //! pack data from user to message buffer
	  template<class MessageBuffer, class EntityType>
	  void gather (MessageBuffer& buff, const EntityType& e) const
	  {
		gfs.dataHandleGlobalIndices(e,global);
		for (size_t i=0; i<global.size(); ++i)
		  t.gather(buff,e,B::access(v,global[i]));
	  }

	  /*! unpack data from message buffer to user

		n is the number of objects sent by the sender
	  */
	  template<class MessageBuffer, class EntityType>
	  void scatter (MessageBuffer& buff, const EntityType& e, size_t n)
	  {
		gfs.dataHandleGlobalIndices(e,global);
		if (global.size()!=n)
          DUNE_THROW(Exception,"size mismatch in generic data handle");
		for (size_t i=0; i<global.size(); ++i)
		  t.scatter(buff,e,B::access(v,global[i]));
	  }

	private:
	  const GFS& gfs;
	  V& v;
	  mutable T t;
	  mutable std::vector<typename GFS::Traits::SizeType> global;
	};



	class AddGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
	  void gather (MessageBuffer& buff, DataType& data)
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
	  void scatter (MessageBuffer& buff, DataType& data)
	  {
		DataType x; 
		buff.read(x);
		data += x;
	  }
	};
	
	template<class GFS, class V>
	class AddDataHandle
	  : public GenericDataHandle<GFS,V,AddGatherScatter>
	{
	  typedef GenericDataHandle<GFS,V,AddGatherScatter> BaseT;

	public:

	  AddDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,AddGatherScatter())
	  {}
	};

	class CopyGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
	  void gather (MessageBuffer& buff, DataType& data)
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
	  void scatter (MessageBuffer& buff, DataType& data)
	  {
		DataType x; 
		buff.read(x);
		data = x;
	  }
	};
	
	template<class GFS, class V>
	class CopyDataHandle
	  : public GenericDataHandle<GFS,V,CopyGatherScatter>
	{
	  typedef GenericDataHandle<GFS,V,CopyGatherScatter> BaseT;

	public:

	  CopyDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,CopyGatherScatter())
	  {}
	};

	class MinGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
	  void gather (MessageBuffer& buff, DataType& data)
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
	  void scatter (MessageBuffer& buff, DataType& data)
	  {
		DataType x; 
		buff.read(x);
		data = std::min(data,x);
	  }
	};
	
	template<class GFS, class V>
	class MinDataHandle
	  : public GenericDataHandle<GFS,V,MinGatherScatter>
	{
	  typedef GenericDataHandle<GFS,V,MinGatherScatter> BaseT;

	public:

	  MinDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,MinGatherScatter())
	  {}
	};

	class MaxGatherScatter
	{
	public:
	  template<class MessageBuffer, class DataType>
	  void gather (MessageBuffer& buff, DataType& data)
	  {
		buff.write(data);
	  }
	  
	  template<class MessageBuffer, class DataType>
	  void scatter (MessageBuffer& buff, DataType& data)
	  {
		DataType x; 
		buff.read(x);
		data = std::max(data,x);
	  }
	};
	
	template<class GFS, class V>
	class MaxDataHandle
	  : public GenericDataHandle<GFS,V,MaxGatherScatter>
	{
	  typedef GenericDataHandle<GFS,V,MaxGatherScatter> BaseT;

	public:

	  MaxDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,MaxGatherScatter())
	  {}
	};

    // assign degrees of freedoms to processors
    // owner is never a ghost
    class PartitionGatherScatter
    {
    public:
      template<class MessageBuffer, class EntityType, class DataType>
      void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
      {
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = (1<<24);
        buff.write(data);
      }
  
      template<class MessageBuffer, class EntityType, class DataType>
      void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
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
	  : public GenericDataHandle2<GFS,V,PartitionGatherScatter>
	{
	  typedef GenericDataHandle2<GFS,V,PartitionGatherScatter> BaseT;

	public:

	  PartitionDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,PartitionGatherScatter())
	  {
        v_ = gfs_.gridview().comm().rank();
      }
	};

    // compute dofs assigned to ghost entities
	class GhostGatherScatter
	{
	public:
      template<class MessageBuffer, class EntityType, class DataType>
      void gather (MessageBuffer& buff, const EntityType& e, DataType& data)
	  {
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = 1;
        buff.write(data);
	  }
	  
      template<class MessageBuffer, class EntityType, class DataType>
      void scatter (MessageBuffer& buff, const EntityType& e, DataType& data)
	  {
		DataType x; 
		buff.read(x);
        if (e.partitionType()!=Dune::InteriorEntity && e.partitionType()!=Dune::BorderEntity)
          data = 1;
	  }
	};
	
	template<class GFS, class V>
	class GhostDataHandle
	  : public Dune::PDELab::GenericDataHandle2<GFS,V,GhostGatherScatter>
	{
	  typedef Dune::PDELab::GenericDataHandle2<GFS,V,GhostGatherScatter> BaseT;

	public:

	  GhostDataHandle (const GFS& gfs_, V& v_) 
		: BaseT(gfs_,v_,GhostGatherScatter())
	  {
        v_ = 0;
      }
	};


  }
}
#endif
