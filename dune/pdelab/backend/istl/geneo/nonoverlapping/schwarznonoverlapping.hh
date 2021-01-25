// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_ISTL_SCHWARZNONOVERLAPPING_GRIDS_HH
#define DUNE_ISTL_SCHWARZNONOVERLAPPING_GRIDS_HH

// dune-istl includes
#include<dune/istl/bvector.hh>
#include<dune/istl/bcrsmatrix.hh>
#include<dune/istl/operators.hh>
#include<dune/istl/preconditioners.hh>
#include<dune/istl/solvers.hh>
#include<dune/istl/io.hh>
#if HAVE_SUITESPARSE_UMFPACK
#include<dune/istl/umfpack.hh>
#endif

// dune-common includes
#include"overlaptools.hh"


namespace Dune {

#if HAVE_SUITESPARSE_UMFPACK
    template<typename Matrix, typename Vector, typename GridView>
    class NonoverlappingSchwarzPreconditioner : public Dune::Preconditioner<Vector,Vector>
    {
      // some types we need
      using GlobalId = typename GridView::Grid::GlobalIdSet::IdType;
      using CollectiveCommunication = typename GridView::CollectiveCommunication;
      using EPIS = Dune::ExtendedParallelIndexSet<CollectiveCommunication,GlobalId,Matrix>;
      using Attribute = EPISAttribute;
      using AttributedLocalIndex = Dune::ParallelLocalIndex<Attribute>;
      using ParallelIndexSet = Dune::ParallelIndexSet<GlobalId,AttributedLocalIndex,256>;
      using RemoteIndices = Dune::RemoteIndices<ParallelIndexSet>;
      using ExactSolver = Dune::UMFPack<Matrix>;
      using LocalIndex = typename Vector::size_type;
      using FieldType = typename Vector::field_type;
      using VectorBlockType = typename Vector::block_type;
      using MatrixBlockType = typename Matrix::block_type;
      using ScalarVector = Dune::BlockVector<Dune::FieldVector<FieldType,1>>;

      enum {components = Vector::block_type::dimension};

      // data
      const CollectiveCommunication& comm; // collective communication object
      int rank;     // my rank
      int p;        // number of subdomains in total
      int overlapsize;// overlap in terms of graph distance
      size_t N_orig; // number of degrees of freedom in the input
      size_t N_prec; // number of degrees of freedom of vectors used in preconditioner
      std::shared_ptr<ScalarVector> pu;  // improved partition of unity
      std::shared_ptr<ExactSolver> subdomainsolver; // umfpack solver for subdomain problem
      std::shared_ptr<ExactSolver> coarsesolver;    // umfpack solver for coarse problem (solved by every rank
      typename Vector::field_type* floating; // vector with size p storing 1 if subdomain is floating and 0 else
      bool coarsespace; // true if we are working with a coarse space
      //std::vector<LocalIndex> new2old_localindex; // maps new local index to a local index in the original input
      std::shared_ptr<ParallelIndexSet> pis; // newly built up index set
      std::shared_ptr<RemoteIndices> si; // shared index information
      std::shared_ptr<Dune::Interface> allinterface; // communication interface
      std::shared_ptr<Dune::BufferedCommunicator> communicator; // build all to all communicator once

      template<typename V>
      struct AddGatherScatter
      {
        static typename V::value_type gather(const V& a, int i)
        {
          return a[i]; // I am sending my value
        }
        static void scatter (V& a, typename V::value_type v, int i)
        {
          a[i]+=v; // add what I receive to my value
        }
      };

      std::shared_ptr<NonoverlappingOverlapAdapter<GridView, Vector, Matrix>> adapter;

    public:
      //! \brief The domain type of the preconditioner.
      typedef Vector domain_type;
      //! \brief The range type of the preconditioner.
      typedef Vector range_type;
      //! \brief The field type of the preconditioner.
      typedef typename Vector::field_type field_type;

      // define the category
      virtual Dune::SolverCategory::Category category() const
      {
        return Dune::SolverCategory::nonoverlapping;
      }

      //! \brief Constructor.
      NonoverlappingSchwarzPreconditioner (const CollectiveCommunication& comm_,    // collective communication object
                                           const GridView& gv,                         // matrix assembled on interior elements and with dirichlet rows replaced
                                           const std::vector<int>& neighbors,       // ranks with whom we share degrees of freedom
                                           const Matrix& A,
                                           bool floatingSubdomain,                  // true if this subdomain has no connection to the Dirichlet boundary
                                           const std::vector<EPISAttribute>& partitiontype,  // vector giving partitiontype for each degree of freedom
                                           const std::vector<GlobalId>& globalid,   // maps local index to globally unique id
                                           int avg_,                                // average number of matrix entries
                                           int overlap_,                            // layers of overlap to add
                                           bool coarsespace_)                       // set to true if a coarse space should be built and used
      : //gv(gv_),
        comm(comm_),
        rank(comm.rank()),
        p(comm.size()),
        overlapsize(overlap_),
        N_orig(A.N()),
        coarsespace(coarsespace_)
      {
        // handle the case p=1 so we do not have to care anymore
        if (p==1)
          {
            N_prec = N_orig;
            subdomainsolver = std::shared_ptr<ExactSolver>(new ExactSolver(A,0));
            pu = std::shared_ptr<ScalarVector>(new ScalarVector(N_prec));
            for (typename ScalarVector::size_type i=0; i<pu->N(); i++)
              (*pu)[i] = 1.0;
            return;
          }

        //adapter = std::make_shared<NonoverlappingOverlapAdapter<Vector, CollectiveCommunication, GlobalId, Matrix>>(comm_, neighbors, A, partitiontype, globalid, avg_, overlap_);
        adapter = std::make_shared<NonoverlappingOverlapAdapter<GridView, Vector, Matrix>>(gv, A, avg_, overlap_);

        N_prec = adapter->getExtendedSize(); // this now the new size after adding overlap

        // construct stiffness matrix on overlapping set
        std::shared_ptr<Matrix> M = adapter->extendMatrix(A);

        // compute LU decomposition
        subdomainsolver = std::shared_ptr<ExactSolver>(new ExactSolver(*M,0));

        // construct partition of unity in a scalar vector
        //pu = Dune::makePartitionOfUnity(epis,*M,overlapsize);
        pu = Dune::makePartitionOfUnity(*adapter, *M);
        for (typename ScalarVector::size_type i=0; i<pu->N(); i++)
          (*pu)[i] = std::sqrt((*pu)[i]); // we expect to store the square root

        // construct owner partition scalar vector
        auto owner = Dune::makeOwner(adapter->getEpis(),*M,overlapsize);


        // build up communication interface using all atributes (thankfully ghosts are out of the way :-)) for later use
        Dune::AllSet<Attribute> allAttribute;
        allinterface = std::shared_ptr<Dune::Interface>(new Dune::Interface());
        allinterface->build(*adapter->getRemoteIndices(),allAttribute,allAttribute); // all to all communication

        // build up buffered communicator allowing communication over a dof vector
        communicator = std::shared_ptr<Dune::BufferedCommunicator>(new Dune::BufferedCommunicator());
        communicator->build<Vector>(*allinterface);

        // now build the Nicolaides coarse space, but only if really needed
        if (!coarsespace) return;

        // build up rows of coarse grid matrix in O(p) algorithm
        MatrixBlockType* myrow = new MatrixBlockType[p]; // my row of coarse grid matrix
        for (int j=0; j<p; j++) // loop over columns of the coarse grid operator
          for (int J=0; J<components; J++)
            {
              // prepare column j, component J
              Vector r(N_prec);
              r = 0.0;
              if (rank==j)
                for (typename Vector::size_type i=0; i<r.N(); i++) r[i][J] = (*pu)[i] * (*pu)[i]; // initialize with partition of unity
              communicator->forward<AddGatherScatter<Vector>>(r,r); // make function known in other subdomains

              // multiply with M
              Vector z(N_prec);
              M->mv(r,z);
              for (typename Vector::size_type i=0; i<z.N(); i++) z[i] *= (*owner)[i]; // the product is correct in all owners
              communicator->forward<AddGatherScatter<Vector>>(z,z); // now everybody has correct values of A*r_j

              // now we have the whole column (j,J) of A in all processors and every one can compute the scalar product with its pu
              for (int I=0; I<components; I++)
                {
                  r = 0.0;
                  for (typename Vector::size_type i=0; i<r.N(); i++) r[i][I] = (*pu)[i] * (*pu)[i];
                  myrow[j][I][J] = r*z;
                }
            }

        // make a BCRSMatrix out of the rows
        int avg = 2.5*avg_;
        Matrix A0(p,p,avg,0.1,Matrix::implicit);
        MatrixBlockType* row = new MatrixBlockType[p];
        for (int i=0; i<p; i++)
          {
            if (rank==i) for (int j=0; j<p; j++) row[j] = myrow[j]; // send my row
            comm.broadcast(row,p,i);
            for (int j=0; j<p; j++)
              if (std::abs(row[j].frobenius_norm())>1e-15)
                A0.entry(i,j) = row[j];
          }
        auto stats0 = A0.compress();

        // make floating subdomains known to everyone
        floating = new field_type[p];
        for (int j=0; j<p; j++) floating[j] = 0.0;
        if (floatingSubdomain)
          floating[rank] = 1.0;
        comm.sum(floating,p);

        // eliminate non-floating subdomains from coarse system
        for (size_t i=0; i<A0.N(); i++)
          {
            auto cIt = A0[i].begin();
            auto cEndIt = A0[i].end();
            for (; cIt!=cEndIt; ++cIt)
              {
                auto j = cIt.index();
                if (i==j)
                  {
                    // this is a diagonal block
                    if (floating[i]==0.0)
                      for (int k=0; k<components; k++)
                        for (int l=0; l<components; l++)
                          (*cIt)[k][l] = (k==l) ? 1.0 : 0.0;
                  }
                else
                  {
                    // this is an offdiagonal block
                    if (floating[i]==0.0 || floating[j]==0.0) *cIt = 0.0;
                  }
              }
          }

        // compute LU decomposition of coarse matrix
        coarsesolver = std::shared_ptr<ExactSolver>(new ExactSolver(A0,0));

        // delete temporary vectors
        delete [] row;
        delete [] myrow;
      }

      template<typename T>
      void getPartitionOfUnity (T& vec)
      {
        adapter->restrictVector(*pu, vec);
        for (auto& val : vec)
          vec *= vec;
      }

      template<typename T>
      void getPartitionOfUnity (T& vec, int j)
      {
        // fill with PU in rank j
        ScalarVector x(N_prec);
        if (pu->N()!=x.size()) { std::cout << "size in getPartitionOfUnity not matching" << std::endl; exit(1); }
        if (rank==j)
          for (typename ScalarVector::size_type i=0; i<x.size(); i++)
            x[i] = (*pu)[i]*(*pu)[i];
        else
          x = 0.0;

        // add up to communicate to all others
        communicator->forward<AddGatherScatter<ScalarVector>>(x,x);

        // pick out on interior and border
        adapter->restrictVector(x, vec);
      }

      ~NonoverlappingSchwarzPreconditioner ()
      {
        if (coarsespace && p>1) delete [] floating;
      }

      /*!
        \brief Prepare the preconditioner.
      */
      virtual void pre (Vector& x, Vector& b)
      {
      }

      /*!
        \brief Apply the precondioner.
      */
      virtual void apply (Vector& v, const Vector& d)
      {
        // copy defect to a vector with the new local index set
        Vector d2(N_prec);
        adapter->extendVector(d, d2);
        if (p>1) communicator->forward<AddGatherScatter<Vector>>(d2,d2); // make it consistent

        // coarse grid correction
        Vector z(N_prec);
        z = 0.0;
        if (coarsespace && p>1)
          {
            // assemble rhs of coarse problem
            VectorBlockType* d0 = new VectorBlockType[p];
            for (int i=0; i<p; i++) d0[i] = 0.0;
            for (typename Vector::size_type i=0; i<N_prec; i++)
              d0[rank].axpy((*pu)[i]*(*pu)[i],d2[i]);
            comm.sum(d0,p);
            Vector R0d(p);
            for (int i=0; i<p; i++)
              {
                R0d[i] = d0[i];
                R0d[i] *= floating[i]; // zero rhs for non-floating domains
              }
            delete [] d0;

            // solve coarse problem
            Vector  v0(p);
            v0 = 0.0;
            Dune::InverseOperatorResult stat;
            coarsesolver->apply(v0,R0d,stat);
            // assemble correction from subdomains
            for (typename Vector::size_type i=0; i<z.N(); i++)
              z[i].axpy((*pu)[i]*(*pu)[i],v0[rank]);
          }

        Dune::InverseOperatorResult stat;
        Vector v2(N_prec);
        v2 = 0.0;
        subdomainsolver->apply(v2,d2,stat);
        if (coarsespace && p>1) v2 += z;

        // add up corrections
        if (p>1) communicator->forward<AddGatherScatter<Vector>>(v2,v2);

        // write back result to old index set
        std::cout << "size v: " << v.N() << " size v2: " << v2.N() << std::endl;
        adapter->restrictVector(v2, v);
      }

      /*!
        \brief Clean up.
      */
      virtual void post (Vector& x)
      {
      }

    };
#endif

} // namespace Dune

#endif // DUNE_ISTL_SCHWARZNONOVERLAPPING_GRIDS_HH
