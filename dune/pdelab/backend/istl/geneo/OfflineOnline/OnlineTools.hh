#ifndef DUNE_PDELAB_ONLINE_TOOLS_HH
#define DUNE_PDELAB_ONLINE_TOOLS_HH

template<typename T>
T FromOffline(std::string path_to_storage, std::string object, int subdomain_number=-1) {
  T output;

  std::string filename = path_to_storage + std::to_string(subdomain_number) + "_" + object + ".mm";
  if(subdomain_number<0)
    filename = path_to_storage + object + ".mm";
  std::ifstream file;
  file.open(filename.c_str(), std::ios::in);
  Dune::readMatrixMarket(output,file);
  file.close();

  return output;
}

template<typename Vector>
std::vector<int> offlineDoF2GI2gmsh2onlineDoF(int subdomain_ID, std::vector<int>& gmsh2dune, Vector& offlineDoF2GI, const std::string path_to_storage, int verbose=0){

  Vector gmsh2GI = FromOffline<Vector>(path_to_storage, "LocalToGlobalNode", subdomain_ID);

  assert(gmsh2GI.size()==offlineDoF2GI.size());

  int v_size = gmsh2GI.size();

  // /* First method to deal with online2GI and offline2GI */
    // std::vector<int> GI2gmsh;
    // auto result = std::max_element(online2GI.begin(), online2GI.end());
    // GI2gmsh.resize(online2GI[std::distance(online2GI.begin(), result)]);
    // for(int i=0; i<GI2gmsh.size(); i++)
    //   GI2gmsh[i]=-1;
    // for(int i=0; i<v_size; i++){
    //   GI2gmsh[online2GI[i]] = i;
    // }
  // /* End first method */

  /* Second method to deal with online2GI and offline2GI */
  std::vector<std::pair<int, int>> gmsh(v_size), offlineDoF(v_size);
  for(int i=0; i<v_size; i++){
    gmsh[i] = std::make_pair(gmsh2GI[i],i);
    offlineDoF[i] = std::make_pair(offlineDoF2GI[i],i);
  }
  std::sort(gmsh.begin(), gmsh.end());
  std::sort(offlineDoF.begin(), offlineDoF.end());
  /* End second method */

  std::vector<int> indices_change(v_size);
  // Change order from Dof offline to offline global
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    // indices_change[i] = GI2gmsh[offlineDoF[i]]; // First method
    indices_change[offlineDoF[i].second] = gmsh[i].second; // second method
  }

  auto tmp = indices_change;
  // We know that a re-ordering of nodes has been done in the GmshReader function
  for(int i=0; i<v_size; i++){ // go through dof offline ordering
    indices_change[i] = gmsh2dune[tmp[i]];
  }

  if(verbose>0){
    /* Visualise changes in terminal */
    for(int i=0; i<v_size; i++){ // go through dof offline ordering
      std::cout << i << " -> " << indices_change[i] << std::endl;
    }
  }

  return indices_change;
}

template<typename Vector, typename Matrix, typename CoarseMatrix, typename vector1i>
void UpdateAH(CoarseMatrix& AH, Matrix A, std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>> online_subdomainbasis, std::vector<std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>> neighbour_subdomainbasis, std::vector<int> local_basis_sizes, std::vector<int> local_offset, vector1i NR, int subdomain_number) {

  typedef typename Matrix::field_type field_type;

  int my_offset = local_offset[subdomain_number];
  int max_local_basis_size = *std::max_element(local_basis_sizes.begin(),local_basis_sizes.end());

  // ~~~~~~~~~~~~~~~~~~
  //  Create a vector of AH entries modification
  // ~~~~~~~~~~~~~~~~~~

  // Set up container for storing rows of coarse matrix associated with current rank
  std::vector<std::vector<std::vector<field_type> > > local_rows;
  local_rows.resize(local_basis_sizes[subdomain_number]);
  for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
    local_rows[basis_index].resize(NR.size()+1);
  }

  for (int basis_index_remote = 0; basis_index_remote < max_local_basis_size; basis_index_remote++) {

    // Compute local products of basis functions with discretization matrix
    if (basis_index_remote < local_basis_sizes[subdomain_number]) {
      auto basis_vector = *online_subdomainbasis->get_basis_vector(basis_index_remote);
      Vector Atimesv(A.N());
      A.mv(basis_vector, Atimesv);
      for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
        field_type entry = *online_subdomainbasis->get_basis_vector(basis_index)*Atimesv;
        local_rows[basis_index][NR.size()].push_back(entry);
      }
    }

    // Compute products of discretization matrix with local and remote vectors
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      if (basis_index_remote >= local_basis_sizes[NR[neighbor_id]])
        continue;
      auto basis_vector = *neighbour_subdomainbasis[neighbor_id]->get_basis_vector(basis_index_remote);
      Vector Atimesv(A.N());
      A.mv(basis_vector, Atimesv);
      for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
        field_type entry = *online_subdomainbasis->get_basis_vector(basis_index)*Atimesv;
        local_rows[basis_index][neighbor_id].push_back(entry);
      }
    }
  }

  // ~~~~~~~~~~~~~~~~~~
  //  Modify AH entries
  // ~~~~~~~~~~~~~~~~~~

  int row_id = my_offset;
  // Modify AH entries with just computed local_rows
  for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
    // Communicate number of entries in this row
    int couplings = local_basis_sizes[subdomain_number];
    for (int neighbor_id : NR) {
      couplings += local_basis_sizes[neighbor_id];
    }

    // Communicate row's pattern
    int entries_pos[couplings];
    int cnt = 0;
    for (int basis_index2 = 0; basis_index2 < local_basis_sizes[subdomain_number]; basis_index2++) {
      entries_pos[cnt] = my_offset + basis_index2;
      cnt++;
    }
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      int neighbor_offset = local_offset[NR[neighbor_id]];
      for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
        entries_pos[cnt] = neighbor_offset + basis_index2;
        cnt++;
      }
    }

    // Communicate actual entries
    field_type entries[couplings];
    cnt = 0;
    for (int basis_index2 = 0; basis_index2 < local_basis_sizes[subdomain_number]; basis_index2++) {
      entries[cnt] = local_rows[basis_index][NR.size()][basis_index2];
      cnt++;
    }
    for (std::size_t neighbor_id = 0; neighbor_id < NR.size(); neighbor_id++) {
      for (int basis_index2 = 0; basis_index2 < local_basis_sizes[NR[neighbor_id]]; basis_index2++) {
        entries[cnt] = local_rows[basis_index][neighbor_id][basis_index2];
        cnt++;
      }
    }

    // Set matrix entries
    for (int i = 0; i < couplings; i++){
      AH[row_id][entries_pos[i]] = entries[i];
      AH[entries_pos[i]][row_id] = entries[i];
    }
    row_id++;
  }
}

template<typename Vector, typename CoarseVector>
void UpdatebH(CoarseVector& bH, Vector& b, std::shared_ptr<Dune::PDELab::SubdomainBasis<Vector>>& subdomainbasis, std::vector<int> local_basis_sizes, std::vector<int> local_offset, int subdomain_number) {

  typedef typename Vector::field_type field_type;

  int my_offset = local_offset[subdomain_number];
  field_type buf_defect_local[local_basis_sizes[subdomain_number]];

  for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
    buf_defect_local[basis_index] = 0.0;
    for (std::size_t i = 0; i < b.N(); i++)
      buf_defect_local[basis_index] += (*subdomainbasis->get_basis_vector(basis_index))[i] * b[i];
  }

  for (int basis_index = 0; basis_index < local_basis_sizes[subdomain_number]; basis_index++) {
    bH[basis_index+my_offset] = buf_defect_local[basis_index];
  }
}

#endif // DUNE_PDELAB_ONLINE_TOOLS_HH
