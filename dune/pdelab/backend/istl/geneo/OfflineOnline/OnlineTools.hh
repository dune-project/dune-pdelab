#ifndef DUNE_PDELAB_ONLINE_TOOLS_HH
#define DUNE_PDELAB_ONLINE_TOOLS_HH

template<typename Vector>
std::vector<int> offlineDoF2GI2gmsh2onlineDoF(int subdomain_ID, std::vector<int>& gmsh2dune, Vector& offlineDoF2GI, const std::string path_to_storage, int verbose=0){

  Vector gmsh2GI;
  std::string filename_gmsh2GI = path_to_storage + std::to_string(subdomain_ID) + "_LocalToGlobalNode.mm";
  std::ifstream file_gmsh2GI;
  file_gmsh2GI.open(filename_gmsh2GI.c_str(), std::ios::in);
  Dune::readMatrixMarket(gmsh2GI,file_gmsh2GI);
  file_gmsh2GI.close();

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


#endif // DUNE_PDELAB_ONLINE_TOOLS_HH
