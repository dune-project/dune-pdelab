#ifndef DUNE_PDELAB_BACKEND_ISTL_GENEO_PARTITIONOFUNITY_HH
#define DUNE_PDELAB_BACKEND_ISTL_GENEO_PARTITIONOFUNITY_HH

/*!
 * \brief Compute a simple partition of unity.
 *
 * The resulting partition of unity is defined as 1/k on each DOF, where k is the
 * number of subdomains covering that DOF. Subdomain boundaries are not counted here
 * provided the constraints container has Dirichlet constraints on subdomain boundaries.
 *
 * \tparam X Vector type
 * \param gfs Grid function space.
 * \param cc Constraints container of the problem to be solved.
 * \return A vector representing the partition of unity.
 */
template<class X, class GFS, class CC>
X standardPartitionOfUnity(const GFS& gfs, const CC& cc) {

  X part_unity(gfs, 1);

  Dune::PDELab::set_constrained_dofs(cc,0.0,part_unity); // Zero on subdomain boundary

  Dune::PDELab::AddDataHandle<GFS,X> parth(gfs,part_unity);
  gfs.gridView().communicate(parth,Dune::All_All_Interface,Dune::ForwardCommunication);

  Dune::PDELab::set_constrained_dofs(cc,0.0,part_unity); // Zero on subdomain boundary (Need that a 2nd time due to add comm before!)

  for (auto iter = part_unity.begin(); iter != part_unity.end(); iter++) {
    if (*iter > 0)
      *iter = 1.0 / *iter;
  }
  return part_unity;
}

/*!
 * \brief Compute a partition of unity according to Sarkis.
 *
 * The resulting partition of unity interpolates linearly where two subdomains
 * overlap and extends consistently to more overlapping subdomains.
 * It is strictly bound to the subdomain arrangement as provided by YaspGrid
 * and only supports the 2D case.
 * The number of subdomains is assumed to divide the number of cells in each dimension.
 *
 * \tparam X Vector type
 * \param gfs Grid function space.
 * \param lsf Local function space.
 * \param cc Constraints container of the problem to be solved.
 * \param cells_x Number of cells in x direction.
 * \param cells_y Number of cells in y direction.
 * \param partition_x Number of partitions in x direction.
 * \param partition_y Number of partitions in y direction.
 * \return A vector representing the partition of unity.
 */
template<class X, class GFS,  class LFS, class CC>
X sarkisPartitionOfUnity(const GFS& gfs, LFS& lfs, const CC& cc, int cells_x, int cells_y, int overlap, int partition_x, int partition_y) {
  using Dune::PDELab::Backend::native;

  int my_rank = gfs.gridView().comm().rank();
  const int dim = 2;

  X part_unity(gfs, 1);

  for (auto it = gfs.gridView().template begin<0>(); it != gfs.gridView().template end<0>(); ++it) {

    lfs.bind(*it);

    auto geo  = it->geometry();
    const auto gt = geo.type();
    const auto& ref_el = Dune::ReferenceElements<double, dim>::general(gt);

    auto& coeffs = lfs.finiteElement().localCoefficients();

    for (std::size_t i = 0; i < coeffs.size(); ++i) {

      auto local_pos = ref_el.position (coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());

      auto global_pos = geo.global(local_pos);

      auto subindex = gfs.entitySet().indexSet().subIndex(*it, coeffs.localKey(i).subEntity(), coeffs.localKey(i).codim());

      double Hx = 1.0 / (double)partition_x;
      double Hy = 1.0 / (double)partition_y;
      double hx = (double)overlap / cells_x;
      double hy = (double)overlap / cells_y;

      int row = std::floor(my_rank / partition_x);
      int col = my_rank - partition_x * row;

      double dx1 = (col + 1) * Hx + hx - global_pos[0];
      double dx2 = global_pos[0] - (col * Hx - hx);

      double dy1 = (row + 1) * Hy + hy - global_pos[1];
      double dy2 = global_pos[1] - (row * Hy - hy);

      if (row == 0) dy2 = 2*Hy;
      if (row == partition_y - 1) dy1 = 2*Hy;
      if (col == 0) dx2 = 2*Hx;
      if (col == partition_x - 1) dx1 = 2*Hx;

      native(part_unity)[subindex] = std::min(std::min(std::min(dx1, dx2), dy1), dy2);
    }
  }

  X sum_dists(part_unity);
  Dune::PDELab::AddDataHandle<GFS,X> addh_dists(gfs,sum_dists);
  gfs.gridView().communicate(addh_dists,Dune::All_All_Interface,Dune::ForwardCommunication);

  auto iter_sum = sum_dists.begin();
  for (auto iter = part_unity.begin(); iter != part_unity.end(); iter++) {
    if (*iter > 0)
      *iter *= 1.0 / *iter_sum;
    iter_sum++;
  }

  Dune::PDELab::set_constrained_dofs(cc,0.0,part_unity); // Zero on Dirichlet domain boundary

  return part_unity;
}


#endif //DUNE_PDELAB_BACKEND_ISTL_GENEO_PARTITIONOFUNITY_HH
