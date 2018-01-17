#ifndef DUNE_GENEO_COARSESPACE_HH
#define DUNE_GENEO_COARSESPACE_HH

/*! \brief Representation of a coarse space intended for two-level Schwarz preconditioners.
 * \tparam M Matrix type on the subdomain
 * \tparam X Vector type on the subdomain
 */
template <class M, class X>
class CoarseSpace {

public:
  typedef Dune::BlockVector<Dune::FieldVector<double,1> > COARSE_V;
  typedef Dune::BCRSMatrix<Dune::FieldMatrix<double,1,1> > COARSE_M;

  /*! \brief Restricts a vector defined on a subdomain to the coarse space
   * \param[in] d The subdomain space vector to be restricted
   * \return The restriction in coarse space
   */
  virtual std::shared_ptr<COARSE_V> restrict_defect (const X& d) const = 0;

  /*! \brief Prolongates a vector defined on the coarse space to the subdomain
   * \param[in] v The coarse space vector to be prolongated
   * \return The prolongation in subdomain space
   */
  virtual std::shared_ptr<X> prolongate_defect (const COARSE_V& v) const = 0;

  /*! \brief Returns the matrix representing the coarse basis
   * \return The coarse matrix
   */
  virtual std::shared_ptr<COARSE_M> get_coarse_system () = 0;

  /*! \brief Returns the size of the coarse basis
   * \return Size of the basis
   */
  virtual int basis_size() = 0;
};

#endif //DUNE_GENEO_COARSESPACE_HH
