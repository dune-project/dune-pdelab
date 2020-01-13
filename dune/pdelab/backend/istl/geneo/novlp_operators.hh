 /**
 * \brief communication data handle for adding data
 *
 * V vector container. Elements of this vector are sent around. Shoul work for std::vector and ISTL::BlockVector
 */
template<typename GV, typename Vector>
class VectorAddDataHandle
  : public Dune::CommDataHandleIF<VectorAddDataHandle<GV,Vector>,typename Vector::block_type>
{
  const typename GV::IndexSet& indexset;
  Vector& v;

public:
  typedef typename Vector::value_type DataType;

  VectorAddDataHandle (const GV& gv, Vector& v_)
    : indexset(gv.indexSet()), v(v_)
  {}

  bool contains (int dim, int codim) const
  {
    return (codim==dim);
  }

  bool fixedSize (int dim, int codim) const
  {
    return true;
  }

  template<class EntityType>
  size_t size (const EntityType& e) const
  {
    return 1;
  }

  template<class MessageBufferImp, class EntityType>
  void gather (MessageBufferImp& buff, const EntityType& e) const
  {
    buff.write(v[indexset.index(e)]);
  }

  template<class MessageBufferImp, class EntityType>
  void scatter (MessageBufferImp& buff, const EntityType& e, size_t n)
  {
    DataType x;
    buff.read(x);
    v[indexset.index(e)] += x;
  }
};

//! Operator for the non-overlapping parallel case
/**
 * Calculate \f$y:=Ax\f$.
 *
 * \tparam GV  Grid View
 * \tparam Matrix   Type of the matrix.  Should be one of the ISTL matrix types.
 * \tparam Vector   Type of the vectors the matrix is applied to.
 */
template<typename GV, typename Matrix, typename Vector>
class NonoverlappingOperator
  : public Dune::LinearOperator<Vector,Vector>
{
  const GV& gv;

public:
  //! export type of matrix
  using matrix_type = Matrix;
  //! export type of vectors the matrix is applied to
  using domain_type = Vector;
  //! export type of result vectors
  using range_type = Vector;
  //! export type of the entries for x
  typedef typename Vector::field_type field_type;

  //! need to have this method now
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::nonoverlapping;
  }

  //! Construct a non-overlapping operator
  /**
   * \param gfs_ GridFunctionsSpace for the vectors.
   * \param A    Matrix for this operator.  This should be the locally
   *             assembled matrix.
   *
   * \note The constructed object stores references to all the objects
   *       given as parameters here.  They should be valid for as long as
   *       the constructed object is used.  They are not needed to
   *       destruct the constructed object.
   */
  NonoverlappingOperator (const GV& gv_, const Matrix& A_)
    : gv(gv_), A(A_)
  { }

  //! apply operator
  /**
   * Compute \f$y:=A(x)\f$ on this process.
   * It is assumed that x is consistent, A is additive
   * then y is additive after multiplication and is immediately made consistent
   */
  virtual void apply (const Vector& x, Vector& y) const
  {
    A.mv(x,y); // A is additive and x is consistent; produces additive result
  }

  //! apply operator to x, scale and add:  \f$ y = y + \alpha A(x) \f$
  /**
   * Compute \f$y += \alpha A(x)\f$ on this process
   * It is assumed that x is consistent, A is additive and y is consistent
   */
  virtual void applyscaleadd (field_type alpha, const Vector& x, Vector& y) const
  {
    A.usmv(alpha,x,y); // A is additive and x is consistent; produces additive result
  }

private:
  const Matrix& A;
};


// parallel scalar product assuming no overlap
template<class GV, class Vector>
class NonoverlappingScalarProduct : public Dune::ScalarProduct<Vector>
{
public:
  //! export types
  typedef Vector domain_type;
  typedef typename Vector::field_type field_type;

  //! define the category
  virtual Dune::SolverCategory::Category category() const
  {
    return Dune::SolverCategory::nonoverlapping;
  }

  /*! \brief Constructor needs to know the grid view
   */
  NonoverlappingScalarProduct (const GV& gv_, const Vector& x)
    : gv(gv_)
  {
  }

  /*! \brief Dot product of two vectors.
   * It is assumed that one vector is consistent and the other is additive
   */
  virtual field_type dot (const Vector& x, const Vector& y)
  {
    field_type sum = 0.0;
    for (typename Vector::size_type i=0; i<x.N(); i++) sum += x[i]*y[i];
    auto sumsum = gv.comm().sum(sum);
    return sumsum;
  }

  /*! \brief Norm of a right-hand side vector.
   * It is assumed that x is additive.
   * This operation requires a local communication in addition
   * to make one argument consistent.
   */
  virtual double norm (const Vector& x)
  {
    // if (!check_vector_isfinite(x))
    //   std::cout << gv.comm().rank() << ": NaN in x detected in NonoverlappingScalarProduct.norm" << std::endl;
    Vector y(x);
    VectorAddDataHandle<GV,Vector> adddh(gv,y);
    if (gv.comm().size()>1)
      gv.communicate(adddh,Dune::All_All_Interface,Dune::ForwardCommunication);
    auto sp = static_cast<double>(this->dot(x,y));
    auto rv = std::sqrt(std::abs(sp)); // due to roundoff this may become negative for close to zero norms
    return rv;
  }

private:
  const GV gv;
};
