/** \brief A function that defines Dirichlet boundary conditions AND
    its extension to the interior */
template<typename GV, typename RF>
class BCExtension
  : public Dune::PDELab::GridFunctionBase<Dune::PDELab::
           GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> >, BCExtension<GV,RF> > {
  const GV& gv;
public:
  typedef Dune::PDELab::GridFunctionTraits<GV,RF,1,Dune::FieldVector<RF,1> > Traits;

  //! construct from grid view
  BCExtension (const GV& gv_) : gv(gv_) {}

  //! evaluate extended function on element
  inline void evaluate (const typename Traits::ElementType& e,
                        const typename Traits::DomainType& xlocal,
                        typename Traits::RangeType& y) const
  {
    const int dim = Traits::GridViewType::Grid::dimension;
    typedef typename Traits::GridViewType::Grid::ctype ctype;
    Dune::FieldVector<ctype,dim> x = e.geometry().global(xlocal);
    if (x[0]<1E-6 && x[1]>0.25 && x[1]<0.5) y = 1.0; else y = 0.0;
    return;
  }

  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};
