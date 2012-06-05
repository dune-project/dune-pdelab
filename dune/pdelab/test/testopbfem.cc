// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>

#include <dune/pdelab/finiteelementmap/opbfem.hh>


namespace Dune {

  namespace PB {

    //! \brief A class to test orthogonality by computing the mass matrix
    template<typename FieldType, int k, int d, Dune::GeometryType::BasicType bt, typename ComputationFieldType=FieldType>
    class MassMatrixTest
    {
    public:
      static void compute ()
      {
        // make basis
        typedef Dune::PB::OrthonormalPolynomialBasis<FieldType,k,d,bt,ComputationFieldType> PolynomialBasis;
        PolynomialBasis polynomialbasis;
        Dune::GeometryType gt(bt,d);
        std::cout << "mass matrix test"
                  << " k=" << k
                  << " on " << gt
                  << " n=" << PolynomialBasis::n;

        // select quadrature rule
        const Dune::QuadratureRule<FieldType,d>& rule = Dune::QuadratureRules<FieldType,d>::rule(gt,2*k);

        // make mass matrix
        Dune::FieldMatrix<FieldType,PolynomialBasis::n,PolynomialBasis::n> massmatrix;
        for (int i=0; i<PolynomialBasis::n; i++)
          for (int j=0; j<PolynomialBasis::n; j++)
            massmatrix[i][j] = 0.0;

        // loop over quadrature points
        for (typename Dune::QuadratureRule<FieldType,d>::const_iterator
               it=rule.begin(); it!=rule.end(); ++it)
          {
            typedef Dune::FieldVector<FieldType,1> RangeType;
            std::vector<RangeType> phi(PolynomialBasis::n);
            polynomialbasis.evaluateFunction(it->position(),phi);

            for (int i=0; i<PolynomialBasis::n; i++)
              for (int j=0; j<PolynomialBasis::n; j++)
                massmatrix[i][j] += phi[j]*phi[i]*it->weight();
          }

        // check mass matrix
        FieldType error=0.0;
        for (int i=0; i<PolynomialBasis::n; i++)
          for (int j=0; j<PolynomialBasis::n; j++)
            if (i==j)
              error = std::max(error,std::abs(massmatrix[i][j]-1.0));
            else
              error = std::max(error,std::abs(massmatrix[i][j]));
        std::cout << " maxerror=" << error << std::endl;
      }
    };

    void testmassmatrix()
    {
      MassMatrixTest<double,1,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,2,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,3,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,4,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,5,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,6,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,7,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,8,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,9,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,10,1,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();

      MassMatrixTest<double,1,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,2,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,3,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,4,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,5,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,6,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,7,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,8,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,9,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,10,2,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();

      MassMatrixTest<double,1,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,2,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,3,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,4,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,5,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,6,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,7,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,8,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,9,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,10,3,Dune::GeometryType::cube,Dune::GMPField<512> >::compute();

      MassMatrixTest<double,1,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,2,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,3,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,4,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,5,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,6,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,7,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,8,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,9,2,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();

      MassMatrixTest<double,1,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,2,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,3,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,4,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,5,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,6,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
      MassMatrixTest<double,7,3,Dune::GeometryType::simplex,Dune::GMPField<512> >::compute();
    }

  } // namespace PB

} // namespace Dune


int main(int argc, char** argv)
{
  try{
    // run tests
    Dune::PB::testmassmatrix();

    // test passed
    return 0;
  }
  catch (Dune::Exception &e){
    std::cerr << "Dune reported error: " << e << std::endl;
    return 1;
  }
  catch (...){
    std::cerr << "Unknown exception thrown!" << std::endl;
    return 1;
  }
}
