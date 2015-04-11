#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <iostream>

#include <dune/common/math.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/yaspgrid.hh>

#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/vtkexport.hh>

#include <dune/pdelab/backend/istlvectorbackend.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>

// an analytic scalar function
template<typename T>
class F : public Dune::PDELab::FunctionInterface<
  Dune::PDELab::FunctionTraits<T,2,Dune::FieldVector<T,2>,
                               T,1,Dune::FieldVector<T,1> >,
  F<T> >
{
public:
  inline void evaluate (const Dune::FieldVector<T,2>& x,
                        Dune::FieldVector<T,1>& y) const
  {
    y = sin(3*Dune::StandardMathematicalConstants<T>::pi()*x[0])
      * cos(7*Dune::StandardMathematicalConstants<T>::pi()*x[1]);
  }
};

// an analytic vector-valued function
template<typename T>
class G : public Dune::PDELab::FunctionInterface<
  Dune::PDELab::FunctionTraits<T,2,Dune::FieldVector<T,2>,
                               T,2,Dune::FieldVector<T,2> >,
  G<T> >
{
public:
  inline void evaluate (const Dune::FieldVector<T,2>& x,
                        Dune::FieldVector<T,2>& y) const
  {
    y[0] =  x.two_norm()*x[1];
    y[1] = -x.two_norm()*x[0];
  }
};


// iterate over grid view and use analytic function as grid function through adapter
template<class GV, class T>
void testgridfunction (const GV& gv, const T& t)
{
  // make a grid function from the analytic function
  typedef Dune::PDELab::FunctionToGridFunctionAdapter<GV,T> GF;
  GF gf(gv,t);

  typedef typename GV::Traits::template Codim<0>::Iterator ElementIterator;
  for (ElementIterator it = gv.template begin<0>();
       it!=gv.template end<0>(); ++it)
    {
      typename T::Traits::DomainType x(0.0);
      typename T::Traits::RangeType y;
      gf.evaluate(*it,x,y);

      // make a function in local coordinates from the grid function
      Dune::PDELab::GridFunctionToLocalFunctionAdapter<GF>
        lf(gf,*it);
      lf.evaluate(x,y);

      // make a function in local coordinates from the global function
      Dune::PDELab::GlobalFunctionToLocalFunctionAdapter<T,typename ElementIterator::Entity>
        lg(t,*it);
      lg.evaluate(x,y);
    }
}

// iterate over grid view and use analytic function as grid function through adapter
template<class GV, class T>
void testvtkexport (const GV& gv, const T& t)
{
  // make a grid function from the analytic function
  typedef Dune::PDELab::FunctionToGridFunctionAdapter<GV,T> GF;
  GF gf(gv,t);

  // make a VTKFunction from grid function
  Dune::PDELab::VTKGridFunctionAdapter<GF> vtkf(gf,"blub");

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  vtkwriter.addVertexData(Dune::make_shared< Dune::PDELab::VTKGridFunctionAdapter<GF> >(gf,"blub")); // VTKWriter takes control
  vtkwriter.write("single",Dune::VTK::ascii);
}

// a grid function
template<typename G, typename T>
class L : public Dune::PDELab::GridFunctionBase<
  Dune::PDELab::GridFunctionTraits<G,T,1,Dune::FieldVector<T,1> >,
  L<G,T> >
{
  typedef typename G::Traits::template Codim<0>::Entity ElementType;
public:
  L (const G& g_) : g(g_) {}

  inline void evaluate (const Dune::FieldVector<T,2>& x,
                        Dune::FieldVector<T,1>& y) const
  {
    y = sin(3*3.1415*x[0])*cos(7*3.1415*x[1]);
  }

  inline void evaluate (const ElementType& e,
                        const Dune::FieldVector<T,2>& x,
                        Dune::FieldVector<T,1>& y) const
  {
    evaluate(e.geometry().global(x),y);
  }

  inline const G& getGridView () const
  {
    return g;
  }

private:
  const G& g;
};


// test function trees
template<class GV>
void testfunctiontree (const GV& gv)
{
  // a leaf function
  typedef L<GV,typename GV::Grid::ctype> A;
  A a(gv);

  // test power
  typedef Dune::PDELab::PowerGridFunction<A,10> B10;
  B10 b10(a);
  typedef Dune::PDELab::PowerGridFunction<A,2> B2;
  B2 b2(a,a);
  typedef Dune::PDELab::PowerGridFunction<A,3> B3;
  B3 b3(a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,4> B4;
  B4 b4(a,a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,5> B5;
  B5 b5(a,a,a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,6> B6;
  B6 b6(a,a,a,a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,7> B7;
  B7 b7(a,a,a,a,a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,8> B8;
  B8 b8(a,a,a,a,a,a,a,a);
  typedef Dune::PDELab::PowerGridFunction<A,9> B9;
  B9 b9(a,a,a,a,a,a,a,a,a);

  // test composite
  typedef Dune::PDELab::CompositeGridFunction<A,A,A,A,A,A,A,A,A> C9;
  C9 c9(a,a,a,a,a,a,a,a,a);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2> C2;
  C2 c2(b10,b2);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3> C3;
  C3 c3(b10,b2,b3);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3,B4> C4;
  C4 c4(b10,b2,b3,b4);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3,B4,B5> C5;
  C5 c5(b10,b2,b3,b4,b5);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3,B4,B5,B6> C6;
  C6 c6(b10,b2,b3,b4,b5,b6);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3,B4,B5,B6,B7> C7;
  C7 c7(b10,b2,b3,b4,b5,b6,b7);
  typedef Dune::PDELab::CompositeGridFunction<B10,B2,B3,B4,B5,B6,B7,B8> C8;
  C8 c8(b10,b2,b3,b4,b5,b6,b7,b8);

  typedef Dune::PDELab::CompositeGridFunction<C2,C9> T;
  T t(c2,c9);

  std::cout << "depth of T is " << Dune::TypeTree::TreeInfo<T>::depth << std::endl;
  std::cout << "number of nodes in T is " << Dune::TypeTree::TreeInfo<T>::nodeCount << std::endl;
  std::cout << "number of leaves in T is " << Dune::TypeTree::TreeInfo<T>::leafCount << std::endl;

  Dune::VTKWriter<GV> vtkwriter(gv,Dune::VTK::conforming);
  Dune::PDELab::vtkwriter_tree_addvertexdata(vtkwriter,t);
  vtkwriter.write("multi",Dune::VTK::ascii);
}

template<class GV>
void testgridviewfunction (const GV& gv)
{
    enum { dim = GV::dimension };
    typedef Dune::PDELab::QkLocalFiniteElementMap<GV,float,double,1> Q1FEM;
    Q1FEM q1fem(gv);
    // make a grid function space
    typedef Dune::PDELab::GridFunctionSpace<GV,Q1FEM> Q1GFS;
    Q1GFS q1gfs(gv,q1fem);
    // make vector
    typedef typename Dune::PDELab::BackendVectorSelector<Q1GFS,double>::Type Vector;
    Vector x(q1gfs);
    // make functions
    typedef Dune::PDELab::DiscreteGridViewFunction<Q1GFS,double> DiscreteFunction;
    DiscreteFunction dgvf(q1gfs,x);
    // make local functions
    typedef typename DiscreteFunction::LocalFunction LocalFunction;
    Dune::shared_ptr<LocalFunction>
        localf = Dune::Functions::localFunction(dgvf);
    // iterate grid and evaluate local function
    static const int maxDiffOrder = LocalFunction::Traits::maxDiffOrder;
    std::cout << "max diff order: " << maxDiffOrder << std::endl;
    std::cout << "checking for:\n";
    std::cout << "\tevaluate\n";
    if (maxDiffOrder >= 1)
        std::cout << "\tjacobian\n";
    if (maxDiffOrder >= 2)
        std::cout << "\thessian\n";
    if (maxDiffOrder >= 3)
        std::cout << "\tdiff(3)\n";
    for (auto it=gv.template begin<0>(); it!=gv.template end<0>(); ++it)
    {
        localf->bind(*it);
        Dune::FieldVector<double,1> value;
        Dune::FieldMatrix<double,1,dim> jacobian;
        Dune::FieldMatrix<double,dim,dim> hessian;
        localf->evaluate(it->geometry().center(), value);
        if (maxDiffOrder >= 1)
            Dune::Functions::derivative(localf)->
                evaluate(it->geometry().center(), jacobian);
        if (maxDiffOrder >= 2)
            Dune::Functions::derivative(
                Dune::Functions::derivative(localf))->
                evaluate(it->geometry().center(), hessian);
        if (maxDiffOrder >= 3)
            Dune::Functions::derivative(
                Dune::Functions::derivative(
                    Dune::Functions::derivative(localf)));
        localf->unbind();
    }
}

int main(int argc, char** argv)
{
  try{
    //Maybe initialize Mpi
    Dune::MPIHelper::instance(argc, argv);

    std::cout << "instantiate and evaluate some functions" << std::endl;

    // instantiate F and evaluate
    F<float> f;
    Dune::FieldVector<float,2> x;
    x[0] = 1.0; x[1] = 2.0;
    Dune::FieldVector<float,1> y;
    f.evaluate(x,y);

    // instantiate G and evaluate
    G<float> g;
    Dune::FieldVector<float,2> u;
    g.evaluate(x,u);

    // need a grid in order to test grid functions
    Dune::FieldVector<double,2> L(1.0);
    Dune::array<int,2> N(Dune::fill_array<int,2>(1));
    Dune::YaspGrid<2> grid(L,N);
    grid.globalRefine(6);

    // run algorithm on a grid
    std::cout << "instantiate grid functions on a grid" << std::endl;
    testgridfunction(grid.leafGridView(),F<Dune::YaspGrid<2>::ctype>());

    // run algorithm on a grid
    std::cout << "testing vtk output" << std::endl;
    testvtkexport(grid.leafGridView(),F<Dune::YaspGrid<2>::ctype>());
    testfunctiontree(grid.leafGridView());

    testgridviewfunction(grid.leafGridView());

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
