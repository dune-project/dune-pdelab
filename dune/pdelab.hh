#ifndef DUNE_PDELAB_HH
#define DUNE_PDELAB_HH

#include <dune/pdelab/gridfunctionspace/dunefunctionslfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspaceutilities.hh>
#include <dune/pdelab/gridfunctionspace/interpolate.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powercompositegridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/lfsindexcache.hh>
#include <dune/pdelab/gridfunctionspace/vectorgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/powergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/dynamicpowergridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/vtk.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/compositegridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/entityindexcache.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionsgridfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/localfunctionspacetags.hh>
#include <dune/pdelab/gridfunctionspace/dunefunctionslocalfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/subspace.hh>
#include <dune/pdelab/gridfunctionspace/loadbalance.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionadapter.hh>
#include <dune/pdelab/gridfunctionspace/subspacelocalfunctionspace.hh>
#include <dune/pdelab/gridfunctionspace/tags.hh>
#include <dune/pdelab/gridfunctionspace/gridfunctionspacebase.hh>
#include <dune/pdelab/gridfunctionspace/localvector.hh>
#include <dune/pdelab/gridfunctionspace/datahandleprovider.hh>
#include <dune/pdelab/gridfunctionspace/genericdatahandle.hh>
#include <dune/pdelab/newton/newton.hh>
#include <dune/pdelab/finiteelement/l2orthonormal.hh>
#include <dune/pdelab/finiteelement/localbasiscache.hh>
#include <dune/pdelab/finiteelement/pk1d.hh>
#include <dune/pdelab/finiteelement/qkdglobatto.hh>
#include <dune/pdelab/finiteelement/qkdglagrange.hh>
#include <dune/pdelab/finiteelement/qkdglegendre.hh>
#include <dune/pdelab/localoperator/numericaljacobian.hh>
#include <dune/pdelab/localoperator/darcyccfv.hh>
#include <dune/pdelab/localoperator/maxwellparameter.hh>
#include <dune/pdelab/localoperator/variablefactories.hh>
#include <dune/pdelab/localoperator/nonlinearconvectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/numericalnonlinearjacobianapply.hh>
#include <dune/pdelab/localoperator/numericaljacobianapply.hh>
#include <dune/pdelab/localoperator/dgnavierstokesvelvecfem.hh>
#include <dune/pdelab/localoperator/navierstokesmass.hh>
#include <dune/pdelab/localoperator/convectiondiffusionparameter.hh>
#include <dune/pdelab/localoperator/interface.hh>
#include <dune/pdelab/localoperator/diffusionmixed.hh>
#include <dune/pdelab/localoperator/zero.hh>
#include <dune/pdelab/localoperator/scaled.hh>
#include <dune/pdelab/localoperator/numericalresidual.hh>
#include <dune/pdelab/localoperator/blockdiagonal.hh>
#include <dune/pdelab/localoperator/dgnavierstokesparameter.hh>
#include <dune/pdelab/localoperator/taylorhoodnavierstokes.hh>
#include <dune/pdelab/localoperator/linearacousticsdg.hh>
#include <dune/pdelab/localoperator/flags.hh>
#include <dune/pdelab/localoperator/permeability_adapter.hh>
#include <dune/pdelab/localoperator/idefault.hh>
#include <dune/pdelab/localoperator/darcyfem.hh>
#include <dune/pdelab/localoperator/l2.hh>
#include <dune/pdelab/localoperator/dgnavierstokes.hh>
#include <dune/pdelab/localoperator/l2volumefunctional.hh>
#include <dune/pdelab/localoperator/electrodynamic.hh>
#include <dune/pdelab/localoperator/pattern.hh>
#include <dune/pdelab/localoperator/callswitch.hh>
#include <dune/pdelab/localoperator/convectiondiffusionfem.hh>
#include <dune/pdelab/localoperator/eval.hh>
#include <dune/pdelab/localoperator/dginteriorpenaltyparameter.hh>
#include <dune/pdelab/localoperator/convectiondiffusiondg.hh>
#include <dune/pdelab/localoperator/linearelasticityparameter.hh>
#include <dune/pdelab/localoperator/defaultimp.hh>
#include <dune/pdelab/localoperator/errorindicatordg.hh>
#include <dune/pdelab/localoperator/linearacousticsparameter.hh>
#include <dune/pdelab/localoperator/twophaseccfv.hh>
#include <dune/pdelab/localoperator/convectiondiffusionccfv.hh>
#include <dune/pdelab/localoperator/maxwelldg.hh>
#include <dune/pdelab/localoperator/linearelasticity.hh>
#include <dune/pdelab/localoperator/stokesparameter.hh>
#include <dune/pdelab/adaptivity/adaptivity.hh>
#include <dune/pdelab/function/product.hh>
#include <dune/pdelab/function/sqr.hh>
#include <dune/pdelab/function/memberadaptor.hh>
#include <dune/pdelab/function/bindtime.hh>
#include <dune/pdelab/function/division.hh>
#include <dune/pdelab/function/sqrt.hh>
#include <dune/pdelab/function/minus.hh>
#include <dune/pdelab/function/localfunctionhelper.hh>
#include <dune/pdelab/function/discretegridviewfunction.hh>
#include <dune/pdelab/function/inverse.hh>
#include <dune/pdelab/function/tags.hh>
#include <dune/pdelab/function/const.hh>
#include <dune/pdelab/function/scalarscaled.hh>
#include <dune/pdelab/function/localfunction.hh>
#include <dune/pdelab/function/callableadapter.hh>
#include <dune/pdelab/function/oldinterfaceadapter.hh>
#include <dune/pdelab/function/selectcomponent.hh>
#include <dune/pdelab/common/logtag.hh>
#include <dune/pdelab/common/utility.hh>
#include <dune/pdelab/common/vtkexport.hh>
#include <dune/pdelab/common/partitionviewentityset.hh>
#include <dune/pdelab/common/elementmapper.hh>
#include <dune/pdelab/common/quadraturerules.hh>
#include <dune/pdelab/common/dofindex.hh>
#include <dune/pdelab/common/benchmarkhelper.hh>
#include <dune/pdelab/common/topologyutility.hh>
#include <dune/pdelab/common/intersectiontype.hh>
#include <dune/pdelab/common/range.hh>
#include <dune/pdelab/common/instationaryfilenamehelper.hh>
#include <dune/pdelab/common/functionwrappers.hh>
#include <dune/pdelab/common/exceptions.hh>
#include <dune/pdelab/common/simpledofindex.hh>
#include <dune/pdelab/common/hostname.hh>
#include <dune/pdelab/common/functionutilities.hh>
#include <dune/pdelab/common/borderindexidcache.hh>
#include <dune/pdelab/common/geometrywrapper.hh>
#include <dune/pdelab/common/clock.hh>
#include <dune/pdelab/common/crossproduct.hh>
#include <dune/pdelab/common/polymorphicbufferwrapper.hh>
#include <dune/pdelab/common/typetraits.hh>
#include <dune/pdelab/common/function.hh>
#include <dune/pdelab/common/globaldofindex.hh>
#include <dune/pdelab/common/multiindex.hh>
#include <dune/pdelab/common/jacobiantocurl.hh>
#include <dune/pdelab/stationary/linearproblem.hh>
#include <dune/pdelab/constraints/noconstraints.hh>
#include <dune/pdelab/constraints/hangingnodemanager.hh>
#include <dune/pdelab/constraints/p0ghost.hh>
#include <dune/pdelab/constraints/common/constraintstransformation.hh>
#include <dune/pdelab/constraints/common/constraints.hh>
#include <dune/pdelab/constraints/common/constraintsparameters.hh>
#include <dune/pdelab/constraints/raviartthomas0.hh>
#include <dune/pdelab/constraints/p0.hh>
#include <dune/pdelab/constraints/hangingnode.hh>
#include <dune/pdelab/constraints/conforming.hh>
#include <dune/pdelab/constraints/interiornode.hh>
#include <dune/pdelab/backend/interface.hh>
#include <dune/pdelab/backend/istl.hh>
#include <dune/pdelab/backend/common/aliasedmatrixview.hh>
#include <dune/pdelab/backend/common/tags.hh>
#include <dune/pdelab/backend/common/uncachedvectorview.hh>
#include <dune/pdelab/backend/common/uncachedmatrixview.hh>
#include <dune/pdelab/backend/common/aliasedvectorview.hh>
#include <dune/pdelab/backend/solver.hh>
#include <dune/pdelab/backend/eigen.hh>
#include <dune/pdelab/backend/eigen/solvers.hh>
#include <dune/pdelab/backend/eigen/descriptors.hh>
#include <dune/pdelab/backend/eigen/vector.hh>
#include <dune/pdelab/backend/eigen/matrix.hh>
#include <dune/pdelab/backend/istl/utility.hh>
#include <dune/pdelab/backend/istl/bcrsmatrix.hh>
#include <dune/pdelab/backend/istl/ovlp_amg_dg_backend.hh>
#include <dune/pdelab/backend/istl/seq_amg_dg_backend.hh>
#include <dune/pdelab/backend/istl/ovlpistlsolverbackend.hh>
#include <dune/pdelab/backend/istl/vectoriterator.hh>
#include <dune/pdelab/backend/istl/vectorhelpers.hh>
#include <dune/pdelab/backend/istl/descriptors.hh>
#include <dune/pdelab/backend/istl/bcrspattern.hh>
#include <dune/pdelab/backend/istl/cg_to_dg_prolongation.hh>
#include <dune/pdelab/backend/istl/forwarddeclarations.hh>
#include <dune/pdelab/backend/istl/novlpistlsolverbackend.hh>
#include <dune/pdelab/backend/istl/seqistlsolverbackend.hh>
#include <dune/pdelab/backend/istl/parallelhelper.hh>
#include <dune/pdelab/backend/istl/vector.hh>
#include <dune/pdelab/backend/istl/bcrsmatrixbackend.hh>
#include <dune/pdelab/backend/istl/tags.hh>
#include <dune/pdelab/backend/istl/patternstatistics.hh>
#include <dune/pdelab/backend/istl/blockmatrixdiagonal.hh>
#include <dune/pdelab/backend/istl/dunefunctions.hh>
#include <dune/pdelab/backend/istl/istlsolverbackend.hh>
#include <dune/pdelab/backend/istl/matrixhelpers.hh>
#include <dune/pdelab/backend/istl/geneo/subdomainbasis.hh>
#include <dune/pdelab/backend/istl/geneo/localoperator_ovlp_region.hh>
#include <dune/pdelab/backend/istl/geneo/subdomainprojectedcoarsespace.hh>
#include <dune/pdelab/backend/istl/geneo/two_level_schwarz.hh>
#include <dune/pdelab/backend/istl/geneo/geneobasis.hh>
#include <dune/pdelab/backend/istl/geneo/multicommdatahandle.hh>
#include <dune/pdelab/backend/istl/geneo/liptonbabuskabasis.hh>
#include <dune/pdelab/backend/istl/geneo/geneo.hh>
#include <dune/pdelab/backend/istl/geneo/partitionofunity.hh>
#include <dune/pdelab/backend/istl/geneo/coarsespace.hh>
#include <dune/pdelab/backend/istl/geneo/arpackpp_geneo.hh>
#include <dune/pdelab/backend/simple/sparse.hh>
#include <dune/pdelab/backend/simple/descriptors.hh>
#include <dune/pdelab/backend/simple/vector.hh>
#include <dune/pdelab/backend/simple/matrix.hh>
#include <dune/pdelab/backend/simple.hh>
#include <dune/pdelab/ordering/utility.hh>
#include <dune/pdelab/ordering/leaflocalordering.hh>
#include <dune/pdelab/ordering/localorderingbase.hh>
#include <dune/pdelab/ordering/directleaflocalordering.hh>
#include <dune/pdelab/ordering/leafgridviewordering.hh>
#include <dune/pdelab/ordering/lexicographicordering.hh>
#include <dune/pdelab/ordering/permutedordering.hh>
#include <dune/pdelab/ordering/decorator.hh>
#include <dune/pdelab/ordering/subordering.hh>
#include <dune/pdelab/ordering/entityblockedlocalordering.hh>
#include <dune/pdelab/ordering/leaforderingbase.hh>
#include <dune/pdelab/ordering/orderingbase.hh>
#include <dune/pdelab/ordering/transformations.hh>
#include <dune/pdelab/ordering/chunkedblockordering.hh>
#include <dune/pdelab/ordering/singlecodimleafordering.hh>
#include <dune/pdelab/ordering/interleavedordering.hh>
#include <dune/pdelab/ordering/gridviewordering.hh>
#include <dune/pdelab/instationary/implicitonestep.hh>
#include <dune/pdelab/instationary/onestep.hh>
#include <dune/pdelab/instationary/explicitonestep.hh>
#include <dune/pdelab/instationary/onestepparameter.hh>
#include <dune/pdelab/finiteelementmap/qkfem.hh>
#include <dune/pdelab/finiteelementmap/utility.hh>
#include <dune/pdelab/finiteelementmap/rt0cube3dfem.hh>
#include <dune/pdelab/finiteelementmap/mimeticfem.hh>
#include <dune/pdelab/finiteelementmap/variableopbfem.hh>
#include <dune/pdelab/finiteelementmap/bdm1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/edges0.5fem.hh>
#include <dune/pdelab/finiteelementmap/p0fem.hh>
#include <dune/pdelab/finiteelementmap/rannacherturekfem.hh>
#include <dune/pdelab/finiteelementmap/rt1cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/pkqkfem.hh>
#include <dune/pdelab/finiteelementmap/pkfem.hh>
#include <dune/pdelab/finiteelementmap/bdm1cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt0simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt1simplex2dfem.hh>
#include <dune/pdelab/finiteelementmap/rt0cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/pk1d.hh>
#include <dune/pdelab/finiteelementmap/rt2cube2dfem.hh>
#include <dune/pdelab/finiteelementmap/raviartthomasfem.hh>
#include <dune/pdelab/finiteelementmap/variableqkdgfem.hh>
#include <dune/pdelab/finiteelementmap/powerfem.hh>
#include <dune/pdelab/finiteelementmap/qkdg.hh>
#include <dune/pdelab/finiteelementmap/global.hh>
#include <dune/pdelab/finiteelementmap/monomfem.hh>
#include <dune/pdelab/finiteelementmap/rt1cube3dfem.hh>
#include <dune/pdelab/finiteelementmap/brezzidouglasmarinifem.hh>
#include <dune/pdelab/finiteelementmap/opbfem.hh>
#include <dune/pdelab/finiteelementmap/finiteelementmap.hh>
#include <dune/pdelab/finiteelementmap/variablemonomfem.hh>
#include <dune/pdelab/boilerplate/pdelab.hh>
#include <dune/pdelab/gridoperator/gridoperator.hh>
#include <dune/pdelab/gridoperator/onestep/residualengine.hh>
#include <dune/pdelab/gridoperator/onestep/localassembler.hh>
#include <dune/pdelab/gridoperator/onestep/jacobianengine.hh>
#include <dune/pdelab/gridoperator/onestep/patternengine.hh>
#include <dune/pdelab/gridoperator/onestep/prestageengine.hh>
#include <dune/pdelab/gridoperator/onestep/jacobianresidualengine.hh>
#include <dune/pdelab/gridoperator/onestep/enginebase.hh>
#include <dune/pdelab/gridoperator/fastdg.hh>
#include <dune/pdelab/gridoperator/common/assemblerutilities.hh>
#include <dune/pdelab/gridoperator/common/diagonallocalmatrix.hh>
#include <dune/pdelab/gridoperator/common/localmatrix.hh>
#include <dune/pdelab/gridoperator/common/borderdofexchanger.hh>
#include <dune/pdelab/gridoperator/common/gridoperatorutilities.hh>
#include <dune/pdelab/gridoperator/common/localassemblerenginebase.hh>
#include <dune/pdelab/gridoperator/fastdg/residualengine.hh>
#include <dune/pdelab/gridoperator/fastdg/localassembler.hh>
#include <dune/pdelab/gridoperator/fastdg/jacobianengine.hh>
#include <dune/pdelab/gridoperator/fastdg/assembler.hh>
#include <dune/pdelab/gridoperator/fastdg/patternengine.hh>
#include <dune/pdelab/gridoperator/fastdg/nonlinearjacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/fastdg/jacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/onestep.hh>
#include <dune/pdelab/gridoperator/default/residualengine.hh>
#include <dune/pdelab/gridoperator/default/localassembler.hh>
#include <dune/pdelab/gridoperator/default/jacobianengine.hh>
#include <dune/pdelab/gridoperator/default/assembler.hh>
#include <dune/pdelab/gridoperator/default/patternengine.hh>
#include <dune/pdelab/gridoperator/default/nonlinearjacobianapplyengine.hh>
#include <dune/pdelab/gridoperator/default/jacobianapplyengine.hh>

#endif // DUNE_PDELAB_HH
