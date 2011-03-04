#include <gridoperatorutilities.hh>
#include <assemblerutilties.hh>


class AssemblerInterface{
public:
  template<class LocalAssembler>
  void assemble(LocalAssembler & local_assembler);
};


class LocalAssemblerInterface : public LocalAssemblerBase
{
public:

  template<class TT>
  void setTime(TT time);

  template<class RF>
  void setWeight(RF weight);

  template<typename TT>
  void preStep (TT time, TT dt, std::size_t stages);

  void postStep ();

  template<typename TT>
  void preStage (TT time, std::size_t stage);

  void postStage ();

  template<typename TT>
  TT suggestTimestep (TT dt) const;

  LocalPatternAssemblerEngine & localPatternAssemblerEngine(R & r, const X & x);
  LocalResidualAssemblerEngine & localResidualAssemblerEngine(R & r, const X & x);
  LocalJacobianAssemblerEngine & localJacobianAssemblerEngine(A & a, const X & x);
  LocalResidualJacobianAssemblerEngine & localResidualJacobianAssemblerEngine(R & r, A & a, const X & x);

  class LocalAssemblerEngine{
  public:
    typedef LocalAssemblerInterface LocalAssembler;
    const LocalAssembler & localAssembler();

    bool requireSkeleton() const;
    bool requireSkeletonTwoSided() const;
    bool requireUVVolume() const;
    bool requireVVolume() const;
    bool requireUVSkeleton() const;
    bool requireVSkeleton() const;
    bool requireUVBoundary() const;
    bool requireVBoundary() const;
    bool requireUVEnrichedCoupling() const;
    bool requireVEnrichedCoupling() const;
    bool requireUVVolumePostSkeleton() const;
    bool requireVVolumePostSkeleton() const;

    template<>
    void assembleUVVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);
    template<>
    void assembleVVolume(const EG & eg, const LFSV & lfsv);
    template<>
    void assembleUVSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                               const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
    template<>
    void assembleVSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n);
    template<>
    void assembleUVBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    template<>
    void assembleVBoundary(const IG & ig, const LFSV_S & lfsv_s);
    template<>
    void assembleUVEnrichedCoupling(const IG & ig,
                                       const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                                       const LFSU_N & lfsu_n, const LFSV_N & lfsv_n,
                                       const LFSU_Coupling & lfsu_coupling, const LFSV_Coupling & lfsv_coupling);
    template<>
    void assembleVEnrichedCoupling(const IG & ig,
                                        const LFSV_S & lfsv_s,
                                        const LFSV_N & lfsv_n,
                                        const LFSV_Coupling & lfsv_coupling);
    template<>
    void assembleUVVolumePostSkeleton(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);
    template<>
    void assembleVVolumePostSkeleton(const EG & eg, const LFSV & lfsv);

    void preAssembly();
    void postAssembly();

    void onBindLFSUV(const EG & eg, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    void onBindLFSV(const EG & eg, const LFSV_S & lfsv_s);
    void onBindLFSUVInside(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    void onBindLFSVInside(const IG & ig, const LFSV_S & lfsv_s);
    void onBindLFSUVOutside(const IG & ig, const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
    void onBindLFSVOutside(const IG & ig, const LFSV_N & lfsv_n);
    void onBindLFSUVCoupling(const IG & ig, const LFSU_Coupling & lfsu_coupling, const LFSV_Coupling & lfsv_coupling);
    void onBindLFSVCoupling(const IG & ig, const LFSV_Coupling & lfsv_coupling);

    void loadCoefficientsLFSUInside(const LFSU_S & lfsu_s);
    void loadCoefficientsLFSUOutside(const LFSU_N & lfsu_n);
    void loadCoefficientsLFSUCoupling(const LFSU_Coupling & lfsu_coupling);

    void onUnbindLFSUV(const EG & eg, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    void onUnbindLFSV(const EG & eg, const LFSV_S & lfsv_s);
    void onUnbindLFSUVInside(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    void onUnbindLFSVInside(const IG & ig, const LFSV_S & lfsv_s);
    void onUnbindLFSUVOutside(const IG & ig, const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
    void onUnbindLFSVOutside(const IG & ig, const LFSV_N & lfsv_n);
    void onUnbindLFSUVCoupling(const IG & ig, const LFSU_Coupling & lfsu_coupling, const LFSV_Coupling & lfsv_coupling);
    void onUnbindLFSVCoupling(const IG & ig, const LFSV_Coupling & lfsv_coupling);

    // these methods are optional - not necessary for things like the PatternEngine
    void setSolution(const X& x);
    void setResidual(const R& r);
  };

  class LocalPatternAssemblerEngine : public LocalAssemblerEngine {};
  class LocalResidualAssemblerEngine : public LocalAssemblerEngine {};
  class LocalJacobianAssemblerEngine : public LocalAssemblerEngine {};
  class LocalResidualJacobianAssemblerEngine : public LocalAssemblerEngine {};
};

template<typename GFSU, typename GFSV,
         typename MB, typename DF, typename RF, typename JF>
class GridOperator{
public:

  typedef GridOperatorTraits<GFSU,
                             GFSV,
                             MB,
                             DF,
                             RF,
                             JF,
                             Assembler,
                             LocalAssembler> Traits;

  template<typename P>
  void fill_pattern (P& globalpattern) const;

  template<typename X, typename R>
  void residual (const X& x, R& r) const;

  template<typename X, typename A>
  void jacobian (const X& x, A& a) const;

  Assembler & assembler();

  LocalAssemblerInterface & localAssembler();

  const GFSU& trialGridFunctionSpace() const;

  const GFSV& testGridFunctionSpace() const;

  typename GFSU::Traits::SizeType globalSizeU () const;

  typename GFSV::Traits::SizeType globalSizeV () const;

  //! Interpolate xnew from f, taking unconstrained values from xold.
  /**
   * \note the exact type of F will depend on the GridOperator and
   *       may be a more complicated object than a simple GridFunction
   *       for scenarios like MultiDomain or grid-glue.
   */
  template<typename F>
  void interpolate(const typename Traits::Domain& xold,
                   const F& f,
                   typename Traits::Domain& xnew);

};
