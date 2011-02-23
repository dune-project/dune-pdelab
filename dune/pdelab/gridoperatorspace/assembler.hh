class Assembler{
public:
  template<class LocalAssembler>
  void assemble(LocalAssembler & local_assembler);
};


class LocalAssemblerInterface : public LocalAssemblerBase 
{
public:
  bool requireIntersections();
  bool requireIntersectionsTwoSided();
  bool requireAlphaVolume();
  bool requireLambdaVolume();
  bool requireAlphaSkeleton();
  bool requireLambdaSkeleton();
  bool requireAlphaBoundary();
  bool requireLambdaBoundary();
  bool requirePostSkeletonAlphaVolume();
  
  template<class TT>
  void setTime(TT time);

  template<class RF>
  void setWeight(RF weight);

  LocalResidualAssembler & localResidualAssembler(R & r);
  LocalJacobianAssembler & localJacobianAssembler(A & a);
  LocalResidualJacobianAssembler & localResidualJacobianAssembler(R & r, A & a);
  
  class LocalAssemblerEngine{
  public:
    typedef LocalAssemblerInterface LocalAssembler; 
    const LocalAssembler & localAssembler();
    
    template<>
    void assembleAlphaVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);
    template<>
    void assembleLambdaVolume(const EG & eg, const LFSV & lfsv);
    template<>
    void assembleAlphaSkeleton(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s,
                               const LFSU_N & lfsu_n, const LFSV_N & lfsv_n);
    template<>
    void assembleLambdaSkeleton(const IG & ig, const LFSV_S & lfsv_s, const LFSV_N & lfsv_n);
    template<>
    void assembleAlphaBoundary(const IG & ig, const LFSU_S & lfsu_s, const LFSV_S & lfsv_s);
    template<>
    void assembleLambdaBoundary(const IG & ig, const LFSV_S & lfsv_s);
    template<>
    void assemblePostSkeletonVolume(const EG & eg, const LFSU & lfsu, const LFSV & lfsv);

    void preAssembly();
    void postAssembly();

    void onBindInside(const LFSU & lfsu_s, const LFSV & lfsv_s);
    void onBindOutside(const LFSU & lfsu_n, const LFSV & lfsv_n);
    void onUnbindInside(const LFSU & lfsu_s, const LFSV & lfsv_s);
    void onUnbindOutside(const LFSU & lfsu_n, const LFSV & lfsv_n);
  };

  class LocalResidualAssembler : public LocalAssemblerEngine {};
  class LocalJacobianAssembler : public LocalAssemblerEngine {};
  class LocalResidualJacobianAssembler : public LocalAssemblerEngine {};
};


class StandardGridAssembler{
public:
  template<typename P>
  void fill_pattern (P& globalpattern) const;

  template<typename X, typename R> 
  void residual (const X& x, R& r) const;

  template<typename X, typename A> 
  void jacobian (const X& x, A& a) const;

};
