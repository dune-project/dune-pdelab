#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <dune/pdelab.hh>

// This tests if the function Dune::PDELab::impl::hasOldLOPInterface correctly
// detects if a local operator only implements the old nonlinear jacobian apply
// interface. See comments in checklopinterface.hh for more details.

// jacobian_apply_volume with new interface
struct A1{
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void jacobian_apply_volume(const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5, T6& t6) const {}
};
// jacobian_apply_volume with old interface
struct A2{
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void jacobian_apply_volume(const T1& t1, const T2& t2, const T3& t3_1, const T3& t3_2, const T4& t4, T5& t5) const {}
};


// jacobian_apply_volume_post_skeleton with new interface
struct B1{
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void jacobian_apply_volume_post_skeleton(const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5, T6& t6) const {}
};
// jacobian_apply_volume_post_skeleton with old interface
struct B2{
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void jacobian_apply_volume_post_skeleton(const T1& t1, const T2& t2, const T3& t3_1, const T3& t3_2, const T4& t4, T5& t5) const {}
};


// jacobian_apply_skeleton with new interface
struct C1{
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void jacobian_apply_skeleton(const T1& t1,
                               const T2& t2_1, const T3& t3_1, const T4& t4_1, const T5& t5_1,
                               const T2& t2_2, const T3& t3_2, const T4& t4_2, const T5& t5_2,
                               T6& t6_1, T6& t6_2) const {}
};
// jacobian_apply_skeleton with old interface
struct C2{
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void jacobian_apply_skeleton(const T1& t1,
                               const T2& t2_1, const T3& t3_1, const T3& t3_2, const T4& t4_1,
                               const T2& t2_2, const T3& t3_3, const T3& t3_4, const T4& t4_2,
                               T5& t5_1, T5& t5_2) const {}
};


// jacobian_apply_volume_boundary with new interface
struct D1{
  template <typename T1, typename T2, typename T3, typename T4, typename T5, typename T6>
  void jacobian_apply_boundary(const T1& t1, const T2& t2, const T3& t3, const T4& t4, const T5& t5, T6& t6) const {}
};
// jacobian_apply_volume_boundary with old interface
struct D2{
  template <typename T1, typename T2, typename T3, typename T4, typename T5>
  void jacobian_apply_boundary(const T1& t1, const T2& t2, const T3& t3_1, const T3& t3_2, const T4& t4, T5& t5) const {}
};


int main(){
  A1 a1;
  std::cout << "A1: " << Dune::PDELab::impl::hasOldLOPInterface(a1) << std::endl;

  A2 a2;
  std::cout << "A2: " << Dune::PDELab::impl::hasOldLOPInterface(a2) << std::endl;

  B1 b1;
  std::cout << "B1: " << Dune::PDELab::impl::hasOldLOPInterface(b1) << std::endl;

  B2 b2;
  std::cout << "B2: " << Dune::PDELab::impl::hasOldLOPInterface(b2) << std::endl;

  C1 c1;
  std::cout << "C1: " << Dune::PDELab::impl::hasOldLOPInterface(c1) << std::endl;

  C2 c2;
  std::cout << "C2: " << Dune::PDELab::impl::hasOldLOPInterface(c2) << std::endl;

  D1 d1;
  std::cout << "D1: " << Dune::PDELab::impl::hasOldLOPInterface(d1) << std::endl;

  D2 d2;
  std::cout << "D2: " << Dune::PDELab::impl::hasOldLOPInterface(d2) << std::endl;

  return 0;
}
