# Module that checks whether the compiler supports
# the GCC extension __typeof__, which can serve as
# a limited fallback for decltype on some older compilers.
#
# Sets the following variable:
# HAVE_GCC___TYPEOF__
#
# perform tests
include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
  template<typename A, typename B>
  struct check_equal;

  template<typename A>
  struct check_equal<A,A>
  {
    static const int result = 0;
  };

  struct A {};

  A foo();

  int main(void){
      return check_equal<__typeof__(foo()),A>::result;
  }"
  HAVE_GCC___TYPEOF__)
