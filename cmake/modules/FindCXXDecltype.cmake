# Module that checks whether the compiler supports
# the C++11 keyword decltype.
#
# Sets the following variable:
# HAVE_STD_DECLTYPE
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

  struct A;

  const A& foo();

  int main(void){
      return check_equal<decltype(foo()),const A&>::result;
  }"
  HAVE_STD_DECLTYPE)
