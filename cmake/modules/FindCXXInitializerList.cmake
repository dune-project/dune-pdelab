# Module that checks whether the compiler supports
# std::initializer_list, both at the compiler level and
# at the library level.
#
# Sets the following variable:
# HAVE_INITIALIZER_LIST
#
# perform tests
include(CheckCXXSourceCompiles)

check_cxx_source_compiles("
  #include <initializer_list>
  #include <vector>

  struct A
  {

    A(std::initializer_list<int> il)
      : vec(il)
    {}

    std::vector<int> vec;
  };

  int main(void){
      A a{1,3,4,5};
      return 0;
  }"
  HAVE_INITIALIZER_LIST)
