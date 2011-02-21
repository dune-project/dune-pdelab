#ifndef DUNE_PDELAB_TEST_TYPETREESWITCH_HH
#define DUNE_PDELAB_TEST_TYPETREESWITCH_HH


#ifdef TEST_TYPETREE

#if !HAVE_VARIADIC_TEMPLATES || !HAVE_VARIADIC_CONSTRUCTOR_SFINAE || !HAVE_RVALUE_REFERENCES
#define TEST_TYPETREE_INVALID 1
#warning This test (full C++0x support) will not perform any checks due to missing compiler support
#endif

#elif TEST_TYPETREE_NO_SFINAE

#if !HAVE_VARIADIC_TEMPLATES || !HAVE_RVALUE_REFERENCES
#define TEST_TYPETREE_INVALID 1
#warning This test (C++0x support without variadic constructor SFINAE) will not perform any checks due to missing compiler support
#else

#ifdef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#undef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#endif

#endif

#elif TEST_TYPETREE_NO_VARIADIC

#if !HAVE_RVALUE_REFERENCES
#define TEST_TYPETREE_INVALID 1
#warning This test (C++0x support without variadic templates) will not perform any checks due to missing compiler support
#else

#ifdef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#undef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#endif

#ifdef HAVE_VARIADIC_TEMPLATES
#undef HAVE_VARIADIC_TEMPLATES
#endif

#endif

#elif TEST_TYPETREE_NO_RVALUE_REFERENCES

#if !HAVE_VARIADIC_TEMPLATES
#define TEST_TYPETREE_INVALID 1
#warning This test (C++0x support without rvalue references) will not perform any checks due to missing compiler support
#else

#ifdef HAVE_RVALUE_REFERENCES
#undef HAVE_RVALUE_REFERENCES
#endif

#endif

#elif TEST_TYPETREE_LEGACY

#ifdef HAVE_RVALUE_REFERENCES
#undef HAVE_RVALUE_REFERENCES
#endif

#ifdef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#undef HAVE_VARIADIC_CONSTRUCTOR_SFINAE
#endif

#ifdef HAVE_VARIADIC_TEMPLATES
#undef HAVE_VARIADIC_TEMPLATES
#endif

#else
#error You need to specify a test case via preprocessor define!
#endif


#endif // DUNE_PDELAB_TEST_TYPETREESWITCH_HH
