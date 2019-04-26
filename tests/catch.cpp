#define CATCH_CONFIG_MAIN

// On Visual Studio compilers we enable CRTDBG macros to detect
// memory leaks.
#ifdef _MSC_VER
#define CATCH_CONFIG_WINDOWS_CRTDBG
#endif

#include "catch.hpp"