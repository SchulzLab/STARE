cmake_minimum_required(VERSION 3.9)
project(STARE)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

file( GLOB LIB_SOURCES STARE_MiscFunctions.cpp )
file( GLOB LIB_HEADERS STARE_MiscFunctions.h )
add_library( STARELib ${LIB_SOURCES} ${LIB_HEADERS} )

add_executable( ReplaceInvalidChars ReplaceInvalidChars.cpp )
target_link_libraries( ReplaceInvalidChars PUBLIC STARELib )

add_executable( Promoter_Windows Promoter_Windows.cpp )
target_link_libraries( Promoter_Windows PUBLIC STARELib )

add_executable( PSCM_to_PSEM PSCM_to_PSEM.cpp )

add_executable( TRAPmulti TRAPmulti.cpp )

add_executable( STARE_ABCpp STARE_ABCpp.cpp )
target_link_libraries( STARE_ABCpp PUBLIC STARELib )

add_executable( TF_Gene_Scorer_Reshape TF_Gene_Scorer_Reshape.cpp )
target_link_libraries( TF_Gene_Scorer_Reshape PUBLIC STARELib )

find_package(Llvm)
find_package(OpenMP REQUIRED)
if(OpenMP_CXX_FOUND)
    target_link_libraries(TRAPmulti PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(STARE_ABCpp PUBLIC OpenMP::OpenMP_CXX)
    target_link_libraries(TF_Gene_Scorer_Reshape PUBLIC OpenMP::OpenMP_CXX)
endif()
find_package(Boost REQUIRED)





