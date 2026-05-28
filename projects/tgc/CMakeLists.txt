cmake_minimum_required(VERSION 3.20)

project(tgc_sim LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(TGC_SHARED_GARFIELD_HINT "${CMAKE_SOURCE_DIR}/../../local/garfield")

if(DEFINED ENV{ROOTSYS})
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{ROOTSYS}")
endif()
if(DEFINED ENV{GARFIELD_INSTALL} AND NOT "$ENV{GARFIELD_INSTALL}" STREQUAL "")
  list(PREPEND CMAKE_PREFIX_PATH "$ENV{GARFIELD_INSTALL}")
elseif(EXISTS "${TGC_SHARED_GARFIELD_HINT}/lib/cmake/Garfield/GarfieldConfig.cmake")
  list(PREPEND CMAKE_PREFIX_PATH "${TGC_SHARED_GARFIELD_HINT}")
endif()

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Gpad)
find_package(Garfield CONFIG QUIET)

if(NOT Garfield_FOUND)
  message(FATAL_ERROR
    "Garfield++ was not found. Set GARFIELD_INSTALL or CMAKE_PREFIX_PATH to a Garfield installation, "
    "or source ../../local/garfield/share/Garfield/setupGarfield.sh when using the shared workspace install."
  )
endif()

add_executable(tgc_sim
  src/tgc_sim.cc
)

target_include_directories(tgc_sim PRIVATE
  ${CMAKE_SOURCE_DIR}/third_party
)

target_link_libraries(tgc_sim PRIVATE
  Garfield::Garfield
  ROOT::Core
  ROOT::RIO
  ROOT::Hist
  ROOT::Graf
  ROOT::Gpad
)

target_compile_options(tgc_sim PRIVATE
  $<$<OR:$<CXX_COMPILER_ID:Clang>,$<CXX_COMPILER_ID:AppleClang>,$<CXX_COMPILER_ID:GNU>>:
    -Wall
    -Wextra
    -Wpedantic
    -Wshadow
  >
)

set(_runtime_paths "")
if(DEFINED ENV{GARFIELD_INSTALL} AND NOT "$ENV{GARFIELD_INSTALL}" STREQUAL "")
  list(APPEND _runtime_paths "$ENV{GARFIELD_INSTALL}/lib")
elseif(EXISTS "${TGC_SHARED_GARFIELD_HINT}/lib")
  list(APPEND _runtime_paths "${TGC_SHARED_GARFIELD_HINT}/lib")
endif()
if(DEFINED ENV{ROOTSYS})
  list(APPEND _runtime_paths "$ENV{ROOTSYS}/lib")
endif()

set_target_properties(tgc_sim PROPERTIES
  BUILD_RPATH "${_runtime_paths}"
  INSTALL_RPATH "${_runtime_paths}"
)

include(CTest)

if(BUILD_TESTING)
  add_test(
    NAME smoke_tgc_sim
    COMMAND tgc_sim
      --config ${CMAKE_SOURCE_DIR}/config/smoke_tgc.json
      --out ${CMAKE_BINARY_DIR}/smoke-output
  )
  set_tests_properties(smoke_tgc_sim PROPERTIES
    ENVIRONMENT
      "GARFIELD_INSTALL=${TGC_SHARED_GARFIELD_HINT};HEED_DATABASE=${TGC_SHARED_GARFIELD_HINT}/share/Heed/database"
    TIMEOUT 300
  )
endif()
