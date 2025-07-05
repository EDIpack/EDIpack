##################################################
# DMFTED PROJECT
##################################################
CMAKE_MINIMUM_REQUIRED(VERSION 3.10.0)
PROJECT(edipack Fortran)
SET(EDI ${PROJECT_NAME})
SET(EDI2INEQ edipack2ineq)
SET(EDI_C edipack_cbindings)
SET(EDI2INEQ_C edipack2ineq_cbindings)
SET(VERSION 5.3.0)


MESSAGE(STATUS "OS: ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_VERSION}")
SET(BUILD_TYPE "RELEASE" CACHE STRING "Build type")
STRING(TOUPPER "${BUILD_TYPE}" BUILD_TYPE)
MESSAGE(STATUS "BUILD_TYPE=${BUILD_TYPE}")
SET(CMAKE_BUILD_TYPE ${BUILD_TYPE} CACHE INTERNAL "Build type" FORCE)


# CMakeGraphVizOptions.cmake
set(GRAPHVIZ_CUSTOM_TARGETS TRUE)


# Add our local CMAKE modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")
INCLUDE(${CMAKE_MODULE_PATH}/ColorsMsg.cmake)

# Set some options the user may choose
OPTION(USE_MPI "Set MPI pre-compilation flag" ON)
if(USE_MPI)
  message(STATUS "USE_MPI is ON")
else()
  message(STATUS "USE_MPI is OFF")
endif()
OPTION(WITH_INEQ "Includes Support for inequivalent impurities" ON)
if(WITH_INEQ)
  message(STATUS "WITH_INEQ is ON")
else()
  message(STATUS "WITH_INEQ is OFF")
endif()
OPTION(VERBOSE "Build a verbose Makefile" OFF)
IF(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE ON)
ELSE(VERBOSE)
  SET(CMAKE_VERBOSE_MAKEFILE OFF)
ENDIF(VERBOSE)




############################################################
# Look for MPI if needed
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MpiConfig.cmake)



############################################################
# Look for required libraries: scifor
# (includes Blas/Lapack/Scalapack/Quadpack/Minpack/FFTpack/Arpack/...) 
############################################################
# Find Scifor
FIND_PACKAGE(PkgConfig)
PKG_SEARCH_MODULE(SCIFOR REQUIRED scifor)
IF(${SCIFOR_FOUND})
  MESSAGE(STATUS "SCIFOR found at: ${SCIFOR_PREFIX}")
ELSE()
  MESSAGE(FATAL_ERROR "
Can not find SciFortran in this system. Check SciFor installation. Try adding the provided PKG-CONFIG file to the environment:
$ export PKG_CONFIG_PATH=<path_to_scifor_dir>/etc:$PKG_CONFIG_PATH")
ENDIF()
INCLUDE_DIRECTORIES(BEFORE ${SCIFOR_INCLUDE_DIRS})

#message(STATUS "SCIFOR_FOUND: ${SCIFOR_FOUND}")
#message(STATUS "SCIFOR_LIBRARIES: ${SCIFOR_LIBRARIES}")
#message(STATUS "SCIFOR_LDFLAGS: ${SCIFOR_LDFLAGS}")
#message(STATUS "SCIFOR_LIBS: ${SCIFOR_LIBS}")
#message(STATUS "SCIFOR_LINK_LIBRARIES: ${SCIFOR_LINK_LIBRARIES}")


############################################################
# Configure Fortran Compiler Main options
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/FortranCompilerConfig.cmake)




############################################################
# Check for .git and extract the current SHA1 revision.
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/GitConfig.cmake)




############################################################
# Define the actual files and folders that make up the build
############################################################
INCLUDE(${CMAKE_MODULE_PATH}/MainConfig.cmake)



############################################################
############################################################
############################################################
# DEFINE HERE THE DETAILS OF THE LIBRARY CONFIGURATION 

# Define some directories
SET(LIB_SRC ${CMAKE_SOURCE_DIR}/src/singlesite)
SET(LIB_ETC ${CMAKE_SOURCE_DIR}/etc)
SET(LIB_TEST ${CMAKE_SOURCE_DIR}/test)
SET(LIB_CBIND ${CMAKE_SOURCE_DIR}/src/c_bindings)
SET(LIB_INEQ ${CMAKE_SOURCE_DIR}/src/ineq)
SET(LIB_ENV ${CMAKE_SOURCE_DIR}/etc/environment_modules)
#
SET(LIB_TMP_INC ${CMAKE_BINARY_DIR}/include)
SET(LIB_TMP_ETC ${CMAKE_BINARY_DIR}/etc)
SET(LIB_TMP_VER ${CMAKE_BINARY_DIR}/version)
#
SET(LIB_TARGET_LIB ${CMAKE_INSTALL_PREFIX}/lib)
SET(LIB_TARGET_INC ${CMAKE_INSTALL_PREFIX}/include)
SET(LIB_TARGET_ETC ${CMAKE_INSTALL_PREFIX}/etc)
SET(LIB_TARGET_BIN ${CMAKE_INSTALL_PREFIX}/bin)
SET(LIB_TARGET_DIR ${CMAKE_INSTALL_PREFIX})

# Have the .mod files placed in the INC folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB_TMP_INC})
INCLUDE_DIRECTORIES(BEFORE ${CMAKE_Fortran_MODULE_DIRECTORY})






####################################################
#       EDIpack library
####################################################
ADD_SUBDIRECTORY(${LIB_SRC}/ED_BATH)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_NORMAL)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_SUPERC)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_NONSU2)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_IO)

ADD_SUBDIRECTORY(${LIB_SRC}/ED_FIT)

ADD_SUBDIRECTORY(${LIB_SRC})



ADD_DEPENDENCIES(ED_BATH_LIB
  ED_INIT_LIB
)


ADD_DEPENDENCIES(ED_NORMAL_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)

ADD_DEPENDENCIES(ED_SUPERC_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)


ADD_DEPENDENCIES(ED_NONSU2_LIB
  ED_INIT_LIB
  ED_BATH_LIB
)



ADD_DEPENDENCIES(ED_FIT_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
)



ADD_DEPENDENCIES(ED_GLOB_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_FIT_LIB
)

ADD_DEPENDENCIES(ED_IO_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_GLOB_LIB
)


ADD_DEPENDENCIES(ED_LIB_LIB
  ED_INIT_LIB
  ED_BATH_LIB
  ED_NORMAL_LIB
  ED_SUPERC_LIB
  ED_NONSU2_LIB
  ED_FIT_LIB
  ED_IO_LIB
)

IF(WITH_INEQ)
  ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_BATH )  
  ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_IO   )   
  ADD_SUBDIRECTORY(${LIB_INEQ}/E2I_FIT  )  
  ADD_SUBDIRECTORY(${LIB_INEQ}          )
  ADD_DEPENDENCIES(E2I_INIT_LIB ED_LIB_LIB)
  ADD_DEPENDENCIES(E2I_BATH_LIB ED_LIB_LIB E2I_INIT_LIB)
  ADD_DEPENDENCIES(E2I_FIT_LIB  ED_LIB_LIB E2I_INIT_LIB E2I_BATH_LIB)
  ADD_DEPENDENCIES(E2I_IO_LIB   ED_LIB_LIB E2I_INIT_LIB E2I_BATH_LIB)
  ADD_DEPENDENCIES(E2I_LIB      ED_LIB_LIB E2I_INIT_LIB E2I_BATH_LIB E2I_FIT_LIB E2I_IO_LIB)
  ADD_LIBRARY(${EDI} STATIC
     $<TARGET_OBJECTS:ED_INIT_LIB>
     $<TARGET_OBJECTS:ED_BATH_LIB>
     $<TARGET_OBJECTS:ED_NORMAL_LIB>
     $<TARGET_OBJECTS:ED_SUPERC_LIB>
     $<TARGET_OBJECTS:ED_NONSU2_LIB>
     $<TARGET_OBJECTS:ED_FIT_LIB>
     $<TARGET_OBJECTS:ED_GLOB_LIB>
     $<TARGET_OBJECTS:ED_IO_LIB>
     $<TARGET_OBJECTS:ED_LIB_LIB>
     $<TARGET_OBJECTS:E2I_INIT_LIB>
     $<TARGET_OBJECTS:E2I_BATH_LIB>
     $<TARGET_OBJECTS:E2I_FIT_LIB>
     $<TARGET_OBJECTS:E2I_IO_LIB>
     $<TARGET_OBJECTS:E2I_LIB>
   )
 ELSE(WITH_INEQ)
  ADD_LIBRARY(${EDI} STATIC
    $<TARGET_OBJECTS:ED_INIT_LIB>
    $<TARGET_OBJECTS:ED_BATH_LIB>
    $<TARGET_OBJECTS:ED_NORMAL_LIB>
    $<TARGET_OBJECTS:ED_SUPERC_LIB>
    $<TARGET_OBJECTS:ED_NONSU2_LIB>
    $<TARGET_OBJECTS:ED_FIT_LIB>
    $<TARGET_OBJECTS:ED_GLOB_LIB>
    $<TARGET_OBJECTS:ED_IO_LIB>
    $<TARGET_OBJECTS:ED_LIB_LIB>
  )
ENDIF(WITH_INEQ)





####################################################
# EDIpack Fortran-C interface
####################################################
ADD_SUBDIRECTORY(${LIB_CBIND} )
ADD_DEPENDENCIES(EDI_C_LIB ${EDI})
ADD_LIBRARY(${EDI_C} SHARED  $<TARGET_OBJECTS:EDI_C_LIB>)
if(CMAKE_SYSTEM_NAME STREQUAL "Linux")
  TARGET_LINK_OPTIONS(${EDI_C} PRIVATE -Wl,-z,noexecstack)
endif()
IF(WITH_INEQ)
  target_compile_definitions(EDI_C_LIB PRIVATE _WINEQ)
  target_compile_definitions(${EDI_C} PRIVATE _WINEQ)
ENDIF()
TARGET_LINK_LIBRARIES(${EDI_C} PRIVATE ${EDI})
TARGET_LINK_LIBRARIES(${EDI_C} PRIVATE ${SCIFOR_LDFLAGS})

#_LINK_LIBRARIES seems to not include absolute paths. This created a bug because
#if scalapack is not found blas and lapack are not linked. Using _LDFLAGS fixes it.
#TARGET_LINK_LIBRARIES(${EDI_C} PRIVATE ${SCIFOR_LINK_LIBRARIES})


############################################################
############################################################
############################################################

SET(PERMISSION_777
  OWNER_WRITE OWNER_READ OWNER_EXECUTE
  GROUP_WRITE GROUP_READ GROUP_EXECUTE
  WORLD_WRITE WORLD_READ WORLD_EXECUTE)





#Final messages + additional custom targets:
# - distclean: cleans the generated cmake files
# - uninstall: provides uninstall of (almost) all the installed files
# - doc: creates documentation using sphinx (if found)
# - test: compiles and launches the test in ${LIB_TEST} directory 
INCLUDE(${CMAKE_MODULE_PATH}/PostConfig.cmake)

