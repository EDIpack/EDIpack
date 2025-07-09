# EDipack C bindings
FILE(GLOB EDI_C_SRC
  EDIPACK_C_COMMON.f90
  EDIPACK_C.f90  
)
ADD_LIBRARY(EDI_C_LIB OBJECT  ${EDI_C_SRC})

set(HEADER_OUTPUT_DIR ${CMAKE_BINARY_DIR}/include/edipack)
configure_file(${LIB_CBIND}/edipack_cbindings.h
               ${HEADER_OUTPUT_DIR}/edipack_cbindings.h)
