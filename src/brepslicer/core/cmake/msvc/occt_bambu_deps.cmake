# Open CASCADE: bundled under third_party/occt (headers+libs) and
# third_party/occt_runtime (DLLs). Falls back to BambuStudio-deps when absent.
set(_occt_bundled "${CMAKE_SOURCE_DIR}/third_party/occt")
if(NOT OCCT_BAMBU_ROOT)
    if(EXISTS "${_occt_bundled}/include/occt/Standard.hxx")
        set(OCCT_BAMBU_ROOT "${_occt_bundled}")
    else()
        set(OCCT_BAMBU_ROOT "E:/src/BambuStudio-deps/usr/local")
    endif()
endif()
set(OCCT_BAMBU_ROOT "${OCCT_BAMBU_ROOT}" CACHE PATH "OCCT headers/libs root (include/occt, lib/occt)")
set(OCCT_RUNTIME_DIR "${CMAKE_SOURCE_DIR}/third_party/occt_runtime" CACHE PATH "Packed OCCT runtime DLLs")

set(OCCT750_INCLUDE_DIR "${OCCT_BAMBU_ROOT}/include/occt")
set(OCCT750_LIBRARY_DIR "${OCCT_BAMBU_ROOT}/lib/occt")

if(EXISTS "${OCCT_RUNTIME_DIR}/TKernel.dll")
    set(OCCT750_BINARY_DIR "${OCCT_RUNTIME_DIR}")
else()
    set(OCCT750_BINARY_DIR "${OCCT_BAMBU_ROOT}/bin/occt")
endif()

if(NOT EXISTS "${OCCT750_INCLUDE_DIR}/Standard.hxx")
    message(FATAL_ERROR "OCCT headers not found at ${OCCT750_INCLUDE_DIR}")
endif()
if(NOT EXISTS "${OCCT750_LIBRARY_DIR}/TKernel.lib")
    message(FATAL_ERROR "OCCT libraries not found at ${OCCT750_LIBRARY_DIR}")
endif()
if(NOT EXISTS "${OCCT750_BINARY_DIR}/TKernel.dll")
    message(FATAL_ERROR "OCCT runtime DLLs not found at ${OCCT750_BINARY_DIR}")
endif()

file(GLOB OCCT750_LIBS "${OCCT750_LIBRARY_DIR}/*.lib")

add_library(OCCT750 INTERFACE)
target_include_directories(OCCT750 INTERFACE ${OCCT750_INCLUDE_DIR})
target_link_directories(OCCT750 INTERFACE ${OCCT750_LIBRARY_DIR})
target_link_libraries(OCCT750 INTERFACE ${OCCT750_LIBS})

list(APPEND BIN_3RD_LIST ${OCCT750_BINARY_DIR})

list(LENGTH OCCT750_LIBS _occt_lib_count)
message(STATUS "OCCT include: ${OCCT750_INCLUDE_DIR}")
message(STATUS "OCCT libs:    ${OCCT750_LIBRARY_DIR} (${_occt_lib_count} .lib)")
message(STATUS "OCCT runtime: ${OCCT750_BINARY_DIR}")
